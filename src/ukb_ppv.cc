#include "common.h"
#include "globalVars.h"
#include "configFile.h"
#include "fileElem.h"
#include "kbGraph.h"
#include "disambGraph.h"
#include "wdict.h"

#include <string>
#include <iostream>
#include <fstream>

// Program options

#include <boost/program_options.hpp>

// timer

#include <boost/timer.hpp>

// bfs

#include <boost/graph/breadth_first_search.hpp>
#include <boost/pending/indirect_cmp.hpp>

#if BOOST_VERSION > 104400
  #include <boost/range/irange.hpp>
#else
  #include <boost/pending/integer_range.hpp>
#endif

using namespace std;
using namespace boost;
using namespace ukb;

const char *kb_default_binfile = "kb_wnet.bin";

static bool opt_normalize_ranks = true;

static bool output_control_line = false;
static bool output_variants_ppv = false;
static float trunc_ppv = 0.0;
static bool nozero = false;
static string ppv_prefix;
static string cmdline("!! -v ");

// - sort all concepts according to their ppv weight, then scan the
// resulting sequence of concetps with a sliding window of length 100,
// and truncate the sequence when the difference in scores between the
// first and last concepts in the window drops below 5% of the
// highest-scoring concept for this word


struct CWSort {

  CWSort(const vector<float> & _v) : v(_v) {}
  int operator () (const int & i, const int & j) {
	// Descending order
	return v[i] > v[j];
  }
  const vector<float> & v;
};

void truncate_ppv(vector<float> & ppv, float thres) {

  size_t n = ppv.size();
  if (n < 100) return;

  vector<int> idx(n);

  for(size_t i=0; i < n; ++i)
	idx[i] = i;
  sort(idx.begin(), idx.end(), CWSort(ppv));

  size_t cut_i = 99;
  float cut_th = ppv[idx[0]]*thres;
  for(; cut_i < n; ++cut_i) {
	if ((ppv[idx[cut_i-99]] - ppv[idx[cut_i]]) < cut_th) break;
  }

  // truncate ppv
  for(; cut_i < n; ++cut_i) {
	ppv[idx[cut_i]] = 0.0;
  }

  // Normalize result

  if(opt_normalize_ranks)
	normalize_pvector(ppv);

}

// Fill with zero's all values except top k

void top_k(vector<float> & ppv, size_t k) {

  size_t n = ppv.size();
  if (k >= n) return;

  vector<int> idx(n);

  for(size_t i=0; i < n; ++i)
	idx[i] = i;
  sort(idx.begin(), idx.end(), CWSort(ppv));

  for(size_t i = k; i < n; ++i) {
	ppv[idx[i]] = 0.0;
  }
  if(opt_normalize_ranks)
	normalize_pvector(ppv);
}


static void output_ppv_stream(const vector<float> & outranks, ostream & os) {

  Kb & kb = Kb::instance();

  vector<float> ranks_trunc;
  const vector<float> * theranks = &outranks;

  if (trunc_ppv > 0.0f) {
	vector<float>(outranks).swap(ranks_trunc);
	theranks = &ranks_trunc;
	if (trunc_ppv < 1.0f)
	  truncate_ppv(ranks_trunc, trunc_ppv);
	else {
	  // For top k calculation
	  //
	  //  - fill with zeros all values except top k
	  //  - set nozero = 1 so only top k are printed
	  //
	  // * could be a problem if top k had zeros in it, as they
	  //   will not be properly printed.

	  top_k(ranks_trunc, lexical_cast<size_t>(trunc_ppv));
	  nozero = true;
	}
  }

  if (output_control_line)
	os << cmdline << "\n";
  for(size_t i = 0; i < theranks->size(); ++i) {
	string sname = kb.get_vertex_name(i);
	if (nozero && (*theranks) [i] == 0.0) continue;
	os << sname << "\t" << (*theranks)[i];
	if (output_variants_ppv) {
	  os << "\t" << WDict::instance().variant(sname);
	}
	os << "\n";
  }
}

static void create_output_ppv(vector<float> & ranks,
							  const string & filename,
							  File_elem & fout) {
  fout.fname = filename;

  ofstream fo(fout.get_fname().c_str(),  ofstream::out);
  if (!fo) {
	cerr << "Error: can't create" << fout.get_fname() << endl;
	exit(-1);
  }
  output_ppv_stream(ranks, fo);
}

static void maybe_postproc_ranks(vector<float> & ranks) {
  if (glVars::csentence::disamb_minus_static) {
	const vector<float> & static_ranks = Kb::instance().static_prank();
	for(size_t s_i = 0, s_m = static_ranks.size();
		s_i != s_m; ++s_i) {
	  ranks[s_i] -= static_ranks[s_i];
	}
  }
}

void compute_sentence_vectors(string & out_dir) {

  File_elem fout("lala", out_dir, ".ppv");

  vector<CSentence> vcs;
  CSentence cs;

  // Read sentences and compute rank vectors
  size_t l_n  = 0;
  while (cs.read_aw(std::cin, l_n)) {
	// Initialize rank vector
	vector<float> ranks;
	bool ok = calculate_kb_ppr(cs,ranks);
	if (!ok) {
	  cerr << "Error when calculating ranks for csentence " << cs.id() << endl;
	  continue;
	}
	maybe_postproc_ranks(ranks);
	create_output_ppv(ranks, ppv_prefix + cs.id(), fout);
	cs = CSentence();
  }
}

void compute_sentence_vectors_w2w(string & out_dir) {

  File_elem fout("lala", out_dir, ".ppv");

  vector<CSentence> vcs;
  CSentence cs;

  // Read sentences and compute rank vectors
  size_t l_n  = 0;
  while (cs.read_aw(std::cin, l_n)) {
	vector<float> ranks;
	int w_n = 1;
	for(vector<CWord>::iterator cw_it = cs.begin(), cw_end = cs.end();
		cw_it != cw_end; ++cw_it, ++w_n) {
	  if(!cw_it->is_tgtword()) continue;
	  bool ok = calculate_kb_ppr_by_word(cs, cw_it, ranks);
	  if (!ok) {
		cerr << "Error when calculating ranks for word " << cw_it->wpos() << " in csentence " << cs.id() << endl;
		continue;
	  }
	  maybe_postproc_ranks(ranks);
	  string ofile = cs.id() + "#";
	  ofile += lexical_cast<string>(w_n);
	  create_output_ppv(ranks, ppv_prefix + ofile, fout);
	}
	cs = CSentence();
  }
}

void compute_static_ppv() {


  // Calculate static (static) pageRank over KB
  const vector<float> & ranks = Kb::instance().static_prank();

  output_ppv_stream(ranks, cout);
}

int main(int argc, char *argv[]) {

  srand(3);

  bool opt_static = false;

  string kb_binfile(kb_default_binfile);

  cmdline += glVars::ukb_version;
  for (int i=0; i < argc; ++i) {
    cmdline += " ";
    cmdline += argv[i];
  }

  string out_dir(".");
  string fullname_in;
  ifstream input_ifs;

  size_t iterations = 0;
  float thresh = 0.0;
  bool check_convergence = false;

  const char desc_header[] = "ukb_ppv: get personalized PageRank vector if a KB\n"
	"Usage examples:\n"
    "ukb_ppv -K kb.bin -D dict.txt -O outdir input.txt\n"
    "  Creates one file per sentence (.ppv extension) with the vector of the PPV vector given the input sentence"
    "Options";

  using namespace boost::program_options;

  options_description po_desc("General options:");

  po_desc.add_options()
    ("help,h", "This help page.")
    ("version", "Show version.")
    ("kb_binfile,K", value<string>(), "Binary file of KB (see compile_kb).")
    ("out_dir,O", value<string>(), "Directory for leaving output PPV files. Default is current directory.")
    ("static,S", "Compute static PageRank ppv. Only -K option is needed. Output to STDOUT.")
    ("dict_file,D", value<string>(), "Word to synset map file.")
    ("verbose,v", "Be verbose.")
    ;

  options_description po_desc_input("Input options");
  po_desc_input.add_options()
	("nopos", "Don't filter words by Part of Speech.")
	("minput", "Do not die when dealing with malformed input.")
	("ctx_noweight", "Do not use weights of input words (defaut is use context weights).")
    ;

  options_description po_desc_prank("pageRank general options");
  po_desc_prank.add_options()
    ("prank_weight,w", "Use weigths in pageRank calculation. Serialized graph edges must have some weight.")
    ("prank_iter", value<size_t>(), "Number of iterations in pageRank. Default is 30.")
    ("prank_threshold", value<float>(), "Threshold for pageRank convergence. Default is 0.0001.")
    ("prank_damping", value<float>(), "Set damping factor in PageRank equation. Default is 0.85.")
    ;

  options_description po_desc_dict("Dictionary options");
  po_desc_dict.add_options()
    ("dict_weight", "Use weights when linking words to concepts (dict file has to have weights).")
    ("dict_weight_smooth", value<float>(), "Smoothing factor to be added to every weight in dictionary concepts. Default is 1.")
    ("dict_strict", "Be strict when reading the dictionary and stop when any error is found.")
    ;

  options_description po_desc_output("Output options");
  po_desc_output.add_options()
    ("nostatic", "Substract static ppv to final ranks.")
    ("trunc_ppv", value<float>(), "Truncate PPV threshold. If arg > 1, return top arg nodes.")
    ("trunc_topK", value<size_t>(), "Return top arg nodes.")
    ("nozero", "Do not return concepts with zero rank.")
    ("variants,r", "Write also concept variants in PPV.")
    ("control_line,l", "First line in PPV files is control.")
	("prefix,p", value<string>(), "Prefix added to all output ppv files.")
    ("ranks_nonorm", "Do not normalize ranks even with topK or threshold cuts.")
    ;

  options_description po_hidden("Hidden");
  po_hidden.add_options()
    ("only_ctx_words,C", "Backward compatibility with -C.")
    ("only_synsets", "Output only (normalized) PPVs for synsets.")
    ("concept_graph,G", "Backward compatibility with -G.")
    ("input-file",value<string>(), "Input file.")
    ("output-file",value<string>(), "Output file.")
    ;

  options_description po_visible(desc_header);
  po_visible.add(po_desc).add(po_desc_input).add(po_desc_prank).add(po_desc_dict).add(po_desc_output);

  options_description po_desc_all("All options");
  po_desc_all.add(po_visible).add(po_hidden);

  positional_options_description po_optdesc;
  po_optdesc.add("input-file", 1);
  po_optdesc.add("output-file", 1);

  try {
    variables_map vm;
    store(command_line_parser(argc, argv).
	  options(po_desc_all).
	  positional(po_optdesc).
	  run(), vm);
    notify(vm);

    // If asked for help, don't do anything more

    if (vm.count("help")) {
      cout << po_visible << endl;
      exit(0);
    }

    if (vm.count("version")) {
      cout << glVars::ukb_version << endl;
      exit(0);
    }

    // verbosity

    if (vm.count("verbose")) {
      glVars::verbose = 1;
      glVars::debug::warning = 1;
    }

    if (vm.count("kb_binfile")) {
      kb_binfile = vm["kb_binfile"].as<string>();
    }

    if (vm.count("nopos")) {
	  glVars::input::filter_pos = false;
    }

    if (vm.count("minput")) {
	  glVars::input::swallow = true;
    }

    if (vm.count("ctx_noweight")) {
	  glVars::input::weight = false;
    }

    if (vm.count("variants")) {
	  output_variants_ppv = true;
    }

    if (vm.count("control_line")) {
	  output_control_line = true;
    }

    if (vm.count("out_dir")) {
      out_dir = vm["out_dir"].as<string>();
    }

    if (vm.count("prefix")) {
      ppv_prefix = vm["prefix"].as<string>();
    }

    if (vm.count("dict_file")) {
      glVars::dict_filename = vm["dict_file"].as<string>();
    }

    if (vm.count("prank_weight")) {
	  glVars::prank::use_weight = true;
    }

    if (vm.count("dict_strict")) {
      glVars::dict::swallow = false;
    }

    if (vm.count("dict_weight")) {
	  glVars::dict::use_weight = true;
      glVars::prank::use_weight = true;
    }

    if (vm.count("dict_weight_smooth")) {
      glVars::dict::weight_smoothfactor = vm["dict_weight_smooth"].as<float>();
    }

    if (vm.count("static")) {
	  opt_static=true;
    }

    if (vm.count("nostatic")) {
	  glVars::csentence::disamb_minus_static = true;
    }

    if (vm.count("prank_iter")) {
	  iterations = vm["prank_iter"].as<size_t>();
	  check_convergence = true;
    }

    if (vm.count("prank_threshold")) {
	  thresh = vm["prank_threshold"].as<float>();
	  check_convergence = true;
    }

    if (vm.count("prank_damping")) {
	  float dp = vm["prank_damping"].as<float>();
	  if (dp <= 0.0 || dp > 1.0) {
		cerr << "Error: invalid prank_damping value " << dp << "\n";
		goto END;
	  }
      glVars::prank::damping = dp;
    }

    if (vm.count("trunc_ppv")) {
      trunc_ppv = vm["trunc_ppv"].as<float>();
    }

    if (vm.count("trunc_topK")) {
	  size_t topK = vm["trunc_topK"].as<size_t>();
	  if(!topK) {
		cerr << "Error: trunc_topK is zero\n.";
		goto END;
	  }
      trunc_ppv = topK;
    }

    if (vm.count("nozero")) {
	  nozero = true;
    }

    if (vm.count("input-file")) {
      fullname_in = vm["input-file"].as<string>();
    }

  }
  catch(std::exception& e) {
    cerr << e.what() << "\n";
	exit(-1);
  }

  if (check_convergence) set_pr_convergence(iterations, thresh);

  if (opt_static) {
	if (glVars::verbose)
	  cerr << "Reading binary kb file " << kb_binfile;
	Kb::create_from_binfile(kb_binfile);
	if (glVars::verbose)
	  Kb::instance().display_info(cerr);
	compute_static_ppv();
	goto END;
  }

  if(!fullname_in.size()) {
    cout << po_visible << endl;
    cout << "Error: No input" << endl;
    exit(-1);
  }

  if (fullname_in == "-" ) {
	// read from <STDIN>
    cmdline += " <STDIN>";
	fullname_in = "<STDIN>";
  } else {
	input_ifs.open(fullname_in.c_str(), ofstream::in);
	if (!input_ifs) {
	  cerr << "Can't open " << fullname_in << endl;
	  exit(-1);
	}
	// redirect std::cin to read from file
	std::cin.rdbuf(input_ifs.rdbuf());
  }


  if (glVars::verbose)
    cerr << "Reading binary kb file " << kb_binfile;
  Kb::create_from_binfile(kb_binfile);
  if (glVars::verbose)
    Kb::instance().display_info(cerr);

  try {
	compute_sentence_vectors(out_dir);
  } catch(std::exception& e) {
    cerr << "Errore reading " << fullname_in << " : " << e.what() << "\n";
	exit(-1);
  }

 END:
  return 0;
}

/*
 * Local Variables:
 * mode: c++
 * compile-command: "make"
 * End:
 */
