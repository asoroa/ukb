#include "common.h"
#include "globalVars.h"
#include "configFile.h"
#include "fileElem.h"
#include "kbGraph.h"
#include "disambGraph.h"

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
#include <boost/pending/integer_range.hpp>

using namespace std;
using namespace boost;
using namespace ukb;

const char *kb_default_binfile = "kb_wnet.bin";

static int filter_nodes = 0; // 0 -> no filter
                             // 1 -> only words
                             // 2 -> only synsets

static bool insert_all_dict = true;

static bool dict_weight = false; // Use W when linking words to concepts

// - sort all concepts according to their ppv weight, then scan the
// resulting sequence of concetps with a sliding window of length 100,
// and truncate the sequence when the difference in scores between the
// first and last concepts in the window drops below 5% of the
// highest-scoring concept for this word


struct CWSort {

  CWSort(const vector<double> & _v) : v(_v) {}
  int operator () (const int & i, const int & j) {
	// Descending order
	return v[i] > v[j];
  }
  const vector<double> & v;
};

void truncate_ppv(vector<double> & ppv, float thres) {

  size_t n = ppv.size();
  if (n < 100) return;

  vector<int> idx(n);

  for(size_t i=0; i < n; ++i)
	idx[i] = i;
  sort(idx.begin(), idx.end(), CWSort(ppv));

  size_t cut_i = 99;
  double cut_th = ppv[idx[0]]*thres;
  for(; cut_i < n; ++cut_i) {
	if ((ppv[idx[cut_i-99]] - ppv[idx[cut_i]]) < cut_th) break;
  }

  // truncate ppv
  for(; cut_i < n; ++cut_i) {
	ppv[idx[cut_i]] = 0.0;
  }

  // Normalize result

  normalize_pvector(ppv);

}

// Fill with zero's all values except top k

void top_k(vector<double> & ppv, size_t k) {

  size_t n = ppv.size();
  if (k >= n) return;

  vector<int> idx(n);

  for(size_t i=0; i < n; ++i)
	idx[i] = i;
  sort(idx.begin(), idx.end(), CWSort(ppv));

  for(size_t i = k; i < n; ++i) {
	ppv[idx[i]] = 0.0;
  }
  normalize_pvector(ppv);
}

void compute_sentence_vectors(string & fullname_in,
							  string & out_dir,
							  float trunc_ppv,
							  bool nozero) {

  Kb & kb = Kb::instance();
  if (glVars::verbose)
    cerr << "Reading words ...\n";

  ifstream fh_in(fullname_in.c_str());
  File_elem fout(fullname_in, out_dir, ".ppv");

  if (!fh_in) {
    cerr << "Can't open " << fullname_in << endl;
    exit(-1);
  }

  vector<CSentence> vcs;
  CSentence cs;

  // add words into Kb

  if (insert_all_dict) {
	if (glVars::verbose)
	  cerr << "Adding words to Kb ...";
	kb.add_dictionary(dict_weight);
  }

  if (glVars::verbose)
    Kb::instance().display_info(cerr);

  // Read sentences and compute rank vectors

  try {
    while (cs.read_aw(fh_in)) {

      // Initialize rank vector
      vector<double> ranks;

	  if(!insert_all_dict) {
		// Add CSentence words to graph
		CSentence::iterator it = cs.begin();
		CSentence::iterator end = cs.end();
 		for(;it != end; ++it) {
		  kb.add_token(it->word(), dict_weight);
		}
	  }
      bool ok = calculate_kb_ppr(cs,ranks);
      if (!ok) {
		cerr << "Error when calculating ranks for csentence " << cs.id() << endl;
		continue;
      }

      fout.fname = cs.id();

      ofstream fo(fout.get_fname().c_str(),  ofstream::out);
      if (!fo) {
		cerr << "Error: can't create" << fout.get_fname() << endl;
		exit(-1);
      }

	  vector<double> outranks;
	  vector<string> vnames;

	  Kb::instance().filter_ranks_vnames(ranks, outranks, vnames, filter_nodes);

	  if (trunc_ppv > 0.0f) {
		if (trunc_ppv < 1.0f)
		  truncate_ppv(outranks, trunc_ppv);
		else {
		  // For top k calculation
		  //
		  //  - fill with zeros all values except top k
		  //  - set nozero = 1 so only top k are printed
		  //
		  // * could be a problem if top k had zeros in it, as they
		  //   will not be properly printed.

		  top_k(outranks, lexical_cast<size_t>(trunc_ppv));
		  nozero = true;
		}
	  }

	  for(size_t i = 0; i < outranks.size(); ++i) {
		if (nozero && outranks [i] == 0.0) continue;
		fo << vnames[i] << "\t" << outranks[i] << "\n";
	  }
      cs = CSentence();
    }
  } catch (string & e) {
    cerr << "Errore reading " << fullname_in << ":" << e << "\n";
    throw(e);
  }
}

void compute_static_ppv() {

  vector<CSentence> vcs;
  CSentence cs;

  // Calculate static (static) pageRank over KB
  const vector<double> & ranks = Kb::instance().static_prank();

  vector<double> outranks;
  vector<string> vnames;

  Kb::instance().filter_ranks_vnames(ranks, outranks, vnames, 2);

  for(size_t i = 0; i < outranks.size(); ++i) {
	cout << vnames[i] << "\t" << outranks[i] << "\n";
  }

}

int main(int argc, char *argv[]) {

  srand(3);

  timer load;

  bool opt_static = false;
  bool opt_nozero = false;
  float opt_trppv = 0.0f;

  string kb_binfile(kb_default_binfile);

  string out_dir(".");
  string fullname_in;

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
    ("kb_binfile,K", value<string>(), "Binary file of KB (see compile_kb). Default is kb_wnet.bin")
    ("only_ctx_words,C", "Insert only words appearing in contexts to the graph (default is insert all dictionary words).")
    ("dict_file,D", value<string>(), "Word to synset map file (default is dict.txt.")
    ("dict_weight", "Use weights when linking words to concepts (dict file has to have weights). Also sets --prank_weight.")
    ("out_dir,O", value<string>(), "Directory for leaving output PPV files. Default is current directory.")
    ("concepts_in", "Let concept ids in input context. Item must have 5 fields, the fourth being 2 and the last one being the weight.")
    ("static,S", "Compute static PageRank ppv. Only -K option is needed. Output to STDOUT.")
    ("verbose,v", "Be verbose.")
	("nopos", "Don't filter words by Part of Speech.")
	("wiki", "Usual options for wikipedia (sets --nopos, --only_ctx_words and --only_synsets).")
    ;

  options_description po_desc_prank("pageRank general options");
  po_desc_prank.add_options()
    ("prank_weight,w", "Use weigths in pageRank calculation. Serialized graph edges must have some weight.")
    ("prank_iter", value<size_t>(), "Number of iterations in pageRank (good value is 30).")
    ("prank_threshold", value<float>(), "Threshold for pageRank convergence. Default is 0.0001.")
    ;

  options_description po_desc_output("Output options");
  po_desc_output.add_options()
    ("only_words", "Output only (normalized) PPVs for words.")
    ("only_synsets", "Output only (normalized) PPVs for synsets.")
    ("trunc_ppv", value<float>(), "Truncate PPV threshold (a la gabrilovich). If arg > 1, return top arg nodes.")
    ("nozero", "Do not return concepts with zero rank.")
    ;

  options_description po_hidden("Hidden");
  po_hidden.add_options()
    ("input-file",value<string>(), "Input file.")
    ("output-file",value<string>(), "Output file.")
    ;

  options_description po_visible(desc_header);
  po_visible.add(po_desc).add(po_desc_prank).add(po_desc_output);

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
    }

    if (vm.count("kb_binfile")) {
      kb_binfile = vm["kb_binfile"].as<string>();
    }

    if (vm.count("nopos")) {
	  glVars::input::filter_pos = false;
    }

    if (vm.count("only_words")) {
      filter_nodes = 1;
    }

    if (vm.count("only_ctx_words")) {
      insert_all_dict = false;
    }

    if (vm.count("only_synsets")) {
      filter_nodes = 2;
    }

    if (vm.count("out_dir")) {
      out_dir = vm["out_dir"].as<string>();
    }

    if (vm.count("dict_file")) {
      glVars::dict_filename = vm["dict_file"].as<string>();
    }

    if (vm.count("prank_weight")) {
	  glVars::prank::use_weight = true;
    }

    if (vm.count("dict_weight")) {
      dict_weight = true;
      glVars::prank::use_weight = true;
    }

    if (vm.count("static")) {
	  opt_static=true;
    }

    if (vm.count("prank_iter")) {
	  size_t iter = vm["prank_iter"].as<size_t>();
	  if (iter == 0) {
		cerr << "Error: prank_iter can not be zero!\n";
		goto END;
	  }
      glVars::prank::num_iterations = iter;
      glVars::prank::threshold = 0.0;
    }

    if (vm.count("prank_threshold")) {
	  float th = vm["prank_threshold"].as<float>();
	  if (th <= 0.0) {
		cerr << "Error: invalid prank_threshold value " << th << "\n";
		goto END;
	  }
      glVars::prank::threshold = th;
    }

    if (vm.count("concepts_in")) {
      glVars::csentence::concepts_in = true;
    }

    if (vm.count("wiki")) {
	  glVars::input::filter_pos = false;
	  insert_all_dict = false;
	  filter_nodes = 2;
    }

    if (vm.count("trunc_ppv")) {
      opt_trppv = vm["trunc_ppv"].as<float>();
    }

    if (vm.count("nozero")) {
      opt_nozero = 1;
    }

    if (vm.count("input-file")) {
      fullname_in = vm["input-file"].as<string>();
    }

  }
  catch(std::exception& e) {
    cerr << e.what() << "\n";
    throw(e);
  }

  if (opt_static) {
	if (glVars::verbose)
	  cerr << "Reading binary kb file " << kb_binfile;
	Kb::create_from_binfile(kb_binfile);
	if (glVars::verbose)
	  Kb::instance().display_info(cerr);
	compute_static_ppv();
	goto END;
  }

  if(fullname_in.size() == 0) {
    cout << po_visible << endl;
    return 1;
  }

  if (glVars::verbose)
    cerr << "Reading binary kb file " << kb_binfile;
  Kb::create_from_binfile(kb_binfile);
  if (glVars::verbose)
    Kb::instance().display_info(cerr);

  compute_sentence_vectors(fullname_in, out_dir, opt_trppv, opt_nozero);

 END:
  return 0;
}

/*
 * Local Variables:
 * mode: c++
 * compile-command: "make"
 * End:
 */
