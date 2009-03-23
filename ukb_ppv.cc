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

void compute_sentence_vectors(string & fullname_in, string & out_dir, bool weight) {

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

  if (glVars::verbose) 
    cerr << "Adding words to Kb ...";

  kb.add_dictionary(false);

  if (glVars::verbose) 
    Kb::instance().display_info(cerr);

  // Read sentences and compute rank vectors

  try {
    while (cs.read_aw(fh_in)) {

      // Initialize rank vector
      vector<double> ranks;

      bool ok = calculate_kb_ppr(cs,ranks, weight);
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

	  for(size_t i = 0; i < outranks.size(); ++i) {
		fo << vnames[i] << "\t" << outranks[i] << "\n";
	  }
      cs = CSentence();
    }
  } catch (string & e) {
    cerr << "Errore reading " << fullname_in << ":" << e << "\n";
    throw(e);
  }
}

void compute_static_ppv(bool with_w) {

  vector<CSentence> vcs;
  CSentence cs;

  // Calculate static (static) pageRank over KB
  size_t N = Kb::instance().size();
  vector<double> ppv(N, 1.0/static_cast<double>(N));
  vector<double> ranks;
  Kb::instance().pageRank_ppv(ppv, ranks, with_w);

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

  bool opt_with_w = false;
  bool opt_static = false;

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
    ("dict_file,D", value<string>(), "Word to synset map file (default is dict.txt.")
    ("out_dir,O", value<string>(), "Directory for leaving output PPV files. Default is current directory.")
    ("concepts_in", "Let concept ids in input context. Item must have 5 fields, the fourth being 2 and the last one being the weight.")
    ("static,S", "Compute static PageRank ppv. Only -K option is needed. Output to STDOUT.")
    ("verbose,v", "Be verbose.")
	("nopos", "Don't filter words by Part of Speech.")
    ;

  options_description po_desc_prank("pageRank general options");
  po_desc_prank.add_options()
    ("with_weight,w", "Use weigths in pageRank calculation. Serialized graph edges must have some weight.")
    ("prank_iter", value<size_t>(), "Number of iterations in pageRank. Default is 30.")
    ;

  options_description po_desc_output("Output options");
  po_desc_output.add_options()
    ("only_words", "Output only (normalized) PPVs for words.")
    ("only_synsets", "Output only (normalized) PPVs for synsets.")
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

    if (vm.count("only_synsets")) {
      filter_nodes = 2;
    }

    if (vm.count("out_dir")) {
      out_dir = vm["out_dir"].as<string>();
    }

    if (vm.count("dict_file")) {
      glVars::dict_filename = vm["dict_file"].as<string>();
    }

    if (vm.count("with_weight")) {
      opt_with_w = true;
    }

    if (vm.count("static")) {
	  opt_static=true;
    }

    if (vm.count("prank_iter")) {
      glVars::prank::num_iterations = vm["prank_iter"].as<size_t>();
    }

    if (vm.count("concepts_in")) {
      glVars::csentence::concepts_in = true;
    }

    if (vm.count("input-file")) {
      fullname_in = vm["input-file"].as<string>();
    }

  }
  catch(exception& e) {
    cerr << e.what() << "\n";
    throw(e);
  }

  if (glVars::verbose)
    cerr << "Reading binary kb file " << kb_binfile;
  Kb::create_from_binfile(kb_binfile);
  if (glVars::verbose)
    Kb::instance().display_info(cerr);

  if (opt_static) {
	compute_static_ppv(opt_with_w);
	goto END;
  }

  if(fullname_in.size() == 0) {
    cout << po_visible << endl;
    return 1;
  }

  compute_sentence_vectors(fullname_in, out_dir, opt_with_w);

 END:
  return 0;
}

/*
 * Local Variables:
 * mode: c++
 * compile-command: "make"
 * End:
 */
