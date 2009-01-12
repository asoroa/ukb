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
      vector<float> ranks;

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

      vector<float> filter_ranks;
      float filter_sum = 0.0;
      double norm;
      switch(filter_nodes) {
      case 1: // only words
		if (glVars::verbose) 
		  cerr << ranks.size() << "\n";
		for(size_t i=0; i < ranks.size(); ++i) {
		  if (kb.vertex_is_word(i)) {
			filter_ranks.push_back(ranks[i]);
			filter_sum+=ranks[i];
		  }
		}
		if (glVars::verbose) 
		  cerr << filter_ranks.size() << " words\n";
		// normalize vector
		norm = 1.0/filter_sum;
		for(vector<float>::iterator it=filter_ranks.begin(); it != filter_ranks.end(); ++it) {
		  fo << *it * norm << "\n";
		}
		break;
      case 2: // only synsets
		if (glVars::verbose) 
		  cerr << ranks.size() << " synsets\n";
		for(size_t i=0; i < ranks.size(); ++i) {
		  if (kb.vertex_is_synset(i)) {
			filter_ranks.push_back(ranks[i]);
			filter_sum+=ranks[i];
		  }
		}
		// normalize vector
		if (glVars::verbose) 
		  cerr << filter_ranks.size() << " synsets\n";
		norm = 1.0/filter_sum;
		for(vector<float>::iterator it=filter_ranks.begin(); it != filter_ranks.end(); ++it) {
		  fo << *it * norm << "\n";
		}
		break;
		
      default:
		copy(ranks.begin(), ranks.end(), ostream_iterator<float>(fo, "\n"));
		break;
      };
      cs = CSentence();
    }
  } catch (string & e) {
    cerr << "Errore reading " << fullname_in << ":" << e << "\n";
    throw(e);    
  }
}

int main(int argc, char *argv[]) {

  srand(3);

  timer load;

  bool opt_with_w = false;

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
    ("kb_binfile,K", value<string>(), "Binary file of KB (see compile_kb). Default is kb_wnet.bin")
    ("dict_file,D", value<string>(), "Word to synset map file (default is dict.txt.")
    ("out_dir,O", value<string>(), "Directory for leaving output PPV files.")
    ("word_weight", "Use weights on words (init values of PPV). Input must have 5 fields, the last one being the weight.")
    ("verbose,v", "Be verbose.")
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

    // verbosity

    if (vm.count("verbose")) {
      glVars::verbose = 1;
    }

    if (vm.count("kb_binfile")) {
      kb_binfile = vm["kb_binfile"].as<string>();
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

    if (vm.count("prank_iter")) {
      glVars::prank::num_iterations = vm["prank_iter"].as<size_t>();
    }

    if (vm.count("word_weight")) {
      glVars::csentence::word_weight = true;
    }

    if (vm.count("input-file")) {
      fullname_in = vm["input-file"].as<string>();
    }

  }
  catch(exception& e) {
    cerr << e.what() << "\n";
    throw(e);
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

  compute_sentence_vectors(fullname_in, out_dir, opt_with_w);
  return 0;
}

/*
 * Local Variables:
 * mode: c++
 * compile-command: "make"
 * End:
 */
