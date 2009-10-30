#include "wdict.h"
#include "common.h"
#include "fileElem.h"
#include "globalVars.h"
#include "kbGraph.h"
#include "disambGraph.h"
#include <string>
#include <iostream>
#include <fstream>

// Basename & friends
#include <boost/filesystem/operations.hpp>
#include "boost/filesystem/path.hpp"

// Program options

#include <boost/program_options.hpp>

using namespace ukb;
using namespace std;
using namespace boost;

/////////////////////////////////////////////////////////////
// Global variables


const char *kb_default_binfile = "kb_wnet.bin";

static bool dict_weight = false; // Use W when linking words to concepts


// Program options stuff

/* Auxiliary functions for checking input for validity. */

/* Function used to check that 'opt1' and 'opt2' are not specified
   at the same time. */
void conflicting_options(const boost::program_options::variables_map& vm,
                         const char* opt1, const char* opt2)
{
  if (vm.count(opt1) && !vm[opt1].defaulted()
      && vm.count(opt2) && !vm[opt2].defaulted())
    throw logic_error(string("Conflicting options '")
					  + opt1 + "' and '" + opt2 + "'.");
}

/* Function used to check that of 'for_what' is specified, then
   'required_option' is specified too. */
void option_dependency(const boost::program_options::variables_map& vm,
					   const char* for_what, const char* required_option)
{
  if (vm.count(for_what) && !vm[for_what].defaulted())
    if (vm.count(required_option) == 0 || vm[required_option].defaulted())
      throw logic_error(string("Option '") + for_what
						+ "' requires option '" + required_option + "'.");
}

///////////////////////////////////////

// Main functions

// Disambiguate using disambiguation graph (dgraph) method

void disamb_dgraph_from_corpus(string & fullname_in,
							   bool out_semcor) {

  ifstream fh_in(fullname_in.c_str());

  if (!fh_in) {
    cerr << "Can't open " << fullname_in << endl;
    exit(-1);
  }

  CSentence cs;
  size_t l_n = 0;

  try {
    while (cs.read_aw(fh_in, l_n)) {
      DisambGraph dgraph;
      fill_disamb_graph(cs, dgraph);
	  pageRank_disg(dgraph.graph());
	  disamb_csentence(cs, dgraph);
      if (out_semcor) cs.print_csent_semcor_aw(cout);
      else cs.print_csent_simple(cout);
      cs = CSentence();
    }
  }
  catch (std::exception & e) {
    cerr << "Errore reading " << fullname_in << " : " << e.what() << "\n";
    throw(e);
  }
}

void dis_csent_ppr(const string & input_file,
				   bool out_semcor) {

  ifstream fh_in(input_file.c_str());

  if (!fh_in) {
    cerr << "Can't open " << input_file << endl;
    exit(-1);
  }

  CSentence cs;

  if (glVars::verbose)
    cerr << "Adding words to Mcr ...\n";

  Kb::instance().add_dictionary(dict_weight);
  size_t l_n = 0;

  try {
    while (cs.read_aw(fh_in, l_n)) {

      vector<float> ranks;
      bool ok = calculate_kb_ppr(cs,ranks);
      if (!ok) {
		cerr << "Error when calculating ranks for sentence " << cs.id() << "\n";
		cerr << "(No word links to KB ?)\n";
		continue;
      }

      disamb_csentence_kb(cs, ranks);
      if (out_semcor) cs.print_csent_semcor_aw(cout);
      else cs.print_csent_simple(cout);
      cs = CSentence();
    }
  }
  catch (std::exception & e) {
    cerr << "Errore reading " << input_file << " : " << e.what() << "\n";
    throw(e);
  }
}


void dis_csent_ppr_by_word(const string & input_file,
						  bool out_semcor) {

  ifstream fh_in(input_file.c_str());

  if (!fh_in) {
    cerr << "Can't open " << input_file << endl;
    exit(-1);
  }

  CSentence cs;

  if (glVars::verbose)
    cerr << "Adding words to Kb ...\n";

  Kb::instance().add_dictionary(dict_weight);
  size_t l_n = 0;

  try {
    while (cs.read_aw(fh_in, l_n)) {

      calculate_kb_ppr_by_word_and_disamb(cs);
      if (out_semcor) cs.print_csent_semcor_aw(cout);
      else cs.print_csent_simple(cout);

      //cout << cs << '\n';
      cs = CSentence();
    }
  }
  catch (std::exception & e) {
    cerr << "Errore reading " << input_file << " : " << e.what() << "\n";
    throw(e);
  }
}

void dis_csent_classic_prank(const string & input_file,
							 bool out_semcor) {

  ifstream fh_in(input_file.c_str());

  if (!fh_in) {
    cerr << "Can't open " << input_file << endl;
    exit(-1);
  }

  CSentence cs;
  size_t l_n = 0;
  const vector<float> ranks = Kb::instance().static_prank();
  try {
    while (cs.read_aw(fh_in, l_n)) {
      disamb_csentence_kb(cs, ranks);
      if (out_semcor) cs.print_csent_semcor_aw(cout);
      else cs.print_csent_simple(cout);
      cs = CSentence();
    }
  }
  catch (std::exception & e) {
    cerr << "Errore reading " << input_file << " : " << e.what() << "\n";
    throw(e);
  }
}


void test(const string & input_file,
		  bool out_semcor) {

  ifstream fh_in(input_file.c_str());

  if (!fh_in) {
    cerr << "Can't open " << input_file << endl;
    exit(-1);
  }

  CSentence cs;

  vector<float> ranks;
  //Kb::instance().indegree_rank(ranks);
  size_t l_n = 0;
  try {
    while (cs.read_aw(fh_in, l_n)) {
      disamb_csentence_kb(cs, ranks);
      if (out_semcor) cs.print_csent_semcor_aw(cout);
      else cs.print_csent_simple(cout);
      cs = CSentence();
    }
  }
  catch (std::exception & e) {
    cerr << "Errore reading " << input_file << " : " << e.what() << "\n";
    throw(e);
  }
}

int main(int argc, char *argv[]) {

  string kb_binfile(kb_default_binfile);

  enum dis_method {
	dgraph,
	ppr,
	ppr_w2w,
	ppr_static
  };

  dis_method dmethod = ppr;

  bool opt_do_test = false;
  bool opt_out_semcor = false;

  string cmdline("!! -v ");
  cmdline += glVars::ukb_version;
  for (int i=0; i < argc; ++i) {
    cmdline += " ";
    cmdline += argv[i];
  }

  vector<string> input_files;
  string fullname_in;

  using namespace boost::program_options;

  const char desc_header[] = "ukb_wsd: perform WSD with KB based algorithm\n"
    "Usage examples:\n"
	"ukb_wsd -D dict.txt -K kb.bin --ppr input.txt    -> Disambiguate input.txt using PPR technique according to graph kb.bin and dictionary dict.txt\n"
    "ukb_wsd -D dict.txt -K kb.bin --dis_dgraph input.txt -> Disambiguate input.txt using Disambiguation graph technique, according to graph kb.bin and dictionary dict.txt\n"
	"Options";

  //options_description po_desc(desc_header);

  options_description po_desc("General");

  po_desc.add_options()
    ("help,h", "This page")
    ("version", "Show version.")
    ("kb_binfile,K", value<string>(), "Binary file of KB (see compile_kb). Default is kb_wnet.bin.")
    ("dict_file,D", value<string>(), "Dictionary text file. Default is dict.txt")
    ("dict_weight", "Use weights when linking words to concepts (dict file has to have weights). Also sets --prank_weight.")
	("nopos", "Don't filter words by Part of Speech.")
    ;

  options_description po_desc_wsd("WSD methods");
  po_desc_wsd.add_options()
    ("ppr", "Given a text input file, disambiguate context using Personalized PageRank method.")
    ("ppr_w2w", "Given a text input file, disambiguate context using Personalized PageRank method word by word (see README).")
    ("static", "Given a text input file, disambiguate context using static pageRank over kb.")
    ("dis_dgraph", "Given a text input file, disambiguate context using disambiguation graph mehod.")
    ("nostatic", "Substract static ppv to final ranks.")
    ;

  options_description po_desc_prank("pageRank general options");
  po_desc_prank.add_options()
    ("prank_weight,w", "Use weigths in pageRank calculation. Serialized graph edges must have some weight.")
    ("prank_iter", value<size_t>(), "Number of iterations in pageRank. Default is 30.")
    ("prank_threshold", value<float>(), "Threshold for stopping PageRank. Default is zero. Good value is 0.0001.")
    ("prank_damping", value<float>(), "Set damping factor in PageRank equation. Default is 0.85.")
    ;

  options_description po_desc_output("Output options");
  po_desc_output.add_options()
    ("semcor", "Output Semcor key file.")
    ("allranks", "Write key file with all synsets associated with ranks.")
    ("verbose,v", "Be verbose.")
    ("no-monosemous", "Don't output anything for monosemous words.")
    ;

  options_description po_visible(desc_header);
  po_visible.add(po_desc).add(po_desc_wsd).add(po_desc_prank).add(po_desc_output);

  options_description po_hidden("Hidden");
  po_hidden.add_options()
    ("bcomp_kb_binfile,M", value<string>(), "Backward compatibility with -K.")
    ("bcomp_dictfile,W", value<string>(), "Backward compatibility with -D.")
    ("test,t", "(Internal) Do a test.")
    ("test,t", "(Internal) Do a test.")
    ("input-file",value<string>(), "Input file.")
    ;
  options_description po_all("All options");
  po_all.add(po_visible).add(po_hidden);

  positional_options_description po_optdesc;
  po_optdesc.add("input-file", 1);
  //    po_optdesc.add("output-file", 1);

  try {

    variables_map vm;
    store(command_line_parser(argc, argv).
		  options(po_all).
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

    if (vm.count("nopos")) {
	  glVars::input::filter_pos = false;
    }

    if (vm.count("ppr")) {
	  dmethod = ppr;
    }

    if (vm.count("ppr_w2w")) {
	  dmethod = ppr_w2w;
    }

    if (vm.count("static")) {
	  dmethod = ppr_static;
    }

    if (vm.count("dis_dgraph")) {
	  dmethod = dgraph;
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
	  if (th <= 0.0 || th > 1.0) {
		cerr << "Error: invalid prank_threshold value " << th << "\n";
		goto END;
	  }
      glVars::prank::threshold = th;
      glVars::prank::num_iterations = 0;
    }

    if (vm.count("prank_damping")) {
	  float dp = vm["prank_damping"].as<float>();
	  if (dp <= 0.0 || dp > 1.0) {
		cerr << "Error: invalid prank_damping value " << dp << "\n";
		goto END;
	  }
      glVars::prank::damping = dp;
    }

    if (vm.count("nostatic")) {
	  glVars::csentence::disamb_minus_static = true;
    }

    if (vm.count("bcomp_dictfile")) {
      glVars::dict_filename = vm["bcomp_dictfile"].as<string>();
    }

    if (vm.count("dict_file")) {
      glVars::dict_filename = vm["dict_file"].as<string>();
    }

    if (vm.count("dict_weight")) {
      dict_weight = true;
      glVars::prank::use_weight = true;
    }

    if (vm.count("bcomp_kb_binfile")) {
      kb_binfile = vm["bcomp_kb_binfile"].as<string>();
    }

    if (vm.count("kb_binfile")) {
      kb_binfile = vm["kb_binfile"].as<string>();
    }

    if (vm.count("rank_alg")) {
      glVars::RankAlg alg = glVars::get_algEnum(vm["rank_alg"].as<string>());
      if (alg == glVars::no_alg) {
		cerr << "Error: Undefined rank algorithm " << vm["rank_alg"].as<string>() << endl;
		exit(-1);
      }
      glVars::rAlg = alg;
    }

    if (vm.count("prank_weight")) {
	  glVars::prank::use_weight = true;
    }


    if (vm.count("input-file")) {
      fullname_in = vm["input-file"].as<string>();
    }

    if (vm.count("verbose")) {
      glVars::verbose = 1;
      glVars::debug::warning = 1;
    }

    if (vm.count("allranks")) {
      glVars::output::allranks = true;
    }

    if (vm.count("semcor")) {
      opt_out_semcor = 1;
    }

    if (vm.count("test")) {
      opt_do_test = true;
    }

    if (vm.count("no-monosemous")) {
      glVars::output::monosemous = false;
    }
  }
  catch(std::exception& e) {
    cerr << e.what() << "\n";
    throw(e);
  }

  if (!fullname_in.size()) {
    cout << po_visible << endl;
    cout << "Error: No input" << endl;
    exit(-1);
  }

  if (opt_do_test) {
    test(fullname_in, false);
	goto END;
  }


  Kb::create_from_binfile(kb_binfile);
  cout << cmdline << "\n";

  switch(dmethod) {

  case dgraph:
    disamb_dgraph_from_corpus(fullname_in, opt_out_semcor);
    goto END;
	break;
  case ppr:
    dis_csent_ppr(fullname_in, opt_out_semcor);
    goto END;
	break;

  case ppr_w2w:
    dis_csent_ppr_by_word(fullname_in, opt_out_semcor);
    goto END;
	break;

  case ppr_static:
    dis_csent_classic_prank(fullname_in, opt_out_semcor);
    goto END;
	break;

  };

 END:
  return 0;
}
