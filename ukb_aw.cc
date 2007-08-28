#include "w2syn.h"
#include "common.h"
#include "fileElem.h"
#include "globalVars.h"
#include "mcrGraph.h"
#include "kGraph.h"
#include "disambGraph.h"
#include <string>
#include <iostream>
#include <fstream>

// Basename & friends
#include <boost/filesystem/operations.hpp>
#include "boost/filesystem/path.hpp"

// Program options

#include <boost/program_options.hpp>

// timer

#include <boost/timer.hpp>

using namespace std;
using namespace boost;

/////////////////////////////////////////////////////////////
// Global variables


const char *mcr_default_binfile = "mcr_wnet.bin";

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



bool extract_input_files(const string & fullname,
			 vector<string> & input_files,
			 const string & extension) {

  namespace fs = boost::filesystem;

  fs::path full_path( fs::initial_path() );
  
  full_path = fs::system_complete( fs::path( fullname, fs::native ) );

  if ( !fs::exists( full_path ) )
  {
    std::cerr << "\nNot found: " << full_path.native_file_string() << std::endl;
    return false;
  }

  if ( fs::is_directory( full_path ) ) {

    fs::directory_iterator end_iter;
    for ( fs::directory_iterator dir_itr( full_path );
          dir_itr != end_iter;
          ++dir_itr ) {
      string dfile = dir_itr->native_file_string();
      size_t ext_i = dfile.find_last_of('.');
      if (ext_i == string::npos) continue;
      string dext(dfile.begin() + ext_i + 1, dfile.end());
      if (dext != extension) continue;
      input_files.push_back(dfile);
    }    
  } else {
    input_files.push_back(full_path.native_file_string());
  }
  return true;
}

///////////////////////////////////////

//Main program


void create_dgraphs_from_corpus(string & fullname_in,
				const string & out_dir,
				bool opt_create_kgraph = false) {

  File_elem fout(fullname_in, out_dir, ".dgraph");
  ifstream fh_in(fullname_in.c_str());

  if (!fh_in) {
    cerr << "Can't open " << fullname_in << endl;
    exit(-1);
  }

  vector<CSentence> vcs;
  CSentence cs;

  try {
    while (cs.read_aw(fh_in)) {
      DisambGraph dgraph;
      fill_disamb_graph(cs, dgraph);
      fout.fname=cs.id();
      fout.ext = ".dgraph";
      dgraph.write_to_binfile(fout.get_fname());
      fout.ext = ".csent";
      cs.write_to_binfile(fout.get_fname());
      if (opt_create_kgraph) {
	fout.ext="kgraph";
	KGraph kg(cs, dgraph);
	kg.write_to_binfile(fout.get_fname());
      }
      cs = CSentence();
    }
  } 
  catch (string & e) {
    cerr << "Errore reading " << fullname_in << ":" << e << "\n";
    throw(e);    
  }
}

void create_kgraph(string & cs_fname,
		   const string & out_dir) {

  File_elem cs_fe(cs_fname);
  File_elem dg_fe(cs_fe.fname, cs_fe.path, ".dgraph");
  CSentence cs;
  cs.read_from_binfile(cs_fe.get_fname());
  DisambGraph dg;
  dg.read_from_binfile(dg_fe.get_fname());

  File_elem kg_fe(cs_fe.fname, cs_fe.path, ".kgraph");
  kg_fe.set_path(out_dir);
  KGraph kg(cs, dg);
  kg.write_to_binfile(kg_fe.get_fname());
  pageRank(kg.graph());
  disamb_csentence(cs, kg);
  print_disamb_csent(cout, cs);  

  kg_fe.ext = ".dot";
  write_dgraph_graphviz(kg_fe.get_fname(), kg.graph());
}


template<class G>
void dis_csent(const vector<string> & input_files, const string & ext,
	       bool edge_weights) {
  

  map<string, size_t> counts;

  if (glVars::mcr_with_freqs) {
    bool ok = W2Syn::instance().syn_counts(counts);
    if (!ok) {
      cerr << "Error! There are no freqs. Check file " << glVars::w2s_filename << endl;
      exit(-1);
    }
  }

  for(vector<string>::const_iterator fname=input_files.begin(); fname != input_files.end(); ++fname) {

    File_elem cs_finfo(*fname);

    CSentence cs;
    cs.read_from_binfile(cs_finfo.get_fname());

    File_elem dg_finfo(cs.id(), cs_finfo.path, ext);
    G dg;
    dg.read_from_binfile(dg_finfo.get_fname());
    if (glVars::mcr_with_freqs) {
      pageRank_ppv(dg.graph(), counts, edge_weights);
    } else {
      pageRank(dg.graph(), edge_weights);
    }
    disamb_csentence(cs, dg);
    print_disamb_csent(cout, cs);
  }
  //print_complete_csent(cout, cs, dg);
}

void do_dgraph_gviz(const vector<string> & input_files,
		    const string & out_dir) {


  for(vector<string>::const_iterator fname=input_files.begin(); fname != input_files.end(); ++fname) {

    File_elem dg_finfo(*fname);

    DisambGraph dg;
    dg.read_from_binfile(dg_finfo.get_fname());
    dg.reset_edge_weigths();
    dg_finfo.set_path(out_dir);
    dg_finfo.ext = ".dot";

    if (glVars::mcr_with_freqs) {
      map<string, size_t> counts;
      bool ok = W2Syn::instance().syn_counts(counts);
      if (!ok) {
	cerr << "Error! There are no freqs. Check file " << glVars::w2s_filename << endl;
	exit(-1);
      }
      
      pageRank_ppv(dg.graph(), counts);
    } else {
      pageRank(dg.graph());
    }
    write_dgraph_graphviz(dg_finfo.get_fname(), dg.graph());
  }
}

void test (string & fullname_in) {

  

//   CSentence cs;
//   cs.read_from_binfile(cs_finfo.get_fname());

//   File_elem dg_finfo(cs.id(), cs_finfo.path, ".dgraph");
  //DisambGraph dg;
  //dg.read_from_binfile(fullname_in);

//   DisambGraph dg;
//   CSentence cs;

//   dg.read_from_binfile(fullname_in);




//   CSentence cs;

//   cs.read_from_binfile(fullname_in);
//   cout << cs;

//   ifstream fh_in(fullname_in.c_str());

//   if (!fh_in) {
//     cerr << "Can't open " << fullname_in << endl;
//     exit(-1);
//   }

//   vector<CSentence> vcs;
//   CSentence cs;

//   cs.read_aw(fh_in);
//   cs.write_to_binfile("/tmp/cskk.bin");
//   CSentence cs2;
//   cs2.read_from_binfile("/tmp/cskk.bin");
//   cout << cs;
//   cout << endl;
//   cout << cs2;
//   exit(0);
//   try {
//     while (cs.read_aw(fh_in)) {
//       vcs.push_back(cs);
//       cs = CSentence();
//     }
//   } 
//   catch (string & e) {
//     cerr << "Errore reading " << fullname_in << ":" << e << "\n";
//     throw(e);    
//   }
 
  
//  for(vector<CSentence>::const_iterator it = vcs.begin(); it != vcs.end(); ++it) {
//     cout << *it;
//     cout << endl;
//   }
}



int main(int argc, char *argv[]) {

  string out_dir;
  string mcr_binfile(mcr_default_binfile);

  bool opt_create_dgraph = false;
  bool opt_create_kgraph = false;
  //bool opt_eval_dgraph = false;
  bool opt_disamb_csent_dgraph = false;
  bool opt_disamb_csent_kgraph = false;
  bool opt_do_gviz = false;
  bool opt_do_test = false;
  bool opt_with_w = false;
 

  vector<string> input_files;
  string fullname_in;

  using namespace boost::program_options;

  const char desc_header[] = "ukb_aw: perform AW with MCR KB based algorithm\n"
    "Usage:\n"
    "ukb_aw --create_dgraph file.txt -> Create a disgraph from the file (.dgraph extension).\n"
    "ukb_aw --dis_csent csentence.csent -> Disambiguate words of a csentence. A dgraph with same id must be in the directory.\n"
    "                                      If a directory name is specified, disambiguate all found csentences.\n"
    "Options:";

  
  options_description po_desc(desc_header);

  po_desc.add_options()
    ("create_dgraph,c", "Create dgraph binary file(s), one per context, with extension .dgraph")
    ("dis_csent,d", "Disambiguate csentence and output result. A dgraph with same id must be in the directory.")
    ("with_freqs,f", "Use freqs when disambiguation (pageRank PPVs).")
    ("with_weights,w", "Use weigths in pageRank.")
    ("help,h", "This page")
    ("graphviz,G", "Dump disambGraph to a graphviz format. Output file has same name and extension .dot")
    ("mcr_binfile,M", value<string>(), "Binary file of MCR (see create_mcrbin). Default is mcr_wnet.bin")
    ("w2syn_file,W", value<string>(), "Word to synset map file. Default is ../Data/Preproc/wn1.6_index.sense_freq")
    ("out_dir,O", value<string>(), "Directory for leaving output files.")
    ("test,t", "(Internal) Do a test.")
    ("verbose,v", "Be verbose.")
    ;

  options_description po_desc_kgraph("KGraph options");
  po_desc_kgraph.add_options()
    ("create_kgraph", "Create kgraph binary file given one csentence (dgraph must have same name and be in same directory)")
    ("dis_csent_kgraph", "Disambiguate csentence and output result. A kgraph with same id must be in the directory.")
    ;
  
  options_description po_visible;
  po_visible.add(po_desc).add(po_desc_kgraph);
  
  options_description po_hidden("Hidden");
  po_hidden.add_options()
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

    if (vm.count("create_dgraph")) {
      opt_create_dgraph = true;
    }

    if (vm.count("create_kgraph")) {
      opt_create_kgraph = true;
    }

    if (vm.count("dis_csent")) {
      opt_disamb_csent_dgraph = true;
    }


    if (vm.count("graphviz")) {
      opt_do_gviz = true;
    }

    if (vm.count("dis_csent_kgraph")) {
      opt_disamb_csent_kgraph = true;
    }

    if (vm.count("w2syn_file")) {
      glVars::w2s_filename = vm["w2syn_file"].as<string>();
    }

    if (vm.count("mcr_binfile")) {
      mcr_binfile = vm["mcr_binfile"].as<string>();
    }

    if (vm.count("out_dir")) {
      out_dir = vm["out_dir"].as<string>();
    }

    if (vm.count("with_freqs")) {
      glVars::mcr_with_freqs = true;
    }

    if (vm.count("with_weights")) {
      opt_with_w = true;
    }


    if (vm.count("input-file")) {
      fullname_in = vm["input-file"].as<string>();
    }

    if (vm.count("test")) {
      opt_do_test = true;
    }
    conflicting_options(vm, "create_dgraph", "dis_csent");
  }
  catch(exception& e) {
    cerr << e.what() << "\n";
    throw(e);
  }

//   writeV(cout, input_files);
//   cout << endl;
//   return 0;

  if(opt_do_test) {
    test(fullname_in);
    return 0;
  }

  if(opt_create_kgraph && !opt_create_dgraph) {
    create_kgraph(fullname_in, out_dir);
  }

  if(opt_create_dgraph) {
    Mcr::create_from_binfile(mcr_binfile);
    create_dgraphs_from_corpus(fullname_in, out_dir, opt_create_kgraph);
    return 0;
  }

  if(opt_do_gviz) {
    extract_input_files(fullname_in, input_files, "dgraph");
    
    if(input_files.empty()) {
      cout << po_visible << endl;
      cerr << "Error: No input files." << endl;
      exit(0);      
    }
    do_dgraph_gviz(input_files, out_dir);
    return 0;
  }


  if (opt_disamb_csent_dgraph || opt_disamb_csent_kgraph) {
    
    extract_input_files(fullname_in, input_files, "csent");

    if(input_files.empty()) {
      cout << po_visible << endl;
      cerr << "Error: No input files." << endl;
      return -1;
    }
    if(opt_disamb_csent_dgraph) {
      dis_csent<DisambGraph>(input_files, ".dgraph", opt_with_w);
    }
    
    if(opt_disamb_csent_kgraph) {
      dis_csent<KGraph>(input_files, ".kgraph", opt_with_w);
    }
  }
  return 0;
}
