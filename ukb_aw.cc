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


/////////////////////////////////////////////////////////////
// Filesystem stuff (paths, extensions, etc)

string::const_iterator find_last(string::const_iterator sit,
				 string::const_iterator sit_end, 
				 char delim) {
  string::const_iterator sit_found = sit_end;
  for(;sit != sit_end;++sit) 
    if (*sit == delim) sit_found = sit;
  return sit_found;
}

struct File_elem {

  File_elem(const string & fname) {
    fill_finfo(fname);
  }

  File_elem(const string & fullname_in,
	    const string & out_dir,
	    const string & new_ext) {

    fill_finfo(fullname_in);

    size_t m = out_dir.size();
    if (m) {
      if(!boost::filesystem::exists(out_dir)) {
	// Better die
	//boost::filesystem::create_directory( out_dir_dir );
	cerr << "Error: " << out_dir << " doesn't exist.\n" << endl;
	exit(-1);
      }
      if (out_dir[m-1] == '/') --m;
      path.assign(out_dir, 0, m);
    }
    
    if (new_ext.size())
      ext = new_ext;
  }

  void fill_finfo(const string & str) {

    boost::filesystem::path p(str);
    p.normalize();

    path = p.branch_path().string();

    string file_fname = p.leaf(); // name + extension

    string::const_iterator beg = file_fname.begin();
    string::const_iterator end = file_fname.end();
    string::const_iterator it = find_last(file_fname.begin(), file_fname.end(), '.');

    ext.assign(it, end);
    fname.assign(beg, it); // without extension
  }

  string get_fname() const {
    string res(path);
    res.append("/");
    res.append(fname);
    res.append(ext);
    return res;
  }

  string path;
  string fname;
  string ext;
};

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

//Main program


void create_dgraph(string & fullname_in,
		   const string & out_dir) {

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
      cs = CSentence();
    }
  } 
  catch (string & e) {
    cerr << "Errore reading " << fullname_in << ":" << e << "\n";
    throw(e);    
  }
}


void dis_csent(string & csentence_fname) {

  File_elem cs_finfo(csentence_fname);

  CSentence cs;
  cs.read_from_binfile(cs_finfo.get_fname());

  File_elem dg_finfo(cs.id(), cs_finfo.path, ".dgraph");
  DisambGraph dg;
  dg.read_from_binfile(dg_finfo.get_fname());
  pageRank(dg.graph());
  disamb_csentence(cs, dg);
  print_csent(cout, cs, dg);
  print_complete_csent(cout, cs, dg);
}

void test (string & fullname_in) {

//   File_elem cs_finfo(fullname_in);

//   CSentence cs;
//   cs.read_from_binfile(cs_finfo.get_fname());

//   File_elem dg_finfo(cs.id(), cs_finfo.path, ".dgraph");
  DisambGraph dg;
  dg.read_from_binfile(fullname_in);
  dg.kk();

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

  string fullname_in;
  string fullname_out;
  string out_dir;
  string mcr_binfile(mcr_default_binfile);

  bool opt_create_dgraph = false;
  //bool opt_eval_dgraph = false;
  bool opt_disamb_csent = false;
  bool opt_do_test = false;

  const char desc_header[] = "ukb_aw: perform AW with MCR KB based algorithm\n"
    "Usage:\n"
    "ukb_aw file.txt -> Disambiguate words. Output key file to STDOUT.\n"
    "ukb_aw --create_dgraph file.txt -> Create a disgraph from the file (.dgraph extension).\n"
    "ukb_aw --dis_csent csentence.csent -> Disambiguate words of a csentence. A dgraph with same id must be in the directory.\n"
    "Options:";

  try {
    using namespace boost::program_options;

    options_description po_desc(desc_header);

    po_desc.add_options()
      ("create_dgraph,c", "Create dgraph binary file(s), one per context, with extension .dgraph")
      ("dis_csent,d", "Disambiguate csentence and output result. A dgraph with same id must be in the directory.")
      //      ("eval_dgraph,e", "Disambiguate dgraph and output result. Input file is a dgraph file.")
      ("help,h", "This page")
      ("mcr_binfile,M", value<string>(), "Binary file of MCR. Create with create_mcrbin application.")
      ("out_dir,O", value<string>(), "Directory for leaving output files.")
      ("test,t", "(Internal) Do a test.")
      ("verbose,v", "Be verbose.")
      ;
    options_description po_desc_hide("Hidden");
    po_desc_hide.add_options()
      ("input-file",value<string>(), "Input file.")
      ("output-file",value<string>(), "Output file.")
      ;
    options_description po_desc_all("All options");
    po_desc_all.add(po_desc).add(po_desc_hide);

    positional_options_description po_optdesc;
    po_optdesc.add("input-file", 1);
    po_optdesc.add("output-file", 1);

    variables_map vm;
    store(command_line_parser(argc, argv).
	  options(po_desc_all).
	  positional(po_optdesc).
	  run(), vm);
    notify(vm);


    // If asked for help, don't do anything more

    if (vm.count("help")) {
      cout << po_desc << endl;
      exit(0);
    }

    if (vm.count("create_dgraph")) {
      opt_create_dgraph = true;
    }

    if (vm.count("dis_csent")) {
      opt_disamb_csent = true;
    }


    if (vm.count("dis_dgraph")) {
      opt_disamb_csent = true;
    }

    if (vm.count("mcr_binfile")) {
      mcr_binfile = vm["mcr_binfile"].as<string>();
    }

    if (vm.count("out_dir")) {
      out_dir = vm["out_dir"].as<string>();
    }

    if (vm.count("input-file")) {
      fullname_in = vm["input-file"].as<string>();
    }

    if (vm.count("output-file")) {
      fullname_out = vm["output-file"].as<string>();
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

  if(opt_do_test) {
    test(fullname_in);
    return 0;
  }

  if(opt_create_dgraph) {
    Mcr::create_from_binfile(mcr_binfile);
    create_dgraph(fullname_in, out_dir);
  }

  if(opt_disamb_csent) {
    dis_csent(fullname_in);
  }

  return 0;
}
