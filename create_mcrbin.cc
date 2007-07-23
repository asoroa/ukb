#include "common.h"
#include "globalVars.h"
#include "configFile.h"
#include "fileElem.h"
#include "mcrGraph.h"
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


void read_skip_relations (const string & file,
			  set<string> & skip_rels) {

  const string delims(" \t");
  ifstream fi(file.c_str(), ifstream::binary|ifstream::in);
  if (!fi) {
    cerr << "Error: can't open " << file << endl;
    exit(-1);
  }
  string line;
  while(fi) {
    getline(fi, line);
    vector<string> rels = split(line, delims);
    for(vector<string>::iterator it = rels.begin(); it != rels.end(); ++it)
      skip_rels.insert(*it);
  }
}

int main(int argc, char *argv[]) {

  srand(3);

  timer load;

  bool opt_info = false;

  string fullname_out("mcr_wnet.bin");
  string relations_file("mcr_source/wei_relations.txt");
  string mcr_file("mcr_source/MCR+TSSemcor.all");
  string out_dir;

  const size_t source_rels_N = 5;
  const char *source_rels_default[source_rels_N] = {"16", "20", "xg", "xn", "xs"};

  const char desc_header[] = "create_mcrbin: create a serialized image of the MCR\n"
    "Usage:\n"
    "create_mcrbin mcr_file.txt [output.bin] -> Create a MCR image.\n"
    "Options:";
  
  try {
    using namespace boost::program_options;

    options_description po_desc(desc_header);

    po_desc.add_options()
      ("help,h", "This help page.")
      ("info,i", "Give info about some Mcr binfile.")
      ("out_dir,O", value<string>(), "Directory for leaving output files.")
      ("relations_file,r", value<string>(), "Specify file about relations (default mcr_source/wei_relations.txt).")
      ("w2syn_file,W", value<string>(), "Word to synset map file (default is ../Data/Preproc/wn1.6_index.sense_freq).")
      ("param,p", value<string>(), "Specify parameter file.")
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

    // verbosity

    if (vm.count("verbose")) {
      glVars::verbose = 1;
    }

    if (vm.count("info")) {
      opt_info = true;
    }

    // Config params first, so that they can be overriden by switch options

    if (vm.count("param")) {
      parse_config(vm["param"].as<string>());
    }

    if (vm.count("out_dir")) {
      out_dir = vm["out_dir"].as<string>();
    }

    if (vm.count("w2syn_file")) {
      glVars::w2s_filename = vm["w2syn_file"].as<string>();
    }

    if (vm.count("input-file")) {
      mcr_file = vm["input-file"].as<string>();
    }

    if (vm.count("output-file")) {
      fullname_out = vm["output-file"].as<string>();
    }
  }
  catch(exception& e) {
    cerr << e.what() << "\n";
    throw(e);
  }

  if (opt_info) {
    Mcr::create_from_binfile(mcr_file);
    Mcr::instance().display_info(cout);
    return 0;
  }

  if (!glVars::rel_source.size()) {
    // Default relations
    copy(source_rels_default, source_rels_default + source_rels_N, 
	 back_inserter(glVars::rel_source));
  }

  if (glVars::verbose) {
    show_global_variables(cerr);
    cerr << "MCR file:" << mcr_file << "\t";
    cerr << "Relations file:" << relations_file << endl;
  }

  set<string> source_rels(glVars::rel_source.begin(), glVars::rel_source.end());

  if (glVars::verbose) 
    cerr << "Reading relations"<< endl;
  Mcr::create_from_txt(relations_file, mcr_file, source_rels);
  File_elem mcr_fe(fullname_out);
  mcr_fe.set_path(out_dir);
  if (glVars::verbose) 
    cerr << "Writing binary file: "<< mcr_fe.get_fname()<< endl;
  Mcr::instance().write_to_binfile(mcr_fe.get_fname());
  if (glVars::verbose) 
    cerr << "Wrote " << num_vertices(Mcr::instance().graph()) << " vertices and " << num_edges(Mcr::instance().graph()) << " edges" << endl;

  return 1;
}

/*
 * Local Variables:
 * mode: c++
 * compile-command: "make -k create_mcrbin"
 * End:
 */
