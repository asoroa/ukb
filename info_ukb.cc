#include "mcrGraph.h"
#include "kGraph.h"
#include "disambGraph.h"
#include <string>
#include <iostream>
#include <fstream>

// Program options

#include <boost/program_options.hpp>

using namespace std;
using namespace boost;
using namespace ukb;

void display_csent(string & csentence_fname) {

  CSentence cs;
  cs.read_from_binfile(csentence_fname);

  cout << cs << endl;
  //print_complete_csent(cout, cs, dg);
}

int main(int argc, char *argv[]) {

  string fullname_in;
  string fullname_out;
  string out_dir;
  //  string mcr_binfile(mcr_default_binfile);

  bool opt_display_csent = false;

  const char desc_header[] = "info_ukb: display various information\n"
    "Usage:\n"
    "info_ukb [options] file\n"
    "Options:";

  try {
    using namespace boost::program_options;

    options_description po_desc(desc_header);

    po_desc.add_options()
      ("display_csentence,c", "Create dgraph binary file(s), one per context, with extension .dgraph")
      ("help,h", "This page")
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

    if (vm.count("display_csentence")) {
      opt_display_csent = true;
    }

    if (vm.count("input-file")) {
      fullname_in = vm["input-file"].as<string>();
    }

    if (vm.count("output-file")) {
      fullname_out = vm["output-file"].as<string>();
    }
  }
  catch(exception& e) {
    cerr << e.what() << "\n";
    throw(e);
  }

  if(opt_display_csent) {
    display_csent(fullname_in);
  }

  return 0;
}

/*
 * Local Variables:
 * mode: c++
 * compile-command: "make -k info_ukb"
 * End:
 */
