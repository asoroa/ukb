#include "globalVars.h"
#include "kbGraph.h"
#include "kbGraph_v16.h"
#include <string>
#include <iostream>
#include <fstream>

// Program options

#include <boost/program_options.hpp>

using namespace ukb;
using namespace std;
using namespace boost;


int main(int argc, char *argv[]) {

	srand(3);

	string fullname_out;
	string kb_file;

	glVars::kb::v1_kb = false; // Use v2 format
	glVars::kb::filter_src = false; // by default, don't filter relations by src

	const char desc_header[] = "convert2.0: convert a 1.6 serialized graph to 2.0\n"
		"Usage:\n"
		"compile_kb -o output.bin kb_v16.bin -> Convert kb_v16.bin to output.bin.\n"
		"Options:";

	using namespace boost::program_options;

	options_description po_desc(desc_header);

	po_desc.add_options()
		("help,h", "This help page.")
		("version,V", "Show version.")
		("verbose,v", "Be verbose.")
		("output,o", value<string>(), "Output file name.")
		;
	options_description po_desc_hide("Hidden");
	po_desc_hide.add_options()
		("input-file",value<string>(), "Input file.")
		;
	options_description po_desc_all("All options");
	po_desc_all.add(po_desc).add(po_desc_hide);

	positional_options_description po_optdesc;
	po_optdesc.add("input-file", -1);

	variables_map vm;

	try {
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

		if (vm.count("version")) {
			cout << glVars::ukb_version << endl;
			exit(0);
		}

		if (vm.count("verbose")) {
			glVars::verbose = 1;
		}

		if (vm.count("input-file")) {
			kb_file = vm["input-file"].as<string>();
		}

		if (vm.count("output")) {
			fullname_out = vm["output"].as<string>();
		}

	}
	catch(std::exception& e) {
		cerr << e.what() << "\n";
		exit(-1);
	}

	if (!kb_file.size()) {
		cerr << po_desc << "\n";
		cerr << "Error: no input files\n" << endl;
		exit(1);
	}

	if (!fullname_out.size()) {
		cerr << po_desc << "\n";
		cerr << "Error: no output files\n" << endl;
		exit(1);
	}

	Kb16::create_from_binfile(kb_file);

	Kb16 & kb16 = Kb16::instance();

	Kb::create_from_kbgraph16(kb16);

	Kb::instance().write_to_binfile(fullname_out);
	if (glVars::verbose)
		cerr << "Wrote " << num_vertices(Kb::instance().graph()) << " vertices and " << num_edges(Kb::instance().graph()) << " edges" << endl;

	return 0;
}
