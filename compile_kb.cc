#include "common.h"
#include "wdict.h"
#include "globalVars.h"
#include "configFile.h"
#include "fileElem.h"
#include "kbGraph.h"
#include <string>
#include <iostream>
#include <fstream>

// Boost libraries

#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>


// timer

#include <boost/timer.hpp>


using namespace std;
using namespace boost;
using namespace ukb;


void query (const string & str) {

  Kb & kb = Kb::instance();
  KbGraph & g = kb.graph();

  bool aux;
  Kb_vertex_t u;

  tie(u, aux) = kb.get_vertex_by_name(str);
  if (aux) {
    cout << get(vertex_name, g, u);
    cout << "\n";
    graph_traits<Kb::boost_graph_t>::out_edge_iterator it , end;
    tie(it, end) = out_edges(u, g);
    for(;it != end; ++it) {
      cout << "  ";
      cout << get(vertex_name, g, target(*it, g));
      cout << ":" << get(edge_weight, g, *it) << "\n";
    }
  }
}

bool parse_csv(const string & str,
			   vector<string> & V) {

  typedef tokenizer<char_separator<char> > tokenizer_t;

  char_separator<char> sep(" \t,");
  tokenizer_t tok(str.begin(), str.end(), sep);
  copy(tok.begin(), tok.end(), back_inserter(V));
  return (V.size() > 0);
}

void set_source_rels(const string & str,
					 set<string> & S) {

  vector<string> V;
  set<string>().swap(S); // empty set

  parse_csv(str, V);

  copy(V.begin(), V.end(),
	   inserter(S, S.end())); // Remove duplicates

}

void iquery() {

  Kb & kb = Kb::instance();
  KbGraph & g = kb.graph();

  bool aux;
  Kb_vertex_t u;

  while	(cin) {
	string str;
	bool inv = false;
	cout << "\nQ:";
	cin >> str;
	if (str == "q") break;
	if (str.substr(0, 2) == "i:") {
	  inv=true;
	  str = str.substr(2);
	}
	tie(u, aux) = kb.get_vertex_by_name(str);
	if (aux) {
	  cout << get(vertex_name, g, u);
	  cout << "\n";
	  if (!inv) {
		graph_traits<Kb::boost_graph_t>::out_edge_iterator it , end;
		tie(it, end) = out_edges(u, g);
		for(;it != end; ++it) {
		  cout << "  ";
		  cout << get(vertex_name, g, target(*it, g));
		  cout << ":" << get(edge_weight, g, *it) << "\n";
		}
	  } else {
		graph_traits<Kb::boost_graph_t>::in_edge_iterator iit , iend;
		tie(iit, iend) = in_edges(u, g);
		for(;iit != iend; ++iit) {
		  cout << "  ";
		  cout << get(vertex_name, g, source(*iit, g));
		  cout << ":" << get(edge_weight, g, *iit) << "\n";
		}
	  }
	} else {
	  cout << "\n"<< str << " not present!";
	}
  }
}

void remove_dangling(string & fname) {

  Kb::create_from_binfile(fname);

  size_t i = 0;
  do {
	i = Kb::instance().unlink_dangling();
	cout << i << " dangling\n";
  } while (i);

  Kb::instance().write_to_binfile("kk.bin");
}


int main(int argc, char *argv[]) {

  srand(3);

  timer load;

  bool opt_info = false;
  bool opt_Info = false;
  bool opt_query = false;
  bool opt_iquery = false;
  bool opt_dump = false;
  bool opt_dang = false;

  string fullname_out("kb_wnet.bin");
  vector<string> kb_files;
  string query_vertex;


  glVars::kb::v1_kb = false; // Use v2 format
  glVars::kb::filter_src = false; // by default, don't filter relations by src

  set<string> src_allowed;

  string cmdline("cmd: !! -v ");
  cmdline += glVars::ukb_version;
  for (int i=0; i < argc; ++i) {
    cmdline += " ";
    cmdline += argv[i];
  }

  const char desc_header[] = "compile_kb: create a serialized image of the KB\n"
    "Usage:\n"
    "compile_kb [-o output.bin] [-f \"src1, src2\"] kb_file.txt kb_file.txt ... -> Create a KB image reading relations textfiles.\n"
    "compile_kb -i kb_file.bin -> Get info of a previously compiled KB.\n"
    "compile_kb -q concept-id kb_file.bin -> Query a node on a previously compiled KB.\n"
    "Options:";

  using namespace boost::program_options;

  options_description po_desc(desc_header);

  po_desc.add_options()
    ("help,h", "This help page.")
    ("version", "Show version.")
    ("filter_src,f", value<string>(), "Filter relations according to their sources.")
    ("info,i", "Give info about some Kb binfile.")
    ("Info,I", "Give more info about Kb binfile. This option can be computationally expensive.")
    ("dump", "Dump a serialized graph. Warning: very verbose!.")
    ("nodangling", "Recursively remove all dangling nodes from graph.")
    ("output,o", value<string>(), "Output file name.")
    ("query,q", value<string>(), "Given a vertex name, display its coocurrences.")
    ("iquery,Q", "Interactively query graph.")
    ("undirected,U", "Force undirected graph.")
    ("verbose,v", "Be verbose.")
    ("rtypes,r", "Keep relation types on edges.")
    ;
  options_description po_desc_hide("Hidden");
  po_desc_hide.add_options()
    ("input-file",value<vector<string> >(), "Input files.")
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

    if (vm.count("info")) {
      opt_info = true;
    }

    if (vm.count("Info")) {
      opt_Info = true;
    }

    if (vm.count("iquery")) {
      opt_iquery = true;
    }

    if (vm.count("nodangling")) {
      opt_dang = true;
    }

    if (vm.count("query")) {
      opt_query = true;
      query_vertex = vm["query"].as<string>();
    }

    if (vm.count("filter_src")) {
	  glVars::kb::filter_src = true;
	  set_source_rels(vm["filter_src"].as<string>(), src_allowed);
    }

    if (vm.count("dump")) {
      opt_dump = true;
    }

    if (vm.count("undirected")) {
	  glVars::kb::keep_directed = false;
    }

    if (vm.count("rtypes")) {
	  glVars::kb::keep_reltypes = true;
    }

    if (vm.count("input-file")) {
      kb_files = vm["input-file"].as<vector<string> >();
    }

     if (vm.count("output")) {
       fullname_out = vm["output"].as<string>();
     }
  }
  catch(std::exception& e) {
    cerr << e.what() << "\n";
	exit(-1);
  }

  if (kb_files.size()==0) {
    cout << po_desc << endl;
    cout << "No input files" << endl;
    exit(1);
  }

  if (opt_info) {
    Kb::create_from_binfile(kb_files[0]);
    Kb::instance().display_info(cout);
    return 0;
  }

  if (opt_Info) {
    Kb::create_from_binfile(kb_files[0]);

	int id_m, id_M;
	int od_m, od_M;
	int comp;

    Kb::instance().display_info(cout);

	tie(id_m, id_M) = Kb::instance().indeg_maxmin();
	tie(od_m, od_M) = Kb::instance().outdeg_maxmin();
	comp = Kb::instance().components();

	cout << "In degree (max, min):  (" << id_m << ", " << id_M << ")\n";
	cout << "Out degree (max, min): (" << od_m << ", " << od_M << ")\n";
	cout << "Number of (strong) components: " << comp << endl;

    return 0;
  }

  if (opt_dang) {
	remove_dangling(kb_files[0]);
    return 0;
  }

  if (opt_iquery) {
    Kb::create_from_binfile(kb_files[0]);
    iquery();
    return 0;
  }

  if (opt_dump) {
    Kb::create_from_binfile(kb_files[0]);
    Kb::instance().dump_graph(cout);
    return 0;
  }

  if (opt_query) {
    Kb::create_from_binfile(kb_files[0]);
    query(query_vertex);
    return 1;
  }

  if (glVars::verbose) {
    show_global_variables(cerr);
  }

  if (glVars::verbose)
    cerr << "Reading relations"<< endl;

  Kb::create_from_txt(kb_files[0], src_allowed );
  for(size_t i=1; i < kb_files.size(); ++i) {
    Kb::instance().add_from_txt(kb_files[i], src_allowed);
  }

  if (glVars::verbose)
    cerr << "Writing binary file: "<< fullname_out<< endl;
  Kb::instance().add_comment(cmdline);
  Kb::instance().write_to_binfile(fullname_out);
  if (glVars::verbose)
    cerr << "Wrote " << num_vertices(Kb::instance().graph()) << " vertices and " << num_edges(Kb::instance().graph()) << " edges" << endl;

  return 0;
}

/*
 * Local Variables:
 * mode: c++
 * End:
 */
