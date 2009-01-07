#include "common.h"
#include "wdict.h"
#include "globalVars.h"
#include "configFile.h"
#include "fileElem.h"
#include "kbGraph.h"
//#include "disambGraph.h"
#include "coocGraph.h"
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

void merge_cooc(string & fName, bool store_w) {

  Kb & kb = Kb::instance();

  CoocGraph coog;
  coog.read_from_binfile(fName);
  map<CoocGraph::vertex_descriptor, Kb_vertex_t> uMap;

  CoocGraph::vertex_iterator v_it, v_end;
  for(tie(v_it, v_end) = vertices(coog.graph()); v_it != v_end; ++v_it) {
    Kb_vertex_t u = kb.find_or_insert_word(get(vertex_name, coog.graph(), *v_it));
    uMap[*v_it] = u;
  }

  CoocGraph::edge_iterator e_it, e_end;
  for(tie(e_it, e_end) = edges(coog.graph()); e_it != e_end; ++e_it) {
    Kb_vertex_t u = uMap[source(*e_it, coog.graph())];
    Kb_vertex_t v = uMap[target(*e_it, coog.graph())];
    float w = store_w ? get(edge_freq, coog.graph(), *e_it) : 1.0;
    kb.find_or_insert_edge(u, v, w);
    kb.find_or_insert_edge(v, u, w);
  }
  kb.add_relSource("Cooc: " + fName);
}


void merge_ts(bool store_w) {

  Kb::instance().add_dictionary(store_w); // with weights!
  Kb::instance().add_relSource("TS: " + glVars::dict_filename);
}

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



int main(int argc, char *argv[]) {

  srand(3);

  timer load;

  bool opt_info = false;
  bool opt_cooc = false;
  bool opt_ts = false;
  bool opt_param = false;
  bool opt_force_param = false;
  bool opt_query = false;
  bool opt_dump = false;

  bool opt_weight_ts = true;    // use weights for TS
  bool opt_weight_cooc = false; // don't use weights for cooc

  string cograph_filename;
  string fullname_out("kb_wnet.bin");
  vector<string> kb_files;
  string out_dir;
  string query_vertex;

  const size_t source_rels_N = 5;
  const char *source_rels_default[source_rels_N] = {"16", "20", "xg", "xn", "xs"};

  string cmdline("cmd:");
  for (int i=0; i < argc; ++i) {
    cmdline += " ";
    cmdline += argv[i];
  }

  const char desc_header[] = "create_mcrbin: create a serialized image of the KB\n"
    "Usage:\n"
    "create_mcrbin [-o output.bin] kb_file.txt kb_file.txt ... -> Create a KB image reading relations textfiles.\n"
    "create_mcrbin [-o output.bin] -c coocgraph.bin kb.bin -> Merge a kb serialization with a coocurrence graph.\n"
    "create_mcrbin [-o output.bin] --ts ts_dict.txt kb.bin -> Merge TS from textfile.\n"
    "Options:";

  using namespace boost::program_options;

  options_description po_desc(desc_header);

  po_desc.add_options()
    ("help,h", "This help page.")
    ("force-default-values,f", "Use default relations.")
    ("info,i", "Give info about some Kb binfile.")
    ("dump", "Dump a serialized graph. Warning: very verbose!.")
    ("cooc,c", value<string>(), "Merge a coocurrence graph to a serialization graph.")
    ("ts,t", value<string>(), "Merge topic signatures in a serialization graph. Asks for the textfile with ts info.")
    ("out_dir,O", value<string>(), "Directory for leaving output files.")
    ("output,o", value<string>(), "Output file name.")
    //    ("dict_file,W", value<string>(), "Word to synset map file (default is ../Data/Preproc/wn1.6_index.sense_freq).")
    ("param,p", value<string>(), "Specify parameter file.")
    ("query,q", value<string>(), "Given a vertex name, display its coocurrences.")
    ("ts_now", "Don't use weights when linking TS words to synsets (default is yes).")
    ("cooc_w", "Use weights when linking cooc words (default is don't).")
    ("verbose,v", "Be verbose.")
    ("kyoto,k", "Use kyoto format for input txt.")
    ("rtypes,r", "Keep relation types on edges. Implies -k.")
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
      opt_param = true;
    }

    if (vm.count("query")) {
      opt_query = true;
      query_vertex = vm["query"].as<string>();
    }

    if (vm.count("dump")) {
      opt_dump = true;
    }

    if (vm.count("cooc")) {
      opt_cooc = true;
      cograph_filename = vm["cooc"].as<string>();
    }

    if (vm.count("kyoto")) {
	  glVars::kb::v1_kb = false;
    }

    if (vm.count("rtypes")) {
	  if (glVars::kb::v1_kb) throw runtime_error("can't use --rtypes without --kyoto");
	  glVars::kb::keep_reltypes = true;
    }

    if (vm.count("ts")) {
      opt_ts = true;
      glVars::dict_filename = vm["ts"].as<string>(); // use w2s to read ts file
    }

    if (vm.count("force-default-values")) {
      opt_force_param = true;
    }

    if (vm.count("out_dir")) {
      out_dir = vm["out_dir"].as<string>();
    }

//     if (vm.count("dict_file")) {
//       if (opt_ts) {
// 	cerr << "Error, --ts and --dict_file options conflict!\n";
// 	exit(-1);
//       }
//       glVars::dict_filename = vm["dict_file"].as<string>();
//     }

    // weights

    if (vm.count("ts_now")) {
      opt_weight_ts = false;
    }

    if (vm.count("cooc_w")) {
      opt_weight_cooc = true;
    }


    if (vm.count("input-file")) {
      kb_files = vm["input-file"].as<vector<string> >();
    }

     if (vm.count("output")) {
       fullname_out = vm["output"].as<string>();
     }
  }
  catch(exception& e) {
    cerr << e.what() << "\n";
    throw(e);
  }

  if (kb_files.size()==0) {
    cout << po_desc << endl;
    cout << "No input files" << endl;
    exit(0);
  }

  if (opt_info) {
    Kb::create_from_binfile(kb_files[0]);
    Kb::instance().display_info(cout);
    return 0;
  }

  if (opt_dump) {
    Kb::create_from_binfile(kb_files[0]);
    Kb::instance().dump_graph(cout);
    return 0;
  }

  if (opt_cooc || opt_ts) {
    Kb::create_from_binfile(kb_files[0]);
    if (opt_cooc) merge_cooc(cograph_filename, opt_weight_cooc);
    if (opt_ts) merge_ts(opt_weight_ts);
    Kb::instance().add_comment(cmdline);
    Kb::instance().write_to_binfile(fullname_out);
    return 1;
  }


  if (opt_query) {
    Kb::create_from_binfile(kb_files[0]);
    query(query_vertex);
    return 1;
  }

  if (!opt_param && !opt_force_param) {
    cerr << "Error: no param file. For using default values, use --force-default-values" << endl;
    return -1;
  }

  if (!glVars::rel_source.size()) {
    // Default relations
    if (!opt_force_param) {
      cerr << "Error: no relations. For using default relations, use --force-default-values" << endl;
      return -1;
    }
    copy(source_rels_default, source_rels_default + source_rels_N,
	 back_inserter(glVars::rel_source));
  }

  if (glVars::verbose) {
    show_global_variables(cerr);
  }

  set<string> source_rels(glVars::rel_source.begin(), glVars::rel_source.end());

  if (glVars::verbose)
    cerr << "Reading relations"<< endl;

  Kb::create_from_txt(kb_files[0], source_rels );
  for(size_t i=1; i < kb_files.size(); ++i) {
    Kb::instance().add_from_txt(kb_files[i]);
  }

  File_elem kb_fe(fullname_out);
  kb_fe.set_path(out_dir);
  if (glVars::verbose)
    cerr << "Writing binary file: "<< kb_fe.get_fname()<< endl;
  Kb::instance().add_comment(cmdline);
  Kb::instance().write_to_binfile(kb_fe.get_fname());
  if (glVars::verbose)
    cerr << "Wrote " << num_vertices(Kb::instance().graph()) << " vertices and " << num_edges(Kb::instance().graph()) << " edges" << endl;

  return 0;
}

/*
 * Local Variables:
 * mode: c++
 * End:
 */
