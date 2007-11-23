#include "common.h"
#include "w2syn.h"
#include "globalVars.h"
#include "configFile.h"
#include "fileElem.h"
#include "mcrGraph.h"
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

void merge_cooc(string & fName) {
  
  Mcr & mcr = Mcr::instance();

  CoocGraph coog;
  coog.read_from_binfile(fName);
  map<CoocGraph::vertex_descriptor, Mcr_vertex_t> uMap;

  CoocGraph::vertex_iterator v_it, v_end;
  for(tie(v_it, v_end) = vertices(coog.graph()); v_it != v_end; ++v_it) {
    Mcr_vertex_t u = mcr.findOrInsertWord(get(vertex_name, coog.graph(), *v_it));
    uMap[*v_it] = u;
  }

  CoocGraph::edge_iterator e_it, e_end;
  for(tie(e_it, e_end) = edges(coog.graph()); e_it != e_end; ++e_it) {
    Mcr_vertex_t u = uMap[source(*e_it, coog.graph())];
    Mcr_vertex_t v = uMap[target(*e_it, coog.graph())];
    float w = get(edge_freq, coog.graph(), *e_it);
    mcr.findOrInsertEdge(u, v, w);
    mcr.findOrInsertEdge(v, u, w);			 
  }
}

void query (const string & str) {

  Mcr & mcr = Mcr::instance();
  McrGraph & g = mcr.graph();

  bool aux;
  Mcr_vertex_t u;
  
  tie(u, aux) = mcr.getVertexByName(str);
  if (aux) {
    cout << get(vertex_name, g, u);
    cout << "\n";
    graph_traits<Mcr::boost_graph_t>::out_edge_iterator it , end;
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
  bool opt_merge = false;
  bool opt_words = false;
  bool opt_param = false;
  bool opt_force_param = false;
  bool opt_query = false;

  string cograph_filename;
  string fullname_out("mcr_wnet.bin");
  string relations_file("mcr_16_source/wei_relations.txt");
  string mcr_file;
  string out_dir;
  string query_vertex;

  const size_t source_rels_N = 5;
  const char *source_rels_default[source_rels_N] = {"16", "20", "xg", "xn", "xs"};

  const char desc_header[] = "create_mcrbin: create a serialized image of the MCR\n"
    "Usage:\n"
    "create_mcrbin mcr_file.txt [output.bin] -> Create a MCR image.\n"
    "create_mcrbin -m coocgraph.bin mcr.bin output.bin -> Merge a mcr serialization with a coocurrence graph.\n"
    "Options:";
  
  using namespace boost::program_options;

  options_description po_desc(desc_header);

  po_desc.add_options()
    ("help,h", "This help page.")
    ("force-default-values,f", "Use default relations.")
    ("info,i", "Give info about some Mcr binfile.")
    ("words", "Insert word and word#pos nodes in the MCR.")
    ("merge,m", value<string>(), "Merge a coocurrence graph to a serialization graph.")
    ("out_dir,O", value<string>(), "Directory for leaving output files.")
    ("relations_file,r", value<string>(), "Specify file about relations (default mcr_16_source/wei_relations.txt).")
    ("w2syn_file,W", value<string>(), "Word to synset map file (default is ../Data/Preproc/wn1.6_index.sense_freq).")
    ("param,p", value<string>(), "Specify parameter file.")
    ("query,q", value<string>(), "Given a vertex name, display its coocurrences.")
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

    if (vm.count("merge")) {
      opt_merge = true;
      cograph_filename = vm["merge"].as<string>();
    }


    if (vm.count("force-default-values")) {
      opt_force_param = true;
    }

    if (vm.count("words")) {
      opt_words = true;
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

  if (mcr_file.size()==0) {
    cout << po_desc << endl;
    cout << "No input files" << endl;
    exit(0);
  }

  if (opt_info) {
    Mcr::create_from_binfile(mcr_file);
    Mcr::instance().display_info(cout);
    return 0;
  }

  if (opt_merge) {
    Mcr::create_from_binfile(mcr_file);
    merge_cooc(cograph_filename);
    Mcr::instance().write_to_binfile(fullname_out);
    return 1;
  }

  if (opt_query) {
    Mcr::create_from_binfile(mcr_file);
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
    cerr << "MCR file:" << mcr_file << "\t";
    cerr << "Relations file:" << relations_file << endl;
  }

  set<string> source_rels(glVars::rel_source.begin(), glVars::rel_source.end());

  if (glVars::verbose) 
    cerr << "Reading relations"<< endl;

  Mcr::create_from_txt(relations_file, mcr_file, source_rels);

  if (!opt_words) {
    if (glVars::verbose) 
      cerr << "Skipping words"<< endl;
  } else {
    if (glVars::verbose) 
      cerr << "Adding word nodes to MCR"<< endl;

    W2Syn & w2syn = W2Syn::instance();
    vector<string>::const_iterator word_it = w2syn.get_wordlist().begin();
    vector<string>::const_iterator word_end = w2syn.get_wordlist().end();
    for(; word_it != word_end; ++word_it) {
      vector<string>::const_iterator syn_it, syn_end;
      tie(syn_it, syn_end) = w2syn.get_wsyns(*word_it);
      Mcr::instance().add_words();
    }
  }

  File_elem mcr_fe(fullname_out);
  mcr_fe.set_path(out_dir);
  if (glVars::verbose) 
    cerr << "Writing binary file: "<< mcr_fe.get_fname()<< endl;
  Mcr::instance().write_to_binfile(mcr_fe.get_fname());
  if (glVars::verbose) 
    cerr << "Wrote " << num_vertices(Mcr::instance().graph()) << " vertices and " << num_edges(Mcr::instance().graph()) << " edges" << endl;

  return 0;
}

/*
 * Local Variables:
 * mode: c++
 * compile-command: "make -k create_mcrbin"
 * End:
 */
