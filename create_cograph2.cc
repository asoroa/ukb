#include "common.h"
#include "globalVars.h"
#include "fileElem.h"
#include "coocGraph2.h"
#include "hlex_agtree.h"

#include <string>
#include <iostream>
#include <fstream>


// Basename & friends
#include <boost/filesystem/operations.hpp>
#include "boost/filesystem/path.hpp"

// Program options

#include <boost/program_options.hpp>

using namespace std;
using namespace boost;

string basename(const string & fname) {
  
  namespace fs = boost::filesystem;

  fs::path full_path( fs::initial_path() );
  
  full_path = fs::system_complete( fs::path( fname, fs::native ) );

  return full_path.leaf();
}

bool exists_file(const string & fname) {

  namespace fs = boost::filesystem;

  fs::path full_path( fs::initial_path() );
  
  full_path = fs::system_complete( fs::path( fname, fs::native ) );
  return exists(full_path);
}

bool extract_input_files(const string & fullname,
			 vector<string> & input_files) {

  namespace fs = boost::filesystem;

  fs::path full_path( fs::initial_path() );
  
  full_path = fs::system_complete( fs::path( fullname, fs::native ) );

  if ( !fs::exists( full_path ) )
    {
      std::cerr << "\nNot found: " << full_path.native_file_string() << std::endl;
      return false;
    }

  if ( fs::is_directory( full_path) ) {

    fs::directory_iterator end_iter;
    for ( fs::directory_iterator dir_itr( full_path );
          dir_itr != end_iter;
          ++dir_itr ) {
      if (fs::is_directory(*dir_itr)) continue;
      input_files.push_back(dir_itr->native_file_string());
    }    
  } else {
    input_files.push_back(full_path.native_file_string());
  }
  return true;
}

void chsq (const string & input_name,
	   const string & output_name,
	   bool normalize) {

  CoocGraph2 coog;

  coog.read_from_binfile(input_name);
  coog.calculate_chisq();
  if(normalize) {
    if (!coog.normalize_edge_freqs()) {
      cerr << "Error: can't normalize freqs." << endl;
    }
  }
  coog.prune_zero_edges();
  coog.write_to_binfile(output_name);
}


void query (const string & coName, const string & str) {

  CoocGraph2 coog;
  coog.read_from_binfile(coName);

  bool aux;
  CoocGraph2::vertex_descriptor u;
  
  tie(u, aux) = coog.get_vertex_by_name(str);
  if (aux) {
    write_vertex(cout, u, coog.graph());
    cout << "\n";
    CoocGraph2::out_edge_iterator it , end;
    tie(it, end) = out_edges(u, coog.graph());
    for(;it != end; ++it) {
      cout << "  ";
      write_vertex(cout, target(*it, coog.graph()), coog.graph());
      cout << " " << get(edge_freq, coog.graph(), *it) << "\n";
    }
  }
}

void info (string & coName) {

  CoocGraph2 coog;
  float min, max;

  coog.read_from_binfile(coName);
  size_t edge_n = coog.num_edges();
  tie(min, max) = coog.minmax();
  if (edge_n & 2) edge_n++; 
  cout << "V:" << coog.num_vertices() << " E:" << edge_n/2 << " Docs:" << coog.num_docs() << endl;
  cout << "min freq:" << min << " max freq:" << max << endl;
}

void histogram (string & coName) {

  CoocGraph2 coog;

  coog.read_from_binfile(coName);
  CoocGraph2::boost_graph_t & g = coog.graph();

  CoocGraph2::vertex_iterator it, end;
  tie(it, end) = vertices(g);
  for(;it != end; ++it) {
    cout << out_degree(*it, g) << "\n";
  }
}

void edge_histogram (string & coName) {

  CoocGraph2 coog;

  coog.read_from_binfile(coName);
  CoocGraph2::boost_graph_t & g = coog.graph();

  CoocGraph2::edge_iterator it, end;
  tie(it, end) = edges(g);
  for(;it != end; ++it) {
    cout << get(edge_freq, g, *it) << "\n";
  }
}

void test(string & coName) {

  CoocGraph2 coog;

  coog.read_from_binfile(coName);

  CoocGraph2::vertex_iterator it, end;
  for(tie(it, end) = vertices(coog.graph()); it != end; ++it) {
    if(out_degree(*it, coog.graph()) == 0) {
      cout << get(vertex_name, coog.graph(), *it) << '\n';
    }
  }
}

void create_from_hlex(const vector<string> & input_files, string & graph_binfile) {

  CoocGraph2 coog;
  map<AgZuhaitza::vertex_descriptor, CoocGraph2::vertex_descriptor> uMap;

  for(vector<string>::const_iterator file=input_files.begin();
      file != input_files.end(); ++file) {

    // See if file ends with .agt extension
    string word = basename(*file);
    string::size_type dot_i = word.find_last_of(".");
    if (word.find(".agt", dot_i) == string::npos) continue;

    // Read agt file
    AgZuhaitza agt;
    agt.read_from_file(*file);
    // Merge

    AgZuhaitza::vertex_iterator v_it, v_end;
    for(tie(v_it, v_end) = vertices(agt.graph()); v_it != v_end; ++v_it) {
      CoocGraph2::vertex_descriptor u = coog.findOrInsertNode(get(vertex_name, agt.graph(), *v_it));
      uMap[*v_it] = u;
    }
    AgZuhaitza::edge_iterator e_it, e_end;
    for(tie(e_it, e_end) = edges(agt.graph()); e_it != e_end; ++e_it) {
      size_t hub_dist = get(vertex_hubDist, agt.graph(), source(*e_it, agt.graph())); // hDist of father
      CoocGraph2::vertex_descriptor u = uMap[source(*e_it, agt.graph())];
      CoocGraph2::vertex_descriptor v = uMap[target(*e_it, agt.graph())];
      CoocGraph2::edge_descriptor cooc_e = coog.find_or_insert_edge(u, v);
      float old_freq = get(edge_freq, coog.graph(), cooc_e);
      float delta = hub_dist ? 1.0 / static_cast<float>(hub_dist) : 1.0;
      put(edge_freq, coog.graph(), cooc_e, old_freq + delta);
    }
  }  
  float min, max;
  if (glVars::verbose) {
    cout << "After processing graph:\n";
    cout << "V:" << coog.num_vertices() << " E:" << coog.num_edges() << " Docs:" << coog.num_docs() << endl;
  }
  tie(min, max) =  coog.minmax();
  cerr <<  min << "," << max << endl;
  coog.normalize_edge_freqs();
  coog.prune_zero_edges(0.0001); // Prec. of 1e-4
  coog.write_to_binfile(graph_binfile);
  if (glVars::verbose) {
    cout << "After pruning:\n";
    cout << "V:" << coog.num_vertices() << " E:" << coog.num_edges() << " Docs:" << coog.num_docs() << endl;
  }
  tie(min, max) =  coog.minmax();
  cerr <<  min << "," << max << endl;
}


int main(int argc, char *argv[]) {

  srand(3);

  bool opt_force = false;
  bool opt_verbose = false;

  bool opt_query = false;
  bool opt_info = false;
  bool opt_histo = false;
  bool opt_edge_histo = false;
  bool opt_test = false;
  bool opt_chsq = false;
  bool opt_normalize = false;
  bool opt_hlex = false;

  string query_vertex;

  string fullname_in;
  string graph_binfile("coocgraph.bin");


  //float threshold = 0.0; // 95.0% confidence
  //float threshold = 2.706; // 80.0% confidence
  //float threshold = 3.84146; // 95.0% confidence
  //float threshold = 5.41189; // 98.0% confidence
  //float threshold = 6.6349;  // 99.0% confidence
  //float threshold = 10.827;  // 99.99% confidence



  const char desc_header[] = "create_cograph: create a coocurrence graph\n"
    "Usage:\n"
    "create_cograph [-f] file_or_dir [coocgraph.bin] \n"
    "\t file_or_dir: input textfile or directory. If directory given, all input files are read.\n"
    "\t coocgraph.bin: Name of serialization graph. Default is coocgraph.bin.\n"
    "create_cograph -c pruned_cograph.bin [-t threshold] cograph.bin\n"
    "create_cograph --hlex hlex_dir \n"
    "Options:";
  
  using namespace boost::program_options;

  options_description po_desc(desc_header);

  po_desc.add_options()
    ("chsq,c", value<string>(), "Given a serialized graph, create another with chisq information. Prune it according threshold (see -t).")
    ("edge_histogram,E", "Output edge freq. values.")
    ("help,h", "This help page.")
    ("hlex", "Create cooc. graph by joining hlex graphs.")
    ("histogram,H", "Output vertex degrees.")
    ("normalize,N", "Normalize chisq. values, so they fit in [0,1] interval.")
    ("info,i", "Get info of a serielized cograph.")
    ("force,f", "Don't use previous coocgraph.bin if exists.")
    ("query,q", value<string>(), "Given a vertex name, display its coocurrences.")
    ("verbose,v", "Be verbose.")
    ("threshold,t", value<float>(), 
     "Set a threshold for prunning the graph (use with -c option).\n"
     "\tTypical threshold values:\n"
     "\t  0.708   -- 60.0%  confidence\n"
     "\t  1.642   -- 80.0%  confidence\n"
     "\t  2.706   -- 90.0%  confidence\n"
     "\t  3.84146 -- 95.0%  confidence\n"
     "\t  5.41189 -- 98.0%  confidence\n "
     "\t  6.6349  -- 99.0%  confidence\n"
     "\t  7.883   -- 99.5%  confidence\n"
     "\t  10.827  -- 99.99% confidence\n")
    ("test,T", "Internal test.")
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

    if (vm.count("hlex")) {
      opt_hlex = true;
    }

    if (vm.count("chsq")) {
      opt_chsq = true;
      graph_binfile= vm["chsq"].as<string>();
    }

    if (vm.count("info")) {
      opt_info = true;
    }

    if (vm.count("histogram")) {
      opt_histo = true;
    }

    if (vm.count("normalize")) {
      opt_normalize = true;
    }

    if (vm.count("edge_histogram")) {
      opt_edge_histo = true;
    }

    if (vm.count("force")) {
      opt_force = true;
    }

    // verbosity
    if (vm.count("verbose")) {
      glVars::verbose = 1;
    }

    if (vm.count("test")) {
      opt_test = 1;
    }

    if (vm.count("query")) {
      opt_query = true;
      query_vertex = vm["query"].as<string>();
    }

    if (vm.count("threshold")) {
      glVars::chsq::threshold = vm["threshold"].as<float>();
    }

    if (vm.count("input-file")) {
      fullname_in = vm["input-file"].as<string>();
    }

    if (vm.count("output-file")) {
      graph_binfile = vm["output-file"].as<string>();
    }

  }
  catch(exception& e) {
    cerr << e.what() << "\n";
    throw(e);
  }


//______________________________________________________________________________

  CoocGraph2 coog;
  vector<string> input_files;

  extract_input_files(fullname_in, input_files);
  if(input_files.empty()) {
    cout << po_desc << endl;
    cerr << "Error: No input files." << endl;
    return -1;
  }

  if (opt_hlex) {
    create_from_hlex(input_files, graph_binfile);
    return 0;
  }

  if (opt_test) {
    test(input_files[0]);
    return 0;
  }

  if (opt_info) {
    info(input_files[0]);
    return 0;
  }

  if (opt_histo) {
    histogram(input_files[0]);
    return 0;
  }

  if (opt_edge_histo) {
    edge_histogram(input_files[0]);
    return 0;
  }

  if (opt_query) {
    query(input_files[0], query_vertex);
    return 0;
  }

  if (opt_chsq) {
    if (glVars::chsq::threshold < 0.0) {
      cerr << "Error: negative threshold\n";
      return -1;
    }
    if (glVars::verbose) {
      cout << "Cooc min: " << glVars::chsq::cooc_min << " Threshold: " << glVars::chsq::threshold << '\n';
    }
    chsq(input_files[0], graph_binfile, opt_normalize);
    return 0;
  }

//   if (!opt_force && exists_file(graph_binfile)) {
//     // Read previous graph
//     coog.read_from_binfile(graph_binfile);
//   }
  

//   if (opt_verbose) {
//     cout << "Before processing graph:\n";
//     cout << "V:" << coog.num_vertices() << " E:" << coog.num_edges() << " Docs:" << coog.num_docs() << endl;
//   }

  for(vector<string>::iterator it=input_files.begin();
      it != input_files.end(); ++it) {
    string word = basename(*it);
    ifstream fh(it->c_str());
    if (!fh) {
      cerr << "Error: can't read " << *it << endl;
      exit(-1);
    }
    if (opt_verbose) {
      cout << "Processing word " << word << " file " << *it << endl;
    }
    coog.fill_cograph(fh);
  }

  coog.write_to_binfile(graph_binfile);
  if (opt_verbose) {
    cout << "After processing graph:\n";
    cout << "V:" << coog.num_vertices() << " E:" << coog.num_edges() << " Docs:" << coog.num_docs() << endl;
  }
  return 0;
}


/*
 * Local Variables:
 * mode: c++
 * End:
 */
