// nedwsd.cc

#include "common.h"
#include "wdict.h"
#include "globalVars.h"
#include "configFile.h"
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

static bool opt_variants = false;

void query (const string & str) {

  Kb & kb = Kb::instance();
  KbGraph & g = kb.graph();

  bool aux;
  Kb_vertex_t u;

  tie(u, aux) = kb.get_vertex_by_name(str);
  if (aux) {
    cout << kb.get_vertex_name(u);
    cout << "\n";
    graph_traits<Kb::boost_graph_t>::out_edge_iterator it , end;
    tie(it, end) = out_edges(u, g);
    for(;it != end; ++it) {
      cout << "  ";
      cout << kb.get_vertex_name(kb.edge_target(*it));
      cout << ":" << g[*it].weight << "\n";
    }
  }
}

/*void show_bfs_path(string & str) {

  Kb & kb = Kb::instance();
  KbGraph & g = kb.graph();
  bool aux;
  Kb_vertex_t u, v;
  vector<Kb_vertex_t> p;
  vector<string> path_str;

  path_str = split(str, " ");
  if(path_str.size() < 2) {
        cout << "\n"<< " p needs two concepts!";
        return;
  }
  tie(u, aux) = kb.get_vertex_by_name(path_str[0]);
  if(!aux) {
        cout << path_str[0] << " not present!";
        return;
  }

  tie(v, aux) = kb.get_vertex_by_name(path_str[1]);
  if(!aux) {
        cout << path_str[1] << " not present!";
        return;
  }
  kb.bfs(u, p);
  Kb_vertex_t w = v;
  while(w != u) {
        print_iquery_v(g, w, 0, 1);
        w = p[w];
  }
  print_iquery_v(g, u, 0, 1);
}*/

/*void show_neighbours(string & str) {

  Kb & kb = Kb::instance();
  KbGraph & g = kb.graph();

  bool inv;
  bool aux;
  Kb_vertex_t u;

  if (str.substr(0, 2) == "i:") {
        inv=true;
        str = str.substr(2);
  }
  tie(u, aux) = kb.get_vertex_by_name(str);
  if (aux) {
        print_iquery_v(g, u);
        if (!inv) {
          graph_traits<Kb::boost_graph_t>::out_edge_iterator it , end;
          tie(it, end) = out_edges(u, g);
          for(;it != end; ++it) {
                print_iquery_v(g, target(*it, g), g[*it].weight, 2);
          }
        } else {
          graph_traits<Kb::boost_graph_t>::in_edge_iterator iit , iend;
          tie(iit, iend) = in_edges(u, g);
          for(;iit != iend; ++iit) {
                print_iquery_v(g, source(*iit, g), g[*iit].weight, 2);
          }
        }
  } else {
        cout << "\n"<< str << " not present!\n";
  }
}*/

void show_dict_entries(const string & str) {

  WDict_entries entries = WDict::instance().get_entries(str);

  if (!entries.size()) {
        cout << str << " not in dictionary.\n";
        return;
  }
  cout << str << " ->";
  for(size_t i = 0; i < entries.size(); i++) {
        cout << " " << entries.get_entry_str(i) << ":" << entries.get_freq(i);
  }
  cout << "\n";
}

void iquery() {

  while        (cin) {
        string str;
        size_t l;
        cout << "\nQ:";
        read_line_noblank(cin, str, l);;
        if (str == "q") break;
        if (str.substr(0, 2) == "p:") {
          str = str.substr(2);
          //show_bfs_path(str);
          continue;
        }
        if (str.substr(0, 2) == "d:") {
          str = str.substr(2);
          show_dict_entries(str);
          continue;
        }
        //show_neighbours(str);
  }
}

void sPath(const string & sPathV) {

  string source;
  vector<string> aux,targets;
  std::vector<std::vector<std::string> > paths;
  aux = split(sPathV, "#");
  if (aux.size() < 2) {
        cerr << "Spath error: you must at least specify two nodes\n.";
        exit(-1);
  }
  source = aux[0];
  set<string> S(aux.begin() + 1, aux.end());
  copy(S.begin(), S.end(), back_inserter(targets));
  Kb::instance().get_shortest_paths(source, targets, paths);
  for(std::vector<std::vector<std::string> >::iterator it = paths.begin(), end = paths.end();
          it != end; ++it) {
        writeV(cout, *it);
        cout << "\n";
  }
}




  //START OF THE MAIN FUNCTION
int main(int argc, char *argv[]) {


  srand(3);

  timer load;

  bool opt_semSign = false;

  // subgraph options
  string subg_init;


  string kb_file;
  string kb_binfile;
  string query_vertex;
  string sPathV;

  glVars::kb::v1_kb = false; // Use v2 format
  glVars::kb::filter_src = false; // by default, don't filter relations by src

  set<string> src_allowed;

  string cmdline("cmd: !! -v ");
  cmdline += glVars::ukb_version;
  for (int i=0; i < argc; ++i) {
    cmdline += " ";
    cmdline += argv[i];
  }

  const char desc_header[] = "nedwsd: create a serialized image of the KB\n"
    "Usage:\n"
    "nedwsd -o output.bin [-f \"src1, src2\"] kb_file.txt kb_file.txt ... -> Create a KB image reading relations textfiles.\n"
    "nedwsd -o dict.bin -D dict_textfile --serialize_dict kb_file.bin -> Create a dictionary image reading text dictionary.\n"
    "nedwsd -i kb_file.bin -> Get info of a previously compiled KB.\n"
    "nedwsd -q concept-id kb_file.bin -> Query a node on a previously compiled KB.\n"
    "Options:";

  using namespace boost::program_options;

  options_description po_desc("General options");
  po_desc.add_options()
    ("help,h", "This help page.")
    ("version", "Show version.")
    ("kb_binfile,K", value<string>(), "Binary file of KB (see compile_kb).")
    ("dict_file,D", value<string>(), "Dictionary text file.")
    ("verbose,v", "Be verbose.")
    ("sem_signatures,s", "Builds semantic signatures.")
    ;

  options_description po_desc_create("Options for creating binary graphs");
  po_desc_create.add_options()
        ("minput", "Do not die when dealing with malformed input.")
        ;


  options_description po_hidden("Hidden");
  po_hidden.add_options()
    ("input-file",value<string>(), "Input files.")
    ;
  options_description po_visible(desc_header);
  //po_visible.add(po_desc).add(po_desc_create).add(po_desc_query);
  po_visible.add(po_desc).add(po_desc_create);

  options_description po_desc_all("All options");
  po_desc_all.add(po_visible).add(po_hidden);

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
      cout << po_visible << endl;
      exit(0);
    }

    if (vm.count("version")) {
      cout << glVars::ukb_version << endl;
      exit(0);
    }

    if (vm.count("verbose")) {
      glVars::verbose = 1;
    }

    if (vm.count("kb_binfile")) {
      kb_binfile = vm["kb_binfile"].as<string>();
    }

    if (vm.count("sem_signatures")) {
      opt_semSign = true;
      kb_binfile = vm["sem_signatures"].as<string>();
    }

    if (vm.count("dict_file")) {
      glVars::dict::text_fname = vm["dict_file"].as<string>();
          opt_variants = true;
    }

    if (vm.count("minput")) {
      glVars::input::swallow = true;
    }

    if (vm.count("input-file")) {
      kb_binfile = vm["input-file"].as<string>();
    }
  }
  catch(std::exception& e) {
    cerr << e.what() << "\n";
        exit(-1);
  }

  if (opt_semSign) {
    Kb::create_from_binfile(kb_binfile);
    Kb & kb = ukb::Kb::instance();
    kb.structural_weighting();
    //kb.dump_graph(cout)
    return 0;
  }



  return 0;
}

