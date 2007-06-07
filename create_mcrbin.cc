#include "mcrGraph.h"
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

// bfs 

#include <boost/graph/breadth_first_search.hpp>
#include <boost/pending/indirect_cmp.hpp>
#include <boost/pending/integer_range.hpp>


using namespace std;
using namespace boost;

std::pair<string, string> basename(const string & str) {

  boost::filesystem::path p(str);
  p.normalize();
  return std::make_pair<string, string>(p.branch_path().string(), p.leaf());
}

string::const_iterator find_last(string::const_iterator sit,
				 string::const_iterator sit_end, 
				 char delim) {
  string::const_iterator sit_found = sit_end;
  for(;sit != sit_end;++sit) 
    if (*sit == delim) sit_found = sit;
  return sit_found;
}


//

// template < typename TimeMap > class bfs_time_visitor : public default_bfs_visitor {
//   typedef typename property_traits < TimeMap >::value_type T;
// public:
//   bfs_time_visitor(TimeMap tmap, T & t):m_timemap(tmap), m_time(t) { }
//   template < typename Vertex, typename Graph >
//   void discover_vertex(Vertex u, const Graph & g) const {
//     put(m_timemap, u, m_time++);
//   }
//   TimeMap m_timemap;
//   T & m_time;
// };

// template <typename DistanceMap>
// class bacon_number_recorder : public default_bfs_visitor
// {
// public:
//   bacon_number_recorder(DistanceMap dist) : d(dist) { }
  
//   template <typename Edge, typename Graph>
//   void tree_edge(Edge e, const Graph& g) const
//   {
//     typename graph_traits<Graph>::vertex_descriptor
//       u = source(e, g), v = target(e, g);
//     d[v] = d[u] + 1;
//   }
// private:
//   DistanceMap d;
// };
// // Convenience function

// template <typename DistanceMap>
// bacon_number_recorder<DistanceMap>
// record_bacon_number(DistanceMap d)
// {
//   return bacon_number_recorder<DistanceMap>(d);
// }


void create_mcr_bin(const string &name) {

  //Mcr wnet("../mcr_source/wei_relations.txt", "../mcr_source/MCR+TSSemcor.all");

  Mcr::create_from_txt("mcr_source/wei_relations.txt", "mcr_source/MCR+TSSemcor.all");
  Mcr::instance().write_to_binfile(name);
  cerr << "Wrote " << num_vertices(Mcr::instance().graph()) << " vertices and " << num_edges(Mcr::instance().graph()) << " edges" << endl;

}

void do_bfs() {

  Mcr & mcr2 = Mcr::instance();

  timer bfs;

  McrGraph & g = mcr2.graph();

  std::vector < Mcr_vertex_size_t  > dtime;

  string src("00370521-n");

  cerr << src << " " << mcr2.getVertexByName(src).first << endl;
  mcr2.bfs(src, dtime);
  cerr << bfs.elapsed() << endl;
  for(unsigned int i = 0; i < num_vertices(g); ++i) 
    cout << mcr2.getVertexName(i) << " " << mcr2.getVertexName(dtime[i]) << endl;
    //cout << mcr2.getVertexName(i) << " " << dtime[i] << endl;
    //cout << i << " " << dtime[i] << endl;
}


int main() {

  srand(3);

  timer load;

  //  create_mcr_bin("mcr_wnet.bin");
  //return 1;

  Mcr::create_from_binfile("mcr_wnet.bin");
  cerr << "Loaded:" << load.elapsed() << endl;
  //do_bfs();
  //return 1;

  char *Sent[] = {"national","gallery","choice","picture","example","problem","artist","colour","video",
		  "talk","problem","connection","practice","painter","composition","work","artist","Poussin",
		  "subject","master","painting"};

  char *Sent2[] = {
    "acquiring",
    "acquisition",
    "blooper",
    "blunder",
    "boner",
    "boo-boo",
    "botch",
    "bungle",
    "flub",
    "foul-up",
    "fuckup",
    "getting",
    "misdoing"};

  vector<string> sentence(Sent, Sent + 21);
  //vector<string> sentence(Sent2, Sent2 + 13);

  DisambGraph dgraph;
  timer dgr_timer;
  CSentence csent(sentence);
  cerr << csent << endl;
  fill_disamb_graph(csent, dgraph);
  cerr << "DisambG created " << dgr_timer.elapsed() << endl;
  dgraph.write_to_binfile("w.dgraph");
  dgraph.prune();
  write_dgraph_graphviz("w.dot", dgraph.graph());
  
  timer dis_pr;

  dgraph.transform_csentence(csent);
  hits(dgraph.graph());
  //pageRank(dgraph.graph());
  cerr << "PageRank " << dis_pr.elapsed() << endl;
  disamb_csentence(csent, dgraph);
  print_csent_dgraph(cerr, csent, dgraph);
  return 1;
}
