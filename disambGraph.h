// -*-C++-*-

#ifndef DISAMBGRAPH_H
#define DISAMBGRAPH_H

#include "mcrGraph.h"
#include "csentence.h"

using boost::adjacency_list;
using boost::graph_traits;
using boost::property;
using boost::property_map;
using boost::vertex_name_t;
using boost::vertex_name;
using boost::vertex_rank_t;
using boost::vertex_rank;

enum edge_freq_t    { edge_freq };      // relation id

namespace boost {
  BOOST_INSTALL_PROPERTY(edge, freq);
}

typedef adjacency_list <
  boost::listS,
  boost::vecS,
  boost::undirectedS,
  property<vertex_name_t, std::string,          // the synset name (WN1.6)
	   property<vertex_rank_t, float> >,    // vertex rank
  property<edge_freq_t, float>
  > DisambG;

typedef graph_traits<DisambG>::vertex_descriptor Dis_vertex_t;
typedef graph_traits<DisambG>::edge_descriptor Dis_edge_t;
typedef graph_traits<DisambG >::vertices_size_type Dis_vertex_size_t;

class DisambGraph {

 public:

  typedef DisambG boost_graph_type; // the underlying graph type

  DisambGraph();

  std::pair<Dis_vertex_t, bool> getVertexByName(const std::string & str) const;

  void fill_graph(Mcr_vertex_t src,
		  Mcr_vertex_t tgt,
		  const std::vector<Mcr_vertex_t> & parents);
  void add_disamb_edge(Dis_vertex_t u, Dis_vertex_t v, size_t w = 1);
  void write_to_binfile (const std::string & fName) const;
  void read_from_binfile (const std::string & fName);

  DisambG & graph() {return g;}
  void prune() {}
private:

  std::vector<Dis_vertex_t> add_vertices_mcr_path(std::vector<std::string>::iterator v_it, 
						  std::vector<std::string>::iterator v_end);

  void read_from_stream (std::ifstream & is);
  std::ofstream & write_to_stream(std::ofstream & o) const;

  //typedef std::vector<Dis_vertex_t> VertexV;
  //std::map<std::string, VertexV> w2syns;
  std::map<std::string, Dis_vertex_t> synsetMap;

  DisambG g;
};

//////////////////////////////////////////////////////////////7
// global functions

// fill dgraph with a sentence

void fill_disamb_graph(const CSentence & sentence,
		       DisambGraph & dgraph);

void disamb_csentence(CSentence & cs, DisambGraph & dgraph);

void hits(DisambG & g);

void pageRank(DisambG & g);

std::ostream & print_disamb_csent(std::ostream & o, CSentence & cs);
std::ostream & print_complete_csent(std::ostream & o, CSentence & cs, DisambGraph & dgraph);

// export to dot format (graphviz)

void write_dgraph_graphviz(const std::string & fname, const DisambG & g);

// streaming functions for disambG type graphs
Dis_vertex_t read_vertex_from_stream(std::ifstream & is, 
				     DisambG & g);
Dis_edge_t read_edge_from_stream(std::ifstream & is, 
				 DisambG & g);
std::ofstream & write_vertex_to_stream(std::ofstream & o,
				       const DisambG & g,
				       const Dis_vertex_t & v);
std::ofstream & write_edge_to_stream(std::ofstream & o,
				     const DisambG & g,
				     const Dis_edge_t & e);
#endif
