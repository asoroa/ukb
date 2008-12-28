// -*-C++-*-

#ifndef KGRAPH_H
#define KGRAPH_H

#include "disambGraph.h"
#include "csentence.h"

// using boost::adjacency_list;
// using boost::graph_traits;
// using boost::property;
// using boost::property_map;
// using boost::vertex_name_t;
// using boost::vertex_name;
// using boost::vertex_rank_t;
// using boost::vertex_rank;

namespace ukb {

  class KGraph {

  public:

	typedef DisambG boost_graph_t; // the underlying graph type

	KGraph() {};
	KGraph(CSentence & cs, DisambGraph & disg);

	std::pair<Dis_vertex_t, bool> get_vertex_by_name(const std::string & str) const;

	void write_to_binfile (const std::string & fName) const;
	void read_from_binfile (const std::string & fName);

	DisambG & graph() {return g;}

    void reset_edge_weigths() {};

  private:

	std::map<std::string, Dis_vertex_t> synsetMap;


	void read_from_stream (std::ifstream & is);
	std::ofstream & write_to_stream(std::ofstream & o) const;

	DisambG g;
  };

  //////////////////////////////////////////////////////////////7
  // global functions

  // fill KGraph with a sentence

  void disamb_csentence(CSentence & cs, KGraph & dgraph);
}
#endif
