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

struct Kg_dg_; // forward declaration

typedef std::vector<std::vector<Kg_dg_> > KG_SynsV;

class KGraph {

 public:

  typedef DisambG boost_graph_type; // the underlying graph type

  KGraph() {};
  KGraph(CSentence & cs, DisambGraph & disg);

  std::pair<Dis_vertex_t, bool> getVertexByName(const std::string & str) const;

  void write_to_binfile (const std::string & fName) const;
  void read_from_binfile (const std::string & fName);

  DisambG & graph() {return g;}

private:

  std::map<std::string, Dis_vertex_t> synsetMap;

  void create_kgSynsV(CSentence & cs, KG_SynsV & theSyns,
		      DisambGraph & dg);
  void fill_kg(Kg_dg_ & src, 
	       KG_SynsV::iterator cw_it,
	       KG_SynsV::iterator cw_end,
	       DisambG & disg);

  void read_from_stream (std::ifstream & is);
  std::ofstream & write_to_stream(std::ofstream & o) const;

  DisambG g;
};

//////////////////////////////////////////////////////////////7
// global functions

// fill KGraph with a sentence

void disamb_csentence(CSentence & cs, KGraph & dgraph);

#endif
