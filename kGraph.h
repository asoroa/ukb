// -*-C++-*-

#ifndef KGRAPH_H
#define KGRAPH_H

#include "disambGraph.h"
#include "csentence.h"

using boost::adjacency_list;
using boost::graph_traits;
using boost::property;
using boost::property_map;
using boost::vertex_name_t;
using boost::vertex_name;
using boost::vertex_rank_t;
using boost::vertex_rank;

class KGraph {

 public:

  KGraph(DisambG & disg, CSentence & cs);

  void write_to_binfile (const std::string & fName) const;
  void read_from_binfile (const std::string & fName);

  DisambG & graph() {return g;}

private:

  void read_from_stream (std::ifstream & is);
  std::ofstream & write_to_stream(std::ofstream & o) const;

  DisambG g;
};

//////////////////////////////////////////////////////////////7
// global functions

// fill KGraph with a sentence

void disamb_csentence(CSentence & cs, KGraph & dgraph);

#endif
