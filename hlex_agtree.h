// -*-C++-*-

#ifndef HLEX_AGTREE_H
#define HLEX_AGTREE_H

#include "kbGraph.h"
#include <map>
#include <string>

// graph

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map.hpp>
#include <boost/graph/properties.hpp>

// AgTree



enum vertex_hubDist_t  { vertex_hubDist };  // Distance to associated hub
enum vertex_hubN_t     { vertex_hubN };     // Associated hub
enum vertex_vOrig_t    { vertex_vOrig };    // Unused

namespace boost {
  BOOST_INSTALL_PROPERTY(vertex, hubDist);
  BOOST_INSTALL_PROPERTY(vertex, hubN);
  BOOST_INSTALL_PROPERTY(vertex, vOrig);  
}

namespace ukb {

  typedef boost::adjacency_list <  
	boost::listS,
	boost::vecS,
	boost::bidirectionalS,
	boost::property<vertex_name_t, std::string,
					boost::property<vertex_hubDist_t, size_t, 
									boost::property<vertex_hubN_t, size_t,
													boost::property<vertex_vOrig_t, Kb_vertex_t> // not used
													> > >
	> AgTree;


class AgZuhaitza {

public:
  // Typedefs
  typedef graph_traits<AgTree>::vertex_descriptor vertex_descriptor;
  typedef graph_traits<AgTree>::vertex_iterator vertex_iterator;
  typedef graph_traits<AgTree>::edge_descriptor edge_descriptor;
  typedef graph_traits<AgTree>::edge_iterator edge_iterator;

private:

  static long magic_id;
  static long magic_idG;

  AgTree g;
  std::map<std::string, vertex_descriptor> hitzMapa;
  std::string targetWord;
  vertex_descriptor targetVertex;
  size_t hubsN;
  std::vector<std::string> hubs;
  
  std::ifstream & read_from_stream (std::ifstream & o);
  std::ofstream & write_to_stream(std::ofstream & o) const;
  
public:

  AgZuhaitza() {};
  void read_from_file (const std::string & o);

  const AgTree & graph() const {return g;}
  const vertex_descriptor & rootVertex() const {return targetVertex;}
  size_t hub_number() const { return hubs.size(); }
};
}


#endif
