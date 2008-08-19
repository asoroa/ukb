// -*-C++-*-

// This is a hack. Same as coocGraph but directed
// Use for hlex cooc graphs.
// @@todo use only 1 coogGraph, but then recreate
//   Data/CoogGraph/coograph.bin et al.

#ifndef COOCGRAPH2_H
#define COOCGRAPH2_H

#include <iosfwd>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map.hpp>
#include <boost/graph/properties.hpp>

#include"coocGraph.h"

using boost::vertex_name_t;
using boost::vertex_name;
using boost::edge_index_t;
using boost::edge_index;


using std::string;

using boost::property;
using boost::adjacency_list;
using boost::graph_traits;

typedef adjacency_list <
  boost::listS,
  boost::vecS,
  boost::bidirectionalS,
  property<vertex_name_t, string, 
	   property<vertex_cfreq_t, size_t> >,
  property<edge_index_t, size_t,
	   property<edge_freq_t, float> > > coGraph2;

class CoocGraph2 {

 public:
  typedef coGraph2 boost_graph_t;
  typedef graph_traits<boost_graph_t>::vertex_descriptor vertex_descriptor;
  typedef graph_traits<boost_graph_t>::vertex_iterator vertex_iterator;
  typedef graph_traits<boost_graph_t>::adjacency_iterator adjacency_iterator;

  typedef graph_traits<boost_graph_t>::edge_descriptor edge_descriptor;
  typedef graph_traits<boost_graph_t>::edge_iterator edge_iterator;
  typedef graph_traits<boost_graph_t>::out_edge_iterator out_edge_iterator;

  // Constructor and destructor
  CoocGraph2() : _docN(0) {};

  CoocGraph2(CoocGraph2 & cooc); // Copy vertices removing isolated ones

  ~CoocGraph2() {};

  // Accessors

  boost_graph_t & graph() { return g; }
  std::pair<vertex_descriptor, bool> getVertexByName(const std::string & str) const;

  vertex_descriptor findOrInsertNode(const string & str);
  edge_descriptor findOrInsertEdge(vertex_descriptor u, vertex_descriptor v );


  // Info

  size_t num_vertices() { return boost::num_vertices(g); }
  size_t num_edges() { return boost::num_edges(g); }
  size_t num_docs() { return _docN; }
  std::pair<float, float> minmax() const;

  // fill graph

  void fill_cograph(std::ifstream & fh);

  // Chi square and prunning

  void calculate_chisq();

  // Normalize freqs. so they fit in [0,1] interval
  bool normalize_edge_freqs(  float cutValue = 190.00);


  void remove_isolated_vertices();

  // Remove edges below theshold
  void prune_zero_edges(float threshold = 0.0);

  void write_to_binfile (const std::string & fName) const;
  void read_from_binfile (const std::string & fName);

private:

  void insert_doc(std::vector<vertex_descriptor> & doc);
  void insert_doc_workaround(const std::string & w1, 
			     std::vector<vertex_descriptor> & doc);

  void read_from_stream (std::ifstream & is);
  std::ofstream & write_to_stream(std::ofstream & o) const;

  std::map<std::string, vertex_descriptor> _nodeMap;
  boost_graph_t g;

  std::set<std::string> _docSet;
  size_t _docN;
};

void write_vertex(std::ostream & o, const CoocGraph2::vertex_descriptor & v,
		  const CoocGraph2::boost_graph_t & g);

#endif

/*
 * Local Variables:
 * mode: c++
 * compile-command: "make -f makefile_cooc"
 * End:
 */
