// -*-C++-*-

#ifndef COOCGRAPH_H
#define COOCGRAPH_H

#include <iosfwd>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map.hpp>
#include <boost/graph/properties.hpp>

using boost::vertex_name_t;
using boost::vertex_name;
using boost::edge_index_t;
using boost::edge_index;



enum vertex_freq_t     { vertex_freq };   // vertex freq
enum vertex_cfreq_t    { vertex_cfreq };  // vertex context (document) freq
enum edge_freq_t       { edge_freq };     // edge freq


namespace boost {
  BOOST_INSTALL_PROPERTY(vertex, freq);
  BOOST_INSTALL_PROPERTY(vertex, cfreq);
  BOOST_INSTALL_PROPERTY(edge, freq);
}

using std::string;

using boost::property;
using boost::adjacency_list;
using boost::graph_traits;

typedef adjacency_list <
  boost::listS,
  boost::vecS,
  boost::undirectedS,
  property<vertex_name_t, string, 
	   property<vertex_cfreq_t, size_t> >,
  property<edge_index_t, size_t,
	   property<edge_freq_t, float> > > coGraph;

class CoocGraph {

 public:
  typedef coGraph boost_graph_t;
  typedef graph_traits<boost_graph_t>::vertex_descriptor vertex_descriptor;
  typedef graph_traits<boost_graph_t>::vertex_iterator vertex_iterator;
  typedef graph_traits<boost_graph_t>::adjacency_iterator adjacency_iterator;

  typedef graph_traits<boost_graph_t>::edge_descriptor edge_descriptor;
  typedef graph_traits<boost_graph_t>::edge_iterator edge_iterator;
  typedef graph_traits<boost_graph_t>::out_edge_iterator out_edge_iterator;

  // Constructor and destructor
  CoocGraph() : _docN(0) {};

  CoocGraph(CoocGraph & cooc); // Copy vertices removing isolated ones

  ~CoocGraph() {};

  void remove_isolated_vertices();

  // Accessors

  boost_graph_t & graph() { return g; }
  std::pair<vertex_descriptor, bool> getVertexByName(const std::string & str) const;

  vertex_descriptor findOrInsertNode(const string & str);
  edge_descriptor findOrInsertEdge(vertex_descriptor u, vertex_descriptor v );


  // Info

  size_t num_vertices() { return boost::num_vertices(g); }
  size_t num_edges() { return boost::num_edges(g); }
  size_t num_docs() { return _docN; }

  // fill graph

  void fill_cograph(std::ifstream & fh);

  // Chi square and prunning

  void CoocGraph::chisq_prune();

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

void write_vertex(std::ostream & o, const CoocGraph::vertex_descriptor & v,
		  const CoocGraph::boost_graph_t & g);

#endif

/*
 * Local Variables:
 * mode: c++
 * compile-command: "make -f makefile_cooc"
 * End:
 */
