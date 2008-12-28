// -*-C++-*-

#ifndef COOCGRAPH_H
#define COOCGRAPH_H

#include <iosfwd>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/property_map.hpp>
#include <boost/graph/properties.hpp>

enum vertex_freq_t     { vertex_freq };   // vertex freq
enum vertex_cfreq_t    { vertex_cfreq };  // vertex context (document) freq
enum edge_freq_t       { edge_freq };     // edge freq

namespace boost {
  BOOST_INSTALL_PROPERTY(vertex, freq);
  BOOST_INSTALL_PROPERTY(vertex, cfreq);
  BOOST_INSTALL_PROPERTY(edge, freq);
}

namespace ukb {

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

  // Accessors

  boost_graph_t & graph() { return g; }
  std::pair<vertex_descriptor, bool> get_vertex_by_name(const std::string & str) const;

  vertex_descriptor findOrInsertNode(const string & str);
  edge_descriptor find_or_insert_edge(vertex_descriptor u, vertex_descriptor v );


  // Info

  size_t num_vertices() { return boost::num_vertices(g); }
  size_t num_edges() { return boost::num_edges(g); }
  size_t num_docs() { return _docN; }
  std::pair<float, float> minmax() const;

  // fill graph

  void fill_cograph(std::ifstream & fh);

  // fill with dling thesaurus

  void fill_dling_th(std::ifstream & fh);

  // Chi square and prunning

  void calculate_chisq();

  // Normalize freqs. so they fit in [0,1] interval
  bool normalize_edge_freqs(  float cutValue = 190.00);


  // Remove all vertices whose degree is <= m
  void remove_isolated_vertices(size_t m = 0);

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

void write_vertex(std::ostream & o, const CoocGraph::vertex_descriptor & v,
				  const CoocGraph::boost_graph_t & g);

}
#endif

/*
 * Local Variables:
 * mode: c++
 * compile-command: "make -f makefile_cooc"
 * End:
 */
