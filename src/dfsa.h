// -*-C++-*-

#ifndef DFSA_H
#define DFSA_H

#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/properties.hpp>

using namespace boost;

/*!
 * adaptor for limited DFS
 */
template <typename Graph> class dfsa {

	typedef graph_traits<Graph> Traits;

public:
	// Graph requirements
	typedef typename Traits::vertex_descriptor         vertex_descriptor;
	typedef typename Traits::edge_descriptor           edge_descriptor;
	typedef typename Traits::directed_category         directed_category;
	typedef typename Traits::edge_parallel_category    edge_parallel_category;
	typedef typename Traits::traversal_category        traversal_category;

	// IncidenceGraph requirements
	typedef typename Traits::out_edge_iterator         out_edge_iterator;
	typedef typename Traits::degree_size_type          degree_size_type;

	// AdjacencyGraph requirements
	typedef typename Traits::adjacency_iterator        adjacency_iterator;

	// VertexListGraph requirements
	typedef typename Traits::vertex_iterator           vertex_iterator;
	typedef typename Traits::vertices_size_type        vertices_size_type;

	// EdgeListGraph requirements
	typedef typename Traits::edge_iterator             edge_iterator;
	typedef typename Traits::edges_size_type           edges_size_type;

	typedef typename Traits::in_edge_iterator          in_edge_iterator;

	typedef typename Graph::edge_property_type         edge_property_type;
	typedef typename Graph::vertex_property_type       vertex_property_type;
	typedef Graph                                      graph_type;
	typedef typename Graph::graph_property_type        graph_property_type;

	dfsa(const Graph& g, int md = (std::numeric_limits<int>::max)()) : m_g(g), m_md(md), m_d(0) {}
	int depth() const { return m_d; }
	int max_depth() const { return m_md;}
	void inc_depth() { m_d++; }
	void dec_depth() { m_d--; }
	void set_max_depth(int md) { m_md = md; }
	bool out_of_depth() const { return m_d >= m_md; }
	const Graph& get_graph() const { return m_g;}
private:
	const Graph& m_g;
	int m_md; //Maximal depth
	int m_d; //Current depth
};


///IncidenceGraph requirements
namespace boost
{

	template <typename Graph>
	inline typename std::pair<
		typename    graph_traits<dfsa<Graph>  >::vertex_iterator,
		typename     graph_traits< dfsa<Graph> >::vertex_iterator >
	vertices(const dfsa<Graph>& g) {
		return vertices(g.get_graph());
	}

	template <typename Graph>
	inline typename graph_traits<dfsa<Graph> >::vertex_descriptor
	target(typename graph_traits<dfsa<Graph> >::edge_descriptor e, const dfsa<Graph>& g) {
		return target(e, g.get_graph());
	}

	template <typename Graph>
	inline typename graph_traits<dfsa<Graph> >::vertex_descriptor
	source(typename graph_traits<dfsa<Graph> >::edge_descriptor e, const dfsa<Graph>& g) {
		return source(e, g.get_graph());
	}

	template <typename Graph>
	inline typename std::pair<
		typename    graph_traits<dfsa<Graph> >::out_edge_iterator,
		typename    graph_traits<dfsa<Graph> >::out_edge_iterator >
	out_edges( typename graph_traits< dfsa<Graph> >::vertex_descriptor v, const  dfsa<Graph>& g) {
		typedef typename graph_traits< dfsa<Graph> >::out_edge_iterator   myIter;
		myIter oei, oeie;
		tie(oei, oeie) =  out_edges(v, g.get_graph());
		if (g.out_of_depth()) return std::make_pair(oeie, oeie);
		return  std::make_pair(oei, oeie);
	}

	template <typename Graph>
	inline typename graph_traits<dfsa<Graph> >::degree_size_type
	out_degree(typename graph_traits<dfsa<Graph> >::vertex_descriptor v, const  dfsa<Graph>&  g) {
		return out_degree(v, g.get_graph());
	}
};

#endif

