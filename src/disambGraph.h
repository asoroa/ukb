// -*-C++-*-

#ifndef DISAMBGRAPH_H
#define DISAMBGRAPH_H

#include "kbGraph.h"
#include "csentence.h"

#include <boost/graph/adjacency_list.hpp>

using boost::adjacency_list;
using boost::graph_traits;
using boost::property;
using boost::property_map;
using boost::vertex_name_t;
using boost::vertex_name;
using boost::edge_weight_t;
using boost::edge_weight;

enum edge_freq_t    { edge_freq };

namespace boost {
	BOOST_INSTALL_PROPERTY(edge, freq);
}

namespace ukb {

	typedef adjacency_list <
		boost::listS,
		boost::vecS,
		boost::undirectedS,
		//  boost::bidirectionalS,
		property<vertex_name_t, std::string>,          // the synset name
		property<edge_weight_t, float>
		> DisambG;

	typedef graph_traits<DisambG>::vertex_descriptor Dis_vertex_t;
	typedef graph_traits<DisambG>::edge_descriptor Dis_edge_t;
	typedef graph_traits<DisambG >::vertices_size_type Dis_vertex_size_t;

	class DisambGraph {

	public:

		typedef DisambG boost_graph_t; // the underlying graph type
		typedef graph_traits<DisambG>::vertex_descriptor vertex_t;
		typedef graph_traits<DisambG>::edge_descriptor edge_t;
		typedef graph_traits<DisambG >::vertices_size_type vertex_size_t;


		DisambGraph();

		size_t size() const {return num_vertices(g); }
		std::pair<Dis_vertex_t, bool> get_vertex_by_name(const std::string & str) const;

		void fill_graph(Kb::vertex_descriptor src,
						Kb::vertex_descriptor tgt,
						const std::vector<Kb::vertex_descriptor> & parents);

		void fill_graph(const std::set<Kb::edge_descriptor> & E);

		Dis_vertex_t add_dgraph_vertex(const std::string & str);
		void add_dgraph_edge(Dis_vertex_t u, Dis_vertex_t v, float w = 1.0);

		void write_to_binfile (const std::string & fName) const;
		void read_from_binfile (const std::string & fName);

		DisambG & graph() {return g;}
		void prune() {}
		void reset_edge_weights();


		// prank

		void pageRank_ppv(const std::vector<float> & ppv_map,
						  std::vector<float> & ranks);

	private:

		void read_from_stream (std::ifstream & is);
		std::ofstream & write_to_stream(std::ofstream & o) const;

		//typedef std::vector<Dis_vertex_t> VertexV;
		//std::map<std::string, VertexV> w2syns;
		boost::unordered_map<std::string, Dis_vertex_t> synsetMap;

		DisambG g;
	};

	//////////////////////////////////////////////////////////////7
	// global functions

	// fill dgraph with a sentence

	void build_dgraph_bfs(const CSentence & sentence,
						  DisambGraph & dgraph);
	void build_dgraph_bfs(const CSentence & cs, DisambGraph & dgraph,
						  const std::vector<float> & ppv_ranks);

	// dfs version

	void build_dgraph_dfs(const CSentence &cs, DisambGraph & dgraph);
	void build_dgraph_dfs_nocosenses(const CSentence &cs, DisambGraph & dgraph);


	bool disamb_csentence_dgraph(CSentence & cs, DisambGraph & dgraph,
								 const std::vector<float> & ranks);

	void disamb_cword_dgraph(CSentence::iterator it, DisambGraph & dgraph,
							 const std::vector<float> & ranks);


	// HITS ranking

	void hits(DisambG & g, std::vector<float> & ranks);

	// PageRank ranking

	bool dgraph_ppr(const CSentence & cs, DisambGraph & dgraph,
					std::vector<float> & ranks,
					CSentence::const_iterator exclude_word_it);

	bool dgraph_ppr(const CSentence & cs, DisambGraph & dgraph,
					std::vector<float> & ranks);

	// apply ppr to mention graph
	bool dgraph_mention_ppr(const CSentence & cs, DisambGraph & dgraph,
					std::vector<float> & ranks);

	// Degree ranking

	bool dgraph_degree(DisambGraph & dgraph, std::vector<float> & ranks);

	// Static

	bool dgraph_static(DisambGraph & dgraph,
					   std::vector<float> & ranks);


	std::ostream & print_complete_csent(std::ostream & o, CSentence & cs, DisambGraph & dgraph);

	// export to dot format (graphviz)

	void write_dgraph_graphviz(const std::string & fname, DisambG & g);

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
}
#endif
