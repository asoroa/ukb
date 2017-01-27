// -*-C++-*-

#ifndef MCRGRAPH_H
#define MCRGRAPH_H

#include <string>
#include <iterator>
#include <map>
#include <vector>
#include <set>
#include <ctime>            // std::time
#include <iosfwd>
#include <memory>

#include "kbGraph_common.h"
#include "kbGraph_v16.h"

// graph

#define BOOST_GRAPH_USE_NEW_CSR_INTERFACE
#include <boost/graph/compressed_sparse_row_graph.hpp>
#include <boost/graph/iteration_macros.hpp>

#if BOOST_VERSION > 103800
#include <boost/property_map/property_map.hpp>
#else
#include <boost/property_map.hpp>
#endif

#include <boost/graph/properties.hpp>

using boost::compressed_sparse_row_graph;
using boost::graph_traits;
using boost::property;
using boost::property_map;
using boost::edge_weight_t;
using boost::edge_weight;

// Properties for graphs

namespace ukb {


	class Kb {

	public:

		typedef compressed_sparse_row_graph<boost::bidirectionalS,
											vertex_prop_t,
											edge_prop_t> boost_graph_t;

		typedef graph_traits<boost_graph_t>::vertex_descriptor vertex_descriptor;
		typedef graph_traits<boost_graph_t>::vertex_iterator vertex_iterator;
		typedef graph_traits<boost_graph_t>::edge_descriptor edge_descriptor;
		typedef graph_traits<boost_graph_t>::out_edge_iterator out_edge_iterator;
		typedef graph_traits<boost_graph_t>::in_edge_iterator in_edge_iterator;

		// Singleton
		static Kb & instance();


		// 2 functions for creating Kb graphs
		//
		// 1. create_from_txt
		//    Create graph by reading a textfile with synset relations (synsFile)
		//    If glVars::kb::filter_src is true, exclude input relations not in
		//    rels_source set

		static void create_from_txt(const std::string & synsFile,
									const std::set<std::string> & rels_source);

		// same thing, but taking an istream as argument

		static void create_from_txt(std::istream & is,
									const std::set<std::string> & rels_source);



		// 2. create_from_binfile
		//    Load a binary snapshot of the graph into memory

		static void create_from_binfile(const std::string & o);


		// write_to_binfile
		// Write kb graph to a binary serialization file

		void write_to_binfile (const std::string & str);

		// write_to_textfile
		// Write kb graph to a text file

		void write_to_textfile (const std::string & fName) const;

		// write_to_textstream
		// Write kb graph to a text file
		std::ostream & write_to_textstream(std::ostream & o) const;

		// read_from_txt
		// add relations from synsFile to the graph

		void read_from_txt(const std::string & synsFile,
						   const std::set<std::string> & rels_source);

		// add relations from synsFile to the graph given an input stream

		void read_from_txt(std::istream & is,
						   const std::set<std::string> & rels_source);

		// add_relSource
		// add a new relation source

		void add_relSource(const std::string & str) { if (str.size()) m_relsSource.insert(str); }

		// Add tokens and link them to their synsets, according to the dictionary.
		// Note: the words are linked to nodes by _directed_ edges

		void add_dictionary(); // Adds all words of the current dictionary
		void add_token(const std::string & str); // Add just a word (lemma)

		// graph
		// Get the underlying boost graph

		boost_graph_t & graph() {return *m_g;}

		// Add relation type to edge

		void edge_add_reltype(edge_descriptor e, const std::string & rel);

		// Ask for a node

		std::pair<vertex_descriptor, bool> get_vertex_by_name(const std::string & str) const;

		// ask for node properties

		std::string & get_vertex_name(vertex_descriptor u) {return (*m_g)[u].name;}
		const std::string & get_vertex_name(vertex_descriptor u) const {return (*m_g)[u].name;}
		//std::string  get_vertex_gloss(vertex_descriptor u) const {return get(vertex_gloss, g, u);}

		// Get vertices iterator

		std::pair<vertex_iterator, vertex_iterator> get_vertices() { return boost::vertices(*m_g); }

		// Get out-edges for vertex u

		std::pair<out_edge_iterator, out_edge_iterator> out_neighbors(vertex_descriptor u);
		std::pair<in_edge_iterator, in_edge_iterator> in_neighbors(vertex_descriptor u);

		bool exists_edge(vertex_descriptor u, vertex_descriptor v) const {
			return edge(u, v, *m_g).second;
		}

		vertex_descriptor edge_source(edge_descriptor e) const { return source(e, *m_g); }
		vertex_descriptor edge_target(edge_descriptor e) const { return target(e, *m_g); }

		// ask for edge preperties

		std::vector<std::string> edge_reltypes(edge_descriptor e) const;

		float get_edge_weight(edge_descriptor e) const;
		void set_edge_weight(edge_descriptor e, float w);

		// get static pageRank

		const std::vector<float> & static_prank() const;

		// Given a previously calculated rank vector, output 2 vector, probably
		// filtering the nodes.
		//
		// Input params:
		//
		// ranks: previously calculated rank vector
		// filter_mode:
		//   0 -> no filter
		//   1 -> only words
		//   2 -> only concepts
		//
		// Output:
		//
		// outranks: new rank vector
		// vnames: node names

		void filter_ranks_vnames(const std::vector<float> & ranks,
								 std::vector<float> & outranks,
								 std::vector<std::string> & vnames,
								 int filter_mode) const;


		// Add a comment to graph

		void add_comment(const std::string & str);

		// Some useful info

		void display_info(std::ostream & o) const;
		size_t size() const {return num_vertices(*m_g); }

		std::pair<size_t, size_t> indeg_maxmin() const;
		std::pair<size_t, size_t> outdeg_maxmin() const;

		// get how many (strong) components the graph has
		int components() const;
		// get components vector cv[u] = number of component where u belongs
		int components(std::vector<size_t> & cv) const;

		const std::vector<std::string> & get_comments() const;

		// Get a random vertex

		vertex_descriptor get_random_vertex() const;

		// Graph algorithms

		bool bfs (vertex_descriptor source_synset, std::vector<vertex_descriptor> & synv) const ;

		bool dijkstra (vertex_descriptor src, std::vector<vertex_descriptor> & parents) const;

		void pageRank_ppv(const std::vector<float> & ppv_map,
						  std::vector<float> & ranks);

		void ppv_weights(const std::vector<float> & ppv);

		// given a source node and a limit (100) return a subgraph by performing a
		// bfs over the graph.
		// Output:
		//  V -> a vector of subgraph nodes
		//  E -> a map representing subgraph edges

		void get_subgraph(const std::string & src,
						  std::vector<std::string> & V,
						  std::vector<std::vector<std::string> > & E,
						  size_t limit = 100);

		// given a source node and a set of targets, compute the shortest paths from
		// source to each target
		// Output:
		//  paths -> a vector with the shortest paths (one elemet per source-target pair)

		bool get_shortest_paths(const std::string & src,
								const std::vector<std::string> & targets,
								std::vector<std::vector<std::string> > & paths);

		// Write graph in txt format (source of ukb)
		std::ostream & write_txt(std::ostream & o) const;
		std::ostream & dump_graph(std::ostream & o) const;

		static void create_from_kbgraph16(Kb16 & kbg) ;

	private:
		// Singleton
		static Kb * p_instance;
		static Kb * create();

		// Private methods
		Kb() : m_g(NULL), m_vertexN(0), m_edgeN(0) {};
		Kb(const Kb &) {};
		Kb &operator=(const Kb &);
		~Kb() {};

		vertex_descriptor InsertNode(const std::string & name, unsigned char flags);

		void read_from_stream (std::istream & o);
		std::ostream & write_to_stream(std::ostream & o) const;
		// Private members
		std::auto_ptr<boost_graph_t> m_g;
		std::set<std::string> m_relsSource;              // Relation sources
		std::map<std::string, vertex_descriptor> m_synsetMap; // synset name to vertex id

		// Registered relation types

		etype_t m_rtypes;

		std::vector<std::string> m_notes;        // Command line which created the graph

		// Aux variables

		std::vector<float> m_out_coefs;          // aux. vector of out-degree coefficients
		size_t m_vertexN;                        // Number of vertices
		size_t m_edgeN;                          // Number of edges
		std::vector<float> m_static_ppv;         // aux. vector with static prank computation
	};
}

#endif
