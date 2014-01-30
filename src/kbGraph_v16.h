// -*-C++-*-

#ifndef KBGRAPHV16_H
#define KBGRAPHV16_H

#include <string>
#include <iterator>
#include <fstream>
#include <map>
#include <vector>
#include <set>
#include <ctime>            // std::time
#include <iostream>

// integer types

#include <boost/cstdint.hpp>

// graph

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>

#include <boost/version.hpp>
#if BOOST_VERSION > 103800
#include <boost/property_map/property_map.hpp>
#else
#include <boost/property_map.hpp>
#endif

#include <boost/graph/properties.hpp>

using boost::adjacency_list;
using boost::graph_traits;
using boost::property;
using boost::property_map;
using boost::vertex_name_t;
using boost::vertex_name;
using boost::edge_weight_t;
using boost::edge_weight;

// Properties for graphs

enum vertex_flags_t { vertex_flags};  // flags for vertex (for knowing
                                      // wether a node is word or
                                      // synset)
enum vertex_gloss_t { vertex_gloss};  // gloss of node
enum edge_id_t      { edge_id };      // relation id
enum edge_rtype_t   { edge_rtype };  // relation type

namespace boost {
  BOOST_INSTALL_PROPERTY(vertex, flags);
  BOOST_INSTALL_PROPERTY(vertex, gloss);
  BOOST_INSTALL_PROPERTY(edge, id);
  BOOST_INSTALL_PROPERTY(edge, rtype);
}

namespace ukb {

  typedef adjacency_list <
	boost::listS,
	boost::vecS,
	boost::bidirectionalS,
	property<vertex_name_t, std::string,
			 property<vertex_gloss_t, std::string,
					  property<vertex_flags_t, unsigned char> > >,
	property<edge_weight_t, float,
			 property<edge_rtype_t, boost::uint32_t> >
	> Kb16Graph;


  typedef graph_traits<Kb16Graph>::vertex_descriptor Kb16_vertex_t;
  typedef graph_traits<Kb16Graph>::edge_descriptor Kb16_edge_t;
  typedef graph_traits < Kb16Graph >::vertices_size_type Kb16_vertex_size_t;

  class Kb16 {

  public:

	enum {
	  is_word = 1,
	  is_concept = 2
	} vflags;

	typedef Kb16Graph boost_graph_t; // the underlying graph type

	// Singleton
	static Kb16 & instance();


	// 2. create_from_binfile
	//    Load a binary snapshot of the graph into memory

	static void create_from_binfile(const std::string & o);


	// write_to_binfile
	// Write kb graph to a binary serialization file

	void write_to_binfile (const std::string & str);

	Kb16Graph & graph() {return g;}

	// Add nodes and relations to the graph

	Kb16_vertex_t find_or_insert_synset(const std::string & str);
	Kb16_vertex_t find_or_insert_word(const std::string & str);

	Kb16_edge_t find_or_insert_edge(Kb16_vertex_t u, Kb16_vertex_t v, float w );

	void unlink_vertex(Kb16_vertex_t u);

	// Unlink dangling_nodes (out_degree == 0).
	// Return num. of unlinked nodes

	size_t unlink_dangling();

	// Add relation type to edge

	void edge_add_reltype(Kb16_edge_t e, const std::string & rel);

	// Ask for a node

	std::pair<Kb16_vertex_t, bool> get_vertex_by_name(const std::string & str,
													  unsigned char flags = Kb16::is_concept | Kb16::is_word) const;

	// ask for node properties

	std::string  get_vertex_name(Kb16_vertex_t u) const {return get(vertex_name, g, u);}
	std::string  get_vertex_gloss(Kb16_vertex_t u) const {return get(vertex_gloss, g, u);}

	// ask for edge preperties

	std::vector<std::string> edge_reltypes(Kb16_edge_t e) const;


	// Nodes can be synsets or words

	bool vertex_is_synset(Kb16_vertex_t u) const;
	bool vertex_is_word(Kb16_vertex_t u) const;

	Kb16_vertex_t InsertNode(const std::string & name, unsigned char flags);

	void read_from_stream (std::ifstream & o);
	std::ofstream & write_to_stream(std::ofstream & o);

	// Private members
	Kb16Graph g;
	std::set<std::string> relsSource;
	std::map<std::string, Kb16_vertex_t> synsetMap; // synset name to vertex id
	std::map<std::string, Kb16_vertex_t> wordMap; // synset name to vertex id

	// Registered relation types

	std::vector<std::string> rtypes;

	std::vector<std::string> notes;        // Command line which created the graph

	// Aux variables

	std::vector<float> out_coefs;          // aux. vector of out-degree coefficients
	size_t N_no_isolated;                  // Number of non-isolated vertices
	char coef_status;                      // 0 invalid
	// 1 calculated without weights
	// 2 calculated with weights
	std::vector<float> static_ppv;         // aux. vector with static prank computation

  private:
	// Singleton
	static Kb16 * p_instance;
	static Kb16 * create();

	// Private methods
	Kb16() : N_no_isolated(0), coef_status(0) {};
	Kb16(const Kb16 &) {};
	Kb16 &operator=(const Kb16 &);
	~Kb16() {};

  };
}

#endif

