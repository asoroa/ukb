// -*-C++-*-

#ifndef MCRGRAPH_H
#define MCRGRAPH_H

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
#include <boost/property_map.hpp>
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
	> KbGraph;


typedef graph_traits<KbGraph>::vertex_descriptor Kb_vertex_t;
typedef graph_traits<KbGraph>::edge_descriptor Kb_edge_t;
typedef graph_traits < KbGraph >::vertices_size_type Kb_vertex_size_t;

class Kb {

public:

  enum {
    is_word = 1,
	is_concept = 2
  } vflags;

  typedef KbGraph boost_graph_t; // the underlying graph type

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

  // 2. create_from_binfile
  //    Load a binary snapshot of the graph into memory

  static void create_from_binfile(const std::string & o);


  // write_to_binfile
  // Write kb graph to a binary serialization file

  void write_to_binfile (const std::string & str);

  // write_to_textfile
  // Write kb graph to a text file

  void write_to_textfile (const std::string & fName);

  // add_from_txt
  // add relations from synsFile to the graph

  void add_from_txt(const std::string & synsFile,
					const std::set<std::string> & rels_source);

  // add_relSource
  // add a new relation source

  void add_relSource(const std::string & str) { if (str.size()) relsSource.insert(str); }

  // Add tokens and link them to their synsets, according to the dictionary.
  // Note: the words are linked to nodes by _directed_ edges

  void add_dictionary(bool with_weight); // Adds all words of the current dictionary
  void add_token(const std::string & str, bool with_weight); // Add just a word (lemma)

  // graph
  // Get the underlying boost graph

  KbGraph & graph() {return g;}

  // Add nodes and relations to the graph

  Kb_vertex_t find_or_insert_synset(const std::string & str);
  Kb_vertex_t find_or_insert_word(const std::string & str);

  Kb_edge_t find_or_insert_edge(Kb_vertex_t u, Kb_vertex_t v, float w );

  // Unlink dangling_nodes (out_degree == 0).
  // Return num. of unlinked nodes

  size_t unlink_dangling();

  // Add relation type to edge

  void edge_add_reltype(Kb_edge_t e, const std::string & rel);

  // Ask for a node

  std::pair<Kb_vertex_t, bool> get_vertex_by_name(const std::string & str,
												  unsigned char flags = Kb::is_concept | Kb::is_word) const;

  // ask for node properties

  std::string  get_vertex_name(Kb_vertex_t u) const {return get(vertex_name, g, u);}
  std::string  get_vertex_gloss(Kb_vertex_t u) const {return get(vertex_gloss, g, u);}

  // ask for edge preperties

  std::vector<std::string> get_edge_reltypes(Kb_edge_t e) const;


  // get static pageRank

  const std::vector<double> & static_prank() const;

  // Nodes can be synsets or words

  bool vertex_is_synset(Kb_vertex_t u) const;
  bool vertex_is_word(Kb_vertex_t u) const;

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

  void filter_ranks_vnames(const std::vector<double> & ranks,
						   std::vector<double> & outranks,
						   std::vector<std::string> & vnames,
						   int filter_mode) const;


  // Add a comment to graph

  void add_comment(const std::string & str);

  // Some useful info

  void display_info(std::ostream & o) const;
  size_t size() const {return num_vertices(g); }

  const std::vector<std::string> & get_comments() const;

  // Get a random vertex

  Kb_vertex_t get_random_vertex() const;

  // Graph algorithms

  bool bfs (Kb_vertex_t source_synset, std::vector<Kb_vertex_t> & synv) const ;

  bool dijkstra (Kb_vertex_t src, std::vector<Kb_vertex_t> & parents) const;

  void pageRank_ppv(const std::vector<double> & ppv_map,
					std::vector<double> & ranks);

  void ppv_weights(const std::vector<double> & ppv);

  // given a source node and a limit (100) return a subgraph by performing a
  // bfs over the graph.
  // Output:
  //  V -> a vector of subgraph nodes
  //  E -> a map representing subgraph edges

  void get_subgraph(const std::string & src,
					std::vector<std::string> & V,
					std::vector<std::vector<std::string> > & E,
					size_t limit = 100);

  std::ostream & dump_graph(std::ostream & o) const;

private:
  // Singleton
  static Kb * p_instance;
  static Kb * create();

  // Private methods
  Kb() : coef_status(0) {};
  Kb(const Kb &) {};
  Kb &operator=(const Kb &);
  ~Kb() {};

  Kb_vertex_t InsertNode(const std::string & name, unsigned char flags);

  void read_from_txt(const std::string & relFile,
					 const std::set<std::string> & rels_source);

  void read_from_stream (std::ifstream & o);
  std::ofstream & write_to_stream(std::ofstream & o);

  // Private members
  KbGraph g;
  std::set<std::string> relsSource;
  std::map<std::string, Kb_vertex_t> synsetMap; // synset name to vertex id
  std::map<std::string, Kb_vertex_t> wordMap; // synset name to vertex id

  // Registered relation types

  std::vector<std::string> rtypes;

  std::vector<std::string> notes;        // Command line which created the graph

  // Aux variables

  std::vector<float> out_coefs;          // aux. vector of out-degree coefficients
  char coef_status;                      // 0 invalid
  // 1 calculated without weights
  // 2 calculated with weights
  std::vector<double> static_ranks;       // aux. vector with static prank computation
  };
}

#endif
