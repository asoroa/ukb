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

  typedef compressed_sparse_row_graph<boost::bidirectionalS,
									  vertex_prop_t,
									  edge_prop_t> KbGraph;

  typedef graph_traits<KbGraph>::vertex_descriptor Kb_vertex_t;
  typedef graph_traits<KbGraph>::vertex_iterator Kb_vertex_iter_t;
  typedef graph_traits<KbGraph>::vertices_size_type Kb_vertex_size_t;
  typedef graph_traits<KbGraph>::edge_descriptor Kb_edge_t;
  typedef graph_traits<KbGraph>::out_edge_iterator Kb_out_edge_iter_t;
  typedef graph_traits<KbGraph>::in_edge_iterator Kb_in_edge_iter_t;

class Kb {

public:

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

  void write_to_textfile (const std::string & fName);

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

  KbGraph & graph() {return *m_g;}

  // Add relation type to edge

  void edge_add_reltype(Kb_edge_t e, const std::string & rel);

  // Ask for a node

  std::pair<Kb_vertex_t, bool> get_vertex_by_name(const std::string & str) const;

  // ask for node properties

  std::string & get_vertex_name(Kb_vertex_t u) {return (*m_g)[u].name;}
  const std::string & get_vertex_name(Kb_vertex_t u) const {return (*m_g)[u].name;}
  //std::string  get_vertex_gloss(Kb_vertex_t u) const {return get(vertex_gloss, g, u);}

  // Get vertices iterator

  std::pair<Kb_vertex_iter_t, Kb_vertex_iter_t> get_vertices() { return boost::vertices(*m_g); }

  // Get out-edges for vertex u

  std::pair<Kb_out_edge_iter_t, Kb_out_edge_iter_t> out_neighbors(Kb_vertex_t u);
  std::pair<Kb_in_edge_iter_t, Kb_in_edge_iter_t> in_neighbors(Kb_vertex_t u);

  bool exists_edge(Kb_vertex_t u, Kb_vertex_t v) const {
	return edge(u, v, *m_g).second;
  }

  Kb_vertex_t edge_source(Kb_edge_t e) const { return source(e, *m_g); }
  Kb_vertex_t edge_target(Kb_edge_t e) const { return target(e, *m_g); }

  // ask for edge preperties

  std::vector<std::string> edge_reltypes(Kb_edge_t e) const;

  float get_edge_weight(Kb_edge_t e) const;
  void set_edge_weight(Kb_edge_t e, float w);

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

  const std::vector<std::string> & get_comments() const;

  // Get a random vertex

  Kb_vertex_t get_random_vertex() const;

  // Graph algorithms

  bool bfs (Kb_vertex_t source_synset, std::vector<Kb_vertex_t> & synv) const ;

  bool dijkstra (Kb_vertex_t src, std::vector<Kb_vertex_t> & parents) const;

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

  Kb_vertex_t InsertNode(const std::string & name, unsigned char flags);

  void read_from_stream (std::istream & o);
  std::ostream & write_to_stream(std::ostream & o) const;

  // Private members
  std::auto_ptr<KbGraph> m_g;
  std::set<std::string> m_relsSource;              // Relation sources
  std::map<std::string, Kb_vertex_t> m_synsetMap; // synset name to vertex id

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
