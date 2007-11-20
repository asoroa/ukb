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

enum vertex_wname_t { vertex_wname};  // word name
enum edge_id_t      { edge_id };      // relation id
enum edge_source_t  { edge_source };  // relation source

// typedef std::set<size_t> EdgeId_t;

namespace boost {
  BOOST_INSTALL_PROPERTY(vertex, wname);
  BOOST_INSTALL_PROPERTY(edge, id);
 //  BOOST_INSTALL_PROPERTY(edge, source);
}

typedef adjacency_list <
  boost::listS,
  boost::vecS,
  //  boost::undirectedS,
  boost::bidirectionalS,
  property<vertex_name_t, std::string>,
  property<edge_weight_t, float>
  > McrGraph;
 
//  property<edge_id_t, EdgeId_t> > McrGraph;

//  property<edge_id_t, size_t,
//	   property<edge_source_t, size_t> > > McrGraph;

typedef graph_traits<McrGraph>::vertex_descriptor Mcr_vertex_t;
typedef graph_traits<McrGraph>::edge_descriptor Mcr_edge_t;
typedef graph_traits < McrGraph >::vertices_size_type Mcr_vertex_size_t;

// typedef property_map<McrGraph, vertex_name_t>::type vname_map_t;
// typedef property_map<McrGraph, edge_type_t>::type etype_map_t;
// typedef property_map<McrGraph, edge_source_t>::type esource_map_t;


class Mcr {

public:

  typedef McrGraph boost_graph_t; // the underlying graph type

  //Singleton
  static Mcr & instance();
  static void create_from_txt(const std::string & relFile,
			      const std::string & synsFile,
			      const std::set<std::string> & rels_source);
  static void create_from_binfile(const std::string & o);

  // Member functions

  McrGraph & graph() {return g;}

  size_t size() const {return num_vertices(g); }

  Mcr_vertex_t getRandomVertex() const;

  std::pair<Mcr_vertex_t, bool> getVertexByName(const std::string & str) const;
  std::string  getVertexName(Mcr_vertex_t u) const {return get(vertex_name, g, u);}

  Mcr_vertex_t findOrInsertNode(const std::string & str);
  Mcr_edge_t findOrInsertEdge(Mcr_vertex_t u, Mcr_vertex_t v, float w = 1.0 );


  bool bfs (Mcr_vertex_t source_synset, std::vector<Mcr_vertex_t> & synv) const ;
  //bool bfs (const std::string & source_synset, std::vector<Mcr_vertex_t> & synv) const ;

  bool dijkstra (Mcr_vertex_t src, std::vector<Mcr_vertex_t> & parents) const;

  void pageRank_ppv(const std::vector<float> & ppv_map,
		    std::vector<float> & ranks,
		    bool use_weight);

  void write_to_binfile (const std::string & str) const;

  void display_info(std::ostream & o) const;

  void ppv_weights(const std::vector<float> & ppv);

  // Add tokens a la hughes&ramage

  void add_words(); // Adds all words of w2syn file
  void add_token(const std::string & str); // Add a word

private:
  // Singleton
  static Mcr * p_instance;
  static Mcr * create();

  // Private methods
  Mcr() {};
  Mcr(const Mcr &) {};
  Mcr &operator=(const Mcr &) { return *this; };
  ~Mcr() {};

  void read_from_txt(const std::string & relFile,
		     const std::string & synsFile,
		     const std::set<std::string> & skip_rels);
  void read_from_stream (std::ifstream & o);
  std::ofstream & write_to_stream(std::ofstream & o) const;

  // Private members
  McrGraph g;
  std::set<std::string> relsSource;
  std::map<std::string, Mcr_vertex_t> synsetMap; // synset name to vertex id
  std::map<std::string, int> relMap;     // maps from relation name to relation id
  std::map<int, std::string> relMapInv;  // maps from relation id to relation name

  //std::map<int, int> relType;            // maps from relation id to relation type id
  //std::map<int, int> relInv;             // maps from relation id to inverse id
  //std::map<std::string, int> sourceMap;  // maps from source name to source id's
};


#endif 
