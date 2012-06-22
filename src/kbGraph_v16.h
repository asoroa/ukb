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

  template<class vertex16_prop_t, class edge16_prop_t>
  struct precsr16_t {

  private:

	typedef std::pair<std::size_t, std::size_t> vertex_pair_t;

	struct precsr_edge_comp
	  : std::binary_function<vertex_pair_t, vertex_pair_t, bool> {
	  bool operator()(const vertex_pair_t & e1, const vertex_pair_t & e2) const {
		return e1.first == e2.first && e1.second == e2.second;
	  }
	};

	struct precsr_edge_hash
	  : std::unary_function<vertex_pair_t, std::size_t> {

	  std::size_t operator()(const vertex_pair_t & e) const {

		std::size_t seed = 0;
		boost::hash_combine(seed, e.first);
		boost::hash_combine(seed, e.second);
		return seed;
	  }
	};

	struct precsr_edge_sort_p {
	  precsr_edge_sort_p() {}
	  int operator()(const vertex_pair_t & e1,
					 const vertex_pair_t & e2) {
		return e1.first - e2.first;
	  }
	};

  public:

	std::vector<vertex_pair_t>                     E;
	std::vector<vertex16_prop_t>                   vProp;
	std::vector<edge16_prop_t>                     eProp;

	std::size_t m_vsize;
	std::size_t m_esize;

	precsr16_t() : m_vsize(0), m_esize(0) {};

	// used by read_kb

	typedef typename boost::unordered_map<vertex_pair_t, std::size_t,
										  precsr_edge_hash,
										  precsr_edge_comp> edge_map_t;
	typedef typename std::map<std::string, std::size_t> vertex_map_t;

	edge_map_t m_eMap;
	vertex_map_t m_vMap;

	// void sort_edges() {
	//   sort(E.begin(), E.end(), precsr_edge_sort_p);
	// }


	// Add a relation type to edge. Also, update kbC rtype list.

	void edge_add_reltype(std::size_t e_i,
						  boost::uint32_t rt) {

	  boost::uint32_t m = eProp[e_i].rtype;
	  m |= rt;
	  eProp[e_i].rtype = m;
	}

	std::size_t insert_vertex(const std::string & ustr) {

	  bool insertedP;
	  vertex_map_t::iterator vit;

	  boost::tie(vit, insertedP) = m_vMap.insert(std::make_pair(ustr, m_vsize));
	  if(insertedP) {
		vProp.push_back(vertex16_prop_t(ustr));
		++m_vsize;
	  }
	  return vit->second;
	}

	std::size_t insert_edge(const std::string & ustr,
							const std::string & vstr,
							float w,
							boost::uint32_t rtype) {

	  std::size_t u = insert_vertex(ustr);
	  std::size_t v = insert_vertex(vstr);

	  typename edge_map_t::iterator eit;
	  typename edge_map_t::key_type k = std::make_pair(u, v);
	  bool insertedP;

	  boost::tie(eit, insertedP) = m_eMap.insert(std::make_pair(k, m_esize));
	  if(insertedP) {
		E.push_back(k);
		eProp.push_back(edge16_prop_t(w));
		++m_esize;
	  }
	  // add rtype
	  edge_add_reltype(eit->second, rtype);
	  return eit->second;
	};
  };

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

  std::vector<std::string> get_edge_reltypes(Kb16_edge_t e) const;


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

