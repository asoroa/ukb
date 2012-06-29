// -*-C++-*-

#ifndef KBGRAPHCOMMON_H
#define KBGRAPHCOMMON_H

// integer types

#include <boost/version.hpp>
#if BOOST_VERSION < 104100
#error You need boost version >= 1.41 for compiling this version of ukb.
#endif

#include <utility>
#include <vector>
#include <map>
#include <string>
#include <boost/cstdint.hpp>
#include <boost/unordered_map.hpp>
#include <boost/tuple/tuple.hpp>   // for "tie"

namespace ukb {

  class Kb; // forward decl.

  // A class for dealing with edge types

  class etype_t {
  public:

	typedef boost::uint32_t value_type;
	etype_t() : m_strtypes() {}
	etype_t(const etype_t & o) : m_strtypes(o.m_strtypes) {}
	etype_t &operator=(const etype_t & o) {
	  if (&o != this) {
		m_strtypes = o.m_strtypes;
	  }
	  return *this;
	}
	void swap(etype_t & o) {
	  m_strtypes.swap(o.m_strtypes);
	}

	size_t size() const { return m_strtypes.size(); }

	void add_type(const std::string & tstr, value_type & val);
	bool has_type(const std::string & tstr, value_type val) const;
	std::vector<std::string> tvector(value_type val) const;

	friend class Kb;

  private:

	void read_from_stream (std::istream & o);
	std::ostream & write_to_stream(std::ostream & o) const;

	int strpos(const std::string & tstr) const;
	size_t stradd(const std::string & tstr);
	std::vector<std::string> m_strtypes;

  };

  // Vertex and edge properties

  struct vertex_prop_t {
	std::string name;

	vertex_prop_t() : name(std::string()) {}
	vertex_prop_t(const std::string & str) : name(str) {}
  };

  struct edge_prop_t {
	float weight;
	etype_t::value_type etype;

	edge_prop_t() {}
	edge_prop_t(float w) : weight(w), etype(0) {}
	edge_prop_t(float w, etype_t::value_type et) : weight(w), etype(et) {}
  };


  // temporary class used for creating CSR graphs. When reading (or converting)
  // a graph, we first fill this structure and then initialize the graph.

  struct precsr_t {

  private:

	typedef std::pair<size_t, size_t> vertex_pair_t;

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

	std::vector<vertex_pair_t>             E;
	std::vector<vertex_prop_t>                     vProp;
	std::vector<edge_prop_t>                     eProp;

	etype_t                                m_rtypes;

	size_t m_vsize;
	size_t m_esize;

	precsr_t() : m_vsize(0), m_esize(0) {};

	// used by read_kb

	typedef boost::unordered_map<vertex_pair_t, size_t,
								 precsr_edge_hash,
								 precsr_edge_comp> edge_map_t;
	typedef std::map<std::string, size_t> vertex_map_t;

	edge_map_t m_eMap;
	vertex_map_t m_vMap;

	// void sort_edges() {
	//   sort(E.begin(), E.end(), precsr_edge_sort_p);
	// }


	size_t insert_vertex(const std::string & ustr) {

	  bool insertedP;
	  vertex_map_t::iterator vit;

	  boost::tie(vit, insertedP) = m_vMap.insert(std::make_pair(ustr, m_vsize));
	  if(insertedP) {
		vProp.push_back(vertex_prop_t(ustr));
		++m_vsize;
	  }
	  return vit->second;
	}

	size_t insert_edge(const std::string & ustr,
					   const std::string & vstr,
					   float w,
					   etype_t::value_type etype) {

	  size_t u = insert_vertex(ustr);
	  size_t v = insert_vertex(vstr);

	  edge_map_t::iterator eit;
	  edge_map_t::key_type k = std::make_pair(u, v);
	  bool insertedP;

	  boost::tie(eit, insertedP) = m_eMap.insert(std::make_pair(k, m_esize));
	  if(insertedP) {
		E.push_back(k);
		eProp.push_back(edge_prop_t(w));
		++m_esize;
	  }
	  eProp[eit->second].etype |= etype;

	  return eit->second;
	}

	size_t insert_edge(const std::string & ustr,
					   const std::string & vstr,
					   float w,
					   const std::string & rtype) {
	  size_t eidx = insert_edge(ustr, vstr, w, etype_t::value_type(0));
	  // add edge type
	  if (rtype.size())
		m_rtypes.add_type(rtype,eProp[eidx].etype);
	  return eidx;
	}
  };
}

#endif

