// -*-C++-*-

// common.h
// + global variables
// + constants
// + singletons

#ifndef COMMON_H
#define COMMON_H

#include <string>
#include <iterator>
#include <map>
#include <vector>
#include <set>
#include <iosfwd>
//#include <algorithm>
#include <numeric>

#include <boost/graph/graphviz.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>

#include <boost/version.hpp>
#if BOOST_VERSION > 103800
#include <boost/property_map/property_map.hpp>
#else
#include <boost/property_map.hpp>
#endif

#include <boost/graph/properties.hpp>

// Stuff for generating random numbers

#include <boost/random.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>

using boost::adjacency_list;
using boost::graph_traits;
using boost::property;
using boost::property_map;

/////////////////////
// Common functions

namespace ukb {
	/////////////////////////////////////////////////////////////////////
	// Streaming

	template<class T>
	void read_atom_from_stream(std::istream & is,
							   T & a) {
		is.read(reinterpret_cast<char *>(&a), sizeof(a));
	}

	template<>
	inline void read_atom_from_stream <std::string> (std::istream & is, std::string & str) {
		size_t len;

		is.read(reinterpret_cast<char *>(&len), sizeof(len));
		if (len) {
			char * aux= new char [len];
			is.read(aux, len);
			str.assign(aux, aux + len);
			delete[] aux;
		}
	}

	template<class Map>
	void read_map_from_stream(std::istream & is, Map & map) {

		typedef typename Map::key_type key_type;
		typedef typename Map::mapped_type data_type;
		typedef typename Map::value_type value_type;

		key_type key;
		data_type value;

		size_t mapSize;
		size_t i;

		read_atom_from_stream(is, mapSize);
		if (!is) return;
		for(i=0; i< mapSize; ++i) {
			read_atom_from_stream(is, key);
			read_atom_from_stream(is, value);
			map.insert(value_type(key, value));
		}
	}

	template<class Vector>
	void read_vector_from_stream(std::istream & is, Vector & v) {

		size_t vSize;
		size_t i;
		Vector auxV;
		typename Vector::value_type value;

		read_atom_from_stream(is, vSize);
		if (!is) return;
		for(i=0; i< vSize; ++i) {
			read_atom_from_stream(is, value);
			auxV.push_back(value);
		}
		// swap from temporary, so that it doesn't take more space than neccesary
		Vector (auxV).swap(v);
	}

	template<class Set>
	void read_set_from_stream(std::istream & is, Set & v) {

		size_t vSize;
		size_t i;
		typename Set::value_type value;

		read_atom_from_stream(is, vSize);
		if (!is) return;
		for(i=0; i< vSize; ++i) {
			read_atom_from_stream(is, value);
			v.insert(value);
		}
	}

	template<class T>
	std::ostream & write_atom_to_stream(std::ostream & o,
										const T & a) {

		o.write(reinterpret_cast<const char *>(&a), sizeof(a));

		return o;
	}

	template<>
	inline std::ostream & write_atom_to_stream <std::string> (std::ostream & o,
															  const std::string & str) {

		size_t len;

		len = str.size();
		o.write(reinterpret_cast<const char *>(&len), sizeof(len));
		if (len)
			o.write(reinterpret_cast<const char *>(str.c_str()), len);
		return o;
	}

	template<class Map>
	std::ostream & write_map_to_stream(std::ostream & o,
									   const Map & m) {

		typename Map::const_iterator map_it, map_end;
		size_t mapSize;

		mapSize = m.size();

		o.write(reinterpret_cast<const char *>(&mapSize), sizeof(mapSize));
		if(mapSize) {
			map_end = m.end();
			for(map_it = m.begin(); map_it != map_end; ++map_it) {
				write_atom_to_stream(o, map_it->first);
				write_atom_to_stream(o, map_it->second);
			}
		}
		return o;
	}

	template<class Set>
	std::ostream & write_set_to_stream(std::ostream & o,
									   const Set & s) {

		typename Set::const_iterator set_it, set_end;
		size_t setSize;

		setSize = s.size();

		o.write(reinterpret_cast<const char *>(&setSize), sizeof(setSize));
		if(setSize) {
			set_end = s.end();
			for(set_it = s.begin(); set_it != set_end; ++set_it) {
				typename Set::value_type aux = *set_it;
				write_atom_to_stream(o, aux);
			}
		}
		return o;
	}

	template<class Vector>
	std::ostream & write_vector_to_stream(std::ostream & o,
										  const Vector & v) {

		typename Vector::const_iterator map_it, map_end;
		size_t vSize;

		vSize = v.size();

		write_atom_to_stream(o, vSize);
		if(vSize) {
			map_end = v.end();
			for(map_it = v.begin(); map_it != map_end; ++map_it) {
				write_atom_to_stream(o, *map_it);
			}
		}
		return o;
	}

	/////////////////////////////////////////////////////////////////////
	// Connected vertices of undirected graphs

	template<typename G>
	struct vertex_is_connected {
		typedef typename boost::graph_traits<G>::vertex_descriptor Vertex_G;
		const G & g;
		vertex_is_connected(const G & g_) : g(g_) {}
		bool operator() (Vertex_G v) {
			return (out_degree(v,g) > 0);
		}
	};

	template<typename G>
	struct sum_connected_vertex_ {
		typedef typename boost::graph_traits<G>::vertex_descriptor Vertex_G;
		const G & g;
		sum_connected_vertex_(const G & g_) : g(g_) {}
		size_t operator() (size_t result, const Vertex_G & v) {
			return (boost::out_degree(v,g) > 0) ? result + 1 : result;
		}
	};


	template<typename G>
	size_t num_connected_vertices(const G & g) {

		// verbose mode.
		//
		typename boost::graph_traits<G>::vertex_iterator vIt, vItEnd;
		tie(vIt, vItEnd) = boost::vertices(g);
		size_t erpin_onak = std::accumulate(vIt, vItEnd, size_t(0), sum_connected_vertex_<G>(g));
		return erpin_onak;
	}

	/////////////////////////////////////////////////////////////////////
	// Random

	extern boost::minstd_rand global_mersenne;

	template <typename Generator>
	struct random_functor {
		random_functor(Generator& g) : g(g) { }
		std::size_t operator()(std::size_t n) {
			boost::uniform_int<std::size_t> distrib(0, n-1);
			boost::variate_generator<boost::mt19937&, boost::uniform_int<std::size_t> >
				x(g, distrib);
			return x();
		}
		Generator& g;
	};

	template<typename Generator>
	inline random_functor <Generator>
	make_random_functor(Generator & g) {
		return random_functor<Generator>(g);
	}

	int g_randTarget(int Target);

	////////////////////////////////////////////////////////////////////
	// Remove isolated vertices of a graph

	// template<class G>
	// void remove_isolated_vertices(G & g) {

	//   typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;
	//   typedef typename boost::graph_traits<G>::vertex_iterator vertex_iterator;
	//   typedef typename boost::graph_traits<G>::edge_iterator edge_iterator;

	//   G newg;
	//   std::map<vertex_decriptor, vertex_descriptor> nMap;
	//   vertex_iterator v_it, v_end;
	//   for(boost::tie(v_it, v_end) = vertices(g); v_it != v_end; ++v_it) {
	//     if (boost::out_degree(*v_it, g) || boost::in_degree(*v_it, g)) {

	//     }
	//   }

	// }


	/////////////////////////////////////////////////////////////////////
	// graphviz

	template < class ValType >
	class my_name_writer {
	public:
		my_name_writer(ValType _val, const std::string & _name) : val(_val), name(_name) {}
		template <class VertexOrEdge>
		void operator()(std::ostream& out, const VertexOrEdge& v) const {
			out << "["<< name <<"=\"" << val[v] << "\"]";
		}
	private:
		ValType val;
		const std::string & name;
	};

	template < class ValType1, class ValType2, class ValType3 >
	class my_name_writer3 {
	public:
		my_name_writer3(ValType1 _val1, ValType2 _val2, ValType3 _val3, const std::string & _name1, const std::string & _name2, const std::string & _name3)
			: val1(_val1), val2(_val2), val3(_val3), name1(_name1), name2(_name2), name3(_name3) {}
		template <class VertexOrEdge>
		void operator()(std::ostream& out, const VertexOrEdge& v) const {
			out << "["<< name1 <<"=\"" << val1[v] << "\" " << name2 << "=\"" <<  val2[v] << "\" " << name3 << "=\"" <<  val3[v] << "\"]";
		}
	private:
		ValType1 val1;
		ValType2 val2;
		ValType3 val3;
		const std::string & name1;
		const std::string & name2;
		const std::string & name3;
	};

	template <class ValType>
	inline my_name_writer<ValType>
	make_my_writer(ValType n, const std::string & ize) {
		return my_name_writer<ValType>(n, ize);
	}

	template <class ValType1, class ValType2, class ValType3>
	inline my_name_writer3<ValType1, ValType2, ValType3>
	make_my_writer3(ValType1 n1, ValType2 n2, ValType3 n3, const std::string & name1, const std::string & name2, const std::string & name3) {
		return my_name_writer3<ValType1, ValType2, ValType3> (n1, n2, n3, name1, name2, name3);
	}

	/////////////////////////////////////////////////////////////////////
	// Display vector

	template<class T>
	std::ostream & writeV(std::ostream & o, const std::vector<T> & v) {

		typename std::vector<T>::const_iterator vIt = v.begin();
		typename std::vector<T>::const_iterator vItEnd = v.end();
		if (vIt != vItEnd) {
			--vItEnd;
			o << "(";
			std::copy(vIt, vItEnd, std::ostream_iterator<T>(o, ","));
			o << *vItEnd << ")";
		}
		return o;
	}

	template<class T>
	std::ostream & writeS(std::ostream & o, const std::set<T> & v) {

		typename std::set<T>::const_iterator vIt = v.begin();
		typename std::set<T>::const_iterator vItEnd = v.end();
		if (vIt != vItEnd) {
			--vItEnd;
			o << "(";
			std::copy(vIt, vItEnd, std::ostream_iterator<T>(o, ","));
			o << *vItEnd << ")";
		}
		return o;
	}

	///////////////////////////////////////////////////////////////////////
	// Normalize probability vector
	//

	template<class T>
	bool normalize_pvector(std::vector<T> & v) {
		T sum(0);
		typename std::vector<T>::iterator it = v.begin();
		typename std::vector<T>::iterator end = v.end();

		for(; it != end; ++it)
			sum += *it;

		if (sum == T(0)) return false;
		T factor(T(1) / sum);

		for(it = v.begin(); it != end; ++it)
			*it *= factor;
		return true;
	}


	/////////////////////////////////////////////////////////////////////
	// copy_if

	template <typename InputIt, typename OutputIt, typename Predicate>
	OutputIt copy_if (InputIt first, InputIt last,
					  OutputIt result, Predicate pred)
	{
		while (first != last) {
			if (pred(*first))  *result++ = *first;
			++first;
		}
		return result;
	}

	//////////////////////////////////////////////////////////////////////
	// split function

	std::vector<std::string> split(const std::string & str, const std::string & delims);


	// join function

	std::string join(const std::string &delim,
					 std::vector<std::string>::const_iterator it,
					 std::vector<std::string>::const_iterator end);

	std::string join(const std::string & delim, const std::vector<std::string> & V);


	// read non blank line

	std::istream & read_line_noblank(std::istream & is, std::string & line, size_t & l_n);

	// trim trailing and leading spaces
	void trim_spaces(std::string &l);

	// Common misc. functions
	void set_pr_convergence(size_t iterations, float threshold);

}

#endif
