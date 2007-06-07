// -*-C++-*-

// common.h 
// + global variables
// + constants
// + singletons

#ifndef COMMON_H
#define COMMOM_H

#include <string>
#include <iterator>
#include <map>
#include <vector>
#include <iosfwd>
//#include <algorithm>
#include <numeric>

#include <boost/graph/graphviz.hpp>

/////////////////////
// Common functions

template<class T>
void read_atom_from_stream(std::ifstream & is,
			   T & a) {
  is.read(reinterpret_cast<char *>(&a), sizeof(a));  
}

template<> 
inline void read_atom_from_stream <std::string> (std::ifstream & is, std::string & str) {
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
void read_map_from_stream(std::ifstream & is, Map & map) {
 
  typedef typename Map::key_type key_type;
  typedef typename Map::mapped_type data_type;
  typedef typename Map::value_type value_type;

  key_type key;
  data_type value;

  size_t mapSize;
  size_t i;

  read_atom_from_stream(is, mapSize);

  for(i=0; i< mapSize; ++i) {
    read_atom_from_stream(is, key);
    read_atom_from_stream(is, value);
    map.insert(value_type(key, value));
  }
}

template<class Vector> 
void read_vector_from_stream(std::ifstream & is, Vector & v) {

  size_t vSize;
  size_t i;
  typename Vector::value_type value;

  read_atom_from_stream(is, vSize);

  for(i=0; i< vSize; ++i) {
    read_atom_from_stream(is, value);
    v.push_back(value);
  }
}

template<class T>
std::ofstream & write_atom_to_stream(std::ofstream & o,
				     const T & a) {

  o.write(reinterpret_cast<const char *>(&a), sizeof(a));

  return o;
}

template<>
inline std::ofstream & write_atom_to_stream <std::string> (std::ofstream & o,
							   const std::string & str) {

  size_t len;
  
  len = str.size();
  o.write(reinterpret_cast<const char *>(&len), sizeof(len));
  if (len)
    o.write(reinterpret_cast<const char *>(str.c_str()), len);  
  return o;
}

template<class Map>
std::ofstream & write_map_to_stream(std::ofstream & o,
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

template<class Vector>
std::ofstream & write_vector_to_stream(std::ofstream & o,
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
  const std::string name;
};

template <class ValType>
inline my_name_writer<ValType>
make_my_writer(ValType n, const std::string & ize) {
  return my_name_writer<ValType>(n, ize);
}

#endif
