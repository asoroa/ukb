#include "kGraph.h"

#include <list>
#include <algorithm>

// bfs 

#include <boost/graph/visitors.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/pending/indirect_cmp.hpp>
#include <boost/pending/integer_range.hpp>
#include <boost/graph/graph_utility.hpp> // for boost::make_list

using namespace std;
using namespace boost;

////////////////////////////////////////////////////////////////
// bfs over DisambG

struct dis_bfs_init:public base_visitor<dis_bfs_init> {
public:
  dis_bfs_init(size_t *v):m_v(v) { }
  typedef on_initialize_vertex event_filter;
  inline void operator()(Dis_vertex_t u, const DisambG & g)
  {
    m_v[u] = 0;
  }
  size_t *m_v;
};

// Record the distance of minimum path

struct dis_bfs_pred:public base_visitor<dis_bfs_pred> {
public:
  dis_bfs_pred(size_t *v):m_v(v) { }
  typedef on_tree_edge event_filter;
  inline void operator()(Dis_edge_t e, const DisambG & g) {
    Dis_vertex_t u = source(e, g), v = target(e, g);
    m_v[v] = m_v[u] + 1;
  }
  size_t *m_v;
};

bool disg_bfs (DisambG & g,
	       Dis_vertex_t src, 
	       std::vector<size_t> & dist) {
  vector<size_t>(num_vertices(g)).swap(dist);  // reset parents
  breadth_first_search(g, 
		       src,
		       boost::visitor(boost::make_bfs_visitor
				      (boost::make_list(dis_bfs_init(&dist[0]),
							dis_bfs_pred(&dist[0])))));
  return true;
}

////////////////////////////////////////////////////////////////
// Constructor

KGraph::KGraph(DisambG & disg, CSentence & cs) {

//   vector<string> d_vertices;
//   map<Dis_vertex_t, Dis_vertex_t> dis_to_kg;
//   cs.distinguished_synsets(d_vertices);

//   vector<Dis_vertex_t>::iterator it_i, end;
//   Dis_vertex_t aux;
//   it_i = d_vertices.begin();
//   end = d_vertices.end();

//   for(; it_i != end; ++it_i) {
//     aux = add_vertex(g);
//     dis_to_kg[*it_i] = aux;
//     put(vertex_name, g, aux,
// 	get(vertex_name, disg, *it_i));
//     put(vertex_wname, g, aux,
// 	get(vertex_wname, disg, *it_i));
//   }

//   it_i = d_vertices.begin();

//   for(; it_i != end; ++it_i) {
//     vector<Dis_vertex_t> dist;
//     disg_bfs(disg, *it_i, dist);
//     vector<Dis_vertex_t>::iterator it_j = it_i;
//     ++it_j;
//     for(; it_j != end; ++it_j) {
//       // Add edge (*it_i, *it_j) with weigth dist[*it_j]
//       assert(dist[*it_j]);
//       add_edge(dis_to_kg[*it_i], dis_to_kg[*it_j], dist[*it_j], g);
//     }
//   }
}

