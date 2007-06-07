#include "disambGraph.h"
#include "common.h"
#include "w2syn.h"
#include "csentence.h"

#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <string>
#include <map>
#include <iterator>
#include <ostream>


// bfs 

#include <boost/graph/visitors.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/pending/indirect_cmp.hpp>
#include <boost/pending/integer_range.hpp>
#include <boost/graph/graph_utility.hpp> // for boost::make_list

using namespace std;
using namespace boost;


////////////////////////////////////////////////////////////////////////////////
// Disamb fill functions


//void add_vertex(){}

DisambGraph::DisambGraph() {
}

void DisambGraph::add_disamb_edge(Dis_vertex_t u, Dis_vertex_t v) {

  Dis_edge_t e;
  bool existP;

  map<Mcr_vertex_t, Dis_vertex_t>::iterator it;

  tie(e, existP) = edge(u, v, g);
  if(!existP) {
    // new edge
    tie(e, existP) = add_edge(u, v, g);
    put(edge_freq, g, e, 1);
  } else {
    // edge already there. Increase freq.
    size_t freq = get(edge_freq, g, e);
    ++freq;
    put(edge_freq, g, e, freq);    
  }
}

void DisambGraph::add_vertices_mcr_path(vector<Dis_vertex_t>::iterator v_it, 
					vector<Dis_vertex_t>::iterator v_end) {

  map<Dis_vertex_t, Dis_vertex_t>::iterator map_it;
  map<Dis_vertex_t, Dis_vertex_t>::iterator map_end = mcr2dis.end();
  bool insertedP;

  McrGraph & mcr_g = Mcr::instance().graph();

  for(;v_it != v_end; ++v_it) {
    tie(map_it, insertedP) = mcr2dis.insert(make_pair(*v_it, Dis_vertex_t()));
    if (insertedP) {
      Dis_vertex_t v = add_vertex(g);
      put(vertex_name, g, v, get(vertex_name, mcr_g, *v_it));
      put(vertex_wname, g, v, get(vertex_wname, mcr_g, *v_it));
      //put(vertex_rank, g, v, numeric_limits<float>::min());
      put(vertex_rank, g, v, 0.0f);
      map_it->second = v;
    }
    *v_it = map_it->second; // update source id
  }
}



void DisambGraph::fill_graph(Mcr_vertex_t src,
			     Mcr_vertex_t tgt,
			     const std::vector<Mcr_vertex_t> & parents) {

  vector<Dis_vertex_t> path;
  Mcr_vertex_t pred;

  pred = parents[tgt];
  while(tgt != pred) {
    path.push_back(static_cast<Dis_vertex_t>(tgt));
    tgt = pred;
    pred = parents[tgt];
  }
  if (tgt != static_cast<Dis_vertex_t>(src)) return;
  
  vector<Dis_vertex_t>::iterator path_it = path.begin();
  vector<Dis_vertex_t>::iterator path_end = path.end();
  vector<Dis_vertex_t>::iterator path_prev;

  if (path_it == path_end) return;

  add_vertices_mcr_path(path_it, path_end);

  path_prev = path_it;
  ++path_it;
  while(path_it != path_end) {
    add_disamb_edge(*path_prev, *path_it);
    path_prev = path_it;
    ++path_it;
  }
}

////////////////////////////////////////////////////////////////
// Global functions

void fill_disamb_synset(Mcr_vertex_t src,
			vector<CWord>::const_iterator s_it,
			vector<CWord>::const_iterator s_end,
			DisambGraph & dgraph) {

  //bfs from src
  std::vector<Mcr_vertex_t> parents;
  Mcr::instance().bfs(src, parents);
  
  //fill disamb graph
  for(;s_it != s_end; ++s_it) {
    vector<Mcr_vertex_t>::const_iterator tg_it = s_it->begin();
    vector<Mcr_vertex_t>::const_iterator tg_end = s_it->end();
    for(;tg_it != tg_end; ++tg_it) {
      dgraph.fill_graph(src, *tg_it, parents);
    }
  }

}

void fill_disamb_graph(const CSentence &cs, DisambGraph & dgraph) {

  vector<CWord>::const_iterator cw_it = cs.begin();
  vector<CWord>::const_iterator cw_end = cs.end();

  //if (cw_it == cw_end) return;
  //cw_end--;
  while(cw_it != cw_end) {
    vector<Mcr_vertex_t>::const_iterator sset_it = cw_it->begin();
    vector<Mcr_vertex_t>::const_iterator sset_end = cw_it->end();
    ++cw_it; // point to next word
    for(;sset_it != sset_end; ++sset_it) {
      fill_disamb_synset(*sset_it, cw_it, cw_end, dgraph);
    }
  }
}


// Given a CSentence where synsets are represented as vertices in the MCR,
// transform the vertex-id' so they point to the dgraph


struct MaxSynsLast {
  int operator()(const CWord & a, const CWord & b) {
    return a.size() < b.size();
  }
};

void DisambGraph::transform_csentence(CSentence & cs) const {
  
  vector<CWord>::iterator cw_it = cs.begin();
  vector<CWord>::iterator cw_end = cs.end();
  for(; cw_it != cw_end; ++cw_it) {
    vector<Mcr_vertex_t> & old_syns = cw_it->get_syns_vector();
    vector<Mcr_vertex_t>::iterator syn_it = old_syns.begin();
    vector<Mcr_vertex_t>::iterator syn_end = old_syns.end();
    vector<Dis_vertex_t> new_syns;
    for(; syn_it != syn_end; ++syn_it) {
      map<Mcr_vertex_t, Dis_vertex_t>::const_iterator m2d_it = mcr2dis.find(*syn_it);
      if (m2d_it != mcr2dis.end()) {
 	new_syns.push_back(m2d_it->second);
      }
    }
    old_syns.swap(new_syns);
  }
}

// Warning: Rank scores of dgraph must be previously computed.
//

struct SortByRank {
  DisambG & g;
  SortByRank(DisambG & g_) : g(g_) {};
  int operator() (const Dis_vertex_t & u, const Dis_vertex_t & v) {
    // Descending order !
    return get(vertex_rank, g, v) < get(vertex_rank, g, u);
  }
};

void disamb_csentence(CSentence & cs, DisambGraph & dgraph) {

  vector<CWord>::iterator cw_it = cs.begin();
  vector<CWord>::iterator cw_end = cs.end();
  for(; cw_it != cw_end; ++cw_it) {
    sort(cw_it->begin(), cw_it->end(), SortByRank(dgraph.graph()));
  }  
}

ostream & print_csent(ostream & o, CSentence & cs, DisambGraph & dgraph) {
  DisambG & g = dgraph.graph();
  vector<CWord>::iterator cw_it = cs.begin();
  vector<CWord>::iterator cw_end = cs.end();
  for(; cw_it != cw_end; ++cw_it) {
    if (cw_it->size() == 0) continue;
    if (!cw_it->is_distinguished()) continue;
    o << cw_it->id() << " ";
    o << get(vertex_name, g, *(cw_it->begin()));
    o << endl;
  }
  return o;
}

ostream & print_complete_csent(ostream & o, CSentence & cs, DisambGraph & dgraph) {
  DisambG & g = dgraph.graph();
  vector<CWord>::iterator cw_it = cs.begin();
  vector<CWord>::iterator cw_end = cs.end();
  for(; cw_it != cw_end; ++cw_it) {
    o << cw_it->word() << "{";
    vector<Mcr_vertex_t>::iterator syn_it = cw_it->begin();
    vector<Mcr_vertex_t>::iterator syn_end = cw_it->end();
    if (syn_it != syn_end) {
      --syn_end;
      for(; syn_it != syn_end; ++syn_it) {
	o << get(vertex_name, g, *syn_it) << ":" << get(vertex_rank, g, *syn_it) << " ,";
      }
      o << get(vertex_name, g, *syn_end) << ":" << get(vertex_rank, g, *syn_it);
    }
    o << "}" << endl;
  }
  return o;
}

// HITS

// Hits

void hits_norm(vector<float> & v) {

  vector<float>::iterator it;
  vector<float>::iterator end = v.end();

  float sum = 0.0;
  for(it = v.begin(); it != end; ++it)
    sum += *it;
  assert(sum);
  float coef = 1.0f / sum;
  for(it = v.begin(); it != end; ++it) 
    *it *= coef;
}

template<typename HProp, typename AProp>
void hits_update_auth (DisambG & g,
		       vector<Dis_vertex_t>::iterator vit,
		       vector<Dis_vertex_t>::iterator end,
		       HProp hProp, AProp aProp) {
  // I operation
  // x^{p} = \sum_{q:(q,p) \in E} w_{qp}*y^{q}

  //typename graph_traits<G>::vertex_iterator vit, end;
  //tie(vit, end) =  vertices(g);
  for(; vit != end; ++vit) {
    float auth = 0;
    graph_traits<DisambG>::in_edge_iterator e, e_end;
    tie(e, e_end) = in_edges(*vit, g);
    for(; e != e_end; ++e) {
      auth+= get(hProp, source(*e, g)) * get(edge_freq, g, *e);
    }
    put(aProp, *vit, auth);
  }
}

template<typename HProp, typename AProp>
void hits_update_hub (DisambG & g,
		      vector<Dis_vertex_t>::iterator vit,
		      vector<Dis_vertex_t>::iterator end,
		      HProp hProp, AProp aProp) {
  // O operation
  // y^{p} = \sum_{q:(p,q) \in E} w_{pq}*x^{q}

  //typename graph_traits<G>::vertex_iterator vit, end;
  //tie(vit, end) =  vertices(g);

  for(; vit != end; ++vit) {
    float hub = 0.0;
    graph_traits<DisambG>::out_edge_iterator e, e_end;
    tie(e, e_end) = out_edges(*vit, g);
    for(; e != e_end; ++e) {
      hub+= get(aProp, target(*e, g)) * get(edge_freq, g, *e);
    }
    put(hProp, *vit, hub);
  }
}

template<typename HProp, typename AProp>
void hits_init_ranks(DisambG & g, 
		     vector<Dis_vertex_t>::iterator vit,
		     vector<Dis_vertex_t>::iterator end,
		     HProp hProp, AProp aProp, size_t n) {

  assert(n);
  float init_v = sqrt(n) / n;

  //graph_traits<DisambG>::vertex_iterator vit, end;
  //tie(vit, end) =  vertices(g);
  for(; vit != end; ++vit) {
    put(aProp, *vit, init_v);
    put(hProp, *vit, init_v);
  }
}

void hits_iterate(DisambG & g,
		  vector<Dis_vertex_t>::iterator vit,
		  vector<Dis_vertex_t>::iterator end,
		  vector<float> & hProp,
		  vector<float> & aProp,
		  size_t iterations) {

  hits_init_ranks(g, vit, end, &hProp[0], &aProp[0], num_connected_vertices(g));
  for(size_t i = 0; i<iterations; ++i) {
    hits_update_auth(g, vit, end, &hProp[0], &aProp[0]);
    hits_update_hub(g, vit, end, &hProp[0], &aProp[0]);
    hits_norm(aProp);
    hits_norm(hProp);
  }
}

void hits(DisambG & g) {

  vector<float> aRank(num_vertices(g), 0.0f);
  vector<float> hRank(num_vertices(g), 0.0f);

  vector<Dis_vertex_t> V(num_connected_vertices(g));

  graph_traits<DisambG>::vertex_iterator vIt, vItEnd;
  tie(vIt, vItEnd) = vertices(g);
  copy_if(vIt, vItEnd, V.begin(), vertex_is_connected<DisambG>(g));

  hits_iterate(g, V.begin(), V.end(), hRank, aRank, 50); // 50 iterations

  vector<Dis_vertex_t>::const_iterator it = V.begin();
  vector<Dis_vertex_t>::const_iterator end = V.end();
  for(; it != end; ++it)
    put(vertex_rank, g, *it, hRank[*it]);
}

/////////////////////////////////////////////////////////////////////
// pageRank


//
// PageRank iteration
//
// it solves pageRank with the so called power method see  
// @article{langville04,
// title = {{Deeper Inside PageRank}},
// author = {A.N. Langville and C.D. Meyer},
// journal = {Internet Mathematics},
// number = {3},
// pages = {335--380},
// volume = {1},
// year = {2004},
// } 
//
// Note: 
//    * as G is an undirected graph, there is no dangling node
//    * by now, an uniform v^{T} vector is used, which is just (1/N)*e^{T}


template<typename G, typename map1, typename map2>
void update_pRank_weight(G & g,
			 vector<typename graph_traits<G>::vertex_descriptor> & V,
			 float Ncoef,
			 const vector<float> & out_coef, // 1/N
			 float dfactor,
			 const map1 & rank_map1,
			 map2 & rank_map2) {

  typedef typename graph_traits<G>::vertex_descriptor vertex_descriptor;
  
  typename vector<vertex_descriptor>::iterator v = V.begin();
  typename vector<vertex_descriptor>::iterator end = V.end();
  for (; v != end; ++v) {
    float rank=0.0;
    typename graph_traits<G>::in_edge_iterator e, e_end;
    tie(e, e_end) = in_edges(*v, g);
    for(; e != e_end; ++e) {
      vertex_descriptor u = source(*e, g);
      rank += rank_map1[u] * get(edge_freq, g, *e) * out_coef[u];
    }
    rank_map2[*v] = dfactor*rank + (1-dfactor)*Ncoef;
  }
}

template<typename G, typename map1, typename map2>
void pageRank_iterate(G & g,
		      vector<typename graph_traits<G>::vertex_descriptor> & V,
		      size_t N,
		      const vector<float> & out_coef,
		      map1 & rank_map1,
		      map2 & rank_map2,
		      int iterations) {
  int erpin_kop = V.size();
  float damping = 0.85;

  // Initialize rank_map1 appropriately
  {
    typename graph_traits<G>::vertex_iterator v, end;
    for (tie(v, end) = vertices(g); v != end; ++v) {
      rank_map1[*v]=1.0/(float) erpin_kop;
    }
  }

  // Continue iterating until the termination condition is met
  // or the algorithm converges (a vertex's rank does not change)
  bool to_map_2 = true;
  while(iterations--) {
    //if ((akOrok::verbose)) cerr << iterations << "...";
    // Update to the appropriate rank map
    if (to_map_2)
      update_pRank_weight(g, V, 1.0/static_cast<float>(N), out_coef,damping, rank_map1, rank_map2);
    else
      update_pRank_weight(g, V, 1.0/static_cast<float>(N), out_coef, damping, rank_map2, rank_map1);
    // The next iteration will reverse the update mapping
    to_map_2 = !to_map_2;
    //cerr << sum << endl;
  }

  // We stopped after writing the latest results to rank_map2, so copy the
  // results back to rank_map1 for the caller
  if (!to_map_2) {
    typename graph_traits<G>::vertex_iterator v, end;
    for (tie(v, end) = vertices(g); v != end; ++v) {
      rank_map1[*v] = rank_map2[*v];
    }
  }
}

void init_out_coefs(const DisambG & g,
		    const vector<Dis_vertex_t> & V,
		    vector<float> & W) {
  vector<Dis_vertex_t>::const_iterator v = V.begin();
  vector<Dis_vertex_t>::const_iterator end = V.end();
  for (; v != end; ++v) {
    graph_traits<DisambG>::out_edge_iterator e, e_end;
    tie(e, e_end) = out_edges(*v, g);
    size_t total_w = 0;
    for(; e != e_end; ++e) {
      total_w += get(edge_freq, g, *e);
    }
    assert(total_w);
    W[*v] = 1.0f / static_cast<float>(total_w);
  }
}

void pageRank(DisambG & g) {

  vector<float> map_tmp(num_vertices(g), 0.0f);
  vector<float> out_coefs(num_vertices(g), 0);

  size_t N = num_connected_vertices(g);

  vector<Dis_vertex_t> V(N);

  graph_traits<DisambG>::vertex_iterator vIt, vItEnd;
  tie(vIt, vItEnd) = vertices(g);
  copy_if(vIt, vItEnd, V.begin(), vertex_is_connected<DisambG>(g));

  init_out_coefs(g, V, out_coefs);
  property_map<DisambG, vertex_rank_t>::type rank_map = get(vertex_rank, g);
  pageRank_iterate(g, V, N, out_coefs, rank_map, map_tmp, 30); // 30 iterations
}


////////////////////////////////////////////////////////////////
// Streaming
// Note: uses template functions in common.h

const size_t magic_id = 0x070517;

// read

Dis_vertex_t read_vertex_from_stream(ifstream & is, 
				     DisambG & g) {

  string name;
  string wname;
  float rank;

  read_atom_from_stream(is, name);
  read_atom_from_stream(is, wname);
  read_atom_from_stream(is, rank);
  Dis_vertex_t v = add_vertex(g);
  put(vertex_name, g, v, name);
  put(vertex_wname, g, v, wname);
  put(vertex_rank, g, v, rank);

  return v;
}

Dis_edge_t read_edge_from_stream(ifstream & is, 
				 DisambG & g) {

  size_t sIdx;
  size_t tIdx; 
  size_t freq;
  //size_t source;
  bool insertedP;
  Dis_edge_t e;

  read_atom_from_stream(is, sIdx);
  read_atom_from_stream(is, tIdx);
  read_atom_from_stream(is, freq);
  tie(e, insertedP) = add_edge(sIdx, tIdx, g);
  assert(insertedP);
  put(edge_freq, g, e, freq);

  return e;
}

void DisambGraph::read_from_stream (std::ifstream & is) {

  size_t vertex_n;
  size_t edge_n;
  size_t i;
  size_t id;

  read_atom_from_stream(is, id);
  if(id != magic_id) {
    cerr << "Error: invalid id (filename is a mcrGraph?)" << endl;
  }
  read_map_from_stream(is, mcr2dis);
  read_atom_from_stream(is, id);
  if(id != magic_id) {
    cerr << "Error: invalid id after reading maps" << endl;
  }

  read_atom_from_stream(is, vertex_n);
  for(i=0; i<vertex_n; ++i) {
    read_vertex_from_stream(is, g);
  }

  read_atom_from_stream(is, id);
  if(id != magic_id) {
    cerr << "Error: invalid id after reading vertices" << endl;
  }

  read_atom_from_stream(is, edge_n);
  for(i=0; i<edge_n; ++i) {
    read_edge_from_stream(is, g);
  }

  read_atom_from_stream(is, id);
  if(id != magic_id) {
    cerr << "Error: invalid id after reading edges" << endl;
  }

//   graph_traits<DisambG>::vertex_iterator v_it, v_end;
//   tie(v_it, v_end) = vertices(g);
//   for(; v_it != v_end; ++v_it) {
//     w2syns[get(vertex_name, g, *v_it)] = *v_it;
//   }
}

void DisambGraph::read_from_binfile (const string & fname) {
  ifstream fi(fname.c_str(), ifstream::binary|ifstream::in);
  if (!fi) {
    cerr << "Error: can't open " << fname << endl;
    exit(-1);
  }
  read_from_stream(fi);
}

// write

ofstream & write_vertex_to_stream(ofstream & o,
				  const DisambG & g,
				  const Dis_vertex_t & v) {
  string name;

  write_atom_to_stream(o, get(vertex_name, g, v));
  write_atom_to_stream(o, get(vertex_wname, g, v));
  write_atom_to_stream(o, get(vertex_rank, g, v));
  return o;
}

ofstream & write_edge_to_stream(ofstream & o,
				const DisambG & g,
				const Dis_edge_t & e) {

  size_t uIdx = get(vertex_index, g, source(e,g));
  size_t vIdx = get(vertex_index, g, target(e,g));
  size_t freq = get(edge_freq, g, e);

  o.write(reinterpret_cast<const char *>(&uIdx), sizeof(uIdx));
  o.write(reinterpret_cast<const char *>(&vIdx), sizeof(vIdx));
  o.write(reinterpret_cast<const char *>(&freq), sizeof(freq));

  return o;
}

ofstream & DisambGraph::write_to_stream(ofstream & o) const {

  // First write maps

  write_atom_to_stream(o, magic_id);
  write_map_to_stream(o, mcr2dis);
  write_atom_to_stream(o, magic_id);

  // Then the graph

  size_t vertex_n = num_vertices(g);

  write_atom_to_stream(o, vertex_n);
  graph_traits<DisambG>::vertex_iterator v_it, v_end;
  tie(v_it, v_end) = vertices(g);
  for(; v_it != v_end; ++v_it) {
    write_vertex_to_stream(o, g, *v_it);
  }

  write_atom_to_stream(o, magic_id);

  size_t edge_n = num_edges(g);

  write_atom_to_stream(o, edge_n);
  graph_traits<DisambG>::edge_iterator e_it, e_end;

  tie(e_it, e_end) = edges(g);
  for(; e_it != e_end; ++e_it) {
    write_edge_to_stream(o, g, *e_it);
  }
  return o;
}

void DisambGraph::write_to_binfile (const string & fName) const {

  ofstream fo(fName.c_str(),  ofstream::binary|ofstream::out);
  if (!fo) {
    cerr << "Error: can't create" << fName << endl;
    exit(-1);
  }
  write_to_stream(fo);
}

//////////////////////////////////////////////////////7
// graphviz

void write_dgraph_graphviz(const string & fname, const DisambG & g) {

  ofstream fo(fname.c_str(), ofstream::out);
  if (!fo) {
    cerr << "Can't create " << fname << endl;
    exit(-1);
  }
  boost::write_graphviz(fo,
                        g,
                        //boost::default_writer(),
			make_my_writer(get(vertex_wname, g), "label"),
                        make_my_writer(get(edge_freq, g), "weigth"));
}
