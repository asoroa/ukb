#include "kGraph.h"
#include "common.h"

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



Dis_vertex_t add_kg_vertex(const string & str, 
			   DisambG & g, 
			   map<string, Dis_vertex_t> & synsetMap) {

  map<string, Dis_vertex_t>::iterator map_it;
  map<string, Dis_vertex_t>::iterator map_end = synsetMap.end();
  bool insertedP;

  tie(map_it, insertedP) = synsetMap.insert(make_pair(str, Dis_vertex_t()));  
  if (insertedP) {
    Dis_vertex_t v = add_vertex(g);
    put(vertex_name, g, v, str);
    put(vertex_rank, g, v, 0.0f);
    map_it->second = v;
  }
  return map_it->second;
}



// Auxiliary struct for constructing KGraphs from DisambGraphs

struct Kg_dg_ {
  Dis_vertex_t dg_v; // vertex id in DisambGraph
  Dis_vertex_t kg_v; // vertex id in KGraph
};

typedef std::vector<std::vector<Kg_dg_> > KG_SynsV;

// given a csentence, do two things:
//
// * create a KG_Synsv structure, where the synsets of every word are
//   stored. 
//
// * add the vertices of those synsets in the newly created Kgraph
//

void create_kgSynsV(CSentence & cs,
		    KG_SynsV & theSyns,
		    DisambGraph & dg,
		    DisambG & g,
		    std::map<std::string, Dis_vertex_t> & synsetMap) {

  Dis_vertex_t aux;
  bool existP;

  vector<CWord>::const_iterator cw_it = cs.begin();
  vector<CWord>::const_iterator cw_end = cs.end();
  
  for(;cw_it != cw_end; ++cw_it) {
    theSyns.push_back(vector<Kg_dg_>());
    vector<string>::const_iterator sset_it = cw_it->begin();
    vector<string>::const_iterator sset_end = cw_it->end();
    vector<Kg_dg_> cw_theSyns;
    for(; sset_it != sset_end; ++sset_it) {
      tie(aux, existP) = dg.getVertexByName(*sset_it);
      assert(existP);
      Kg_dg_ new_elem;
      new_elem.dg_v = aux;
      new_elem.kg_v = add_kg_vertex(*sset_it, g, synsetMap);
      cw_theSyns.push_back(new_elem);
    }
    theSyns.back().swap(cw_theSyns);
  }
}


void fill_kg(Kg_dg_ & src, 
	     vector<vector<Kg_dg_> >::iterator cw_it,
	     vector<vector<Kg_dg_> >::iterator cw_end,
	     DisambG & disg,
	     DisambG & g) {
  
  vector<Dis_vertex_t> dist;
  Dis_edge_t e;
  bool P;

  disg_bfs(disg, src.dg_v, dist);
  for(; cw_it != cw_end; ++cw_it) {
    vector<Kg_dg_>::iterator syn_it, syn_end;
    syn_it  = cw_it->begin();
    syn_end = cw_it->end();
    for(;syn_it != syn_end; ++syn_it) {
      tie (e, P) = edge(src.kg_v, syn_it->kg_v, g);
      //cerr << get(vertex_name, disg, src.dg_v) << " " <<
      //	get(vertex_name, disg, syn_it->dg_v) << endl;
//       if (P) {
// 	cerr << get(vertex_name, disg, src.dg_v) << " " <<
// 	  get(vertex_name, disg, syn_it->dg_v) << " repeated\n";
//       }
      tie (e,P) = add_edge(src.kg_v, syn_it->kg_v, g);
      put(edge_freq, g, e, dist[syn_it->dg_v]);
    }
  }
}

KGraph::KGraph(CSentence & cs, DisambGraph & disg) {

  KG_SynsV theSyns;

  if (cs.size() == 0) return;
  
  create_kgSynsV(cs, theSyns, disg, g, synsetMap); // This function also adds the vertices

  // Add edges
  vector<vector<Kg_dg_> >::iterator cw_it, cw_end;
  cw_it =  theSyns.begin();
  cw_end = theSyns.end();
  while(cw_it != cw_end) {
    vector<Kg_dg_>::iterator syn_it, syn_end;
    syn_it = cw_it->begin();
    syn_end = cw_it->end();
    ++cw_it; // point to next word
    for(; syn_it != syn_end; ++syn_it) {
      fill_kg(*syn_it, cw_it, cw_end, disg.graph(), g);
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
// vertex_id <-> strings 


pair<Dis_vertex_t, bool> KGraph::getVertexByName(const std::string & str) const {
  map<string, Dis_vertex_t>::const_iterator it = synsetMap.find(str);
  if (it == synsetMap.end()) return make_pair(Dis_vertex_t(), false);
  return make_pair(it->second, true);
}

void disamb_csentence(CSentence & cs, KGraph & kgraph) {

  cerr << "kgraph.cc:disamb_csentence Not implemented!!" <<endl;
  exit(-1);

//   vector<CWord>::iterator cw_it = cs.begin();
//   vector<CWord>::iterator cw_end = cs.end();
//   for(; cw_it != cw_end; ++cw_it) {
//     Syn2vert_tie<KGraph> scs(cw_it->begin(), cw_it->end(), kgraph);
//     sort(scs.V.begin(), scs.V.end(), 
// 	 make_SortByRank(kgraph.graph(),
// 			 get(vertex_rank, kgraph.graph())));
//     vector<string> new_v;
//     vector<pair<string *, Dis_vertex_t> >::iterator it, end;
//     it  = scs.V.begin();
//     end = scs.V.end();
//     for(;it != end; ++it) {
//       new_v.push_back(*(it->first));
//     }
//     cw_it->get_syns_vector().swap(new_v);
//   }  
}

////////////////////////////////////////////////////////////////
// Streaming
// Note: uses template functions in common.h
//       uses read vertices/edges functions from disambGraph.cc

const size_t magic_id = 0x070706;

void KGraph::read_from_stream (std::ifstream & is) {

  size_t vertex_n;
  size_t edge_n;
  size_t i;
  size_t id;

  read_atom_from_stream(is, id);
  if(id != magic_id) {
    cerr << "Error: invalid id (filename is a kGraph?)" << endl;
  }
  read_map_from_stream(is, synsetMap);
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

void KGraph::read_from_binfile (const string & fname) {
  ifstream fi(fname.c_str(), ifstream::binary|ifstream::in);
  if (!fi) {
    cerr << "Error: can't open " << fname << endl;
    exit(-1);
  }
  read_from_stream(fi);
}

ofstream & KGraph::write_to_stream(ofstream & o) const {

  // First write maps

  write_atom_to_stream(o, magic_id);
  write_map_to_stream(o, synsetMap);
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

void KGraph::write_to_binfile (const string & fName) const {

  ofstream fo(fName.c_str(),  ofstream::binary|ofstream::out);
  if (!fo) {
    cerr << "Error: can't create" << fName << endl;
    exit(-1);
  }
  write_to_stream(fo);
}

//////////////////////////////////////////////////////7
