#include "disambGraph.h"
#include "common.h"
#include "mcrGraph.h"
#include "csentence.h"
#include "prank.h"
#include "globalVars.h"

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


// graphviz

#include <boost/graph/graphviz.hpp>
#include <boost/graph/filtered_graph.hpp>

// bfs 

#include <boost/graph/visitors.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/pending/indirect_cmp.hpp>
#include <boost/pending/integer_range.hpp>
#include <boost/graph/graph_utility.hpp> // for boost::make_list

namespace ukb {

  using namespace std;
  using namespace boost;


  ////////////////////////////////////////////////////////////////////////////////
  // Disamb fill functions


  DisambGraph::DisambGraph() {
  }

  void DisambGraph::add_dgraph_edge(Dis_vertex_t u, Dis_vertex_t v, float w) {

	Dis_edge_t e;
	bool existP;

	map<Mcr_vertex_t, Dis_vertex_t>::iterator it;

	tie(e, existP) = edge(u, v, g);
	if(!existP) {
	  // new edge
	  tie(e, existP) = add_edge(u, v, g);
	  put(edge_freq, g, e, w);
	} else {
	  // edge already there. Increase freq.
	  put(edge_freq, g, e, get(edge_freq, g, e) + w);
	}
  }

  Dis_vertex_t DisambGraph::add_dgraph_vertex(const string & str) {

	map<string, Dis_vertex_t>::iterator map_it;
	bool insertedP;
	tie(map_it, insertedP) = synsetMap.insert(make_pair(str, Dis_vertex_t()));
	if (insertedP) {
	  Dis_vertex_t v = add_vertex(g);
	  put(vertex_name, g, v, str);
	  put(vertex_freq, g, v, 0.0);
	  put(vertex_mcrSource, g, v, ukb::Mcr::instance().get_vertex_by_name(str).first);
	  map_it->second = v;
	} else {
	  put(vertex_freq, g, map_it->second, 
		  get(vertex_freq, g, map_it->second) + 1.0);
	}
	return map_it->second;
  }

  void DisambGraph::fill_graph(Mcr_vertex_t src,
							   Mcr_vertex_t tgt,
							   const std::vector<Mcr_vertex_t> & parents) {

	//   if (tgt == 81369) {
	//     int deb=0;
	//     deb++;
	//   }

	vector<string> path_str;
	Mcr_vertex_t pred;
	McrGraph & mcr_g = ukb::Mcr::instance().graph();

	pred = parents[tgt];
	while(tgt != pred) {
	  path_str.push_back(get(vertex_name, mcr_g, tgt));
	  tgt = pred;
	  pred = parents[tgt];
	}
	if (tgt != src) return;
	path_str.push_back(get(vertex_name, mcr_g, src));
  
	vector<Dis_vertex_t> path_v;
	size_t length = 0;
	for(vector<string>::iterator v_it = path_str.begin(); 
		v_it != path_str.end(); 
		++v_it) {
      Dis_vertex_t u = add_dgraph_vertex(*v_it);
      path_v.push_back(u);
      ++length;
	}

	vector<Dis_vertex_t>::iterator path_it = path_v.begin();
	vector<Dis_vertex_t>::iterator path_end = path_v.end();
	vector<Dis_vertex_t>::iterator path_prev;

	path_prev = path_it;
	++path_it;
	while(path_it != path_end) {
	  add_dgraph_edge(*path_it, *path_prev, 1.0 / static_cast<float>(length));
	  path_prev = path_it;
	  ++path_it;
	}
  }

  ////////////////////////////////////////////////////////////////////////////////
  // vertex_id <-> strings 


  pair<Dis_vertex_t, bool> DisambGraph::get_vertex_by_name(const std::string & str) const {
	map<string, Dis_vertex_t>::const_iterator it = synsetMap.find(str);
	if (it == synsetMap.end()) return make_pair(Dis_vertex_t(), false);
	return make_pair(it->second, true);
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Reset edge weigths

  void DisambGraph::reset_edge_weigths() {
	graph_traits<DisambG>::edge_iterator it, end;
	tie(it, end) = edges(g);
	for(;it != end; ++it)
	  put(edge_freq, g, *it, 1.0);
  }

  ////////////////////////////////////////////////////////////////
  // Global functions

  void fill_disamb_synset_bfs(const string & src_str,
							  vector<CWord>::const_iterator s_it,
							  vector<CWord>::const_iterator s_end,
							  DisambGraph & dgraph) {

	//   if (src_str == "00663525-a") {
	//     int deb;
	//     cerr << "Eooo!" << endl;
	//     deb++;
	//   }

	//bfs from src
	std::vector<Mcr_vertex_t> parents;
	Mcr & mcr = ukb::Mcr::instance();
	bool existP;
	Mcr_vertex_t src, tgt;
  
	tie(src, existP) = mcr.get_vertex_by_name(src_str);
	assert(existP);

	mcr.bfs(src, parents);

	// insert src vertex in dgraph (fixes a bug)
	dgraph.add_dgraph_vertex(src_str);
  
	//fill disamb graph

	for(;s_it != s_end; ++s_it) {
	  vector<string>::const_iterator tg_it = s_it->begin();
	  vector<string>::const_iterator tg_end = s_it->end();
	  for(;tg_it != tg_end; ++tg_it) {
		tie(tgt, existP) = mcr.get_vertex_by_name(*tg_it);
		assert(existP);
		dgraph.fill_graph(src, tgt, parents);
	  }
	}

  }

  void fill_disamb_synset_dijkstra(const string & src_str,
								   vector<CWord>::const_iterator s_it,
								   vector<CWord>::const_iterator s_end,
								   DisambGraph & dgraph) {

	std::vector<Mcr_vertex_t> parents;
	Mcr & mcr = ukb::Mcr::instance();
	bool existP;
	Mcr_vertex_t src, tgt;
  
	tie(src, existP) = mcr.get_vertex_by_name(src_str);
	assert(existP);

	mcr.dijkstra(src, parents);

	// insert src vertex in dgraph (fixes a bug)
	dgraph.add_dgraph_vertex(src_str);
  
	//fill disamb graph

	for(;s_it != s_end; ++s_it) {
	  vector<string>::const_iterator tg_it = s_it->begin();
	  vector<string>::const_iterator tg_end = s_it->end();
	  for(;tg_it != tg_end; ++tg_it) {
		tie(tgt, existP) = mcr.get_vertex_by_name(*tg_it);
		assert(existP);
		dgraph.fill_graph(src, tgt, parents);
	  }
	}

  }


  void fill_disamb_graph(const CSentence &cs, DisambGraph & dgraph) {

	vector<CWord>::const_iterator cw_it = cs.begin();
	vector<CWord>::const_iterator cw_end = cs.end();

	//if (cw_it == cw_end) return;
	//cw_end--;
	while(cw_it != cw_end) {
	  vector<string>::const_iterator sset_it = cw_it->begin();
	  vector<string>::const_iterator sset_end = cw_it->end();
	  ++cw_it; // point to next word
	  for(;sset_it != sset_end; ++sset_it) {
		//tie(src_v, existP) = mcr.get_vertex_by_name(*sset_it);
		//assert(existP);
		fill_disamb_synset_bfs(*sset_it, cw_it, cw_end, dgraph);
	  }
	}
  }

  // fill dgraph with ppv ranks
  // using edge weights in mcr derived from ppv_rank

  void fill_disamb_graph(const CSentence & cs, DisambGraph & dgraph,
						 const vector<float> & ppv_ranks) {

  
	// First, update mcr's edge weights
	ukb::Mcr::instance().ppv_weights(ppv_ranks);  

	vector<CWord>::const_iterator cw_it = cs.begin();
	vector<CWord>::const_iterator cw_end = cs.end();

	//if (cw_it == cw_end) return;
	//cw_end--;
	while(cw_it != cw_end) {
	  vector<string>::const_iterator sset_it = cw_it->begin();
	  vector<string>::const_iterator sset_end = cw_it->end();
	  ++cw_it; // point to next word
	  for(;sset_it != sset_end; ++sset_it) {
		//tie(src_v, existP) = mcr.get_vertex_by_name(*sset_it);
		//assert(existP);
		fill_disamb_synset_dijkstra(*sset_it, cw_it, cw_end, dgraph);
	  }
	}
  }

  // Warning: Rank scores of dgraph must be previously computed.
  //

  void disamb_csentence(CSentence & cs, DisambGraph & dgraph) {

	vector<CWord>::iterator cw_it = cs.begin();
	vector<CWord>::iterator cw_end = cs.end();
	for(; cw_it != cw_end; ++cw_it) {
	  cw_it->rank_synsets(dgraph, get(vertex_rank, dgraph.graph()));
	  cw_it->disamb_cword();
	}
  }

  ostream & print_complete_csent(ostream & o, CSentence & cs, DisambGraph & dgraph) {
	DisambG & g = dgraph.graph();
	vector<CWord>::iterator cw_it = cs.begin();
	vector<CWord>::iterator cw_end = cs.end();
	Dis_vertex_t v;
	bool P;

	for(; cw_it != cw_end; ++cw_it) {
	  o << cw_it->word() << "{";
	  vector<string>::iterator syn_it = cw_it->begin();
	  vector<string>::iterator syn_end = cw_it->end();
	  if (syn_it != syn_end) {
		--syn_end;
		for(; syn_it != syn_end; ++syn_it) {
		  tie(v,P) = dgraph.get_vertex_by_name(*syn_it);
		  assert(P);
		  o << *syn_it << ":" << get(vertex_rank, g, v) << " ,";
		}
		tie(v,P) = dgraph.get_vertex_by_name(*syn_end);
		assert(P);
		o << *syn_end << ":" << get(vertex_rank, g, v);
	  }
	  o << "}" << endl;
	}
	return o;
  }


  ////////////////////////////////////////////////////////////////////////////////
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
	float init_v = sqrt(n) / static_cast<float>(n);

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

  void pageRank_disg(DisambG & g,
					 bool use_weigth) {

	vector<float> map_tmp(num_vertices(g), 0.0f);

	size_t N = num_connected_vertices(g);

	vector<Dis_vertex_t> V(N);

	// property maps
	prank::constant_property_map<Dis_vertex_t, float> ppv(1.0/static_cast<float>(N)); // always return 1/N
	property_map<DisambG, vertex_rank_t>::type rank_map = get(vertex_rank, g);

	graph_traits<DisambG>::vertex_iterator vIt, vItEnd;
	tie(vIt, vItEnd) = vertices(g);
	copy_if(vIt, vItEnd, V.begin(), vertex_is_connected<DisambG>(g));
	//cerr << num_vertices(g) << endl;
	//cerr << num_connected_vertices(g) << endl;
  
	if (use_weigth) {
	  property_map<DisambG, edge_freq_t>::type weight_map = get(edge_freq, g);
	  //init_out_coefs(g, V, &out_coefs[0], weight_map);
	  prank::pageRank_iterate(g, V, ppv, 
							  weight_map, rank_map, &map_tmp[0], glVars::prank::num_iterations);
	} else {
	  prank::pageRank_iterate_now(g, V, ppv, 
								  rank_map, &map_tmp[0], 
								  glVars::prank::num_iterations); 

	}
  }

  void pageRank_ppv_disg(DisambG &g,
						 const map<string, size_t> & syn_n,
						 bool use_weigth) {

	// Fill rank freqs

	map<string, size_t>::const_iterator syn_n_end = syn_n.end();

	size_t total_count = 0;
	graph_traits<DisambG>::vertex_iterator u, end;
	tie(u, end) = vertices(g);
	for(; u != end; ++u) {
	  map<string, size_t>::const_iterator syn_n_it = syn_n.find(get(vertex_name, g, *u));
	  if (syn_n_it != syn_n_end) {
		put(vertex_freq, g, *u, syn_n_it->second);
		total_count += syn_n_it->second;
	  } else {
		cerr << "W: " << syn_n_it->first << " synset not found in dgraph!" << endl;
		put(vertex_freq, g, *u, 0.0);
	  }
	} 

	assert(total_count);
	double factor = double(1.0) / static_cast<double>(total_count);

	// Make freqs a prob. dist
	tie(u, end) = vertices(g);
	for(; u != end; ++u) {
	  put(vertex_freq, g, *u, 
		  get(vertex_freq, g, *u) * factor);
	}

	// usual pagerank

	vector<float> map_tmp(num_vertices(g), 0.0f);

	size_t N = num_connected_vertices(g);

	vector<Dis_vertex_t> V(N);

	tie(u, end) = vertices(g);
	copy_if(u, end, V.begin(), vertex_is_connected<DisambG>(g));

	// property maps

	property_map<DisambG, vertex_freq_t>::type ppv_map = get(vertex_freq, g);
	property_map<DisambG, vertex_rank_t>::type rank_map = get(vertex_rank, g);

	if (use_weigth) {
	  property_map<DisambG, edge_freq_t>::type weight_map = get(edge_freq, g);
	  //init_out_coefs(g, V, &out_coefs[0], weight_map);
	  prank::pageRank_iterate(g, V, ppv_map, weight_map, 
							  rank_map, map_tmp, glVars::prank::num_iterations);
	} else {
	  prank::pageRank_iterate_now(g, V, ppv_map, 
								  rank_map, map_tmp, glVars::prank::num_iterations);
	}
  }




  void degreeRank(DisambG & g) {
	prank::constant_property_map <Dis_edge_t, float> cte_weight(1); // always return 1  
	property_map<DisambG, vertex_rank_t>::type rank_map = get(vertex_rank, g);

	init_degree(g, rank_map, cte_weight);
  }

  ////////////////////////////////////////////////////////////////
  // Streaming
  // Note: uses template functions in common.h

  const size_t magic_id = 0x070517;

  // read

  Dis_vertex_t read_vertex_from_stream(ifstream & is, 
									   DisambG & g) {

	string name;
	Mcr_vertex_t mcr_source;
	float rank;

	read_atom_from_stream(is, name);
	read_atom_from_stream(is, mcr_source);
	read_atom_from_stream(is, rank);
	Dis_vertex_t v = add_vertex(g);
	put(vertex_name, g, v, name);
	//  put(vertex_wname, g, v, wname);
	put(vertex_rank, g, v, rank);
	put(vertex_mcrSource, g, v, mcr_source);
	return v;
  }

  Dis_edge_t read_edge_from_stream(ifstream & is, 
								   DisambG & g) {

	size_t sIdx;
	size_t tIdx; 
	float freq;
	//size_t source;
	bool insertedP;
	Dis_edge_t e;

	read_atom_from_stream(is, tIdx);
	read_atom_from_stream(is, sIdx);
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
	  cerr << "Error: invalid id (filename is a disambGraph?)" << endl;
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
	write_atom_to_stream(o, get(vertex_mcrSource, g, v));
	write_atom_to_stream(o, get(vertex_rank, g, v));
	return o;
  }

  ofstream & write_edge_to_stream(ofstream & o,
								  const DisambG & g,
								  const Dis_edge_t & e) {

	size_t uIdx = get(vertex_index, g, source(e,g));
	size_t vIdx = get(vertex_index, g, target(e,g));
	float freq = get(edge_freq, g, e);

	o.write(reinterpret_cast<const char *>(&vIdx), sizeof(vIdx));
	o.write(reinterpret_cast<const char *>(&uIdx), sizeof(uIdx));
	o.write(reinterpret_cast<const char *>(&freq), sizeof(freq));

	return o;
  }

  ofstream & DisambGraph::write_to_stream(ofstream & o) const {

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


  class myV_writer {
  public:
	myV_writer(const DisambG & g_) : g(g_) {};
	void operator()(std::ostream& out, const Dis_vertex_t & v) const {
	  out << " [label=\"" << get(vertex_name, g, v) << "\"]";
	}
	const DisambG & g;
  };

  class myE_writer {
  public:
	myE_writer(const DisambG & g_) : g(g_) {};
	void operator()(std::ostream& out, const Dis_edge_t & e) const {
	  out << " [weigth=\"" << get(edge_freq, g, e) << "\"]";
	}
	const DisambG & g;
  };


  //
  // Predicate to make a filtered graph
  //

  template<class G>
  struct connected_vertex {
	G * g;
	connected_vertex( G & g_) : g(&g_) {};
	connected_vertex( ) {};
	typedef typename boost::graph_traits<G>::vertex_descriptor vertex_descriptor;
	bool operator ()(vertex_descriptor u) const {
	  return out_degree(u, *g) || in_degree(u, *g);
	}
  };


  void write_dgraph_graphviz(const string & fname, DisambG & g) {

	ofstream fo(fname.c_str(), ofstream::out);
	if (!fo) {
	  cerr << "Can't create " << fname << endl;
	  exit(-1);
	}
	connected_vertex<DisambG> vp(g);
	//   boost::filtered_graph<DisambG, keep_all, keep_all >
	//     fg(g, keep_all(), keep_all());

	boost::filtered_graph<DisambG, keep_all, connected_vertex<DisambG> >
	  fg(g, keep_all(), vp);

	boost::write_graphviz(fo,
						  fg,
						  //boost::default_writer(),
						  //make_my_writer(get(vertex_name, g), "label"),
						  make_my_writer3(get(vertex_name, g), 
										  get(vertex_rank, g), 
										  get(vertex_mcrSource, g),
										  "label", "rank", "mcr"),
						  make_my_writer(get(edge_freq, g), "weigth"));
  }
}
