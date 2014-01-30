#include "disambGraph.h"
#include "common.h"
#include "kbGraph.h"
#include "dfsa.h"
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

#if BOOST_VERSION > 104400
#include <boost/range/irange.hpp>
#else
#include <boost/pending/integer_range.hpp>
#endif

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

	if (u == v)
	  throw runtime_error("Can't insert self loop !");

	map<Kb_vertex_t, Dis_vertex_t>::iterator it;

	tie(e, existP) = edge(u, v, g);
	if(!existP) {
	  // new edge
	  tie(e, existP) = add_edge(u, v, g);
	  put(edge_weight, g, e, w);
	} else {
	  // edge already there. Increase freq.
	  put(edge_weight, g, e, get(edge_weight, g, e) + w);
	}
  }

  Dis_vertex_t DisambGraph::add_dgraph_vertex(const string & str) {

	map<string, Dis_vertex_t>::iterator map_it;
	bool insertedP;
	tie(map_it, insertedP) = synsetMap.insert(make_pair(str, Dis_vertex_t()));
	if (insertedP) {
	  Dis_vertex_t v = add_vertex(g);
	  put(vertex_name, g, v, str);
	  map_it->second = v;
 	}
	return map_it->second;
  }

  void DisambGraph::fill_graph(Kb_vertex_t src,
							   Kb_vertex_t tgt,
							   const std::vector<Kb_vertex_t> & parents) {

	//   if (tgt == 81369) {
	//     int deb=0;
	//     deb++;
	//   }

	vector<string> path_str;
	Kb_vertex_t pred;
	Kb & kb = ukb::Kb::instance();

	pred = parents[tgt];
	while(tgt != pred) {
	  path_str.push_back(kb.get_vertex_name(tgt));
	  tgt = pred;
	  pred = parents[tgt];
	}
	if (tgt != src) return;
	path_str.push_back(kb.get_vertex_name(src));

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


  void DisambGraph::fill_graph(const set<Kb_edge_t> & E) {
	Kb & kb = Kb::instance();
	for(set<Kb_edge_t>::const_iterator it = E.begin(), end = E.end();
		it != end; ++it) {
	  Kb_vertex_t uu = kb.edge_source(*it);
	  Kb_vertex_t vv = kb.edge_target(*it);
	  if (uu == vv) continue;
	  Dis_vertex_t u = add_dgraph_vertex(kb.get_vertex_name(uu));
	  Dis_vertex_t v = add_dgraph_vertex(kb.get_vertex_name(vv));
	  add_dgraph_edge(u, v, 1.0);
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
  // Reset edge weights

  void DisambGraph::reset_edge_weights() {
	graph_traits<DisambG>::edge_iterator it, end;
	tie(it, end) = edges(g);
	for(;it != end; ++it)
	  put(edge_weight, g, *it, 1.0);
  }

  ////////////////////////////////////////////////////////////////
  // Global functions

  void fill_disamb_synset_bfs(const string & src_str,
							  vector<CWord>::const_iterator s_it,
							  vector<CWord>::const_iterator s_end,
							  DisambGraph & dgraph) {

	//bfs from src
	std::vector<Kb_vertex_t> parents;
	Kb & kb = ukb::Kb::instance();
	bool existP;
	Kb_vertex_t src, tgt;

	tie(src, existP) = kb.get_vertex_by_name(src_str);
	assert(existP);

	kb.bfs(src, parents);

	// insert src vertex in dgraph (fixes a bug)
	dgraph.add_dgraph_vertex(src_str);

	//fill disamb graph

	for(;s_it != s_end; ++s_it) {
	  vector<string>::const_iterator tg_it = s_it->begin();
	  vector<string>::const_iterator tg_end = s_it->end();
	  for(;tg_it != tg_end; ++tg_it) {
		tie(tgt, existP) = kb.get_vertex_by_name(*tg_it);
		assert(existP);
		dgraph.fill_graph(src, tgt, parents);
	  }
	}

  }

  void fill_disamb_synset_dijkstra(const string & src_str,
								   vector<CWord>::const_iterator s_it,
								   vector<CWord>::const_iterator s_end,
								   DisambGraph & dgraph) {

	std::vector<Kb_vertex_t> parents;
	Kb & kb = ukb::Kb::instance();
	bool existP;
	Kb_vertex_t src, tgt;

	tie(src, existP) = kb.get_vertex_by_name(src_str);
	assert(existP);

	kb.dijkstra(src, parents);

	// insert src vertex in dgraph (fixes a bug)
	dgraph.add_dgraph_vertex(src_str);

	//fill disamb graph

	for(;s_it != s_end; ++s_it) {
	  vector<string>::const_iterator tg_it = s_it->begin();
	  vector<string>::const_iterator tg_end = s_it->end();
	  for(;tg_it != tg_end; ++tg_it) {
		tie(tgt, existP) = kb.get_vertex_by_name(*tg_it);
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
		//tie(src_v, existP) = kb.get_vertex_by_name(*sset_it);
		//assert(existP);
		fill_disamb_synset_bfs(*sset_it, cw_it, cw_end, dgraph);
	  }
	}
  }

  // fill dgraph with ppv ranks
  // using edge weights in kb derived from ppv_rank

  void fill_disamb_graph(const CSentence & cs, DisambGraph & dgraph,
						 const vector<float> & ppv_ranks) {


	// First, update kb's edge weights
	ukb::Kb::instance().ppv_weights(ppv_ranks);

	vector<CWord>::const_iterator cw_it = cs.begin();
	vector<CWord>::const_iterator cw_end = cs.end();

	//if (cw_it == cw_end) return;
	//cw_end--;
	while(cw_it != cw_end) {
	  vector<string>::const_iterator sset_it = cw_it->begin();
	  vector<string>::const_iterator sset_end = cw_it->end();
	  ++cw_it; // point to next word
	  for(;sset_it != sset_end; ++sset_it) {
		//tie(src_v, existP) = kb.get_vertex_by_name(*sset_it);
		//assert(existP);
		fill_disamb_synset_dijkstra(*sset_it, cw_it, cw_end, dgraph);
	  }
	}
  }

  // dfs visitor (should be in kbgraph but ...)

  // implements as dfs with max depth. Here is the main idea (taken form a mail
  // in boost.users mailing list), due to David Abrahams:

  // It should be pretty easy to see how to a different visitor adaptor, sort of
  // like time_stamper, which increments a depth count when you discover a
  // vertex and decrements when you finish it.

  // Then you create a graph adaptor that refers to the depth count, and returns
  // an empty range for out_edges when the depth reaches your limit.  This
  // relies on knowing a lot about the structure of the DFS implementation, but
  // after all, the library documents it, so it should be OK to take advantage
  // of that.


  class dfsa_visitor : public default_dfs_visitor {
  public:
	dfsa_visitor(Kb_vertex_t s, const set<Kb_vertex_t> & S, const set<Kb_vertex_t> & C,
				 set<Kb_edge_t> & E)
	  : m_s(s), m_S(S), m_C(C), m_E(E), m_P(list<Kb_edge_t> ()), m_S_end(S.end()), m_C_end(C.end()) {}

	void discover_vertex(Kb_vertex_t u, const dfsa<KbGraph>& ag) {
	  if (u == m_s) return;
	  if(m_S.find(u) == m_S_end) return; // if u not in S, return
	  // Insert current path m_P into m_E
	  // Do nothing if there are cosenses in path
	  if(m_C.size()) {
		for(list<Kb_edge_t>::iterator it = m_P.begin(), end = m_P.begin();
			it != end; ++it) {
		  if(m_C.find(target(*it, ag)) != m_C_end) return; // if edge target is co-sense, do nothing
		}
	  }
	  m_E.insert(m_P.begin(), m_P.end());
	}

	void finish_vertex(Kb_vertex_t u, const dfsa<KbGraph>& ag) {
	  if (u == m_s) return;
	  const_cast<dfsa<KbGraph>&>(ag).dec_depth(); // decrease depth
	  m_P.pop_back();
	}

	void tree_edge(Kb_edge_t e, const dfsa<KbGraph> & ag) {
	  const_cast<dfsa<KbGraph>&>(ag).inc_depth(); // increase depth
	  m_P.push_back(e);
	}

  private:
	Kb_vertex_t m_s;              // source vertex
	const set<Kb_vertex_t> & m_S; // set of target synsets
	const set<Kb_vertex_t> & m_C; // the cosenses of source word
	set<Kb_edge_t> & m_E;         // the result set of edges
	list<Kb_edge_t> m_P;          // path of DFS so far
	set<Kb_vertex_t>::const_iterator m_S_end;
	set<Kb_vertex_t>::const_iterator m_C_end;
  };


  void fill_disamb_graph_dfs_nocosenses(const CSentence &cs, DisambGraph & dgraph) {
	set<Kb_vertex_t> S;
	KbGraph g = Kb::instance().graph();
	dfsa<KbGraph> ag(g, glVars::dGraph::max_depth);
	map<string, set<Kb_vertex_t> > coSenses;

	// Init S with all target synsets
	// Also, update coSense set for target words
	for(vector<CWord>::const_iterator cw_it = cs.begin(), cw_end = cs.end();
		cw_it != cw_end; ++cw_it) {
	  map<string, set<Kb_vertex_t> >::iterator coS_it;
	  if(cw_it->is_tgtword())
		coS_it = coSenses.insert(make_pair(cw_it->wpos(), set<Kb_vertex_t> ())).first;
	  for(vector<pair<Kb_vertex_t, float> >::const_iterator v_it = cw_it->V_vector().begin(),
			v_end = cw_it->V_vector().end();
		  v_it != v_end; ++v_it) {
		S.insert((*v_it).first);
		if(cw_it->is_tgtword()) (*coS_it).second.insert((*v_it).first);
	  }
	}

	std::vector<default_color_type> colors(num_vertices(g));
	typedef color_traits<default_color_type> Color;

	for(vector<CWord>::const_iterator cw_it = cs.begin(), cw_end = cs.end();
		cw_it != cw_end; ++cw_it) {
	  if (!cw_it->is_tgtword()) continue;
	  set<Kb_vertex_t> & coS = coSenses[cw_it->wpos()];
	  for(vector<pair<Kb_vertex_t, float> >::const_iterator wit = cw_it->V_vector().begin(), wend = cw_it->V_vector().end();
		  wit != wend; ++wit) {
		set<Kb_edge_t> subg;
		Kb_vertex_t src = (*wit).first;
		fill(colors.begin(), colors.end(), Color::white());
		dfsa_visitor vis(src, S, coS, subg);
		depth_first_visit(ag, src, vis, &colors[0]);
		// Now  populate disambGraph with edges in subg
		dgraph.fill_graph(subg);
	  }
	}
  }

  void fill_disamb_graph_dfs(const CSentence &cs, DisambGraph & dgraph) {
	set<Kb_vertex_t> S;
	set<Kb_vertex_t> TW_S;
	Kb & kb = Kb::instance();
	dfsa<KbGraph> ag(kb.graph(), glVars::dGraph::max_depth);

	// Init S with all target synsets
	for(vector<CWord>::const_iterator cw_it = cs.begin(), cw_end = cs.end();
		cw_it != cw_end; ++cw_it) {
	  for(vector<pair<Kb_vertex_t, float> >::const_iterator v_it = cw_it->V_vector().begin(),
			v_end = cw_it->V_vector().end();
		  v_it != v_end; ++v_it) {
		S.insert((*v_it).first);
		if(cw_it->is_tgtword()) TW_S.insert((*v_it).first);
	  }
	}

	std::vector<default_color_type> colors(kb.size());
	typedef color_traits<default_color_type> Color;

	for(set<Kb_vertex_t>::iterator it = TW_S.begin(), end = TW_S.end(); it != end; ++it) {
	  fill(colors.begin(), colors.end(), Color::white());
	  set<Kb_edge_t> subg;
	  dfsa_visitor vis(*it, S, set<Kb_vertex_t>(), subg);
	  depth_first_visit(ag, *it, vis, &colors[0]);
	  // Now  populate disambGraph with edges in subg
	  dgraph.fill_graph(subg);
	}
  }

  // Convert a pv vector of Kb_vertex_t to the equivalent for Dis_vertex_t

  size_t pv_to_dgraph(DisambGraph & dgraph,
					  const vector<float> & pv,
					  vector<float> & pv_dgraph) {

	Dis_vertex_t u;
	bool P;
	size_t k = 0;
	Kb & kb = ukb::Kb::instance();
	for(size_t i = 0, end = pv.size(); i != end; ++i) {
	  if (pv[i] == 0.0) continue;
	  tie(u, P) = dgraph.get_vertex_by_name(kb.get_vertex_name(i));
	  if (!P) continue;
	  ++k;
	  pv_dgraph[u] = pv[i];
	}
	return k;
  }

  bool csentence_dgraph_ppr(const CSentence & cs, DisambGraph & dgraph,
							vector<float> & ranks,
							CSentence::const_iterator exclude_word_it) {

	// get pv pointing to KbGraph vertex_t
	// transform into Dis_vertex_t

	if (!cs.has_tgtwords()) return false; // no target words

	vector<float> pv;
	size_t  pv_m = pv_from_cs_onlyC(cs, pv, exclude_word_it);
	if (!pv_m) return false;

	// create pv_dgraph (map Kb_vertex_t to Dis_vertex_t)
	vector<float> pv_dgraph(dgraph.size(), 0.0);
	size_t pv_dgraph_m = pv_to_dgraph(dgraph, pv, pv_dgraph);
	if (!pv_dgraph_m) return false;

	if (pv_m != pv_dgraph_m) {
	  // Some synsets from the CSentence are not in the dgraph
	  // Renormalize pv_dgraph
	  normalize_pvector(pv_dgraph);
	}

	// Execute PageRank using pv_dgraph
	dgraph.pageRank_ppv(pv_dgraph, ranks);
	return true;
  }

  bool csentence_dgraph_ppr(const CSentence & cs, DisambGraph & dgraph,
							vector<float> & ranks) {
	return csentence_dgraph_ppr(cs, dgraph, ranks, cs.end());
  }


  bool dgraph_static(DisambGraph & dgraph,
					 vector<float> & ranks) {
	size_t N = dgraph.size();
	if (!N) return false;
	float factor = 1 / N;
	vector<float> pv_dgraph(dgraph.size(), factor);
	dgraph.pageRank_ppv(pv_dgraph, ranks);
	return true;
 }


  bool dgraph_degree(DisambGraph & dgraph, vector<float> & ranks) {

	size_t N = dgraph.size();

	if (N == ranks.size()) {
	  std::fill(ranks.begin(), ranks.end(), 0.0);
	} else {
	  vector<float>(N, 0.0).swap(ranks); // Initialize rank vector
	}

	prank::init_degree(dgraph.graph(), &ranks[0]);
	return true;
  }


  bool disamb_csentence_dgraph(CSentence & cs, DisambGraph & dgraph,
							   const vector<float> & ranks) {

	if (!cs.has_tgtwords()) return false; // no target words

	vector<CWord>::iterator cw_it = cs.begin();
	vector<CWord>::iterator cw_end = cs.end();
	for(; cw_it != cw_end; ++cw_it) {
	  if(!cw_it->is_tgtword()) continue;
	  disamb_cword_dgraph(cw_it, dgraph, ranks);
	}
	return true;
  }

  void disamb_cword_dgraph(CSentence::iterator it, DisambGraph & dgraph,
						   const std::vector<float> & ranks) {
	it->rank_synsets(dgraph, ranks);
	it->disamb_cword();
  }

  ostream & print_complete_csent(ostream & o, CSentence & cs, DisambGraph & dgraph) {

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
		  o << *syn_it << " ,";
		}
		tie(v,P) = dgraph.get_vertex_by_name(*syn_end);
		assert(P);
		o << *syn_end;
	  }
	  o << "}" << endl;
	}
	return o;
  }


  /////////////////////////////////////////////////////////////////////
  // pageRank
  //

  void DisambGraph::pageRank_ppv(const vector<float> & ppv_map,
								 vector<float> & ranks) {

	typedef graph_traits<DisambG>::edge_descriptor edge_descriptor;
	property_map<DisambGraph::boost_graph_t, edge_weight_t>::type weight_map = get(edge_weight, g);
	prank::constant_property_map <edge_descriptor, float> cte_weight(1.0); // always return 1

	size_t N = num_vertices(g);
	size_t N_no_isolated;
	vector<float> out_coefs(N, 0.0);

	if (N == ranks.size()) {
	  std::fill(ranks.begin(), ranks.end(), 0.0);
	} else {
	  vector<float>(N, 0.0).swap(ranks); // Initialize rank vector
	}
	vector<float> rank_tmp(N, 0.0);    // auxiliary rank vector

	if (glVars::prank::use_weight) {
	  N_no_isolated = prank::init_out_coefs(g,  &out_coefs[0], weight_map);
	} else {
	  N_no_isolated = prank::init_out_coefs(g, &out_coefs[0], cte_weight);
	}

	if (glVars::prank::use_weight) {
	  prank::do_pageRank(g, N_no_isolated, &ppv_map[0],
						 weight_map, &ranks[0], &rank_tmp[0],
						 glVars::prank::num_iterations,
						 glVars::prank::threshold,
						 glVars::prank::damping,
						 out_coefs);
	} else {
	  prank::do_pageRank(g, N_no_isolated, &ppv_map[0],
						 cte_weight, &ranks[0], &rank_tmp[0],
						 glVars::prank::num_iterations,
						 glVars::prank::threshold,
						 glVars::prank::damping,
						 out_coefs);
	}

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
		auth+= get(hProp, source(*e, g)) * get(edge_weight, g, *e);
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
		hub+= get(aProp, target(*e, g)) * get(edge_weight, g, *e);
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
	float init_v = sqrt((float)n)/n;

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

  void hits(DisambG & g, vector<float> & rank) {

	vector<float> aRank(num_vertices(g), 0.0f);
	vector<float> hRank(num_vertices(g), 0.0f);

	vector<Dis_vertex_t> V(num_connected_vertices(g));

	graph_traits<DisambG>::vertex_iterator vIt, vItEnd;
	tie(vIt, vItEnd) = vertices(g);
	ukb::copy_if(vIt, vItEnd, V.begin(), vertex_is_connected<DisambG>(g));

	hits_iterate(g, V.begin(), V.end(), hRank, aRank, 50); // 50 iterations

	//vector<Dis_vertex_t>::const_iterator it = V.begin();
	//vector<Dis_vertex_t>::const_iterator end = V.end();
	rank.swap(hRank);
  }

  ////////////////////////////////////////////////////////////////
  // Streaming
  // Note: uses template functions in common.h

  const size_t magic_id = 0x070517;

  // read

  Dis_vertex_t read_vertex_from_stream(ifstream & is,
									   DisambG & g) {

	string name;

	read_atom_from_stream(is, name);
	Dis_vertex_t v = add_vertex(g);
	put(vertex_name, g, v, name);
	//  put(vertex_wname, g, v, wname);
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
	put(edge_weight, g, e, freq);

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
	return o;
  }

  ofstream & write_edge_to_stream(ofstream & o,
								  const DisambG & g,
								  const Dis_edge_t & e) {

	size_t uIdx = get(vertex_index, g, source(e,g));
	size_t vIdx = get(vertex_index, g, target(e,g));
	float freq = get(edge_weight, g, e);

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
	  out << " [weight=\"" << get(edge_weight, g, e) << "\"]";
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
						  make_my_writer(get(vertex_name, g), "label"),
						  make_my_writer(get(edge_weight, g), "weight"));
  }
}
