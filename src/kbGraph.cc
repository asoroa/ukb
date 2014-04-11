#include "kbGraph.h"
#include "common.h"
#include "globalVars.h"
#include "wdict.h"
#include "prank.h"

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <string>
#include <map>
#include <iterator>
#include <algorithm>
#include <ostream>

// Tokenizer
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

// Stuff for generating random numbers

#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>

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

// dijkstra

#include <boost/graph/dijkstra_shortest_paths.hpp>

// strong components

#include <boost/graph/strong_components.hpp>


namespace ukb {

  using namespace std;
  using namespace boost;

  ////////////////////////////////////////////////////////////////////////////////
  // Class Kb


  ////////////////////////////////////////////////////////////////////////////////
  // Singleton stuff

  Kb* Kb::p_instance = 0;

  Kb *Kb::create() {

	static Kb theKb;
	return &theKb;
  }

  Kb & Kb::instance() {
	if (!p_instance) {
	  throw runtime_error("KB not initialized");
	}
	return *p_instance;
  }

  void Kb::create_from_txt(const string & synsFileName,
						   const std::set<std::string> & src_allowed) {
	if (p_instance) return;
	Kb *tenp = create();
	tenp->read_from_txt(synsFileName, src_allowed);
	p_instance = tenp;
  }

  void Kb::create_from_txt(std::istream & is,
						   const std::set<std::string> & src_allowed) {
	if (p_instance) return;
	Kb *tenp = create();
	tenp->read_from_txt(is, src_allowed);
	p_instance = tenp;
  }

  void Kb::create_from_binfile(const std::string & fname) {

	if (p_instance) return;
	Kb *tenp = create();

	ifstream fi(fname.c_str(), ifstream::binary|ifstream::in);
	if (!fi) {
	  cerr << "Error: can't open " << fname << endl;
	  exit(-1);
	}
	try {
	  tenp->read_from_stream(fi);
	} catch(std::exception& e) {
	  cerr << e.what() << "\n";
	  exit(-1);
	}
	p_instance = tenp;
  }


  void Kb::create_from_kbgraph16(Kb16 & kbg) {
	if (p_instance) return;
	Kb *tenp = create();
	precsr_t precsr16;

	Kb16::boost_graph_t oldg = kbg.g;
	graph_traits<Kb16::boost_graph_t>::edge_iterator eit, eend;
	tie(eit, eend) = edges(oldg);
	for(; eit != eend; ++eit) {
	  string ustr(get(vertex_name, oldg, source(*eit, oldg)));
	  string vstr(get(vertex_name, oldg, target(*eit, oldg)));
	  precsr16.insert_edge(ustr, vstr, get(edge_weight, oldg, *eit), get(edge_rtype, oldg, *eit));
	}

	KbGraph *new_g = new KbGraph(boost::edges_are_unsorted_multi_pass,
								 precsr16.E.begin(), precsr16.E.end(),
								 precsr16.eProp.begin(),
								 precsr16.m_vsize);

	tenp->m_g.reset(new_g);

	BGL_FORALL_VERTICES(v, *(tenp->m_g), Kb::boost_graph_t) {
	  (*(tenp->m_g))[v].name = precsr16.vProp[v].name;
	}

	tenp->m_vertexN = num_vertices(*(tenp->m_g));
	tenp->m_edgeN = num_edges(*(tenp->m_g));
	// relation sources
	std::set<std::string>(kbg.relsSource).swap(tenp->m_relsSource);
	// vertex map
	tenp->m_synsetMap.swap(precsr16.m_vMap);
	// relation types
	tenp->m_rtypes.m_strtypes.swap(kbg.rtypes);
	// Notes
	tenp->m_notes = kbg.notes;
	tenp->m_notes.push_back("--");
	tenp->m_notes.push_back("converted_to_2.0");

	p_instance = tenp;
  }


  ////////////////////////////////////////////////////////////////////////////////


  void Kb::add_comment(const string & str) {
	m_notes.push_back(str);
  }

  const vector<string> & Kb::get_comments() const {
	return m_notes;
  }

  ////////////////////////////////////////////////////////////////////////////////
  // bfs

  struct kb_bfs_init:public base_visitor<kb_bfs_init> {
  public:
	kb_bfs_init(Kb_vertex_t *v):m_v(v) { }
	typedef on_initialize_vertex event_filter;
	inline void operator()(Kb_vertex_t u, const KbGraph & g)
	{
	  m_v[u] = u;
	}
	Kb_vertex_t *m_v;
  };

  struct kb_bfs_pred:public base_visitor<kb_bfs_pred> {
  public:
	kb_bfs_pred(Kb_vertex_t *v):m_v(v) { }
	typedef on_tree_edge event_filter;
	inline void operator()(Kb_edge_t e, const KbGraph & g) {
	  m_v[target(e, g)] = source(e, g);
	}
	Kb_vertex_t *m_v;
  };


  // Note:
  //
  // after bfs, if (parents[v] == v) and (v != u), then u and v are not
  // connected in the graph.

  bool Kb::bfs(Kb_vertex_t src,
			   std::vector<Kb_vertex_t> & parents) const {

	size_t m = num_vertices(*m_g);
	if(parents.size() == m) {
	  std::fill(parents.begin(), parents.end(), Kb_vertex_t());
	} else {
	  vector<Kb_vertex_t>(m).swap(parents);  // reset parents
	}

	breadth_first_search(*m_g,
						 src,
						 boost::visitor(boost::make_bfs_visitor
										(boost::make_list(kb_bfs_init(&parents[0]),
														  kb_bfs_pred(&parents[0])))));
	return true;
  }


  bool Kb::dijkstra (Kb_vertex_t src,
					 std::vector<Kb_vertex_t> & parents) const {

	size_t m = num_vertices(*m_g);
	if(parents.size() == m) {
	  std::fill(parents.begin(), parents.end(), Kb_vertex_t());
	} else {
	  vector<Kb_vertex_t>(m).swap(parents);  // reset parents
	}

	// Hack to remove const-ness
    Kb & me = const_cast<Kb &>(*this);
	vector<float> w;
	vector<float> dist(m);
	property_map<Kb::boost_graph_t, boost::vertex_index_t>::type indexmap = get(vertex_index, *m_g);
	property_map<Kb::boost_graph_t, float edge_prop_t::*>::type wmap = get(&edge_prop_t::weight, *(me.m_g));

	dijkstra_shortest_paths(*m_g,
							src,
							predecessor_map(make_iterator_property_map(parents.begin(),
																	   get(vertex_index, *m_g))).
							distance_map(make_iterator_property_map(dist.begin(),
																	get(vertex_index, *m_g))).
							weight_map(wmap).
							vertex_index_map(indexmap));

	return true;
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Get shortest subgraphs

  class bfs_subg_terminate : public std::exception {};

  struct subg {
	vector<Kb_vertex_t> V;
	vector<vector<Kb_vertex_t> > E;
  };


  class bfs_subg_visitor : public default_bfs_visitor {

  public:
	bfs_subg_visitor(subg & s_, Kb_vertex_t u, int limit)
	  : m_sg(s_), m_idx(), m_i(0), m_t(0), m_max(limit) {
	  add_v(u);
	}

    void tree_edge(Kb_edge_t e, const KbGraph & g)
	{

	  Kb_vertex_t u = source(e,g);
	  Kb_vertex_t v = target(e,g);

	  // vertex v is new, but yet undiscovered
	  int v_i = add_v(v);
	  if (v_i == -1) return; // max limit reached.
	  int u_i = get_v(u);
	  add_e(u_i, v_i);

	  Kb_edge_t aux;
	  bool existsP;
	  tie(aux, existsP) = edge(v, u, g);
	  if (existsP) add_e(v_i, u_i); // as this edge is no more traversed.
	}

	void non_tree_edge(Kb_edge_t e, const KbGraph & g)
 	{
 	  // cross edge. source is previously stored for sure. target probably
 	  // is, unless max limit was reached
 	  int v_i = get_v(target(e,g));
	  if (v_i == -1) return; // target vertex not stored because max limit.
 	  int u_i = get_v(source(e,g));

 	  add_e(u_i, v_i);
 	}

	void discover_vertex(Kb_vertex_t u, const KbGraph & g)
	{
	  if (m_t == m_max) throw bfs_subg_terminate();
	  ++m_t;
	}

	int add_v(Kb_vertex_t v)
	{
	  if(m_i == m_max) return -1;
	  m_sg.V.push_back(v);
	  m_idx[v] = m_i;
	  m_sg.E.push_back(vector<Kb_vertex_t>());
	  int res = m_i;
	  ++m_i;
	  return res;
	}

	int get_v(Kb_vertex_t v) {
	  map<Kb_vertex_t, int>::iterator it=m_idx.find(v);
	  if(it == m_idx.end()) return -1;
	  return it->second;
	}

	void add_e(int u_i, int v_i) {
	  m_sg.E[u_i].push_back(m_sg.V[v_i]);
	}

  private:
	subg & m_sg;
	map<Kb_vertex_t, int> m_idx;
	int m_i; // num of inserted vertices
	int m_t; // time
	int m_max;
  };


  void Kb::get_subgraph(const string & src,
						vector<string> & V,
						vector<vector<string> > & E,
						size_t limit) {

	Kb_vertex_t u;
	bool aux;
	tie(u,aux) = get_vertex_by_name(src);
	if(!aux) return;

	subg sg;
	bfs_subg_visitor vis(sg, u, limit);

	try {
	  breadth_first_search(*m_g, u, boost::visitor(vis));
	} catch (bfs_subg_terminate & ) {}

	size_t N = sg.V.size();
	vector<string>(N).swap(V);
	vector<vector<string> >(N).swap(E);

	for(size_t i=0; i < N; ++i) {
	  V[i] = (*m_g)[sg.V[i]].name;
	  size_t m = sg.E[i].size();
	  vector<string> l(m);
	  for(size_t j=0; j < m; ++j) {
		l[j] =  (*m_g)[sg.E[i].at(j)].name;
	  }
	  E[i].swap(l);
	}
  }

  bool Kb::get_shortest_paths(const std::string & src,
							  const std::vector<std::string> & targets,
							  std::vector<std::vector<std::string> > & paths) {
	vector<Kb_vertex_t> parents;
	Kb_vertex_t u;
	bool aux;
	tie(u,aux) = get_vertex_by_name(src);
	if(!aux) return false;
	std::vector<std::vector<std::string> >().swap(paths);
	this->bfs(u, parents);
	for(std::vector<std::string>::const_iterator it = targets.begin(), end = targets.end();
		it != end; ++it) {
	  Kb_vertex_t v;
	  tie(v,aux) = get_vertex_by_name(*it);
	  if (!aux) continue;
	  if (parents[v] == v) continue; // either (u == v) or v is not connected to u.
	  paths.push_back(vector<string>());
	  vector<string> & P = paths.back();
	  // iterate until source is met
	  P.push_back(get_vertex_name(v));
	  while(1) {
		v = parents[v];
		P.push_back(get_vertex_name(v));
		if (v == u) break;
	  }
	  std::reverse(P.begin(), P.end());
	}
	return paths.size();
  }


  ////////////////////////////////////////////////////////////////////////////////
  // strings <-> vertex_id

  pair<Kb_vertex_t, bool> Kb::get_vertex_by_name(const std::string & str) const {
	map<string, Kb_vertex_t>::const_iterator it;

	it = m_synsetMap.find(str);
	if (it != m_synsetMap.end()) return make_pair(it->second, true);
	return make_pair(Kb_vertex_t(), false);
  }

  void Kb::edge_add_reltype(Kb_edge_t e, const string & rel) {
	m_rtypes.add_type(rel, (*m_g)[e].etype);
  }

  std::vector<std::string> Kb::edge_reltypes(Kb_edge_t e) const {
	return m_rtypes.tvector((*m_g)[e].etype);
  }

  float Kb::get_edge_weight(Kb_edge_t e) const {
	return (*m_g)[e].weight;
  }

  void Kb::set_edge_weight(Kb_edge_t e, float w) {
	(*m_g)[e].weight = w;
  }

  std::pair<Kb_out_edge_iter_t, Kb_out_edge_iter_t> Kb::out_neighbors(Kb_vertex_t u) {
  	return out_edges(u, *m_g);
  }

  std::pair<Kb_in_edge_iter_t, Kb_in_edge_iter_t> Kb::in_neighbors(Kb_vertex_t u) {
  	return in_edges(u, *m_g);
  }


  ////////////////////////////////////////////////////////////////////////////////
  // Query and retrieval

  // filter_mode
  //   0 -> no filter
  //   1 -> only words
  //   2 -> only concepts


  void vname_filter(const map<string, Kb_vertex_t> & theMap,
					const vector<float> & ranks,
					const KbGraph & g,
					vector<float> & outranks,
					vector<string> & vnames) {

	size_t v_m = theMap.size();

	vector<Kb_vertex_t> V(v_m);
	// empty output vectors
	vector<float>(v_m).swap(outranks);
	vector<string>(v_m).swap(vnames);

	map<string, Kb_vertex_t>::const_iterator m_it = theMap.begin();
	map<string, Kb_vertex_t>::const_iterator m_end = theMap.end();
	size_t v_i = 0;

	// Fill vertices index vector
	for(; m_it != m_end; ++m_it, ++v_i) {
	  V[v_i] = m_it->second;
	}
	// Sort to guarantee uniqueness
	sort(V);

	for(v_i = 0; v_i < v_m; ++v_i) {
	  outranks[v_i] = ranks[V[v_i]];
	  vnames[v_i] = g[V[v_i]].name;
	}
  }


  void Kb::filter_ranks_vnames(const vector<float> & ranks,
							   vector<float> & outranks,
							   vector<string> & vnames,
							   int filter_mode) const {

	size_t v_i, v_m;

	// No filtering
	v_m = ranks.size();
	outranks.resize(v_m);
	vnames.resize(v_m);
	for(v_i = 0; v_i < v_m; ++v_i) {
	  outranks[v_i] = ranks[v_i];
	  vnames[v_i] = (*m_g)[v_i].name;
	}
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Get static pageRank vector

  const std::vector<float> & Kb::static_prank() const {
	if (m_static_ppv.size()) return m_static_ppv;

	// Hack to remove const-ness
    Kb & me = const_cast<Kb &>(*this);

	if (m_vertexN == 0) return m_static_ppv; // empty graph
	vector<float> pv(m_vertexN, 1.0/static_cast<float>(m_vertexN));
	me.pageRank_ppv(pv, me.m_static_ppv);
	return m_static_ppv;
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Random

  Kb_vertex_t Kb::get_random_vertex() const {

	int r = g_randTarget(num_vertices(*m_g));

	return r;
  }

  ////////////////////////////////////////////////////////////////////////////////
  // read from textfile and create graph


  // Line format:
  //
  // u:synset v:synset t:rel i:rel s:source d:directed w:weight
  //
  // u: source vertex. Mandatory.
  // v: target vertex. Mandatory.
  // t: relation type (hyperonym, meronym, etc) of edge u->v. Optional.
  // i: (inverse) relation type of edge v->u (hyponym, etc). Optional. Useless on undirected graphs.
  // s: source of relation (wn30, kb17, etc). Optional.
  // d: wether the relation is directed. Optional, default is undirected.
  // w: relation weight. Must be positive. Optional.


  struct rel_parse {
	string u;
	string v;
	string rtype;
	string irtype;
	string src;
	float w;
	bool directed;

	rel_parse() : u(), v(), rtype(), irtype(), src(), w(0.0), directed(false) {}

  };

  bool parse_line(const string & line, rel_parse & out) {

	rel_parse res;

	char_separator<char> sep(" \t");
	tokenizer<char_separator<char> > tok(line, sep);
	tokenizer<char_separator<char> >::iterator it = tok.begin();
	tokenizer<char_separator<char> >::iterator end = tok.end();
	if (it == end) return false; // empty line
	for(;it != end; ++it) {

	  string str = *it;
	  if (str.length() < 3 || str[1] != ':') {
		throw runtime_error("parse_line error. Malformed line: " + line);
	  }
	  char f = str[0];
	  string val = str.substr(2);
	  if(!val.size()) continue;
	  switch (f) {
	  case 'u':
		res.u = val;
		break;
	  case 'v':
		res.v = val;
		break;
	  case 't':
		res.rtype = val;
		break;
	  case 'i':
		res.irtype = val;
		break;
	  case 's':
		res.src = val;
		break;
	  case 'w':
		res.w = lexical_cast<float>(val);
		break;
	  case 'd':
		res.directed = glVars::kb::keep_directed && lexical_cast<bool>(val);
		break;
	  default:
		throw runtime_error("parse_line error. Unknown value " + str);
		break;
	  }
	}
	if (!res.u.size()) throw runtime_error("parse_line error. No source vertex.");
	if (!res.v.size()) throw runtime_error("parse_line error. No target vertex.");
	out = res;
	return true;
  }

  void Kb::read_from_txt(istream & kbFile,
						 const set<string> & src_allowed) {
	string line;
	size_t line_number = 0;
	precsr_t csr_pre;

	set<string>::const_iterator srel_end = src_allowed.end();
	while(kbFile) {
	  vector<string> fields;
	  read_line_noblank(kbFile, line, line_number);
	  if(!kbFile) continue;
	  if (line[0] == '#') continue;
	  rel_parse f;
	  try {
		if (!parse_line(line, f)) continue;

		if (glVars::kb::filter_src) {
		  if (src_allowed.find(f.src) == srel_end) continue; // Skip this relation
		}

		if (f.u == f.v) continue; // no self-loops

		if (f.src.size()) {
		  this->add_relSource(f.src);
		}

		float w = f.w ? f.w : 1.0;
		// add edge

		// relation type
		// empty f.rtype unless glVars::kb::keep_reltypes

		if (!glVars::kb::keep_reltypes)
		  string().swap(f.rtype);

		csr_pre.insert_edge(f.u, f.v, w, f.rtype);

		// Insert v->u if undirected relation

		if (!f.directed || !glVars::kb::keep_directed) {
		  csr_pre.insert_edge(f.v, f.u, w, f.rtype);
		}
	  } catch (std::exception & e) {
		string msg(string(e.what()) + " in line " + lexical_cast<string>(line_number));
		if(!glVars::input::swallow) throw std::runtime_error(msg);
		if (glVars::debug::warning) {
		  cerr << msg << " (Skipping)\n";
		}
	  }
	}

	KbGraph *new_g = new KbGraph(boost::edges_are_unsorted_multi_pass,
								 csr_pre.E.begin(), csr_pre.E.end(),
								 csr_pre.eProp.begin(),
								 csr_pre.m_vsize);

	m_g.reset(new_g);
	// add_edges(csr_pre.E.begin(), csr_pre.E.end(),
	// 		  // csr_pre.eProp.begin(),
	// 		  // csr_pre.eProp.end(),
	// 		  g);

	BGL_FORALL_VERTICES(v, *m_g, Kb::boost_graph_t) {
	  (*m_g)[v].name = csr_pre.vProp[v].name;
	}

	m_vertexN = num_vertices(*m_g);
	m_edgeN = num_edges(*m_g);
	// Init vertex map
	m_synsetMap.swap(csr_pre.m_vMap);
  }

  void Kb::read_from_txt(const std::string & synsFileName,
						 const set<string> & src_allowed) {

	std::ifstream input_ifs(synsFileName.c_str(), ofstream::in);
	if (!input_ifs) {
	  throw runtime_error("Kb::read_from_txt error: Can't open " + synsFileName);
	}
	if(glVars::kb::v1_kb) {
	  throw runtime_error(synsFileName + " :sorry, the binary representation has an old format.");
	} else {
	  std::istream input_is(input_ifs.rdbuf());
	  this->read_from_txt(input_is, src_allowed);
	}
  }

  void Kb::display_info(std::ostream & o) const {

	o << "Relation sources: ";
	writeS(o, m_relsSource);
	if (m_notes.size()) {
	  o << "\nM_Notes: ";
	  writeV(o, m_notes);
	}
	size_t edge_n = num_edges(*m_g);

	o << "\n" << num_vertices(*m_g) << " vertices and " << edge_n << " edges.\n(Note that if graph is undirected you should divide the edge number by 2)" << endl;
	if (m_rtypes.size()) {
	  o << "Relations:";
	  writeV(o, m_rtypes.m_strtypes);
	  o << endl;
	}
  }


  std::pair<size_t, size_t> Kb::indeg_maxmin() const {

	size_t m = std::numeric_limits<size_t>::max();
	size_t M = std::numeric_limits<size_t>::min();

	size_t d = 0;

	graph_traits<KbGraph>::vertex_iterator it, end;
	tie(it, end) = vertices(*m_g);
	for(; it != end; ++it) {
	  d = in_degree(*it, *m_g);
	  if (d > M) M = d;
	  if (d < m) m = d;
	}
	return make_pair<size_t, size_t>(M, m);
  }

  std::pair<size_t, size_t> Kb::outdeg_maxmin() const {

	size_t m = std::numeric_limits<size_t>::max();
	size_t M = std::numeric_limits<size_t>::min();

	size_t d;

	graph_traits<KbGraph>::vertex_iterator it, end;
	tie(it, end) = vertices(*m_g);
	for(; it != end; ++it) {
	  d = out_degree(*it, *m_g);
	  if (d > M) M = d;
	  if (d < m) m = d;
	}
	return make_pair<size_t, size_t>(M, m);
  }


  int Kb::components() const {

	vector<size_t> v(num_vertices(*m_g));
	boost::iterator_property_map<
	  std::vector<size_t>::iterator,
	  boost::property_map<KbGraph, boost::vertex_index_t>::type> pm(v.begin(), get(boost::vertex_index, *m_g));

	int i = boost::strong_components(*m_g, pm);

	return i;

  }

  void Kb::ppv_weights(const vector<float> & ppv) {

	graph_traits<KbGraph>::edge_iterator it, end;

	tie(it, end) = edges(*m_g);
	for(; it != end; ++it) {
	  (*m_g)[*it].weight = ppv[target(*it, *m_g)];
	}
  }

  ////////////////////////////////////////////////////////////////////////////////
  // PageRank in KB


  // PPV version

  void Kb::pageRank_ppv(const vector<float> & ppv_map,
						 vector<float> & ranks) {

	typedef graph_traits<KbGraph>::edge_descriptor edge_descriptor;
	property_map<Kb::boost_graph_t, float edge_prop_t::*>::type weight_map = get(&edge_prop_t::weight, *m_g);
	prank::constant_property_map <edge_descriptor, float> cte_weight(1.0); // always return 1

	if (0 == m_out_coefs.size()) {
	  vector<float>(m_vertexN, 0.0f).swap(m_out_coefs);
	  if (glVars::prank::use_weight) {
		prank::init_out_coefs(*m_g,  &m_out_coefs[0], weight_map);
	  } else {
		prank::init_out_coefs(*m_g,  &m_out_coefs[0], cte_weight);
	  }
	}
	if (m_vertexN == ranks.size()) {
	  std::fill(ranks.begin(), ranks.end(), 0.0);
	} else {
	  vector<float>(m_vertexN, 0.0).swap(ranks); // Initialize rank vector
	}
	vector<float> rank_tmp(m_vertexN, 0.0);    // auxiliary rank vector

	if (glVars::prank::use_weight) {
	  prank::do_pageRank(*m_g, m_vertexN, &ppv_map[0],
						 weight_map, &ranks[0], &rank_tmp[0],
						 glVars::prank::num_iterations,
						 glVars::prank::threshold,
						 glVars::prank::damping,
						 m_out_coefs);
	} else {
	  prank::do_pageRank(*m_g, m_vertexN, &ppv_map[0],
						 cte_weight, &ranks[0], &rank_tmp[0],
						 glVars::prank::num_iterations,
						 glVars::prank::threshold,
						 glVars::prank::damping,
						 m_out_coefs);
	}
  }


  ////////////////////////////////////////////////////////////////////////////////
  // Debug

  ostream & Kb::dump_graph(std::ostream & o) const {
	o << "Sources: ";
	writeS(o, m_relsSource);
	o << endl;
	graph_traits<KbGraph>::vertex_iterator it, end;
	tie(it, end) = vertices(*m_g);
	for(;it != end; ++it) {
	  o << (*m_g)[*it].name;
	  graph_traits<KbGraph>::out_edge_iterator e, e_end;
	  tie(e, e_end) = out_edges(*it, *m_g);
	  if (e != e_end)
		o << "\n";
	  for(; e != e_end; ++e) {
		o << "  ";
		vector<string> r = edge_reltypes(*e);
		writeV(o, r);
		o << " " << (*m_g)[target(*e, *m_g)].name;
		o << " (" << (*m_g)[*e].weight << ")\n";
	  }
	}
	return o;
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Streaming

  static const size_t magic_id_v1 = 0x070201;
  static const size_t magic_id = 0x080826;
  static const size_t magic_id_csr = 0x110501;

  // CSR read

  vertex_prop_t read_vertex_prop_from_stream(istream & is) {

	string name;

	read_atom_from_stream(is, name);
	return vertex_prop_t(name);
  }

  edge_prop_t read_edge_prop_from_stream(istream & is) {

	float w;
	etype_t::value_type etype;

	read_atom_from_stream(is, w);
	read_atom_from_stream(is, etype);

	return edge_prop_t(w, etype);

  }

  void  Kb::read_from_stream (std::istream & is) {

	size_t vertex_n;
	size_t edge_n;
	size_t id;
	KbGraph *new_g;

	try {
	  read_atom_from_stream(is, id);
	  if (id != magic_id_csr) {
		if (id == magic_id_v1 || id == magic_id)
		  throw runtime_error("Old (pre 2.0) binary serialization format. Convert the graph to new format using the \"convert2.0\" utility.");
		else
		  throw runtime_error("Invalid id (same platform used to compile the KB?)");
	  }
	  read_set_from_stream(is, m_relsSource);
	  m_rtypes.read_from_stream(is);
	  read_map_from_stream(is, m_synsetMap);

	  read_atom_from_stream(is, id);
	  if (id != magic_id_csr) {
		throw runtime_error("Invalid id after reading maps");
	  }

	  read_atom_from_stream(is, edge_n);
	  read_atom_from_stream(is, vertex_n);
	  read_atom_from_stream(is, id);

	  if (id != magic_id_csr) {
		throw runtime_error("Invalid id after reading graph sizes");
	  }
	  new_g = new KbGraph();

	  new_g->m_forward.m_rowstart.resize(0);
	  read_vector_from_stream(is, new_g->m_forward.m_rowstart);
	  read_vector_from_stream(is, new_g->m_forward.m_column);
	  new_g->m_backward.m_rowstart.resize(0);
	  read_vector_from_stream(is, new_g->m_backward.m_rowstart);
	  read_vector_from_stream(is, new_g->m_backward.m_column);
	  read_vector_from_stream(is, new_g->m_backward.m_edge_properties);

	  for(size_t i = 0; i != vertex_n; ++i) {
		new_g->vertex_properties().m_vertex_properties.push_back(read_vertex_prop_from_stream(is));
	  }

	  for(size_t i = 0; i != edge_n; ++i) {
		new_g->m_forward.m_edge_properties.push_back(read_edge_prop_from_stream(is));
	  }

	  read_atom_from_stream(is, id);
	  if (id != magic_id_csr) {
		throw runtime_error("Invalid id after reading graph");
	  }
	  read_vector_from_stream(is, m_notes);
	} catch (std::exception & e) {
	  throw runtime_error(string("Error when reading serialized graph: ") + e.what());
	}

	m_g.reset(new_g);
	vector<float>().swap(m_static_ppv); // empty static rank vector

	m_vertexN = vertex_n;
	m_edgeN = edge_n;
	assert(num_vertices(*m_g) == m_vertexN);
	assert(num_edges(*m_g) == m_edgeN);
  }

  // write

  ostream & write_vertex_prop_to_stream(ostream & o,
										const vertex_prop_t & p) {
	write_atom_to_stream(o, p.name);
	return o;
  }

  ostream & write_edge_prop_to_stream(ostream & o,
									  const edge_prop_t & ep) {
	write_atom_to_stream(o, ep.weight);
	write_atom_to_stream(o, ep.etype);
	return o;
  }

  ostream & Kb::write_to_stream(ostream & o) const {

	// First write maps

	assert(m_vertexN == num_vertices(*m_g));
	assert(m_edgeN == num_edges(*m_g));

	write_atom_to_stream(o, magic_id_csr);

	write_vector_to_stream(o, m_relsSource);
	m_rtypes.write_to_stream(o);
	write_map_to_stream(o, m_synsetMap);

	write_atom_to_stream(o, magic_id_csr);

	write_atom_to_stream(o, m_edgeN);
	write_atom_to_stream(o, m_vertexN);

	write_atom_to_stream(o, magic_id_csr);

	write_vector_to_stream(o, m_g->m_forward.m_rowstart);
	write_vector_to_stream(o, m_g->m_forward.m_column);
	write_vector_to_stream(o, m_g->m_backward.m_rowstart);
	write_vector_to_stream(o, m_g->m_backward.m_column);
	write_vector_to_stream(o, m_g->m_backward.m_edge_properties);

	//	write_vector_to_stream(o, m_g->vertex_properties().m_vertex_properties);
	//	write_vector_to_stream(o, m_g->m_forward.m_edge_properties);

	size_t vProp_n = m_g->vertex_properties().m_vertex_properties.size();
	assert(vProp_n == m_vertexN);
	for(size_t i = 0; i != vProp_n; ++i) {
	  write_vertex_prop_to_stream(o, m_g->vertex_properties().m_vertex_properties[i]);
	}

	size_t eProp_n = m_g->m_forward.m_edge_properties.size();
	assert(eProp_n == m_edgeN);
	for(size_t i = 0; i != eProp_n; ++i) {
	  write_edge_prop_to_stream(o, m_g->m_forward.m_edge_properties[i]);
	}

	write_atom_to_stream(o, magic_id_csr);

	write_vector_to_stream(o, m_notes);
	return o;
  }


  void Kb::write_to_binfile (const string & fName) {

	ofstream fo(fName.c_str(),  ofstream::binary|ofstream::out);
	if (!fo) {
	  cerr << "Error: can't create" << fName << endl;
	  exit(-1);
	}
	write_to_stream(fo);
  }

  // text write

  ostream & write_to_textstream(const KbGraph & g, ostream & o) {

	graph_traits<KbGraph>::edge_iterator e_it, e_end;

	tie(e_it, e_end) = edges(g);
	for(; e_it != e_end; ++e_it) {
	  o << "u:" << g[source(*e_it, g)].name << " ";
	  o << "v:" << g[target(*e_it, g)].name << " d:1\n";
	}
	return o;
  }

  void Kb::write_to_textfile (const string & fName) {

	ofstream fo(fName.c_str(),  ofstream::out);
	if (!fo) {
	  cerr << "Error: can't create" << fName << endl;
	  exit(-1);
	}
	write_to_textstream(*m_g, fo);
  }


   /**
   *
   * TODO:  Egiteko funtzioa
   *
   */

     std::map<Kb_vertex_t, std::pair<std::map<Kb_vertex_t, Kb_vertex_t>, std::map<Kb_vertex_t, Kb_vertex_t> > > graphToMap(Kb_vertex_t initial, std::map<Kb_vertex_t, std::pair<std::map<Kb_vertex_t, Kb_vertex_t>,std::map<Kb_vertex_t, Kb_vertex_t> > > graph){
 		cout << "Graph is going to be created" << endl;
      	Kb::create_from_binfile("/media/datuak/Proiektua/konpilatutakoProbaGrafoa.bin");
      	Kb & kb = ukb::Kb::instance();
      	Kb_vertex_t u;
      	bool aux = false;
      	if(initial == 0){
      		tie(u, aux) = kb.get_vertex_by_name("bn00000001n");
      	}else{
      		u = initial;
      	}
      	try{ //If the vertex u exists in the graph we will end the execution.  If not, we will get all its IN and OUT vertices and put in the graph
      		graph.at(u);
      		//return graph;
      	}catch(...){
      		std::map<Kb_vertex_t, Kb_vertex_t> outVertexMap = Kb::getVertices(u, false);
      		std::map<Kb_vertex_t, Kb_vertex_t> inVertexMap = Kb::getVertices(u, true);
      		std::pair <std::map<Kb_vertex_t>, std::map<Kb_vertex_t> > maps;  
      		maps = std::make_pair(outVertexMap, inVertexMap);
      		graph.insert(u, maps);

			typedef std::map<Kb_vertex_t, Kb_vertex_t>::iterator it_type;
			for(it_type iterator = outVertexMap.begin(); iterator != outVertexMap.end(); iterator++) {
    		
    				graph = graphToMap(iterator->second, graph);	
			}
			for(it_type iterator = inVertexMap.begin(); iterator != inVertexMap.end(); iterator++) {
    			
    				graph = graphToMap(iterator->second, graph);	
			}
      		//return graph;
      	}
      return graph;
     }

     void buildSemanticSignatures (Kb_vertex_t initial){
     	std::map<Kb_vertex_t, std::pair<std::map<Kb_vertex_t, Kb_vertex_t>, std::map<Kb_vertex_t, Kb_vertex_t> > > graph;
     	typedef std::map<Kb_vertex_t, std::pair<std::map<Kb_vertex_t, Kb_vertex_t>, std::map<Kb_vertex_t, Kb_vertex_t> > >::iterator itGraph;


		std::map<Kb_vertex_t, Kb_vertex_t> vSecondInMap, vSecondOutMap, inMap, outMap;
     	Kb_vertex_t v;
     	bool aux  = true;
     	std::vector<Kb_vertex_t> inVector, outVector;
     	Kb_in_edge_iter_t edge_it;

     		graph = graphToMap(0, graph);
     		if(v == 0){
     			Kb & kb = ukb::Kb::instance();
      			tie(v, aux) = kb.get_vertex_by_name("bn00000001n");
      		}else{
      			v = initial;
      		}

      		if(aux != 0){
    			tie (inMap, outMap) = graph.at(v);  // Get IN and OUT vertices of the first vertex (v)
    		
     				for(itGraph it = outMap.begin(); it != outMap.end(); ++it){ 
     					
     						tie (vSecondInMap, vSecondOutMap) = graph.at(it->second); 

		     				for(itGraph itSecond = vSecondOutMap.begin(); itSecond != vSecondOutMap.end(); ++itSecond){ 
		     					try{

		     						inMap.at(itSecond->second);
		     						//Increase weight of the edge between the first vertex and the out vertex

		     					}catch(...){
		     						//Nothing to do
		     					}
		     				}
     					
     				}

     				//Iterate recursively trough the graph

     		}

     }

   /*
  void Kb::buildSemanticSignatures(string vertexName){
		
		//Creates two vectors (all IN and OUT vertices) of the given vertex.
		//The IN vector will be a simple vector.  
		//Each element of the OUT vector will contain the vertex itself and its all OUT vertices vector
		//Then we will analyze the all elements in IN vector seeing if they are in the each element of OUT vectors OUT vertex.  When this condition happens we add in one unit to the weight of the edge.
		//Finishing, we will do this again for each and every one of the OUT vertices, recursively.

	    Kb_out_edge_iter_t it, end, itOut, itOutEnd;
	    Kb_in_edge_iter_t in_it, in_end, itIn, itInEnd;
	    Kb_vertex_t u; //Source vertex
	    Kb_vertex_t v; //Auxiliar vertex, i.e returning OUT or IN vertex
	    Kb_vertex_t x; //Another auxiliar vertex.
	    bool aux;
	    float w = 1.0;

	    
		std::map<std::string, Kb_vertex_t> inVertexMap; //HashMap for the IN vertices. (UNUSED?!?!) It'll be used to collect the IN vertices of each out vertex of the main vertex.  KEY:  Vertex name.  VALUE:  Vertex object.
	    std::map<std::string, Kb_vertex_t> outVertexMap; //It has the same purpose of the inVertexMap but instead it0ll be used for the OUT vertices.  KEY:  Vertex name.  VALUE:  Vertex object.

	    std::vector<Kb_vertex_t> inVertexVector; //Vector containing IN vertices of the main vertex.  
	    //The main vertex will be considered as the vertes wich is analyzed. This is going to be changing to analyze all vertices of the graph.
	    std::vector<std::pair<Kb_vertex_t, outVertexMap>> outVertexVector; //Vector containing OUT vertices of the main vertex
	    


	    Kb::create_from_binfile("/media/datuak/Proiektua/konpilatutakoProbaGrafoa.bin");
	    Kb & kb = ukb::Kb::instance();
		//tie(u, aux) = kb.get_vertex_by_name("00000001n");
		tie(u, aux) = kb.get_vertex_by_name(vertexName);
	    if(!aux){
	    	tie(itIn, itInEnd) = Kb::in_neighbors(u);
	    	tie(itOut, itOutEnd) = kb.out_neighbors(u);
	    	for(; itOut != itOutEnd; ++itOut) {
	    		 v = kb.edge_source(*itOut);
	    		 outVertexVector.push_back(pair<v, outVertexMap);  //Now we are assigning an empty map to the pair <vertex, map>.  We will fill it up 

	    		if(itIn != itInEnd){
	    		 v = kb.edge_source(*itIn);
	    		 inVertexVector.push_back(v);
	    		 ++itIn;
	    		}
	    	}
	    	//Finished analyzing in and out vertices of the node.  
	    	//Now it's time to identify triangular relations between nodes. 
	    	for(int i=0; i < outVertexVector.length; i++){
				if(outVertexMap.size > 0){ //If the map is not empty clear the content
					outVertexMap.clear();
				}
	    		tie (v, outVertexMap) = outVertexVector[i]; //Fetch one OUT vertex and its OUT vertices
	    		tie(itOut, itOutEnd) = kb.out_neighbors(v); //Take all it's OUT vertices and put in a map
	    		for(; itOut != itOutEnd; ++itOutEnd){
	    			x = kb.edge_source(*itOut);
	    			outVertexMap.insert(VERTEXNAME, x);
	    		}
	    		outVertexVector[i] = std::make_pair(v, outVertexMap); //Reinsert in the ORIGINAL VERTEX OUT vertex vector.  
	    	}

	    	//Now it's prepared to see if it has triangular relations.  Let's analyze
	    	outVertexMap.clear();
	    	for(int i = 0; i < inVertexVector.length; i++){
	    		for(int j = 0; j < outVertexVector.length; j++){
	    			tie (x, outVertexMap) = outVertexVector[j];
	    			outVertexMap.at(inVertexVector[i]); //Catch the exception out_of_range.  This happens when the provided key is not a key of an element in the map.  
	    			//In this case we have to BREAK the outVertexVector loop.  BUT BEFORE, we have to assign weight 1 to the edge.
	    			//If there is not an exception we say that there is a triangular relationship 
	    			w  + = 0.01;

	    			//In all cases, 
	    			kb.set_edge_weight(*in_it, w);
	    		}
	    	}

	    	//For each vertex run this procedure
	    	for( ;itOut != itOutEnd; ++itOut){
				v = kb.edge_source(*itOut);
				buildSemanticSignatures(v);//
	    	}
	    	
	    	
	    }
	    kb.dump_graph(cout);
 
  }*/


	void Kb::rwr(Kb_vertex_t v, float alpha, int n, float p){

	/**
	*v:  The starting vertex
	*alpha:  Restart probability.  Should be between 0 and 1.
	*n: Number of step to be executed
	*p:  The transition probability. Should be between 0 and 1.
	*/

	  srand (static_cast <unsigned> (time(0)));

      //typedef boost::mt19937 RNGType; //The mersenne twister generator.
      //RNGType rng(time(0));  //Instance of the twister generator.  If we don't put time(0)  the random numbers will be the same always.  
      //That's the main reason to put a time dependant generator creator.  In all executions the generated numbers will be different because the generator itself is different.
      //boost::uniform_float<> distr(0, 1);  //Distribution distance.  We want from 0 to 1.
      //boost::variate_generator< RNGType, boost::uniform_float<> > dice(rng, distr);  //Takes the raw numbers and the distribution, and creates the random numbers.
   	  std::vector<Kb_vertex_t> semSignVector;
   	  Kb_vertex_t initial = v;

      
      while(n > 0){
      	
      	float random = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
      	if(random > alpha){
      		typedef std::map<Kb_vertex_t, Kb_vertex_t>::iterator it_type;
      		std::map<Kb_vertex_t, Kb_vertex_t> auxiliar = getVertices(v, false);
      		int i = int ((auxiliar.size() * random)+0.5); //Choose an out neighbor using the random number.  We have to use an integer as an index and now we are rounding it before use.
      		int j = 0;
      		for(it_type iterator = auxiliar.begin(); iterator != auxiliar.end(); ++iterator){
      			if(j == i){
      				break;
      			}
      		}
      		semSignVector.push_back(auxiliar.at(iterator->first)); //Add to the semantic signature group
      		v = auxiliar.at(iterator->first);
       	}else{
       		v = initial;
      	}
      	n--;
    	 }
	}

	std::map<Kb_vertex_t, Kb_vertex_t> Kb::getVertices(Kb_vertex_t source, bool inFlag){
		//Returns a map containing source vertex's all in and out vertices in a array.  
		//PARAMETERS:  
		//source:  The source vertex to analyze.  
		//inFlag:  A boolean flag.  When TRUE this function will return sources all IN vertices.  When FALSE, all OUT vertices are returned.
	 	Kb_in_edge_iter_t itOut, itOutEnd;
	 	std::vector<Kb_vertex_t> vertexMap;
	 	Kb_vertex_t x;
	 	Kb::create_from_binfile("/media/datuak/Proiektua/konpilatutakoProbaGrafoa.bin");
	 	Kb & kb = ukb::Kb::instance();
	 	if(inFlag){
			tie(itOut, itOutEnd) = kb.in_neighbors(source);
			for(; itOut < itOutEnd; ++itOut){
				x = kb.edge_target(*itOut);
	    		vertexMap.insert(x, x);
			}
	 	}else{
			tie(itOut, itOutEnd) = kb.out_neighbors(source);
			for(; itOut < itOutEnd; ++itOut){
				x = kb.edge_source(*itOut);
	    		vertexMap.insert(x, x);
			}
	 	}
		return vertexMap;
	}




}
