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

  ////////////////////////////////////////////////////////////////////////////////


  void Kb::add_comment(const string & str) {
	notes.push_back(str);
  }

  const vector<string> & Kb::get_comments() const {
	return notes;
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
	  m_v[target(e, g)] = source(e,g);
	}
	Kb_vertex_t *m_v;
  };


  bool Kb::bfs (Kb_vertex_t src,
				std::vector<Kb_vertex_t> & parents) const {

	size_t m = num_vertices(g);
	if(parents.size() == m) {
	  std::fill(parents.begin(), parents.end(), Kb_vertex_t());
	} else {
	  vector<Kb_vertex_t>(m).swap(parents);  // reset parents
	}

	breadth_first_search(g,
						 src,
						 boost::visitor(boost::make_bfs_visitor
										(boost::make_list(kb_bfs_init(&parents[0]),
														  kb_bfs_pred(&parents[0])))));
	return true;
  }


  bool Kb::dijkstra (Kb_vertex_t src,
					 std::vector<Kb_vertex_t> & parents) const {

	size_t m = num_vertices(g);
	if(parents.size() == m) {
	  std::fill(parents.begin(), parents.end(), Kb_vertex_t());
	} else {
	  vector<Kb_vertex_t>(m).swap(parents);  // reset parents
	}

	vector<float> dist(m);

	dijkstra_shortest_paths(g,
							src,
							predecessor_map(&parents[0]).
							distance_map(&dist[0]));
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
	  breadth_first_search(g, u, boost::visitor(vis));
	} catch (bfs_subg_terminate & ) {}

	size_t N = sg.V.size();
	vector<string>(N).swap(V);
	vector<vector<string> >(N).swap(E);

	for(size_t i=0; i < N; ++i) {
	  V[i] = get(vertex_name, g, sg.V[i]);
	  size_t m = sg.E[i].size();
	  vector<string> l(m);
	  for(size_t j=0; j < m; ++j) {
		l[j] =  get(vertex_name, g, sg.E[i].at(j));
	  }
	  E[i].swap(l);
	}
  }


  ////////////////////////////////////////////////////////////////////////////////
  // strings <-> vertex_id

  pair<Kb_vertex_t, bool> Kb::get_vertex_by_name(const std::string & str,
												 unsigned char flags) const {
	map<string, Kb_vertex_t>::const_iterator it;

	if(flags & Kb::is_concept) {
	  it = synsetMap.find(str);
	  if (it != synsetMap.end()) return make_pair(it->second, true);
	}
	// is it a word ?
	if(flags & Kb::is_word) {
	  it = wordMap.find(str);
	  if (it != wordMap.end()) return make_pair(it->second, true);
	}
	return make_pair(Kb_vertex_t(), false);
  }

  Kb_vertex_t Kb::InsertNode(const string & name, unsigned char flags) {
	coef_status = 0; // reset out degree coefficients
	if (static_ppv.size()) vector<float>().swap(static_ppv); // empty static rank vector
	Kb_vertex_t u = add_vertex(g);
	put(vertex_name, g, u, name);
	put(vertex_flags, g, u, flags);
	return u;
  }

  Kb_vertex_t Kb::find_or_insert_synset(const string & str) {
	map<string, Kb_vertex_t>::iterator it;
	bool insertedP;
	tie(it, insertedP) = synsetMap.insert(make_pair(str, Kb_vertex_t()));
	if(insertedP) {
	  // new vertex
	  it->second = InsertNode(str, Kb::is_concept);
	}
	return it->second;
  }

  Kb_vertex_t Kb::find_or_insert_word(const string & str) {
	map<string, Kb_vertex_t>::iterator it;
	bool insertedP;
	tie(it, insertedP) = wordMap.insert(make_pair(str, Kb_vertex_t()));
	if(insertedP) {
	  // new vertex
	  it->second = InsertNode(str, Kb::is_word);
	}
	return it->second;
  }

  Kb_edge_t Kb::find_or_insert_edge(Kb_vertex_t u, Kb_vertex_t v,
									float w) {

	Kb_edge_t e;
	bool existsP;

	if (u == v)
	  throw runtime_error("Can't insert self loop !");
	//if (w != 1.0) ++w; // minimum weight is 1
	tie(e, existsP) = edge(u, v, g);
	if(!existsP) {
	  coef_status = 0; // reset out degree coefficients
	  if (static_ppv.size()) vector<float>().swap(static_ppv); // empty static rank vector
	  e = add_edge(u, v, g).first;
	  put(edge_weight, g, e, w);
	  put(edge_rtype, g, e, static_cast<boost::uint32_t>(0));
	}
	return e;
  }

  void Kb::unlink_vertex(Kb_vertex_t u) {
	clear_vertex(u, g);
  }

  size_t Kb::unlink_dangling() {

	size_t n = 0;
	graph_traits<KbGraph>::vertex_iterator it, end;
	tie(it, end) = vertices(g);
	for(; it != end; ++it) {
	  if(out_degree(*it, g) == 0 && in_degree(*it, g) != 0) {
		clear_vertex(*it, g);
		++n;
	  }
	}
	return n;
  }

  vector<string>::size_type get_reltype_idx(const string & rel,
											vector<string> & rtypes) {

	vector<string>::iterator it = rtypes.begin();
	vector<string>::iterator end = rtypes.end();
	vector<string>::size_type idx = 0;

	for(;it != end; ++it) {
	  if (*it == rel) break;
	  ++idx;
	}
	if (it == end) {
	  // new relation type
	  rtypes.push_back(rel);
	}
	if (idx > 31) {
	  throw runtime_error("get_rtype_idx error: too many relation types !");
	}
	return idx;
  }

  void Kb::edge_add_reltype(Kb_edge_t e, const string & rel) {
	boost::uint32_t m = get(edge_rtype, g, e);
	vector<string>::size_type idx = get_reltype_idx(rel, rtypes);
	m |= (1UL << idx);
	put(edge_rtype, g, e, m);
  }

  std::vector<std::string> Kb::get_edge_reltypes(Kb_edge_t e) const {
	vector<string> res;
	boost::uint32_t m = get(edge_rtype, g, e);
	vector<string>::size_type idx = 0;
	boost::uint32_t i = 1;
	while(idx < 32) {
	  if (m & i) {
		res.push_back(rtypes[idx]);
	  }
	  idx++;
	  i <<= 1;
	}
	return res;
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Query and retrieval

  bool Kb::vertex_is_synset(Kb_vertex_t u) const {
	return !vertex_is_word(u);
  }

  bool Kb::vertex_is_word(Kb_vertex_t u) const {
	return (get(vertex_flags, g, u) & Kb::is_word);
  }

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
	  vnames[v_i] = get(vertex_name, g, V[v_i]);
	}
  }


  void Kb::filter_ranks_vnames(const vector<float> & ranks,
							   vector<float> & outranks,
							   vector<string> & vnames,
							   int filter_mode) const {

	size_t v_i, v_m;

	switch(filter_mode) {
	case 0:
	  // No filtering
	  v_m = ranks.size();
	  outranks.resize(v_m);
	  vnames.resize(v_m);
	  for(v_i = 0; v_i < v_m; ++v_i) {
		outranks[v_i] = ranks[v_i];
		vnames[v_i] = get(vertex_name, g, v_i);
	  }
	  break;
	case 1:
	  // Only words
	  vname_filter(wordMap, ranks, g, outranks, vnames);
	  break;
	case 2:
	  // Only concepts
	  vname_filter(synsetMap, ranks, g, outranks, vnames);
	  break;
	}
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Get static pageRank vector

  const std::vector<float> & Kb::static_prank() const {
	if (static_ppv.size()) return static_ppv;
	size_t N = num_vertices(g);

	// Hack to remove const-ness
    Kb & me = const_cast<Kb &>(*this);
	me.out_coefs.resize(N);

	size_t N_no_isolated = prank::init_out_coefs(g, &me.out_coefs[0]);

	if (N_no_isolated == 0) return static_ppv; // empty graph
	vector<float> pv(N, 1.0/static_cast<float>(N_no_isolated));
	me.pageRank_ppv(pv, me.static_ppv);
	return static_ppv;
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Random

  Kb_vertex_t Kb::get_random_vertex() const {

	int r = g_randTarget(num_vertices(g));

	return r;
  }

  ////////////////////////////////////////////////////////////////////////////////
  // read from textfile and create graph

  // Read the actual KB file

  void read_kb_v1(ifstream & kbFile,
				  const set<string> & src_allowed,
				  Kb * kb) {
	string line;
	int line_number = 0;

	set<string>::const_iterator srel_end = src_allowed.end();
	while(kbFile) {
	  vector<string> fields;
	  std::getline(kbFile, line, '\n');
	  trim_spaces(line);
	  line_number++;
	  char_separator<char> sep(" \t");
	  tokenizer<char_separator<char> > tok(line, sep);
	  copy(tok.begin(), tok.end(), back_inserter(fields));
	  if (fields.size() == 0) continue; // blank line
	  if (fields.size() < 4) {
		throw runtime_error("read_kb error: Bad line: " + lexical_cast<string>(line_number));
	  }

	  if (glVars::kb::filter_src) {
		if (src_allowed.find(fields[3]) == srel_end) continue; // Skip this relation
	  }

	  Kb_vertex_t u = kb->find_or_insert_synset(fields[0]);
	  Kb_vertex_t v = kb->find_or_insert_synset(fields[1]);

	  if (u == v) continue; // no self-loops

	  kb->add_relSource(fields[3]);

	  // last element says if relation is directed
	  bool directed = glVars::kb::keep_directed && (fields.size() > 4 && lexical_cast<int>(fields[4]) != 0);
	  // add edge
	  kb->find_or_insert_edge(u, v, 1.0);
	  if (!directed)
		kb->find_or_insert_edge(v, u, 1.0);
	}
  }


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
  // w: relation weigth. Must be positive. Optional.


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

	set<string>::const_iterator srel_end = src_allowed.end();
	try {
	  while(kbFile) {
		vector<string> fields;
		read_line_noblank(kbFile, line, line_number);
		if(!kbFile) continue;
		if (line[0] == '#') continue;
		rel_parse f;
		if (!parse_line(line, f)) continue;

		if (glVars::kb::filter_src) {
		  if (src_allowed.find(f.src) == srel_end) continue; // Skip this relation
		}
		Kb_vertex_t u = this->find_or_insert_synset(f.u);
		Kb_vertex_t v = this->find_or_insert_synset(f.v);

		if (u == v) continue; // no self-loops

		this->add_relSource(f.src);

		float w = f.w ? f.w : 1.0;
		// add edge
		Kb_edge_t e1 = this->find_or_insert_edge(u, v, w);
		// relation type
		if (glVars::kb::keep_reltypes && f.rtype.size()) {
		  this->edge_add_reltype(e1, f.rtype);
		}
		if (!f.directed || !glVars::kb::keep_directed) {
		  Kb_edge_t e2 = this->find_or_insert_edge(v, u, w);
		  if(glVars::kb::keep_reltypes && f.rtype.size()) {
			string aux = f.irtype.size() ? f.irtype : f.rtype;
			if(aux.size()) {
			  this->edge_add_reltype(e2, aux);
			}
		  }
		}
	  }
	} catch (std::exception & e) {
	  throw std::runtime_error(string(e.what()) + " in line " + lexical_cast<string>(line_number));
	}
  }

  void Kb::read_from_txt(const std::string & synsFileName,
						 const set<string> & src_allowed) {

	std::ifstream input_ifs(synsFileName.c_str(), ofstream::in);
	if (!input_ifs) {
	  throw runtime_error("Kb::read_from_txt error: Can't open " + synsFileName);
	}
	if(glVars::kb::v1_kb) {
	  read_kb_v1(input_ifs, src_allowed, this);
	} else {
	  std::istream input_is(input_ifs.rdbuf());
	  this->read_from_txt(input_is, src_allowed);
	}
  }

  void Kb::display_info(std::ostream & o) const {

	o << "Relation sources: ";
	writeS(o, relsSource);
	if (notes.size()) {
	  o << "\nNotes: ";
	  writeV(o, notes);
	}
	size_t edge_n = num_edges(g);

	o << "\n" << num_vertices(g) << " vertices and " << edge_n << " edges.\n(Note that if graph is undirected you should divide the edge number by 2)" << endl;
	if (rtypes.size()) {
	  o << "Relations:";
	  writeV(o, rtypes);
	  o << endl;
	}
  }


  std::pair<size_t, size_t> Kb::indeg_maxmin() const {

	size_t m = std::numeric_limits<size_t>::max();
	size_t M = std::numeric_limits<size_t>::min();

	size_t d;

	graph_traits<KbGraph>::vertex_iterator it, end;
	tie(it, end) = vertices(g);
	for(; it != end; ++it) {
	  d = in_degree(*it, g);
	  if (d > M) M = d;
	  if (d < m) m = d;
	}
	return make_pair<size_t, size_t>(m, M);
  }

  std::pair<size_t, size_t> Kb::outdeg_maxmin() const {

	size_t m = std::numeric_limits<size_t>::max();
	size_t M = std::numeric_limits<size_t>::min();

	size_t d;

	graph_traits<KbGraph>::vertex_iterator it, end;
	tie(it, end) = vertices(g);
	for(; it != end; ++it) {
	  d = out_degree(*it, g);
	  if (d > M) M = d;
	  if (d < m) m = d;
	}
	return make_pair<size_t, size_t>(m, M);
  }


  int Kb::components() const {

	//	std::vector<int> component(num_vertices(g)), discover_time(num_vertices(g));
	//	std::vector<default_color_type> color(num_vertices(g));
	//	std::vector<Vertex> root(num_vertices(g));
	vector<int> v(num_vertices(g));
	int i = boost::strong_components(g,&v[0]);

	return i;

  }


  ////////////////////////////////////////////////////////////////////////////////
  // Add token vertices and link them to synsets


  typedef pair<Kb_vertex_t, float> Syn_elem;

  void create_w2wpos_maps(const string & word,
						  vector<string> & wPosV,
						  map<string, vector<Syn_elem> > & wPos2Syns) {

	bool auxP;
	const Kb & kb = Kb::instance();

	WDict_entries syns = WDict::instance().get_entries(word);

	for(size_t i = 0; i < syns.size(); ++i) {

	  string wpos(word.size() + 2, '#');
	  string::iterator sit = copy(word.begin(), word.end(), wpos.begin());
	  ++sit; // '#' char
	  *sit = syns.get_pos(i); // the pos

	  Kb_vertex_t u;
	  tie(u, auxP) = kb.get_vertex_by_name(syns.get_entry(i), Kb::is_concept);
	  if(!auxP) {
		if (glVars::debug::warning)
		  cerr << "W:Kb::add_tokens: warning: " << syns.get_entry(i) << " is not in KB.\n";
		continue;
	  }

	  map<string, vector<Syn_elem> >::iterator m_it;

	  tie(m_it, auxP) = wPos2Syns.insert(make_pair(wpos, vector<Syn_elem>()));
	  if (auxP) {
		// first appearence of word#pos
		wPosV.push_back(wpos);
	  }
	  m_it->second.push_back(make_pair(u, syns.get_freq(i)));
	}
  }

  void insert_wpos(const string & word,
				   vector<string> & wPosV,
				   map<string, vector<Syn_elem> > & wPos2Syns) {

	Kb & kb = Kb::instance();

	// insert word
	Kb_vertex_t word_v = kb.find_or_insert_word(word);

	vector<string>::iterator wpos_str_it = wPosV.begin();
	vector<string>::iterator wpos_str_end = wPosV.end();
	for (; wpos_str_it != wpos_str_end; ++wpos_str_it) {
	  //insert word#pos
	  Kb_vertex_t wpos_v = kb.find_or_insert_word(*wpos_str_it);
	  // link word to word#pos
	  kb.find_or_insert_edge(word_v, wpos_v, 1.0);
	  // link word#pos to synsets
	  vector<Syn_elem>::const_iterator syns_it = wPos2Syns[*wpos_str_it].begin();
	  vector<Syn_elem>::const_iterator syns_end = wPos2Syns[*wpos_str_it].end();
	  for(;syns_it != syns_end; ++syns_it) {
		kb.find_or_insert_edge(wpos_v, syns_it->first, syns_it->second);
	  }
	}
  }

  void insert_word(const string & word) {

	Kb & kb = Kb::instance();

	WDict_entries syns = WDict::instance().get_entries(word);

	if (!syns.size()) return;

	Kb_vertex_t u = kb.find_or_insert_word(word);

	for(size_t i = 0; i < syns.size(); ++i) {
	  bool auxP;
	  Kb_vertex_t v;
	  tie(v, auxP) = kb.get_vertex_by_name(syns.get_entry(i), Kb::is_concept);
	  if(!auxP) {
		if (glVars::debug::warning)
		  cerr << "W:Kb::add_tokens: warning: " << syns.get_entry(i) << " is not in KB.\n";
		continue;
	  }
	  float w = syns.get_freq(i);
	  // (directed) link word -> concept
	  kb.find_or_insert_edge(u, v, w);
	}
  }

  void Kb::add_dictionary() {

	WDict & w2syn = WDict::instance();

	vector<string>::const_iterator word_it = w2syn.get_wordlist().begin();
	vector<string>::const_iterator word_end = w2syn.get_wordlist().end();

	for(; word_it != word_end; ++word_it) {
	  add_token(*word_it);
	}
  }


  void Kb::add_token(const string & token) {

	vector<string> wPosV;
	map<string, vector<Syn_elem> > wPos2Syns;

	if (glVars::input::filter_pos) {
	  // Create w2wPos and wPos2Syns maps

	  create_w2wpos_maps(token, wPosV, wPos2Syns);

	  // Add vertices and link them in the KB
	  insert_wpos(token, wPosV, wPos2Syns);
	} else {
	  insert_word(token);
	}
  }

  void Kb::ppv_weights(const vector<float> & ppv) {

	graph_traits<KbGraph>::edge_iterator it, end;

	tie(it, end) = edges(g);
	for(; it != end; ++it) {
	  put(edge_weight, g, *it,
		  ppv[target(*it, g)]);
	}
  }

  ////////////////////////////////////////////////////////////////////////////////
  // PageRank in KB


  // PPV version

  void Kb::pageRank_ppv(const vector<float> & ppv_map,
						vector<float> & ranks) {

	typedef graph_traits<KbGraph>::edge_descriptor edge_descriptor;
	property_map<Kb::boost_graph_t, edge_weight_t>::type weight_map = get(edge_weight, g);
	prank::constant_property_map <edge_descriptor, float> cte_weight(1.0); // always return 1

	size_t N = num_vertices(g);

	if (N == ranks.size()) {
	  std::fill(ranks.begin(), ranks.end(), 0.0);
	} else {
	  vector<float>(N, 0.0).swap(ranks); // Initialize rank vector
	}
	vector<float> rank_tmp(N, 0.0);    // auxiliary rank vector

	if (glVars::prank::use_weight) {
	  if(coef_status != 2) {
		out_coefs.resize(N);
		fill(out_coefs.begin(), out_coefs.end(), 0.0);
		N_no_isolated = prank::init_out_coefs(g,  &out_coefs[0], weight_map);
		coef_status = 2;
	  }
	} else {
	  if(coef_status != 1) {
		out_coefs.resize(N);
		fill(out_coefs.begin(), out_coefs.end(), 0.0);
		N_no_isolated = prank::init_out_coefs(g, &out_coefs[0], cte_weight);
		coef_status = 1;
	  }
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
  // Debug

  ostream & Kb::dump_graph(std::ostream & o) const {
	o << "Sources: ";
	writeS(o, relsSource);
	o << endl;
	graph_traits<KbGraph>::vertex_iterator it, end;
	tie(it, end) = vertices(g);
	for(;it != end; ++it) {
	  o << get(vertex_name, g, *it);
	  graph_traits<KbGraph>::out_edge_iterator e, e_end;
	  tie(e, e_end) = out_edges(*it, g);
	  if (e != e_end)
		o << "\n";
	  for(; e != e_end; ++e) {
		o << "  ";
		vector<string> r = get_edge_reltypes(*e);
		writeV(o, r);
		o << " " << get(vertex_name, g, target(*e, g)) << "\n";
	  }
	}
	return o;
  }

  ////////////////////////////////////////////////////////////////////////////////
  // Streaming

  const size_t magic_id_v1 = 0x070201;
  const size_t magic_id = 0x080826;

  // read

  Kb_vertex_t read_vertex_from_stream_v1(ifstream & is,
										 KbGraph & g) {

	string name;

	read_atom_from_stream(is, name);
	Kb_vertex_t v = add_vertex(g);
	put(vertex_name, g, v, name);
	put(vertex_flags, g, v, 0);
	return v;
  }

  Kb_edge_t read_edge_from_stream_v1(ifstream & is,
									 KbGraph & g) {

	size_t sIdx;
	size_t tIdx;
	float w = 0.0;
	//size_t source;
	bool insertedP;
	Kb_edge_t e;

	read_atom_from_stream(is, sIdx);
	read_atom_from_stream(is, tIdx);
	read_atom_from_stream(is, w);
	//read_atom_from_stream(is, id);
	//read_atom_from_stream(is, source);
	tie(e, insertedP) = add_edge(sIdx, tIdx, g);
	assert(insertedP);
	put(edge_weight, g, e, w);
	//put(edge_source, g, e, source);

	return e;
  }

  Kb_vertex_t read_vertex_from_stream(ifstream & is,
									  KbGraph & g) {

	string name;
	string gloss;

	read_atom_from_stream(is, name);
	read_atom_from_stream(is, gloss);
	Kb_vertex_t v = add_vertex(g);
	put(vertex_name, g, v, name);
	put(vertex_flags, g, v, static_cast<unsigned char>(Kb::is_concept));
	return v;
  }

  Kb_edge_t read_edge_from_stream(ifstream & is,
								  KbGraph & g) {

	size_t sIdx;
	size_t tIdx;
	float w = 0.0;
	boost::uint32_t rtype;
	bool insertedP;
	Kb_edge_t e;

	read_atom_from_stream(is, sIdx);
	read_atom_from_stream(is, tIdx);
	read_atom_from_stream(is, w);
	read_atom_from_stream(is, rtype);
	//read_atom_from_stream(is, source);
	tie(e, insertedP) = add_edge(sIdx, tIdx, g);
	assert(insertedP);
	put(edge_weight, g, e, w);
	put(edge_rtype, g, e, rtype);

	return e;
  }

  void  Kb::read_from_stream (std::ifstream & is) {

	size_t vertex_n;
	size_t edge_n;
	size_t i;
	size_t id;

	std::map<std::string, int> relMap_aux;     // Obsolete map from relation name to relation id

	try {
	  coef_status = 0;
	  vector<float>().swap(static_ppv); // empty static rank vector
	  read_atom_from_stream(is, id);
	  if (id == magic_id_v1) {

		// Backward compatibility with binary v1 format

		read_set_from_stream(is, relsSource);
		read_map_from_stream(is, relMap_aux);

		read_map_from_stream(is, synsetMap);
		read_map_from_stream(is, wordMap);
		//read_map_from_stream(is, sourceMap);

		read_atom_from_stream(is, id);
		if(id != magic_id_v1) {
		  cerr << "Error: invalid id after reading maps" << endl;
		  exit(-1);
		}

		read_atom_from_stream(is, vertex_n);
		for(i=0; i<vertex_n; ++i) {
		  read_vertex_from_stream_v1(is, g);
		}

		read_atom_from_stream(is, id);
		if(id != magic_id_v1) {
		  cerr << "Error: invalid id after reading vertices" << endl;
		  exit(-1);
		}

		read_atom_from_stream(is, edge_n);
		for(i=0; i<edge_n; ++i) {
		  read_edge_from_stream_v1(is, g);
		}

		read_atom_from_stream(is, id);
		if(id != magic_id_v1) {
		  cerr << "Error: invalid id after reading edges" << endl;
		  exit(-1);
		}
		read_vector_from_stream(is, notes);
		if(id != magic_id_v1) {
		  cerr << "Error: invalid id (filename is a kbGraph?)" << endl;
		  exit(-1);
		}
	  } else {
		// Normal case
		read_set_from_stream(is, relsSource);
		read_vector_from_stream(is, rtypes);

		read_map_from_stream(is, synsetMap);
		read_map_from_stream(is, wordMap);

		read_atom_from_stream(is, id);
		if(id != magic_id) {
		  cerr << "Error: invalid id after reading maps" << endl;
		  exit(-1);
		}

		read_atom_from_stream(is, vertex_n);
		for(i=0; i<vertex_n; ++i) {
		  read_vertex_from_stream(is, g);
		}

		read_atom_from_stream(is, id);
		if(id != magic_id) {
		  cerr << "Error: invalid id after reading vertices" << endl;
		  exit(-1);
		}

		read_atom_from_stream(is, edge_n);
		for(i=0; i<edge_n; ++i) {
		  read_edge_from_stream(is, g);
		}

		read_atom_from_stream(is, id);
		if(id != magic_id) {
		  cerr << "Error: invalid id after reading edges" << endl;
		  exit(-1);
		}
		read_vector_from_stream(is, notes);
		if(id != magic_id) {
		  cerr << "Error: invalid id (filename is a kbGraph?)" << endl;
		  exit(-1);
		}
	  }
	} catch (...) {
	  throw runtime_error("Error when reading serialized graph (same platform used to compile the KB?)\n");
	}

	map<string, Kb_vertex_t>::iterator m_it(wordMap.begin());
	map<string, Kb_vertex_t>::iterator m_end(wordMap.end());
	for(; m_it != m_end; ++m_it) {
	  put(vertex_flags, g, m_it->second,
		  get(vertex_flags, g, m_it->second) || Kb::is_word);
	}
  }

  // write

  //
  // Auxiliary functions for removing isolated vertices
  //


  static size_t vdelta_isolated = numeric_limits<size_t>::max();

  static size_t get_vdeltas(const KbGraph & g,
							vector<size_t> & vdeltas) {

	size_t d = 0;

	graph_traits<KbGraph>::vertex_iterator vit, vend;
	tie(vit, vend) = vertices(g);
	for(;vit != vend; ++vit) {
	  if (out_degree(*vit, g) + in_degree(*vit, g) == 0) {
		// isolated vertex
		vdeltas[*vit] = vdelta_isolated;
		++d;
	  } else {
		vdeltas[*vit] = d;
	  }
	}
	return d;
  }

  static void map_update(const vector<size_t> & vdelta,
						 map<string, Kb_vertex_t> & theMap) {

	map<string, Kb_vertex_t>::iterator it = theMap.begin();
	map<string, Kb_vertex_t>::iterator end = theMap.end();

	while(it != end) {
	  if (vdelta[it->second] == vdelta_isolated) {
		// erase element
		theMap.erase(it++);
	  } else {
		// update vertex id
		it->second -= vdelta[it->second];
		++it;
	  }
	}
  }


  // write functions

  ofstream & write_vertex_to_stream(ofstream & o,
									const KbGraph & g,
									const vector<size_t> & vdelta,
									const Kb_vertex_t & v) {
	string name;

	if (vdelta[v] != vdelta_isolated) {
	  write_atom_to_stream(o, get(vertex_name, g, v));
	  write_atom_to_stream(o, get(vertex_gloss, g, v));
	}
	return o;
  }


  ofstream & write_edge_to_stream(ofstream & o,
								  const KbGraph & g,
								  const vector<size_t> & vdelta,
								  const Kb_edge_t & e) {

	size_t uIdx = get(vertex_index, g, source(e,g));
	uIdx -= vdelta[uIdx];
	size_t vIdx = get(vertex_index, g, target(e,g));
	vIdx -= vdelta[vIdx];

	float w = get(edge_weight, g, e);
	boost::uint32_t rtype = get(edge_rtype, g, e);

	o.write(reinterpret_cast<const char *>(&uIdx), sizeof(uIdx));
	o.write(reinterpret_cast<const char *>(&vIdx), sizeof(vIdx));
	o.write(reinterpret_cast<const char *>(&w), sizeof(w));
	o.write(reinterpret_cast<const char *>(&rtype), sizeof(rtype));
	return o;
  }

  ofstream & Kb::write_to_stream(ofstream & o) {

	// First remove isolated vertices and
	// - get delta vector
	// - remove from map

	// - get deltas
	vector<size_t> vdelta(num_vertices(g), 0);
	size_t visol_size = get_vdeltas(g, vdelta);

	// - update the maps

	if (visol_size) {
	  map_update(vdelta, synsetMap);
	  map_update(vdelta, wordMap);
	}

	// Write maps

	write_atom_to_stream(o, magic_id);

	write_vector_to_stream(o, relsSource);
	write_vector_to_stream(o, rtypes);

	write_map_to_stream(o, synsetMap);
	write_map_to_stream(o, wordMap);
	//write_map_to_stream(o, sourceMap);

	write_atom_to_stream(o, magic_id);

	size_t vertex_n = num_vertices(g) - visol_size;

	write_atom_to_stream(o, vertex_n);
	graph_traits<KbGraph>::vertex_iterator v_it, v_end;

	tie(v_it, v_end) = vertices(g);
	for(; v_it != v_end; ++v_it) {
	  write_vertex_to_stream(o, g, vdelta, *v_it);
	}

	write_atom_to_stream(o, magic_id);

	size_t edge_n = num_edges(g);

	write_atom_to_stream(o, edge_n);
	graph_traits<KbGraph>::edge_iterator e_it, e_end;

	tie(e_it, e_end) = edges(g);
	for(; e_it != e_end; ++e_it) {
	  write_edge_to_stream(o, g, vdelta, *e_it);
	}
	write_atom_to_stream(o, magic_id);
	if(notes.size()) write_vector_to_stream(o, notes);
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

  ofstream & write_to_textstream(const KbGraph & g, ofstream & o) {

	graph_traits<KbGraph>::edge_iterator e_it, e_end;

	tie(e_it, e_end) = edges(g);
	for(; e_it != e_end; ++e_it) {
	  o << "u:" << get(vertex_name, g, source(*e_it, g)) << " ";
	  o << "v:" << get(vertex_name, g, target(*e_it, g)) << " d:1\n";
	}
	return o;
  }

  void Kb::write_to_textfile (const string & fName) {

	ofstream fo(fName.c_str(),  ofstream::out);
	if (!fo) {
	  cerr << "Error: can't create" << fName << endl;
	  exit(-1);
	}
	write_to_textstream(g, fo);
  }


}
