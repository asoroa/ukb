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
#include <boost/pending/integer_range.hpp>
#include <boost/graph/graph_utility.hpp> // for boost::make_list

// dijkstra

#include <boost/graph/dijkstra_shortest_paths.hpp>

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

  void Kb::create_from_binfile(const std::string & fname) {

	if (p_instance) return;
	Kb *tenp = create();

	ifstream fi(fname.c_str(), ifstream::binary|ifstream::in);
	if (!fi) {
	  cerr << "Error: can't open " << fname << endl;
	  exit(-1);
	}
	tenp->read_from_stream(fi);
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

	vector<Kb_vertex_t>(num_vertices(g)).swap(parents);  // reset parents

	breadth_first_search(g,
						 src,
						 boost::visitor(boost::make_bfs_visitor
										(boost::make_list(kb_bfs_init(&parents[0]),
														  kb_bfs_pred(&parents[0])))));
	return true;
  }


  bool Kb::dijkstra (Kb_vertex_t src,
					  std::vector<Kb_vertex_t> & parents) const {

	vector<Kb_vertex_t>(num_vertices(g)).swap(parents);  // reset parents
	vector<float> dist(num_vertices(g));

	dijkstra_shortest_paths(g,
							src,
							predecessor_map(&parents[0]).
							distance_map(&dist[0]));
	return true;
  }

  ////////////////////////////////////////////////////////////////////////////////
  // strings <-> vertex_id

  pair<Kb_vertex_t, bool> Kb::get_vertex_by_name(const std::string & str) const {
	map<string, Kb_vertex_t>::const_iterator it;

	it = synsetMap.find(str);
	if (it != synsetMap.end()) return make_pair(it->second, true);

	// is it a word ?
	it = wordMap.find(str);
	if (it == wordMap.end()) return make_pair(Kb_vertex_t(), false);
	return make_pair(it->second, true);
  }

  Kb_vertex_t Kb::InsertNode(const string & name, unsigned char flags) {
	coef_status = 0; // reset out degree coefficients
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
	  it->second = InsertNode(str, 0);
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

	//if (w != 1.0) ++w; // minimum weight is 1
	tie(e, existsP) = edge(u, v, g);
	if(!existsP) {
	  coef_status = 0; // reset out degree coefficients
	  e = add_edge(u, v, g).first;
	  put(edge_weight, g, e, w);
	  put(edge_rtype, g, e, static_cast<boost::uint32_t>(0));
	}
	return e;
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
	normalize_pvector(outranks);
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
	  line_number++;
	  char_separator<char> sep(" ");
	  tokenizer<char_separator<char> > tok(line, sep);
	  copy(tok.begin(), tok.end(), back_inserter(fields));
	  if (fields.size() == 0) continue; // blank line
	  if (fields.size() < 4) {
		throw runtime_error("read_kb error: Bad line: " + lexical_cast<string>(line_number));
	  }

	  if (glVars::kb::filter_src) {
		if (src_allowed.find(fields[3]) == srel_end) continue; // Skip this relation
	  }

	  kb->add_relSource(fields[3]);

	  // last element says if relation is directed
	  bool directed = (fields.size() > 4 && lexical_cast<int>(fields[4]) != 0);

	  Kb_vertex_t u = kb->find_or_insert_synset(fields[0]);
	  Kb_vertex_t v = kb->find_or_insert_synset(fields[1]);

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

	char_separator<char> sep(" ");
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
		res.directed = lexical_cast<bool>(val);
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

  void read_kb(ifstream & kbFile,
			   const set<string> & src_allowed,
			   Kb *kb) {
	string line;
	int line_number = 0;

	set<string>::const_iterator srel_end = src_allowed.end();
	try {
	  while(kbFile) {
		vector<string> fields;
		std::getline(kbFile, line, '\n');
		line_number++;
		rel_parse f;
		if (!parse_line(line, f)) continue;

		if (glVars::kb::filter_src) {
		  if (src_allowed.find(f.src) == srel_end) continue; // Skip this relation
		}
		kb->add_relSource(f.src);

		Kb_vertex_t u = kb->find_or_insert_synset(f.u);
		Kb_vertex_t v = kb->find_or_insert_synset(f.v);
		float w = f.w ? f.w : 1.0;
		// add edge
		Kb_edge_t e1 = kb->find_or_insert_edge(u, v, w);
		// relation type
		if (glVars::kb::keep_reltypes && f.rtype.size()) {
		  kb->edge_add_reltype(e1, f.rtype);
		}
		if (!f.directed) {
		  Kb_edge_t e2 = kb->find_or_insert_edge(v, u, w);
		  if(glVars::kb::keep_reltypes && f.rtype.size()) {
			string aux = f.irtype.size() ? f.irtype : f.rtype;
			if(aux.size()) {
			  kb->edge_add_reltype(e2, aux);
			}
		  }
		}
	  }
	} catch (std::exception & e) {
	  throw std::runtime_error(string(e.what()) + " in line " + lexical_cast<string>(line_number));
	}
  }

  void Kb::read_from_txt(const string & synsFileName,
						 const set<string> & src_allowed) {

	// optimize IO
	std::ios::sync_with_stdio(false);

	std::ifstream syns_file(synsFileName.c_str(), ofstream::in);
	if (!syns_file) {
	  throw runtime_error("Kb::read_from_txt error: Can't open " + synsFileName);
	}
	if(glVars::kb::v1_kb) {
	  read_kb_v1(syns_file, src_allowed, this);
	} else {
	  read_kb(syns_file, src_allowed, this);
	}
  }

  ////////////////////////////////////////////////7
  // public function

  void Kb::add_from_txt(const std::string & synsFileName,
						const set<string> & src_allowed) {

	Kb::instance().read_from_txt(synsFileName, src_allowed);
  }

  void Kb::display_info(std::ostream & o) const {

	o << "Relation sources: ";
	writeS(o, relsSource);
	if (notes.size()) {
	  o << "\nNotes: ";
	  writeV(o, notes);
	}
	size_t edge_n = num_edges(g);
	if (edge_n & 2) edge_n++;

	o << "\n" << num_vertices(g) << " vertices and " << edge_n/2 << " edges" << endl;
	if (rtypes.size()) {
	  o << "Relations:";
	  writeV(o, rtypes);
	  o << endl;
	}
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
	  char synpos = syns.get_pos(i);
	  if (glVars::input::filter_pos && synpos) {
		++sit; // '#' char
		*sit = syns.get_pos(i); // the pos
	  }

	  Kb_vertex_t u;
	  tie(u, auxP) = kb.get_vertex_by_name(syns.get_entry(i));
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
				   map<string, vector<Syn_elem> > & wPos2Syns,
				   bool use_weights) {

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
		if (use_weights) {
		  kb.find_or_insert_edge(wpos_v, syns_it->first, syns_it->second);
		} else {
		  kb.find_or_insert_edge(wpos_v, syns_it->first, 1.0);
		}
	  }
	}
  }

  void insert_word(const string & word, bool use_w) {

	Kb & kb = Kb::instance();

	WDict_entries syns = WDict::instance().get_entries(word);

	if (!syns.size()) return;

	Kb_vertex_t u = kb.find_or_insert_word(word);

	for(size_t i = 0; i < syns.size(); ++i) {
	  bool auxP;
	  Kb_vertex_t v;
	  tie(v, auxP) = kb.get_vertex_by_name(syns.get_entry(i));
	  if(!auxP) {
		if (glVars::debug::warning)
		  cerr << "W:Kb::add_tokens: warning: " << syns.get_entry(i) << " is not in KB.\n";
		continue;
	  }
	  // (directed) link word -> concept
	  kb.find_or_insert_edge(u, v, 1.0);
	}
  }

  void Kb::add_dictionary(bool with_weight) {

	WDict & w2syn = WDict::instance();

	vector<string>::const_iterator word_it = w2syn.get_wordlist().begin();
	vector<string>::const_iterator word_end = w2syn.get_wordlist().end();

	for(; word_it != word_end; ++word_it) {
	  add_token(*word_it, with_weight);
	}
  }


  void Kb::add_token(const string & token, bool with_weight) {

	vector<string> wPosV;
	map<string, vector<Syn_elem> > wPos2Syns;

	if (glVars::input::filter_pos) {
	  // Create w2wPos and wPos2Syns maps

	  create_w2wpos_maps(token, wPosV, wPos2Syns);

	  // Add vertices and link them in the KB
	  insert_wpos(token, wPosV, wPos2Syns, with_weight);
	} else {
	  insert_word(token, with_weight);
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

  template<typename G, typename ppvMap_t, typename wMap_t, typename map1_t, typename map2_t>
  void pageRank_dispatch(G & g,
						 vector<Kb_vertex_t> & V,
						 ppvMap_t ppv_map, 
						 wMap_t & wmap,
						 map1_t rank,
						 map2_t rank_tmp,
						 const vector<float> & out_coefs) {

	if (glVars::prank::threshold != 0.0) 
	  prank::do_pageRank_l1(g, V, 
							&ppv_map[0], wmap, 
							&rank[0], &rank_tmp[0],
							glVars::prank::threshold,
							out_coefs);
	else
	  prank::do_pageRank(g, V, 
						 &ppv_map[0], wmap, 
						 &rank[0], &rank_tmp[0],
						 glVars::prank::num_iterations,
						 out_coefs);
  }
						 

  void Kb::pageRank_ppv(const vector<float> & ppv_map,
						 vector<float> & ranks,
						 bool use_weight) {

	size_t N = num_vertices(g);
	vector<float>(N, 0.0).swap(ranks); // Initialize rank vector
	vector<float> rank_tmp(N, 0.0);    // auxiliary rank vector

	// ugly ugly hack @@CHANGE ME !!!

	vector<Kb_vertex_t> V(N);

	graph_traits<KbGraph>::vertex_iterator it, end;
	tie(it, end) = vertices(g);
	copy(it, end, V.begin());

	if (use_weight) {
	  property_map<Kb::boost_graph_t, edge_weight_t>::type weight_map = get(edge_weight, g);
	  if(coef_status != 2) {
		vector<float>(num_vertices(g), 0.0f).swap(out_coefs);
		prank::init_out_coefs(g, V, &out_coefs[0], weight_map);
		coef_status = 2;
	  }
	  pageRank_dispatch(g, V, &ppv_map[0],
						weight_map, &ranks[0], &rank_tmp[0],
						out_coefs);
// 	  prank::do_pageRank(g, V, &ppv_map[0],
// 						 weight_map, &ranks[0], &rank_tmp[0],
// 						 glVars::prank::num_iterations,
// 						 out_coefs);
	} else {
	  typedef graph_traits<KbGraph>::edge_descriptor edge_descriptor;
	  prank::constant_property_map <edge_descriptor, float> cte_weight(1); // always return 1
	  if(coef_status != 1) {
		vector<float>(num_vertices(g), 0.0f).swap(out_coefs);
		prank::init_out_coefs(g, V, &out_coefs[0], cte_weight);
		coef_status = 1;
	  }
	  pageRank_dispatch(g, V, &ppv_map[0],
						cte_weight, &ranks[0], &rank_tmp[0],
						out_coefs);
// 	  prank::do_pageRank(g, V, &ppv_map[0],
// 						 cte_weight, &ranks[0], &rank_tmp[0],
// 						 glVars::prank::num_iterations,
// 						 out_coefs);
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
	put(vertex_flags, g, v, 0);
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

	coef_status = 0;
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

	map<string, Kb_vertex_t>::iterator m_it(wordMap.begin());
	map<string, Kb_vertex_t>::iterator m_end(wordMap.end());
	for(; m_it != m_end; ++m_it) {
	  put(vertex_flags, g, m_it->second,
		  get(vertex_flags, g, m_it->second) || Kb::is_word);
	}
  }

  // write

  ofstream & write_vertex_to_stream(ofstream & o,
									const KbGraph & g,
									const Kb_vertex_t & v) {
	string name;

	write_atom_to_stream(o, get(vertex_name, g, v));
	write_atom_to_stream(o, get(vertex_gloss, g, v));
	return o;
  }

  ofstream & write_edge_to_stream(ofstream & o,
								  const KbGraph & g,
								  const Kb_edge_t & e) {

	size_t uIdx = get(vertex_index, g, source(e,g));
	size_t vIdx = get(vertex_index, g, target(e,g));
	float w = get(edge_weight, g, e);
	boost::uint32_t rtype = get(edge_rtype, g, e);

	o.write(reinterpret_cast<const char *>(&uIdx), sizeof(uIdx));
	o.write(reinterpret_cast<const char *>(&vIdx), sizeof(vIdx));
	o.write(reinterpret_cast<const char *>(&w), sizeof(w));
	o.write(reinterpret_cast<const char *>(&rtype), sizeof(rtype));
	return o;
  }

  ofstream & Kb::write_to_stream(ofstream & o) const {

	// First write maps

	write_atom_to_stream(o, magic_id);

	write_vector_to_stream(o, relsSource);
	write_vector_to_stream(o, rtypes);

	write_map_to_stream(o, synsetMap);
	write_map_to_stream(o, wordMap);
	//write_map_to_stream(o, sourceMap);

	write_atom_to_stream(o, magic_id);

	size_t vertex_n = num_vertices(g);

	write_atom_to_stream(o, vertex_n);
	graph_traits<KbGraph>::vertex_iterator v_it, v_end;

	tie(v_it, v_end) = vertices(g);
	for(; v_it != v_end; ++v_it) {
	  write_vertex_to_stream(o, g, *v_it);
	}

	write_atom_to_stream(o, magic_id);

	size_t edge_n = num_edges(g);

	write_atom_to_stream(o, edge_n);
	graph_traits<KbGraph>::edge_iterator e_it, e_end;

	tie(e_it, e_end) = edges(g);
	for(; e_it != e_end; ++e_it) {
	  write_edge_to_stream(o, g, *e_it);
	}
	write_atom_to_stream(o, magic_id);
	if(notes.size()) write_vector_to_stream(o, notes);
	return o;
  }

  void Kb::write_to_binfile (const string & fName) const {

	ofstream fo(fName.c_str(),  ofstream::binary|ofstream::out);
	if (!fo) {
	  cerr << "Error: can't create" << fName << endl;
	  exit(-1);
	}
	write_to_stream(fo);
  }
}
