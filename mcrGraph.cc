#include "mcrGraph.h"
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

using namespace std;
using namespace boost;



////////////////////////////////////////////////////////////////////////////////
// Auxiliary functions
template<class Map, class InvMap>
void create_inverse_map(const Map & source, InvMap & target) {
 
  typedef typename InvMap::value_type value_type_inv;

  typename Map::const_iterator it = source.begin();
  typename Map::const_iterator it_end = source.end();

  for(; it != it_end; ++it) {
    target.insert(value_type_inv(it->second,it->first));
  }
}

////////////////////////////////////////////////////////////////////////////////
// Class Mcr


////////////////////////////////////////////////////////////////////////////////
// Singleton stuff

Mcr* Mcr::p_instance = 0;

Mcr *Mcr::create() {

  static Mcr theMcr;
  return &theMcr;
}

Mcr & Mcr::instance() {
  if (!p_instance) {
    throw runtime_error("MCR not initialized");
  }
  return *p_instance;
}

void Mcr::create_from_txt(const string & synsFileName,
			  const std::set<std::string> & rels_source,
			  string relFileName) {
  if (p_instance) return;
  Mcr *tenp = create();

  tenp->relsSource = rels_source;
  tenp->read_from_txt(relFileName, synsFileName);
  p_instance = tenp;
}

void Mcr::create_from_binfile(const std::string & fname) {

  if (p_instance) return;
  Mcr *tenp = create();

  ifstream fi(fname.c_str(), ifstream::binary|ifstream::in);
  if (!fi) {
    cerr << "Error: can't open " << fname << endl;
    exit(-1);
  }
  tenp->read_from_stream(fi);
  p_instance = tenp;
}

////////////////////////////////////////////////////////////////////////////////


void Mcr::add_comment(const string & str) {
  notes.push_back(str);
}

const vector<string> & Mcr::get_comments() const {
  return notes;
}

////////////////////////////////////////////////////////////////////////////////
// bfs

struct mcr_bfs_init:public base_visitor<mcr_bfs_init> {
public:
  mcr_bfs_init(Mcr_vertex_t *v):m_v(v) { }
  typedef on_initialize_vertex event_filter;
  inline void operator()(Mcr_vertex_t u, const McrGraph & g)
  {
    m_v[u] = u;
  }
  Mcr_vertex_t *m_v;
};

struct mcr_bfs_pred:public base_visitor<mcr_bfs_pred> {
public:
  mcr_bfs_pred(Mcr_vertex_t *v):m_v(v) { }
  typedef on_tree_edge event_filter;
  inline void operator()(Mcr_edge_t e, const McrGraph & g) {
    m_v[target(e, g)] = source(e,g);
  }
  Mcr_vertex_t *m_v;
};


bool Mcr::bfs (Mcr_vertex_t src, 
	       std::vector<Mcr_vertex_t> & parents) const {

  vector<Mcr_vertex_t>(num_vertices(g)).swap(parents);  // reset parents
  
  breadth_first_search(g, 
		       src,
		       boost::visitor(boost::make_bfs_visitor
				      (boost::make_list(mcr_bfs_init(&parents[0]),
							mcr_bfs_pred(&parents[0])))));
  return true;
}


bool Mcr::dijkstra (Mcr_vertex_t src, 
		    std::vector<Mcr_vertex_t> & parents) const {

  vector<Mcr_vertex_t>(num_vertices(g)).swap(parents);  // reset parents
  vector<float> dist(num_vertices(g));
  
  dijkstra_shortest_paths(g, 
			  src,
			  predecessor_map(&parents[0]).
			  distance_map(&dist[0]));
  return true;
}

////////////////////////////////////////////////////////////////////////////////
// strings <-> vertex_id

pair<Mcr_vertex_t, bool> Mcr::getVertexByName(const std::string & str) const {
  map<string, Mcr_vertex_t>::const_iterator it;

  it = synsetMap.find(str);
  if (it != synsetMap.end()) return make_pair(it->second, true);

  // is it a word ?
  it = wordMap.find(str);
  if (it == wordMap.end()) return make_pair(Mcr_vertex_t(), false);
  return make_pair(it->second, true);
}

Mcr_vertex_t Mcr::InsertNode(const string & name, unsigned char flags) {
  coef_status = 0; // reset out degree coefficients
  Mcr_vertex_t u = add_vertex(g);
  put(vertex_name, g, u, name);
  put(vertex_flags, g, u, flags);
  return u;
}

Mcr_vertex_t Mcr::findOrInsertSynset(const string & str) {
  map<string, Mcr_vertex_t>::iterator it;
  bool insertedP;
  tie(it, insertedP) = synsetMap.insert(make_pair(str, Mcr_vertex_t()));
  if(insertedP) {
    // new vertex
    it->second = InsertNode(str, 0);
  }
  return it->second;
}

Mcr_vertex_t Mcr::findOrInsertWord(const string & str) {
  map<string, Mcr_vertex_t>::iterator it;
  bool insertedP;
  tie(it, insertedP) = wordMap.insert(make_pair(str, Mcr_vertex_t()));
  if(insertedP) {
    // new vertex
    it->second = InsertNode(str, Mcr::is_word);
  }
  return it->second;
}

Mcr_edge_t Mcr::findOrInsertEdge(Mcr_vertex_t u, Mcr_vertex_t v,
				 float w) {

  Mcr_edge_t e;
  bool existsP;

  //if (w != 1.0) ++w; // minimum weight is 1
  tie(e, existsP) = edge(u, v, g);
  if(!existsP) {
    coef_status = 0; // reset out degree coefficients
    e = add_edge(u, v, g).first;
    put(edge_weight, g, e, w);
  }
  return e;
}

vector<string>::size_type get_reltype_idx(const string & rel,
										vector<string> rtypes) {

  vector<string>::iterator it = rtypes.begin();
  vector<string>::iterator end = rtypes.end();
  vector<string>::size_type idx = 0;

  for(;it != end; ++it) {
	if (*it == rel) break;
	++it;
	++idx;
  }
  if (it == end) {
	// new relation type
	rtypes.push_back(rel);
  }
  if (idx > 32) {
	throw runtime_error("get_rtype_idx error: too many relation types !");
  }
  return idx;
}

void Mcr::edge_add_reltype(Mcr_edge_t e, const string & rel) {
  boost::uint32_t m = get(edge_rtype, g, e);
  vector<string>::size_type idx = get_rtype_idx(rel, rtypes);
  m &= (1L << idx);
  put(edge_rtype, g, e, m);
}

////////////////////////////////////////////////////////////////////////////////
// Query

bool Mcr::vertexIsSynset(Mcr_vertex_t u) const {
  return !vertexIsWord(u);
}

bool Mcr::vertexIsWord(Mcr_vertex_t u) const {
  return (get(vertex_flags, g, u) & Mcr::is_word);
}

////////////////////////////////////////////////////////////////////////////////
// Random

Mcr_vertex_t Mcr::getRandomVertex() const {

  int r = g_randTarget(num_vertices(g));

  return r;
}

////////////////////////////////////////////////////////////////////////////////
// read from textfile and create graph


// Read relations name, id and inverse relations

void read_reltypes(ifstream & relFile,
				   map<string, int> & relMap,
				   map<int, string> & relMapInv) {

  string line;
  int line_number = 0;

  while(relFile) {
    // fields: id name code inverse_rel_id type NULL
    vector<string> fields;
    std::getline(relFile, line, '\n');
    line_number++;
    char_separator<char> sep(" ");
    tokenizer<char_separator<char> > tok(line, sep);
    copy(tok.begin(), tok.end(), back_inserter(fields));
    if (fields.size() == 0) continue; // blank line
    if (fields.size() < 2) {
      cerr << "read_reltypes error. Bad line: " << line_number << endl;
      exit(-1);
    }
    int relId = lexical_cast<int>(fields[0]);
    relMap[fields[1]] = relId;
  }
  create_inverse_map(relMap, relMapInv);
}

// Read the actual MCR file


void read_mcr_v1(ifstream & mcrFile, 
				 const set<string> & rels_source,
				 Mcr * mcr) {
  string line;
  int line_number = 0;

  set<string>::const_iterator srel_end = rels_source.end();
  while(mcrFile) {
    vector<string> fields;
    std::getline(mcrFile, line, '\n');
    line_number++;
    char_separator<char> sep(" ");
    tokenizer<char_separator<char> > tok(line, sep);
    copy(tok.begin(), tok.end(), back_inserter(fields));
    if (fields.size() == 0) continue; // blank line
    if (fields.size() < 4) {
      throw runtime_error("read_mcr error: Bad line: " + line_number);
    }

    if (rels_source.find(fields[3]) == srel_end) continue; // Skip this relation

    // last element says if relation is directed
    bool directed = (fields.size() > 4 && lexical_cast<int>(fields[4]) != 0);

    Mcr_vertex_t u = mcr->findOrInsertSynset(fields[0]);
    Mcr_vertex_t v = mcr->findOrInsertSynset(fields[1]);

    // add edge
    mcr->findOrInsertEdge(u, v, 1.0);
    if (!directed)
      mcr->findOrInsertEdge(v, u, 1.0);
  }
}


// Line format:
//
// u:synset v:synset t:rel i:rel s:source d:directed
//
// u: source vertex. Mandatory.
// v: target vertex. Mandatory.
// t: relation type (hyperonym, meronym, etc) of edge u->v. Optional.
// i: (inverse) relation type of edge v->u (hyponym, etc). Optional.
// s: source of relation (wn30, mcr17, etc). Optional.
// d: wether the relation is directed. Optional, default is undirected.


struct rel_parse {
  string u;
  string v;
  string rtype;
  string irtype;
  string src;
  bool directed;
};

rel_parse parse_line(const string & line) {

  rel_parse res = {"","","","","",false};

  char_separator<char> sep(" ");
  tokenizer<char_separator<char> > tok(line, sep);
  tokenizer<char_separator<char> >::iterator it = tok.begin();
  tokenizer<char_separator<char> >::iterator end = tok.end();
  for(;it != end; ++it) {
	// parse field
	vector<string> field;
	char_separator<char> sep_field(":");
	tokenizer<char_separator<char> > tok_field(*it, sep_field);
	copy(tok_field.begin(), tok_field.end(), back_inserter(field));
	if (field.size() != 2) {
      throw runtime_error("parse_line error. Malformed line.");
	}
	char f = field[0].at(0);
	switch (f) {
	case 'u':
	  res.u = field[1];
	  break;
	case 'v':
	  res.v = field[1];
	  break;
	case 't':
	  res.rtype = field[1];
	  break;
	case 'i':
	  res.irtype = field[1];
	  break;
	case 's':
	  res.src = field[1];
	  break;
	case 'd':
	  res.directed = lexical_cast<bool>(field[1]);
	  break;
	}
  }
  if (!res.u.size()) throw runtime_error("parse_line error. No source vertex.");
  if (!res.v.size()) throw runtime_error("parse_line error. No target vertex.");
  return res;
}

void read_mcr(ifstream & mcrFile, 
			  const set<string> & rels_source,
			  Mcr *mcr) {
  string line;
  int line_number = 0;

  set<string>::const_iterator srel_end = rels_source.end();
  while(mcrFile) {
    vector<string> fields;
    std::getline(mcrFile, line, '\n');
    line_number++;
	try {
	  rel_parse f = parse_line(line);
	  if (f.src.size() && rels_source.find(f.src) == srel_end) continue; // Skip this relation
	  Mcr_vertex_t u = mcr->findOrInsertSynset(f.u);
	  Mcr_vertex_t v = mcr->findOrInsertSynset(f.v);

	  // add edge
	  Mcr_edge_t e1 = mcr->findOrInsertEdge(u, v, 1.0);
	  // relation type
	  if (glVars::kb::keep_reltypes && f.rtype.size()) {
		mcr->edge_add_reltype(e1, f.rtype);
	  }
	  if (!f.directed) {
		Mcr_edge_t e2 = mcr->findOrInsertEdge(v, u, 1.0);
		if(glVars::kb::keep_reltypes) {
		  string aux = f.irtype.size() ? f.irtype : f.rtype;
		  if(aux.size()) {
			mcr->edge_add_reltype(e2, aux);
		  }
		}
	  }	  
	} catch (std::runtime_error & e) {
	  throw std::runtime_error(string(e.what()) + " in line " + lexical_cast<string>(line_number));
	}
  }
}

void Mcr::read_from_txt(const string & synsFileName,
						const string & reltypeFileName) {

  // optimize IO
  std::ios::sync_with_stdio(false);

  if (reltypeFileName.length()) {
    std::ifstream rt_file(reltypeFileName.c_str(), ofstream::in);
    if (!rt_file) {
      throw runtime_error("Mcr::read_from_txt error: Can't open " + reltypeFileName);
    }
    read_reltypes(rt_file, relMap, relMapInv);
  }

  std::ifstream syns_file(synsFileName.c_str(), ofstream::in);
  if (!syns_file) {
    throw runtime_error("Mcr::read_from_txt error: Can't open " + synsFileName);
  }
  if(glVars::kb::v1_kb) {
    read_mcr_v1(syns_file, relsSource, this);
  } else {
    read_mcr(syns_file, relsSource, this);  
  }
}

////////////////////////////////////////////////7
// public function

void Mcr::add_from_txt(const std::string & synsFileName) {

  Mcr::instance().read_from_txt(synsFileName, "");
}

void Mcr::display_info(std::ostream & o) const {

  o << "Relation sources: ";
  writeS(o, relsSource);
  if (notes.size()) {
    o << "\nNotes: ";
    writeV(o, notes);
  }
  size_t edge_n = num_edges(g);
  if (edge_n & 2) edge_n++; 

  o << "\n" << num_vertices(g) << " vertices and " << edge_n/2 << " edges" << endl;

}

////////////////////////////////////////////////////////////////////////////////
// Add token vertices and link them to synsets


typedef pair<Mcr_vertex_t, float> Syn_elem;

void create_w2wpos_maps(const string & word,
			vector<string> & wPosV,
			map<string, vector<Syn_elem> > & wPos2Syns) {

  bool auxP;
  const Mcr & mcr = Mcr::instance();

  WDict_entries syns = WDict::instance().get_entries(word);

  for(size_t i = 0; i < syns.size(); ++i) {

    string wpos(word.size() + 2, '#');
    string::iterator sit = copy(word.begin(), word.end(), wpos.begin());
    ++sit; // '#' char
    *sit = syns.get_pos(i); // the pos

    Mcr_vertex_t u;
    tie(u, auxP) = mcr.getVertexByName(syns.get_entry(i));
    if(!auxP) {
      if (glVars::debug::warning) 
	cerr << "W:Mcr::add_tokens: warning: " << syns.get_entry(i) << " is not in MCR.\n";
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

  Mcr & mcr = Mcr::instance();

  // insert word
  Mcr_vertex_t word_v = mcr.findOrInsertWord(word); 

  vector<string>::iterator wpos_str_it = wPosV.begin();
  vector<string>::iterator wpos_str_end = wPosV.end();
  for (; wpos_str_it != wpos_str_end; ++wpos_str_it) {
    //insert word#pos
    Mcr_vertex_t wpos_v = mcr.findOrInsertWord(*wpos_str_it);
    // link word to word#pos
    mcr.findOrInsertEdge(word_v, wpos_v, 1.0);
    // link word#pos to synsets
    vector<Syn_elem>::const_iterator syns_it = wPos2Syns[*wpos_str_it].begin();
    vector<Syn_elem>::const_iterator syns_end = wPos2Syns[*wpos_str_it].end();
    for(;syns_it != syns_end; ++syns_it) {
      if (use_weights) {
	mcr.findOrInsertEdge(wpos_v, syns_it->first, syns_it->second);
      } else {
	mcr.findOrInsertEdge(wpos_v, syns_it->first, 1.0);
      }
    }
  }
}

void Mcr::add_dictionary(bool with_weight) {

  WDict & w2syn = WDict::instance();

  vector<string>::const_iterator word_it = w2syn.get_wordlist().begin();
  vector<string>::const_iterator word_end = w2syn.get_wordlist().end();

  for(; word_it != word_end; ++word_it) {
    add_token(*word_it, with_weight);
  }
}

void Mcr::add_token(const string & token, bool with_weight) {

  vector<string> wPosV;
  map<string, vector<Syn_elem> > wPos2Syns;

    // Create w2wPos and wPos2Syns maps

  create_w2wpos_maps(token, wPosV, wPos2Syns);

    // Add vertices and link them in the MCR
  insert_wpos(token, wPosV, wPos2Syns, with_weight);


}

void Mcr::ppv_weights(const vector<float> & ppv) {

  graph_traits<McrGraph>::edge_iterator it, end;

  tie(it, end) = edges(g);
  for(; it != end; ++it) {
    put(edge_weight, g, *it,
	ppv[target(*it, g)]);
  }
}

////////////////////////////////////////////////////////////////////////////////
// PageRank in MCR


// PPV version

void Mcr::pageRank_ppv(const vector<float> & ppv_map,
		       vector<float> & ranks,
		       bool use_weight) {

  size_t N = num_vertices(g);
  vector<float>(N, 0.0).swap(ranks); // Initialize rank vector
  vector<float> rank_tmp(N, 0.0);    // auxiliary rank vector

  // ugly ugly hack @@CHANGE ME !!!

  vector<Mcr_vertex_t> V(N);

  graph_traits<McrGraph>::vertex_iterator it, end;
  tie(it, end) = vertices(g);
  copy(it, end, V.begin());

  if (use_weight) {
    property_map<Mcr::boost_graph_t, edge_weight_t>::type weight_map = get(edge_weight, g);
    if(coef_status != 2) {
      vector<float>(num_vertices(g), 0.0f).swap(out_coefs);
      prank::init_out_coefs(g, V, &out_coefs[0], weight_map);
      coef_status = 2;
    }
    prank::do_pageRank(g, V, &ppv_map[0],
		       weight_map, &ranks[0], &rank_tmp[0], 
		       glVars::prank::num_iterations,
		       out_coefs);
  } else {
    typedef graph_traits<McrGraph>::edge_descriptor edge_descriptor;
    prank::constant_property_map <edge_descriptor, float> cte_weight(1); // always return 1
    if(coef_status != 1) {
      vector<float>(num_vertices(g), 0.0f).swap(out_coefs);
      prank::init_out_coefs(g, V, &out_coefs[0], cte_weight);
      coef_status = 1;
    }
    prank::do_pageRank(g, V, &ppv_map[0],
		       cte_weight, &ranks[0], &rank_tmp[0], 
		       glVars::prank::num_iterations,
		       out_coefs);
  }
}

////////////////////////////////////////////////////////////////////////////////
// Streaming

const size_t magic_id = 0x070201;

// read

Mcr_vertex_t read_vertex_from_stream(ifstream & is, 
				     McrGraph & g) {

  string name;

  read_atom_from_stream(is, name);
  Mcr_vertex_t v = add_vertex(g);
  put(vertex_name, g, v, name);
  put(vertex_flags, g, v, 0);
  return v;
}

Mcr_edge_t read_edge_from_stream(ifstream & is, 
				 McrGraph & g) {

  size_t sIdx;
  size_t tIdx; 
  float w = 0.0;
  //size_t source;
  bool insertedP;
  Mcr_edge_t e;

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

void  Mcr::read_from_stream (std::ifstream & is) {

  size_t vertex_n;
  size_t edge_n;
  size_t i;
  size_t id;

  coef_status = 0;
  read_atom_from_stream(is, id);
  if(id != magic_id) {
    cerr << "Error: invalid id (filename is a mcrGraph?)" << endl;
    exit(-1);
  }

  read_set_from_stream(is, relsSource);
  read_map_from_stream(is, relMap);
  create_inverse_map(relMap, relMapInv);

  read_map_from_stream(is, synsetMap);
  read_map_from_stream(is, wordMap);
  //read_map_from_stream(is, sourceMap);

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

  map<string, Mcr_vertex_t>::iterator m_it(wordMap.begin());
  map<string, Mcr_vertex_t>::iterator m_end(wordMap.end());
  for(; m_it != m_end; ++m_it) {
    put(vertex_flags, g, m_it->second, 
	get(vertex_flags, g, m_it->second) || Mcr::is_word);
  }

//   graph_traits<McrGraph>::vertex_iterator v_it, v_end;
//   tie(v_it, v_end) = vertices(g);
//   for(; v_it != v_end; ++v_it) {
//     synsetMap[get(vertex_name, g, *v_it)] = *v_it;
//   }
}

// write

ofstream & write_vertex_to_stream(ofstream & o,
				  const McrGraph & g,
				  const Mcr_vertex_t & v) {
  string name;

  write_atom_to_stream(o, get(vertex_name, g, v));
  return o;
}

ofstream & write_edge_to_stream(ofstream & o,
				const McrGraph & g,
				const Mcr_edge_t & e) {

  size_t uIdx = get(vertex_index, g, source(e,g));
  size_t vIdx = get(vertex_index, g, target(e,g));
  float w = get(edge_weight, g, e);
  //size_t source = get(edge_source, g, e);

  o.write(reinterpret_cast<const char *>(&uIdx), sizeof(uIdx));
  o.write(reinterpret_cast<const char *>(&vIdx), sizeof(vIdx));
  o.write(reinterpret_cast<const char *>(&w), sizeof(w));
  //o.write(reinterpret_cast<const char *>(&source), sizeof(source));
  return o;
}

ofstream & Mcr::write_to_stream(ofstream & o) const {

  // First write maps

  write_atom_to_stream(o, magic_id);

  write_vector_to_stream(o, relsSource);
  write_map_to_stream(o, relMap);

  write_map_to_stream(o, synsetMap);
  write_map_to_stream(o, wordMap);
  //write_map_to_stream(o, sourceMap);

  write_atom_to_stream(o, magic_id);

  size_t vertex_n = num_vertices(g);

  write_atom_to_stream(o, vertex_n);
  graph_traits<McrGraph>::vertex_iterator v_it, v_end;

  tie(v_it, v_end) = vertices(g);
  for(; v_it != v_end; ++v_it) {
    write_vertex_to_stream(o, g, *v_it);
  }

  write_atom_to_stream(o, magic_id);

  size_t edge_n = num_edges(g);

  write_atom_to_stream(o, edge_n);
  graph_traits<McrGraph>::edge_iterator e_it, e_end;

  tie(e_it, e_end) = edges(g);
  for(; e_it != e_end; ++e_it) {
    write_edge_to_stream(o, g, *e_it);
  }
  write_atom_to_stream(o, magic_id);
  if(notes.size()) write_vector_to_stream(o, notes);
  return o;
}

void Mcr::write_to_binfile (const string & fName) const {

  ofstream fo(fName.c_str(),  ofstream::binary|ofstream::out);
  if (!fo) {
    cerr << "Error: can't create" << fName << endl;
    exit(-1);
  }
  write_to_stream(fo);
}
