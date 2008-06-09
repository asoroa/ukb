#include "mcrGraph.h"
#include "common.h"
#include "globalVars.h"
//#include "w2syn.h"
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
    //throw "MCR not initialized";
    cerr << "MCR not initialized!" << endl;
    exit(-1);
  }
  return *p_instance;
}

void Mcr::create_from_txt(const string & relFileName,
			  const string & synsFileName,
			  const std::set<std::string> & rels_source) {
  if (p_instance) return;
  Mcr *tenp = create();

  tenp->read_from_txt(relFileName, synsFileName, rels_source);
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
  

  // From Inet
//   boost::dijkstra_shortest_paths(g, s, weight_map(weight).
// 				 visitor(boost::make_dijkstra_visitor(std::make_pair(
// 										     boost::record_distances(node_distance, boost::on_edge_relaxed()),
// 										     update_position))));

  dijkstra_shortest_paths(g, 
			  src,
			  predecessor_map(&parents[0]).
			  distance_map(&dist[0]));

// 			  boost::visitor(boost::make_dijkstra_visitor
// 					 (boost::make_list(mcr_bfs_init(&parents[0]),
// 							   mcr_bfs_pred(&parents[0])))));

			  //			  predecessor_map(&parents[0]).
			  //			  boost::visitor(boost::make_dijkstra_visitor(mcr_bfs_init(&parents[0]))));
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

Mcr_vertex_t Mcr::findOrInsertSynset(const string & str) {
  Mcr_vertex_t u;
  map<string, Mcr_vertex_t>::iterator it;
  bool insertedP;

  tie(it, insertedP) = synsetMap.insert(make_pair(str, Mcr_vertex_t()));
  if(insertedP) {
    // new vertex
    u = add_vertex(g);
    put(vertex_name, g, u, str);
    it->second = u;
  }
  return it->second;
}

Mcr_vertex_t Mcr::findOrInsertWord(const string & str) {
  Mcr_vertex_t u;
  map<string, Mcr_vertex_t>::iterator it;
  bool insertedP;

  tie(it, insertedP) = wordMap.insert(make_pair(str, Mcr_vertex_t()));
  if(insertedP) {
    // new vertex
    u = add_vertex(g);
    put(vertex_name, g, u, str);
    it->second = u;
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
    e = add_edge(u, v, g).first;
    put(edge_weight, g, e, w);
  }
  return e;
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

void read_relations(ifstream & relFile,
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
      cerr << "read_relations error. Bad line: " << line_number << endl;
      exit(-1);
    }
    int relId = lexical_cast<int>(fields[0]);
    relMap[fields[1]] = relId;
  }
  create_inverse_map(relMap, relMapInv);
}

// Read the actual MCR file

void read_mcr(ifstream & mcrFile, 
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
      cerr << "read_mcr error. Bad line: " << line_number << endl;
      exit(-1);
    }
    // last element says if relation is directed
    bool directed = (fields.size() > 4 && lexical_cast<int>(fields[4]) != 0);
    if (rels_source.find(fields[3]) == srel_end) continue; // Skip this relation
//     //if (fields[3] == "16") continue; // Skip WNet 1.6
//     //if (fields[3] == "20") continue; // Skip WNet 2.0
//     if (fields[3] == "ek") continue; // Skip Semcor selectional preferences (Eneko)
//     //if (fields[3] == "xg") continue; // Skip Xwnet gold relations
//     if (fields[3] == "su") continue; // Skip BNC selec. preferences (Diana)
//     if (fields[3] == "sc") continue; // Skip Semcor coocurrence relations
//     //if (fields[3] == "xn") continue; // Skip Xwnet normal relations
//     //if (fields[3] == "xs") continue; // Skip Xwnet silver relations

    Mcr_vertex_t u = mcr->findOrInsertSynset(fields[0]);
    Mcr_vertex_t v = mcr->findOrInsertSynset(fields[1]);
    //Mcr_vertex_t u = insert_synset_vertex(g, synsetMap, fields[0]);
    //Mcr_vertex_t v = insert_synset_vertex(g, synsetMap, fields[1]);

    // add edge
    mcr->findOrInsertEdge(u, v, 1.0);
    if (!directed)
      mcr->findOrInsertEdge(v, u, 1.0);
  }
}

void Mcr::read_from_txt(const string & relFileName,
			const string & synsFileName,
			const set<string> & rels_source) {

  relsSource = rels_source;
  // optimize IO
  std::ios::sync_with_stdio(false);

  std::ifstream relFile(relFileName.c_str(), ofstream::in);
  if (relFile) {
    read_relations(relFile, relMap, relMapInv);
  }

  std::ifstream synsFile(synsFileName.c_str(), ofstream::in);
  if (!synsFile) {
    cerr << "Can't open " << synsFileName << endl;
    throw;
  }

  read_mcr(synsFile, rels_source, this); // g, synsetMap, relsSource); //, relInv, sourceMap, vertexNames);
}

////////////////////////////////////////////////7
// public function

void Mcr::add_from_txt(const std::string & synsFileName) {
  if (!p_instance) {
    cerr << "add_relations_from_txt. MCR not initialized!" << endl;
    exit(-1);
  }
  // optimize IO
  std::ios::sync_with_stdio(false);

  std::ifstream synsFile(synsFileName.c_str(), ofstream::in);
  if (!synsFile) {
    cerr << "Can't open " << synsFileName << endl;
    throw;
  }

  read_mcr(synsFile, relsSource, this); // g, synsetMap, relsSource); //, relInv, sourceMap, vertexNames);  
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

    Mcr_vertex_t u; // = insert_synset_vertex(g, synsetMap, *syns_it);
    tie(u, auxP) = mcr.getVertexByName(syns.get_entry(i));
    if(!auxP) {
      if (glVars::verbose) 
	cerr << "W: Mcr::add_tokens warning: " << syns.get_entry(i) << " is not in MCR.\n";
      // Warning! 
      // do NOT insert node, because it becomes a dangling node and therefore
      // PageRank can't be applied
      // u = insert_synset_vertex(g, synsetMap, *syns_it); // NO
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

// void create_w2wpos_maps(const string & word,
// 			vector<string> & wPosV,
// 			map<string, vector<Syn_elem> > & wPos2Syns) {

//   const Mcr & mcr = Mcr::instance();
//   W2Syn & w2syn = W2Syn::instance();

//   vector<string>::const_iterator syns_it, syns_end;
//   vector<float>::const_iterator weight_it, weight_end;
//   tie(syns_it, syns_end) = w2syn.get_wsyns(word);
//   tie (weight_it, weight_end) = w2syn.get_weights(word);

//   char_separator<char> sep("-");
//   bool auxP;

//   for(;syns_it != syns_end; ++syns_it, ++weight_it) {
//     assert(weight_it != weight_end);
//     vector<string> fields(2);
//     tokenizer<char_separator<char> > tok(*syns_it, sep);
//     copy(tok.begin(), tok.end(), fields.begin());

//     string wpos(word.size() + 2, '#');
//     string::iterator sit = copy(word.begin(), word.end(), wpos.begin());
//     ++sit; // '#' char
//     *sit = fields[1].at(0); // the pos

//     Mcr_vertex_t u; // = insert_synset_vertex(g, synsetMap, *syns_it);
//     tie(u, auxP) = mcr.getVertexByName(*syns_it);
//     if(!auxP) {
//       if (glVars::verbose) 
// 	cerr << "W: Mcr::add_tokens warning: " << *syns_it << " is not in MCR.\n";
//       // Warning! 
//       // do NOT insert node, because it becomes a dangling node and therefore
//       // PageRank can't be applied
//       // u = insert_synset_vertex(g, synsetMap, *syns_it); // NO
//       continue;
//     }
    
//     map<string, vector<Syn_elem> >::iterator m_it;

//     tie(m_it, auxP) = wPos2Syns.insert(make_pair(wpos, vector<Syn_elem>()));
//     if (auxP) {
//       // first appearence of word#pos
//       wPosV.push_back(wpos);
//     }
//     m_it->second.push_back(make_pair(u, *weight_it));
//   }
// }


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
    Mcr_vertex_t wpos_v = mcr.findOrInsertSynset(*wpos_str_it);
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
    //init_out_coefs(g, V, &out_coefs[0], weight_map);
    prank::pageRank_iterate(g, V, &ppv_map[0],
			    weight_map, &ranks[0], &rank_tmp[0], 
			    glVars::prank::num_iterations);
  } else {
    //init_out_coefs(g, V, &out_coefs[0], cte_weight);
    prank::pageRank_iterate_now(g, V, &ppv_map[0],
				&ranks[0], &rank_tmp[0], 
				glVars::prank::num_iterations);
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

  read_atom_from_stream(is, id);
  if(id != magic_id) {
    cerr << "Error: invalid id (filename is a mcrGraph?)" << endl;
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
  read_vector_from_stream(is, notes);

  graph_traits<McrGraph>::vertex_iterator v_it, v_end;
  tie(v_it, v_end) = vertices(g);
  for(; v_it != v_end; ++v_it) {
    synsetMap[get(vertex_name, g, *v_it)] = *v_it;
  }
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
