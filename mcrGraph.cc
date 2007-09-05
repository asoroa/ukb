#include "mcrGraph.h"
#include "common.h"
#include "globalVars.h"
#include "w2syn.h"
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

bool Mcr::bfs (const string & source_synset, 
	       std::vector<Mcr_vertex_t> & parents) const {

  std::map<std::string, Mcr_vertex_t>::const_iterator src_it = synsetMap.find(source_synset);
  if (src_it == synsetMap.end()) return false;
  bfs(src_it->second, parents);
  return true;
}


////////////////////////////////////////////////////////////////////////////////
// strings <-> vertex_id


pair<Mcr_vertex_t, bool> Mcr::getVertexByName(const std::string & str) const {
  map<string, Mcr_vertex_t>::const_iterator it = synsetMap.find(str);
  if (it == synsetMap.end()) return make_pair(Mcr_vertex_t(), false);
  return make_pair(it->second, true);
}


////////////////////////////////////////////////////////////////////////////////
// Random

boost::minstd_rand global_mersenne;

int g_randTarget(int Target) {
   // Return rand up to, but not including, Target
   boost::uniform_int<> local_int_dist(0, Target-1);
   //boost::variate_generator<boost::mt19937&, boost::uniform_int<> > local_uni_int(global_mersenne, local_int_dist);
   boost::variate_generator<boost::minstd_rand, boost::uniform_int<> > local_uni_int(global_mersenne, local_int_dist);
   return local_uni_int();
}
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


template <class G, class Map>
typename graph_traits<G>::vertex_descriptor insert_synset_vertex(G & g, 
								 Map & map,
								 const string & name) {
  typedef typename graph_traits<G>::vertex_descriptor Vertex_Desc;  
  typename Map::iterator pos;
  bool insertedP;
  Vertex_Desc u;

  tie(pos, insertedP) = map.insert(make_pair(name, Vertex_Desc()));
  if (insertedP) {
    // new vertex
    u = add_vertex(g);
    put(vertex_name, g, u, name);
    pos->second = u;
  } else {
    u = pos->second;
  }
  return u;
}

// Read the actual MCR file

void read_mcr(ifstream & mcrFile,
	      McrGraph & g,
	      map<std::string, Mcr_vertex_t> & synsetMap,
	      const set<string> & rels_source) {
  string line;
  int line_number = 0;
  bool aux;

  set<string>::const_iterator srel_end = rels_source.end();
  while(mcrFile) {
    vector<string> fields;
    std::getline(mcrFile, line, '\n');
    line_number++;
    char_separator<char> sep(" ");
    tokenizer<char_separator<char> > tok(line, sep);
    copy(tok.begin(), tok.end(), back_inserter(fields));
    if (fields.size() == 0) continue; // blank line
    if (fields.size() != 4) {
      cerr << "read_mcr error. Bad line: " << line_number << endl;
      exit(-1);
    }
    if (rels_source.find(fields[3]) == srel_end) continue; // Skip this relation
//     //if (fields[3] == "16") continue; // Skip WNet 1.6
//     //if (fields[3] == "20") continue; // Skip WNet 2.0
//     if (fields[3] == "ek") continue; // Skip Semcor selectional preferences (Eneko)
//     //if (fields[3] == "xg") continue; // Skip Xwnet gold relations
//     if (fields[3] == "su") continue; // Skip BNC selec. preferences (Diana)
//     if (fields[3] == "sc") continue; // Skip Semcor coocurrence relations
//     //if (fields[3] == "xn") continue; // Skip Xwnet normal relations
//     //if (fields[3] == "xs") continue; // Skip Xwnet silver relations

    Mcr_vertex_t u = insert_synset_vertex(g, synsetMap, fields[0]);
    Mcr_vertex_t v = insert_synset_vertex(g, synsetMap, fields[1]);

    // add edge
    Mcr_edge_t e;
    tie(e, aux) = edge(u, v, g);
    if(!aux) {
      tie(e, aux) = add_edge(u, v, g);
    }    
    // the other way around
    tie(e, aux) = edge(v, u, g);
    if(!aux) {
      tie(e, aux) = add_edge(v, u, g);
    }    
  }
}

void Mcr::read_from_txt(const string & relFileName,
			const string & synsFileName,
			const set<string> & rels_source) {

  relsSource = rels_source;
  // optimize IO
  std::ios::sync_with_stdio(false);

  std::ifstream relFile(relFileName.c_str(), ofstream::in);
  std::ifstream synsFile(synsFileName.c_str(), ofstream::in);

  read_relations(relFile, relMap, relMapInv);
  read_mcr(synsFile, g, synsetMap, relsSource); //, relInv, sourceMap, vertexNames);
}

void Mcr::display_info(std::ostream & o) const {

  o << "Relation sources: ";
  writeS(o, relsSource);
  o << "\n" << num_vertices(g) << " vertices and " << num_edges(g) << " edges" << endl;

}

////////////////////////////////////////////////////////////////////////////////
// Add token vertices and link them to synsets

void Mcr::add_tokens(const string & word, 
		     vector<string>::const_iterator syns_it, 
		     vector<string>::const_iterator syns_end) {

  map<string, vector<string> > w2wPos;
  map<string, vector<Mcr_vertex_t> > wPos2Syns;

  // Create w2wPos and wPos2Syns maps

  char_separator<char> sep("-");
  bool auxP;

  //vector<string>::const_iterator it = syns.begin();
  //vector<string>::const_iterator end = syns.end();
  for(;syns_it != syns_end; ++syns_it) {

    vector<string> fields(2);
    tokenizer<char_separator<char> > tok(*syns_it, sep);
    copy(tok.begin(), tok.end(), fields.begin());

    string wpos(word.size() + 2, '#');
    string::iterator sit = copy(word.begin(), word.end(), wpos.begin());
    ++sit; // '#' char
    *sit = fields[1].at(0); // the pos

    Mcr_vertex_t u; // = insert_synset_vertex(g, synsetMap, *syns_it);
    tie(u, auxP) = getVertexByName(*syns_it);
    if(!auxP) {
      if (glVars::verbose) 
	cerr << "W: Mcr::add_tokens warning: " << *syns_it << " is not in MCR.\n";
      // Warning! 
      // do NOT insert node, because it becomes a dangling node and therefore
      // PageRank can't be applied
      // u = insert_synset_vertex(g, synsetMap, *syns_it); // NO
      continue;
    }
    
    map<string, vector<Mcr_vertex_t> >::iterator m_it;

    tie(m_it, auxP) = wPos2Syns.insert(make_pair(wpos, vector<Mcr_vertex_t>()));
    if (auxP) {
      // first appearence of word#pos
      w2wPos[word].push_back(wpos);
    }
    m_it->second.push_back(u);
  }

  // Add vertices and link them in the MCR

  map<string, vector<string> >::iterator w2wpos_it = w2wPos.begin();
  map<string, vector<string> >::iterator w2wpos_end = w2wPos.end();
  for(; w2wpos_it != w2wpos_end; ++w2wpos_it) {
    // insert word
    Mcr_vertex_t word_v = insert_synset_vertex(g, synsetMap, w2wpos_it->first); 

    vector<string>::iterator wpos_str_it = w2wpos_it->second.begin();
    vector<string>::iterator wpos_str_end = w2wpos_it->second.end();
    for (; wpos_str_it != wpos_str_end; ++wpos_str_it) {
      //insert word#pos
      Mcr_vertex_t wpos_v = insert_synset_vertex(g, synsetMap, *wpos_str_it);
      // link word to word#pos
      add_edge(word_v, wpos_v, g);
      // link word#pos to synsets
      vector<Mcr_vertex_t> & syns = wPos2Syns[*wpos_str_it];
      vector<Mcr_vertex_t>::iterator syns_it = syns.begin();
      vector<Mcr_vertex_t>::iterator syns_end = syns.end();
      for(;syns_it != syns_end; ++syns_it) {
	add_edge(wpos_v, *syns_it, g);
      }
    }      
  }
}


////////////////////////////////////////////////////////////////////////////////
// PageRank in MCR


// PPV version

void Mcr::pageRank_ppv(const vector<float> & ppv_map,
		       vector<float> & ranks) {

  size_t N = num_vertices(g);
  vector<float> map_tmp(N, 0.0f);
  vector<float> out_coefs(N, 0);
  vector<float>(N, 0.0).swap(ranks); // Initialize rank vector
  vector<float> rank_tmp(N, 0.0);    // auxiliary rank vector

  // ugly ugly hack @@CHANGE ME !!!

  vector<Mcr_vertex_t> V(N);

  graph_traits<McrGraph>::vertex_iterator it, end;
  tie(it, end) = vertices(g);
  copy(it, end, V.begin());

  constant_property_map <Mcr_edge_t, float> cte_weight(1); // always return 1  
  init_out_coefs(g, V, out_coefs, cte_weight);
  pageRank_iterate(g, V, &ppv_map[0], out_coefs, cte_weight, &ranks[0], &rank_tmp[0], 30); // 30 iterations
}



////////////////////////////////////////////////////////////////////////////////
// Streaming

const size_t magic_id = 0x070201;

// read

Mcr_vertex_t read_vertex_from_stream(ifstream & is, 
				     McrGraph & g) {

  string name;
  string wname;

  read_atom_from_stream(is, name);
  read_atom_from_stream(is, wname);
  Mcr_vertex_t v = add_vertex(g);
  put(vertex_name, g, v, name);

  return v;
}

Mcr_edge_t read_edge_from_stream(ifstream & is, 
				 McrGraph & g) {

  size_t sIdx;
  size_t tIdx; 
  //size_t id;
  //size_t source;
  bool insertedP;
  Mcr_edge_t e;

  read_atom_from_stream(is, sIdx);
  read_atom_from_stream(is, tIdx);
  //read_atom_from_stream(is, id);
  //read_atom_from_stream(is, source);
  tie(e, insertedP) = add_edge(sIdx, tIdx, g);
  assert(insertedP);
  //put(edge_id, g, e, id);
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
  //read_map(is, relType);
  //read_map(is, relInv);
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
  string wname("");

  write_atom_to_stream(o, get(vertex_name, g, v));
  write_atom_to_stream(o, wname);
  return o;
}

ofstream & write_edge_to_stream(ofstream & o,
				const McrGraph & g,
				const Mcr_edge_t & e) {

  size_t uIdx = get(vertex_index, g, source(e,g));
  size_t vIdx = get(vertex_index, g, target(e,g));
  //size_t id = get(edge_id, g, e);
  //size_t source = get(edge_source, g, e);

  o.write(reinterpret_cast<const char *>(&uIdx), sizeof(uIdx));
  o.write(reinterpret_cast<const char *>(&vIdx), sizeof(vIdx));
  //o.write(reinterpret_cast<const char *>(&id), sizeof(id));
  //o.write(reinterpret_cast<const char *>(&source), sizeof(source));
  return o;
}

ofstream & Mcr::write_to_stream(ofstream & o) const {

  // First write maps

  write_atom_to_stream(o, magic_id);

  write_vector_to_stream(o, relsSource);
  write_map_to_stream(o, relMap);
  //write_map_to_stream(o, relType);
  //write_map_to_stream(o, relInv);
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
