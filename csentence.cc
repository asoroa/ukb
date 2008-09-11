#include "csentence.h"
#include "common.h"
#include "globalVars.h"
#include "mcrGraph.h"
#include "wdict.h"

// Tokenizer
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

#include<boost/tuple/tuple.hpp> // for "tie"


using namespace std;
using namespace boost;

void fill_syns(const string & w,
			   vector<string> & syns,
			   vector<float> & ranks,
			   char pos = 0) {

  //W2Syn & w2syn = W2Syn::instance();
  Mcr & mcr = Mcr::instance();

  vector<string>::const_iterator str_it;
  vector<string>::const_iterator str_end;

  WDict_entries entries = WDict::instance().get_entries(w);

  for(size_t i= 0; i < entries.size(); ++i) {

    const string & syn_str = entries.get_entry(i);

    bool existP;
    Mcr_vertex_t mcr_v;
    if(pos) {
      // filter synsets by pos
      char synpos = entries.get_pos(i);
      if(!synpos) {
		throw std::runtime_error("CWord: Error reading context. " + syn_str + " has no POS\n");
      }
      if (pos != synpos) continue;
    }
    tie(mcr_v, existP) = mcr.get_vertex_by_name(syn_str);
    if (existP) {
      syns.push_back(syn_str);
      ranks.push_back(0.0f);
    } else {
      if (glVars::debug::warning) {
		cerr << "W:CWord: synset " << syn_str << " of word " << w << " is not in MCR" << endl;
      }
      // debug: synset  which is not in mcr
    }
  }
}

void CWord::shuffle_synsets() {
  boost::random_number_generator<boost::mt19937, long int> rand_dist(glVars::rand_generator);
  //std::random_shuffle(syns.begin(), syns.end(), rand_dist);
  std::random_shuffle(m_syns.begin(), m_syns.end(), rand_dist);
}

// CWord::CWord(const string & w_) :
//   w(w_), m_id(string()), pos(0), m_weight(1.0),
//   m_is_synset(false), m_distinguished(true), m_ranks_equal(true) {
//   fill_syns(w_, syns, ranks, 0);
//   m_disamb = (1 == syns.size()); // monosemous words are disambiguated
//   shuffle_synsets();
// }

CWord::CWord(const string & w_, const string & id_, char pos_, bool dist_)
  : w(w_), m_id(id_), m_pos(pos_), m_weight(1.0), m_is_synset(false),
	m_distinguished(dist_) {
  fill_syns(w_, m_syns, m_ranks, pos_);
  m_disamb = (1 == m_syns.size()); // monosemous words are disambiguated
  shuffle_synsets();
}

CWord & CWord::operator=(const CWord & cw_) {
  if (&cw_ != this) {
    w = cw_.w;
    m_id = cw_.m_id;
    m_weight = cw_.m_weight;
    m_pos = cw_.m_pos;
    m_syns = cw_.m_syns;
    m_ranks = cw_.m_ranks;
    m_is_synset = cw_.m_is_synset;
    m_distinguished = cw_.m_distinguished;
    m_disamb = cw_.m_disamb;
  }
  return *this;
}

// Create a special CWord which is a synset and a weigth.

CWord CWord::create_synset_cword(const string & syn, const string & id_, float w) {
  Mcr_vertex_t aux;
  bool P;
  tie(aux, P) = Mcr::instance().get_vertex_by_name(syn);
  if (!P) {
	throw std::runtime_error("CWord::create_synset_cword " + syn + " not in MCR");
	return CWord();
  }
  CWord cw;
  cw.m_weight = w;
  cw.w = syn;
  cw.m_id = id_;
  cw.m_is_synset = false;
  cw.m_distinguished=false;
  cw.m_syns.push_back(syn);
  cw.m_ranks.push_back(0.0f);
  return cw;
}

string CWord::wpos() const {

  if (m_is_synset) {
	std::runtime_error("CWoord::wpos: can't get wpos of Cword synset " + w);
  }
  string wpos(w);
  wpos.append("#");
  char pos = get_pos();
  wpos.append(1,pos);
  return wpos;
};

struct CWSort {

  CWSort(const vector<float> & _v) : v(_v) {}
  int operator () (const int & i, const int & j) {
    // Descending order
    return v[i] > v[j];
  }
  const vector<float> & v;
};

void CWord::disamb_cword() {

  size_t n = m_syns.size();
  if (!n) return;
  if (n == 1) {
	m_disamb = true; // Monosemous words are disambiguated
	return;
  }

  vector<string> syns(n);
  vector<float> ranks(n);
  vector<int> idx(n);

  for(size_t i=0; i < n; ++i)
    idx[i] = i;
  sort(idx.begin(), idx.end(), CWSort(m_ranks));

  bool ranks_equal = true;
  for(size_t i=0; i < n; ++i) {
    syns[i]  = m_syns[idx[i]];
    ranks[i] = m_ranks[idx[i]];
	if (ranks[i] != ranks[0])
	  ranks_equal = false;
  }
  if (ranks_equal) return; // If all ranks have same value the word is not disambiguated
  syns.swap(m_syns);
  ranks.swap(m_ranks);
  m_disamb = true;
}


std::ostream& operator<<(std::ostream & o, const CWord & cw_) {

  //McrGraph & g = Mcr::instance().graph();

  o << cw_.w;
  if (cw_.m_pos)
    o << "-" << cw_.m_pos;
  if(glVars::csentence::word_weight)
    o << "#" << cw_.m_weight;
  o << "#" << cw_.m_id << "#" << cw_.m_distinguished << " " << cw_.m_disamb;
  o << '\n';
  for(size_t i = 0; i < cw_.m_syns.size(); ++i) {
    assert(i < cw_.m_ranks.size());
    o << cw_.m_syns[i] << ":" << cw_.m_ranks[i] << '\n';
  }
  return o;
}

ostream & cw_aw_print_best(ostream & o,
						   const vector<string> & syns,
						   const vector<float> & ranks) {
  size_t n = syns.size();
  o << " " << syns[0];
  for(size_t i = 1; i != n; ++i) {
    if (ranks[i] != ranks[0]) break;
    o << " " << syns[i];
  }
  return o;
}

ostream & cw_aw_print_all(ostream & o,
						  const vector<string> & syns,
						  const vector<float> & ranks) {
  float rsum = 0.0;
  for(size_t i = 0; i != syns.size(); ++i) {
    rsum += ranks[i];
  }
  float norm_factor = 1.0 / rsum;
  for(size_t i = 0; i != syns.size(); ++i) {
    o << " " << syns[i] << "/" << ranks[i]*norm_factor;
  }
  return o;
}

ostream & CWord::print_cword_aw(ostream & o) const {

  vector<string> id_fields(split(m_id, "."));
  assert(id_fields.size() > 0);
  o << id_fields[0] << " " << m_id;
  if(!glVars::output::allranks) cw_aw_print_best(o, m_syns, m_ranks);
  else cw_aw_print_all(o, m_syns, m_ranks);
  o << " !! " << w << "\n";
  return o;
}

ostream & CWord::print_cword_semcor_aw(ostream & o) const {

  if ((1 == m_syns.size()) && !glVars::output::monosemous) return o; // Don't output monosemous words

  vector<string> id_fields(split(m_id, "."));
  assert(id_fields.size() > 0);
  o << id_fields[0] << "." << id_fields[1] << " " << m_id;
  if(!glVars::output::allranks) cw_aw_print_best(o, m_syns, m_ranks);
  else cw_aw_print_all(o, m_syns, m_ranks);
  o << " !! " << w << "\n";
  return o;
}

ostream & CWord::print_cword_simple(ostream & o) const {

  o << m_id << " ";
  if(!glVars::output::allranks) cw_aw_print_best(o, m_syns, m_ranks);
  else cw_aw_print_all(o, m_syns, m_ranks);
  o << " !! " << w << "\n";
  return o;
}

////////////////////////////////////////////////////////
// Streaming

std::ofstream & CWord::write_to_stream(std::ofstream & o) const {

  write_atom_to_stream(o,w);
  write_atom_to_stream(o,m_id);
  write_atom_to_stream(o,m_pos);
  write_vector_to_stream(o,m_syns);
  write_atom_to_stream(o,m_distinguished);
  return o;

};

void CWord::read_from_stream(std::ifstream & i) {

  read_atom_from_stream(i,w);
  read_atom_from_stream(i,m_id);
  read_atom_from_stream(i,m_pos);
  read_vector_from_stream(i,m_syns);
  vector<float>(m_syns.size()).swap(m_ranks); // Init ranks vector
  read_atom_from_stream(i,m_distinguished);

  m_disamb = (1 == m_syns.size());
  shuffle_synsets();
};

////////////////////////////////////////////////////////////////
// CSentence

// AW file read (create a csentence from a context)

static bool skip_line(const string &l) {

  if (l.size() == 0) return true;
  string::const_iterator sIt = l.begin();
  string::const_iterator sItEnd = l.end();
  while(sIt != sItEnd && isspace(*sIt)) sIt++;
  if (sIt == sItEnd) return true;
  if (*sIt == '#') return true;
  return false;
}

istream & CSentence::read_aw(istream & is) {

  string line;
  int l_n = 0;

  if(is) {
	do {
	  getline(is, line, '\n');
	  l_n++;
	} while(is && skip_line(line));
	if(!is) return is;

    // first line is id
    char_separator<char> sep(" ");
    vector<string> ctx;

    tokenizer<char_separator<char> > tok_id(line, sep);
    copy(tok_id.begin(), tok_id.end(), back_inserter(ctx));
    if (ctx.size() == 0) return is; // blank line or EOF
    cs_id = ctx[0];
    vector<string>().swap(ctx);
	// next comes the context
	do {
	  getline(is, line, '\n');
	  l_n++;
	} while(is && skip_line(line));
	if(!is) return is;

    tokenizer<char_separator<char> > tok_ctx(line, sep);
    copy(tok_ctx.begin(), tok_ctx.end(), back_inserter(ctx));
    if (ctx.size() == 0) return is; // blank line or EOF
    int i = 0;
	try {
	  for(vector<string>::const_iterator it = ctx.begin(); it != ctx.end(); ++i, ++it) {
		vector<string> fields;
		char_separator<char> wsep("#");
		tokenizer<char_separator<char> > wtok(*it, wsep);
		copy(wtok.begin(), wtok.end(), back_inserter(fields));
		if (fields.size() < 4)
		  throw std::runtime_error("Bad word " + *it + " " + lexical_cast<string>(i));
		int dist = lexical_cast<int>(fields[3]);
		CWord new_cw(fields[0], fields[2], fields[1].at(0), dist > 0);
		if (dist == 2) {
		  new_cw = CWord::create_synset_cword(fields[0], fields[2], 1.0);
		}
		if(glVars::csentence::word_weight) {
		  if (fields.size() != 5)
			throw std::runtime_error("Bad word " + *it + " " + lexical_cast<string>(i) + ": no weight");
		  // optional weight
		  new_cw.set_weight(lexical_cast<float>(fields[4]));
		}
		if (new_cw.size()) {
		  v.push_back(new_cw);
		} else {
		  // No synset for that word.
		  cerr << "W: " << fields[0] << "-" << fields[1] << " can't be mapped to MCR." << endl;
		}		
	  }
	} catch (std::exception & e) {
	  throw std::runtime_error("Context error in line " + lexical_cast<string>(l_n) + " : " + e.what() );
	}
  }
  return is;
}


// CSentence::CSentence(const vector<string> & sent) : cs_id(string()) {

//   set<string> wordS;

//   vector<string>::const_iterator s_it = sent.begin();
//   vector<string>::const_iterator s_end = sent.end();
//   for(;s_it != s_end; ++s_it) {
//     bool insertedP;
//     set<string>::iterator it;
//     tie(it, insertedP) = wordS.insert(*s_it);
//     if (insertedP) {
//       v.push_back(CWord(*s_it));
//     }
//   }
// }

CSentence & CSentence::operator=(const CSentence & cs_) {
  if (&cs_ != this) {
    v = cs_.v;
    cs_id = cs_.cs_id;
  }
  return *this;
}

void CSentence::append(const CSentence & cs_) {
  vector<CWord> tenp(v.size() + cs_.v.size());
  vector<CWord>::iterator aux_it;
  aux_it = copy(v.begin(), v.end(), tenp.begin());
  copy(cs_.v.begin(), cs_.v.end(), aux_it);
  v.swap(tenp);
}

void CSentence::distinguished_synsets(vector<string> & res) const {

  vector<CWord>::const_iterator cw_it, cw_end;
  cw_it = v.begin();
  cw_end = v.end();
  for(; cw_it != cw_end; ++cw_it) {
    if (!cw_it->is_distinguished()) continue;
    copy(cw_it->begin(), cw_it->end(), back_inserter(res));
  }
}

std::ostream& operator<<(std::ostream & o, const CSentence & cs_) {
  o << cs_.cs_id << endl;
  copy(cs_.begin(), cs_.end(), ostream_iterator<CWord>(o, "\n"));
  return o;
}

std::ostream & CSentence::print_csent_aw(std::ostream & o) const {

  vector<CWord>::const_iterator cw_it = v.begin();
  vector<CWord>::const_iterator cw_end = v.end();

  for(; cw_it != cw_end; ++cw_it) {
    if (cw_it->size() == 0) continue;
    if (!cw_it->is_distinguished()) continue;
    if (!cw_it->is_disambiguated()) continue;
	if (cw_it->is_monosemous() && !glVars::output::monosemous) return o; // Don't output monosemous words

    cw_it->print_cword_aw(o);
  }
  return o;
}

std::ostream & CSentence::print_csent_semcor_aw(std::ostream & o) const {

  vector<CWord>::const_iterator cw_it = v.begin();
  vector<CWord>::const_iterator cw_end = v.end();

  for(; cw_it != cw_end; ++cw_it) {
    if (cw_it->size() == 0) continue;
    if (!cw_it->is_distinguished()) continue;
    if (!cw_it->is_disambiguated()) continue;
	if (cw_it->is_monosemous() && !glVars::output::monosemous) return o; // Don't output monosemous words

    cw_it->print_cword_semcor_aw(o);
  }
  return o;
}

std::ostream & CSentence::print_csent_simple(std::ostream & o) const {

  vector<CWord>::const_iterator cw_it = v.begin();
  vector<CWord>::const_iterator cw_end = v.end();

  for(; cw_it != cw_end; ++cw_it) {
    if (cw_it->size() == 0) continue;
    if (!cw_it->is_distinguished()) continue;
    if (!cw_it->is_disambiguated()) continue;
	if (cw_it->is_monosemous() && !glVars::output::monosemous) return o; // Don't output monosemous words
	o << cs_id << " ";
    cw_it->print_cword_simple(o);
  }
  return o;
}

////////////////////////////////////////////////////////////////////////////////
//
// PageRank in Mcr


// Given a CSentence obtain it's PageRank vector
// Initial PPV is computed a la hughes&ramage97

bool calculate_mcr_hr(const CSentence & cs,
					  vector<float> & res,
					  bool with_weight) {

  Mcr & mcr = Mcr::instance();
  bool aux;

  // Initialize result vector
  vector<float> (mcr.size(), 0.0).swap(res);

  // Initialize PPV vector
  vector<float> ppv(mcr.size(), 0.0);
  CSentence::const_iterator it = cs.begin();
  CSentence::const_iterator end = cs.end();
  float K = 0.0;
  for(;it != end; ++it) {
    //if(!it->is_distinguished()) continue;

    //string wpos = it->wpos();
    string wpos = it->word();

    Mcr_vertex_t u;
    tie(u, aux) = mcr.get_vertex_by_name(wpos);
    float w = 1.0;
    if(glVars::csentence::word_weight) w = it->get_weight();
    if (aux) {
      ppv[u] = w;
      K +=w;
    }
  }
  if (K == 0.0) return false;
  // Normalize PPV vector
  float div = 1.0 / static_cast<float>(K);
  for(vector<float>::iterator rit = ppv.begin(); rit != ppv.end(); ++rit)
    *rit *= div;

  // Execute PageRank
  mcr.pageRank_ppv(ppv, res, with_weight);
  return true;
}

// Second method for hughes & ramage
// given a word,
// 1. put a ppv in the synsets of the rest of words.
// 2. Pagerank
// 3. use rank for disambiguating word

void calculate_mcr_hr_by_word_and_disamb(CSentence & cs,
										 bool with_weight) {

  Mcr & mcr = Mcr::instance();
  bool aux;

  vector<CWord>::iterator cw_it = cs.begin();
  vector<CWord>::iterator cw_end = cs.end();
  for(; cw_it != cw_end; ++cw_it) {

    // Target word must be distinguished.
    if(!cw_it->is_distinguished()) continue;

    vector<float> ranks (mcr.size(), 0.0);
    // Initialize PPV vector
    vector<float> ppv(mcr.size(), 0.0);
    float K = 0.0;
    // put ppv to the synsets of words except cw_it
    for(CSentence::const_iterator it = cs.begin(); it != cw_end; ++it) {
      if(it == cw_it) continue;

      string wpos = it->word();

      Mcr_vertex_t u;
      tie(u, aux) = mcr.get_vertex_by_name(wpos);
      float w =  1.0;
      if(glVars::csentence::word_weight) w = it->get_weight();
      if (aux) {
		ppv[u] = w;
		K +=w;
      }
    }
    if (K == 0.0) continue;
    // Normalize PPV vector
    float div = 1.0 / static_cast<float>(K);
    for(vector<float>::iterator rit = ppv.begin(); rit != ppv.end(); ++rit)
      *rit *= div;

    // Execute PageRank
    mcr.pageRank_ppv(ppv, ranks, with_weight);
    // disambiguate cw_it
    cw_it->rank_synsets(mcr, ranks);
    cw_it->disamb_cword();
  }
}

//
// Given a previously disambiguated CSentence (all synsets of words
// have a rank), calculate a mcr prgaRank where PPV is formed by
// 'activating' just the synsets of CWords with appropiate
// (normalized) rank
//

bool calculate_mcr_ppv_csentence(CSentence & cs, vector<float> & res) {

  Mcr & mcr = Mcr::instance();
  bool aux;

  // Initialize result vector
  vector<float> (mcr.size(), 0.0).swap(res);

  // Initialize PPV vector
  vector<float> ppv(mcr.size(), 0.0);

  vector<CWord>::iterator cw_it = cs.begin();
  vector<CWord>::iterator cw_end = cs.end();
  float K = 0.0;
  for(; cw_it != cw_end; ++cw_it) {
    for(size_t i = 0; i != cw_it->size(); ++i) {
      Mcr_vertex_t u;
      tie(u, aux) = mcr.get_vertex_by_name(cw_it->syn(i));
      if (aux) {
		ppv[u] = cw_it->rank(i);
		K += cw_it->rank(i);
      }
    }
  }

  if (K == 0.0) return false;
  // Normalize PPV vector
  float div = 1.0 / K;
  for(vector<float>::iterator rit = ppv.begin(); rit != ppv.end(); ++rit)
    *rit *= div;

  // Execute PageRank
  mcr.pageRank_ppv(ppv, res, false);
  return true;
}


//
// Disambiguate a CSentence given a vector of ranks
//

void disamb_csentence_mcr(CSentence & cs,
						  vector<float> & ranks) {

  Mcr & mcr = Mcr::instance();

  vector<CWord>::iterator cw_it = cs.begin();
  vector<CWord>::iterator cw_end = cs.end();
  for(; cw_it != cw_end; ++cw_it) {
    cw_it->rank_synsets(mcr, ranks);
    cw_it->disamb_cword();
  }
}


////////////////////////////////////////////////////////
// Streaming

const size_t magic_id = 0x070704;

ofstream & CSentence::write_to_stream (std::ofstream & o) const {
  size_t vSize = v.size();

  write_atom_to_stream(o, magic_id);
  write_atom_to_stream(o,vSize);
  for(size_t i=0; i< vSize; ++i) {
    v[i].write_to_stream(o);
  }
  write_atom_to_stream(o,cs_id);
  return o;
}

void CSentence::read_from_stream (std::ifstream & is) {

  size_t vSize;
  size_t id;

  // Read a vector of CWords
  read_atom_from_stream(is, id);
  if (id != magic_id) {
    cerr << "Invalid file. Not a csentence." << endl;
    exit(-1);
  }
  read_atom_from_stream(is, vSize);
  v.resize(vSize);
  for(size_t i=0; i<vSize; ++i) {
    v[i].read_from_stream(is);
  }
  read_atom_from_stream(is,cs_id);
}

void CSentence::write_to_binfile (const std::string & fName) const {
  ofstream fo(fName.c_str(),  ofstream::binary|ofstream::out);
  if (!fo) {
    cerr << "Error: can't create" << fName << endl;
    exit(-1);
  }
  write_to_stream(fo);
}


void CSentence::read_from_binfile (const std::string & fName) {
  ifstream fi(fName.c_str(), ifstream::binary|ifstream::in);
  if (!fi) {
    cerr << "Error: can't open " << fName << endl;
    exit(-1);
  }
  read_from_stream(fi);
}
