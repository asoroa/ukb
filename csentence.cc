#include "csentence.h"
#include "common.h"
#include "globalVars.h"
#include "mcrGraph.h"
#include "w2syn.h"

// Tokenizer
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

#include<boost/tuple/tuple.hpp> // for "tie"


using namespace std;
using namespace boost;

CWord & CWord::operator=(const CWord & cw_) {
  if (&cw_ != this) {
    w = cw_.w;
    cw_id = cw_.cw_id;
    pos = cw_.pos;
    syns = cw_.syns;
    ranks = cw_.ranks;
    distinguished = cw_.distinguished;
    distinguished = cw_.ranks_equal;
    disamb = cw_.disamb;
  }
  return *this;
}

void fill_syns(const string & w,
	       vector<string> & syns,
	       vector<float> & ranks,
	       char pos = 0) {

  W2Syn & w2syn = W2Syn::instance();
  Mcr & mcr = Mcr::instance();

  vector<string>::const_iterator str_it;
  vector<string>::const_iterator str_end;
  tie(str_it, str_end) = w2syn.get_wsyns(w);
  for(;str_it != str_end; ++str_it) {
    bool existP;
    Mcr_vertex_t mcr_v;
    if(pos) {
      // filter synsets by pos
      size_t m = str_it->size();
      if (str_it->at(m-1) != pos) continue;
    }
    tie(mcr_v, existP) = mcr.getVertexByName(*str_it);
    if (existP) {
      syns.push_back(*str_it);
      ranks.push_back(0.0f);
    } else {
      cerr << "W: synset " << *str_it << " of word " << w << " is not in MCR" << endl;
      // debug: synset  which is not in mcr
    }
  }
}

void CWord::shuffle_synsets() {
  boost::random_number_generator<boost::mt19937, long int> rand_dist(glVars::rand_generator);
  //std::random_shuffle(syns.begin(), syns.end(), rand_dist);
  std::random_shuffle(syns.begin(), syns.end(), rand_dist);
}

CWord::CWord(const string & w_) : 
  w(w_), cw_id(string()), pos(0), distinguished(true), ranks_equal(true) {
  fill_syns(w_, syns, ranks, 0);
  disamb = (1 == syns.size()); // monosemous words are disambiguated
  shuffle_synsets();
}

CWord::CWord(const string & w_, const string & id_, char pos_, bool dist_) 
  : w(w_), cw_id(id_), pos(pos_), distinguished(dist_), ranks_equal(true) {
  fill_syns(w_, syns, ranks, pos_);
  disamb = (1 == syns.size()); // monosemous words are disambiguated
  shuffle_synsets();
}

string CWord::wpos() const {

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
  if (ranks_equal) return;
  size_t n = syns.size();
  vector<string> n_syns(n);
  vector<float> n_ranks(n);
  vector<int> idx(n);

  for(size_t i=0; i < n; ++i)
    idx[i] = i;
  sort(idx.begin(), idx.end(), CWSort(ranks));

  for(size_t i=0; i < n; ++i) {
    n_syns[i]  = syns[idx[i]];
    n_ranks[i] = ranks[idx[i]];
  }
  n_syns.swap(syns);
  n_ranks.swap(ranks);
  disamb = true;
}


std::ostream& operator<<(std::ostream & o, const CWord & cw_) {

  //McrGraph & g = Mcr::instance().graph();

  o << cw_.w;
  if (cw_.pos) 
    o << "-" << cw_.pos;
  o << "#" << cw_.cw_id << "#" << cw_.distinguished << " " << cw_.ranks_equal << " " << cw_.disamb;
  o << '\n';
  for(size_t i = 0; i < cw_.syns.size(); ++i) {
    assert(i < cw_.ranks.size());
    o << cw_.syns[i] << ":" << cw_.ranks[i] << '\n';
  }
  return o;
}

ostream & CWord::print_cword_aw(ostream & o) const {

  if (!disamb) return o;
  if ((1 == syns.size()) && !glVars::output_monosemous) return o; // Don't output monosemous words

  vector<string> id_fields(split(cw_id, "."));
  assert(id_fields.size() > 0);
  o << id_fields[0] << " " << cw_id << " ";
  size_t n = syns.size();
  o << syns[0];
  for(size_t i = 1; i != n; ++i) {
    if (ranks[i] != ranks[0]) break;
    o << " " << syns[i];
  }
  o << " !! " << w << "\n";
  return o;
}

////////////////////////////////////////////////////////
// Streaming

std::ofstream & CWord::write_to_stream(std::ofstream & o) const {

  write_atom_to_stream(o,w);
  write_atom_to_stream(o,cw_id);
  write_atom_to_stream(o,pos);
  write_vector_to_stream(o,syns);
  write_atom_to_stream(o,distinguished);
  return o;

};

void CWord::read_from_stream(std::ifstream & i) {

  read_atom_from_stream(i,w);
  read_atom_from_stream(i,cw_id);
  read_atom_from_stream(i,pos);
  read_vector_from_stream(i,syns);
  vector<float>(syns.size()).swap(ranks); // Init ranks vector
  read_atom_from_stream(i,distinguished);

  disamb = (1 == syns.size());
  shuffle_synsets();
};

////////////////////////////////////////////////////////////////
// CSentence

// AW file read (create a csentence from a context)

istream & CSentence::read_aw(istream & is) {

  string line;

  if(is) {
    getline(is, line, '\n');
    if (!is) return is; // EOF

    // first line is id
    char_separator<char> sep(" ");
    vector<string> ctx;

    tokenizer<char_separator<char> > tok_id(line, sep);
    copy(tok_id.begin(), tok_id.end(), back_inserter(ctx));
    if (ctx.size() == 0) return is; // blank line or EOF
    cs_id = ctx[0];
    vector<string>().swap(ctx);
    // next comes the context
    getline(is, line, '\n');
    tokenizer<char_separator<char> > tok_ctx(line, sep);
    copy(tok_ctx.begin(), tok_ctx.end(), back_inserter(ctx));
    if (ctx.size() == 0) return is; // blank line or EOF
    int i = 0;
    for(vector<string>::const_iterator it = ctx.begin(); it != ctx.end(); ++i, ++it) {
      vector<string> fields;
      char_separator<char> wsep("#");
      tokenizer<char_separator<char> > wtok(*it, wsep);
      copy(wtok.begin(), wtok.end(), back_inserter(fields));
      if (fields.size() != 4) 
	throw "Bad word " + lexical_cast<string>(i);
      CWord new_cw(fields[0], fields[2], fields[1].at(0), lexical_cast<bool>(fields[3]));
      if (new_cw.size()) {
	v.push_back(new_cw);
      } else {
	// No synset for that word. 
	cerr << "W: " << fields[0] << "-" << fields[1] << " can't be mapped to MCR." << endl;
      }
    }
  }
  return is;
}


CSentence::CSentence(const vector<string> & sent) : cs_id(string()) {

  set<string> wordS;
  
  vector<string>::const_iterator s_it = sent.begin();
  vector<string>::const_iterator s_end = sent.end();
  for(;s_it != s_end; ++s_it) {
    bool insertedP;
    set<string>::iterator it;
    tie(it, insertedP) = wordS.insert(*s_it);
    if (insertedP) {
      v.push_back(CWord(*s_it));
    }
  }
}

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
    cw_it->print_cword_aw(o);
  }
  return o;
}

////////////////////////////////////////////////////////////////////////////////
// 
// PageRank in Mcr


// Given a CSentence obtain it's PageRank vector
// Initial PPV is computed a la hughes&ramage97

bool calculate_mcr_ranks(const CSentence & cs,
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
  size_t K = 0;
  for(;it != end; ++it) {
    //if(!it->is_distinguished()) continue;
    
    //string wpos = it->wpos();
    string wpos = it->word();

    Mcr_vertex_t u;
    tie(u, aux) = mcr.getVertexByName(wpos);
    if (aux) {
      ppv[u] = 1;
      ++K;
    }
  }
  if (!K) return false;
  // Normalize PPV vector
  float div = 1.0 / static_cast<float>(K);
  for(vector<float>::iterator rit = ppv.begin(); rit != ppv.end(); ++rit) 
    *rit *= div;

  // Execute PageRank
  mcr.pageRank_ppv(ppv, res, with_weight);
  return true;
}

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
