#include "csentence.h"
#include "common.h"
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
    syns = cw_.syns;
  }
  return *this;
}

void fill_syns(const string & w,
	       vector<string> & syns,
	       char pos = 0) {

  W2Syn & w2syn = W2Syn::instance();
  Mcr & mcr = Mcr::instance();

  const vector<string> & syns_stringV = w2syn.get_syns(w);

  vector<string>::const_iterator str_it = syns_stringV.begin();
  vector<string>::const_iterator str_end = syns_stringV.end();
  for(;str_it != str_end; ++str_it) {
    bool existP;
    Mcr_vertex_t mcr_v;
    if(pos) {
      size_t m = str_it->size();
      if (str_it->at(m-1) != pos) continue;
    }
    tie(mcr_v, existP) = mcr.getVertexByName(*str_it);
    if (existP) {
      syns.push_back(*str_it);
    } else {
      cerr << "W: synset " << *str_it << " of word " << w << " is not in MCR" << endl;
      // debug: synset  which is not in mcr
    }
  }
}

CWord::CWord(const string & w_) : 
  w(w_), cw_id(string()), pos(0), distinguished(true) {
  fill_syns(w_, syns, 0);
}

CWord::CWord(const string & w_, const string & id_, char pos_, bool dist_) 
  : w(w_), cw_id(id_), pos(pos_), distinguished(dist_) {
  fill_syns(w_, syns, pos_);
}

std::ostream& operator<<(std::ostream & o, const CWord & cw_) {

  //McrGraph & g = Mcr::instance().graph();

  o << cw_.w;
  if (cw_.pos) 
    o << "-" << cw_.pos;
  o << "#" << cw_.cw_id << "#" << cw_.distinguished;
  o << "{";
  std::vector<string>::const_iterator it = cw_.syns.begin();
  std::vector<string>::const_iterator it_end = cw_.syns.end();
  if (it != it_end) {
    --it_end;
    copy(it, it_end, ostream_iterator<string>(o, ","));
    o << *it_end;
//     while(it != it_end) {
//       o << get(vertex_name, g, *it) << ",";
//       ++it;
//     }
//     o << get(vertex_name, g, *it_end);
  }
  o << "}";
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
  read_atom_from_stream(i,distinguished);
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
	cerr << "W: no synset for word " << fields[0] << "-" << fields[1] << endl;
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

////////////////////////////////////////////////////////
// Streaming

ofstream & CSentence::write_to_stream (std::ofstream & o) const {
  size_t vSize = v.size();

  write_atom_to_stream(o,vSize);
  for(size_t i=0; i< vSize; ++i) {
    v[i].write_to_stream(o);
  }
  write_atom_to_stream(o,cs_id);
  return o;
}

void CSentence::read_from_stream (std::ifstream & is) {

  size_t vSize;
  // Read a vector of CWords
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
