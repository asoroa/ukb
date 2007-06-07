
#include "w2syn.h"

#include <fstream>
#include <iostream>

#include<boost/tuple/tuple.hpp> // for "tie"

// Tokenizer
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace boost;

////////////////////////////////////////////////
// Global variables
const std::string w2s_filename = "mcr_source/enWN16";


////////////////////////////////////////////////
// Word2Synset stuff

void read_w2syn_file(const string & fname,
		     map<string, vector<string> > & w2syns) {
    
  // optimize IO
  std::ios::sync_with_stdio(false);
    
  std::ifstream fh(fname.c_str(), ofstream::in);
  if(!fh) {
    cerr << "Can't open " << fname << endl;
    return;
  }

  string line;
  int line_number = 0;
  bool insertedP;

  while(fh) {
    vector<string> fields;
    std::getline(fh, line, '\n');
    line_number++;
    char_separator<char> sep(" ");
    tokenizer<char_separator<char> > tok(line, sep);
    copy(tok.begin(), tok.end(), back_inserter(fields));
    if (fields.size() == 0) continue; // blank line
    if (fields.size() < 2) {
      cerr << "read_w2syn_file error. Bad line: " << line_number << endl;
      exit(-1);
    }
    vector<string>::const_iterator fields_it = fields.begin();
    vector<string>::const_iterator fields_end = fields.end();
    ++fields_it;
    map<string, vector<string> >::iterator map_value_it;
    tie(map_value_it, insertedP) = w2syns.insert(make_pair(fields[0], vector<string>(fields.size() - 1)));
    if (!insertedP) {
      //cerr << "read_w2syn_file warning in line: " << line_number << " " << fields[0] << " element is repeated." <<endl;
      for(;fields_it != fields_end; ++fields_it) 
	map_value_it->second.push_back(*fields_it);
    } else {
      copy(fields_it, fields_end, map_value_it->second.begin());
    }
  }
}
  
W2Syn::W2Syn() {
  read_w2syn_file(w2s_filename, m_w2syns);
}

W2Syn & W2Syn::instance() {
  static W2Syn inst;
  return inst;
}

const vector<std::string> & W2Syn::get_syns(const std::string & word) const {
  static vector<string> nullVector;
  map<string, vector<string> >::const_iterator map_value_it = m_w2syns.find(word);
  if (map_value_it == m_w2syns.end()) return nullVector; // null
  return map_value_it->second;
}

void W2Syn::w2syn_reverse(map<string, string> & rev) const {

  map<string, string>().swap(rev); // empty rev

  map<string, vector<string> >::const_iterator m_it = m_w2syns.begin();
  map<string, vector<string> >::const_iterator m_end = m_w2syns.end();
  for(;m_it != m_end; ++m_it) {
    vector<string>::const_iterator v_it = m_it->second.begin();
    vector<string>::const_iterator v_end =  m_it->second.end();
    for(; v_it != v_end; ++v_it)
      rev.insert(make_pair(*v_it, m_it->first));
  }
}
