
#include "w2syn.h"
#include "globalVars.h"
#include "common.h"

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
//const std::string w2s_filename = "mcr_source/enWN16";

std::ostream & operator<<(std::ostream & o, const W2Syn_item & item) {
  o << "S: ";
  writeV(o, item.wsyns);
  o << "\n";
  if (item.has_freq) {
    o << "F: ";
    writeV(o, item.syns_count);
  } else {
    o << "No freqs";
  }
  o << endl;
  return o;
};


////////////////////////////////////////////////
// Word2Synset stuff

size_t count_lines(const string & fname) {

  size_t N = 0;
  std::ifstream fh(fname.c_str(), ofstream::in);
  if(!fh) {
    cerr << "Can't open " << fname << endl;
    throw;
  }
  // First pass to count total number of words

  while(fh) {
    string line;
    std::getline(fh, line, '\n');
    N++;
  }
  return N;
}

void read_w2syn_file(const string & fname,
		     vector<string> & words,
		     W2Syn::w2syns_t & w2syns) {
    
  // optimize IO
  std::ios::sync_with_stdio(false);

  size_t N = count_lines(fname);
  vector<string>(N).swap(words);
    
  std::ifstream fh(fname.c_str(), ofstream::in);
  if(!fh) {
    cerr << "Can't open " << fname << endl;
    return;
  }

  string line;
  int line_number = 0;
  bool insertedP;
  int words_I = 0;

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
    W2Syn::w2syns_t::iterator map_value_it;

    words[words_I] = fields[0]; // insert word
    tie(map_value_it, insertedP) = w2syns.insert(make_pair(&words[words_I], W2Syn_item()));
    words_I++;

    W2Syn_item & item = map_value_it->second;

    for(; fields_it != fields_end; ++fields_it) {
      char_separator<char> sf_sep(":");
      tokenizer<char_separator<char> > sf_tok(*fields_it, sf_sep);
      vector<string> syn_freq;
      copy(sf_tok.begin(), sf_tok.end(), back_inserter(syn_freq));
      if (syn_freq.size() == 0) {
	cerr << "read_w2syn_file error. Bad line: " << line_number << endl;
	exit(-1);
      }
      item.wsyns.push_back(syn_freq[0]);
      if (syn_freq.size() > 1) {
	item.has_freq = true;
	item.syns_count.push_back(lexical_cast<size_t>(syn_freq[1]));
      } else {
	item.has_freq = false;
      }
    }    
  }
}
  
W2Syn::W2Syn() {
  read_w2syn_file(glVars::w2s_filename, words, m_w2syns);
}

W2Syn & W2Syn::instance() {
  static W2Syn inst;
  return inst;
}

std::pair<vector<std::string>::const_iterator, vector<std::string>::const_iterator>
W2Syn::get_wsyns(const std::string & word) const {
  vector<string>::const_iterator null_it;
  w2syns_t::const_iterator map_value_it = m_w2syns.find(&word);
  if (map_value_it == m_w2syns.end()) return make_pair(null_it, null_it); // null
  return make_pair(map_value_it->second.wsyns.begin(),map_value_it->second.wsyns.end());
}

bool W2Syn::syn_counts(map<string, size_t> & res) const {

  map<string, size_t>().swap(res); // empty res

  w2syns_t::const_iterator m_it = m_w2syns.begin();
  w2syns_t::const_iterator m_end = m_w2syns.end();
  for(;m_it != m_end; ++m_it) {
    const W2Syn_item & item = m_it->second;
    if (!item.has_freq) return false;
    assert(item.wsyns.size() == item.syns_count.size());
    size_t m = item.wsyns.size();
    for(size_t i = 0; i != m; ++i) {
      bool insertedP;
      map<string, size_t>::iterator map_it;
      tie(map_it, insertedP) = res.insert(make_pair(item.wsyns[i], item.syns_count[i]));
      if(!insertedP) {
	map_it->second += item.syns_count[i];
      }
    }
  }
  return true;
}
