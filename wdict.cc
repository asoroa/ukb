
#include "wdict.h"
#include "globalVars.h"
#include "common.h"

#include <fstream>
#include <iostream>

#include<boost/tuple/tuple.hpp> // for "tie"

// Tokenizer
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

namespace ukb {

  using namespace std;
  using namespace boost;

  ////////////////////////////////////////////////
  // Global variables
  //const std::string dict_filename = "kb_source/enWN16";

  std::ostream & operator<<(std::ostream & o, const WDict_item_t & item) {
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

  //////////////////////////////////////////////////////7
  // join

  std::string join(const std::string &delim,
				   std::vector<std::string>::const_iterator it,
				   std::vector<std::string>::const_iterator end) {

	std::string res("");

	if (it == end) return res;
	--end;
	for(;it != end; ++it) {
	  res+=*it;
	  res+=delim;
	}
	res+=*end;
	return res;
  }


  size_t count_lines(const string & fname) {

	size_t N = 0;
	std::ifstream fh(fname.c_str(), ofstream::in);
	if(!fh) {
	  cerr << "Can't open " << fname << endl;
	  throw;
	}
	// First pass to count total number of words
	string line;
	size_t line_number;
	while(read_line_noblank(fh, line, line_number)) {
	  N++;
	}
	return N;
  }


  static pair<string, float> wdict_parse_weight(const string & str) {

	bool has_w = false;
	float weight = 0.0f;
	string concept_id(str);

	char_separator<char> sf_sep("", ":"); // keep delimiters
	tokenizer<char_separator<char> > sf_tok(str, sf_sep);
	vector<string> syn_freq;
	copy(sf_tok.begin(), sf_tok.end(), back_inserter(syn_freq));

	size_t m = syn_freq.size();

	if (m > 2 ) {
	  // Warning. concept-id can have ":" characters in. So, just
	  // take the last field as weight, and join the rest.
	  try {
		weight = lexical_cast<float>(syn_freq[m - 1]);
		has_w = true;
	  } catch (boost::bad_lexical_cast &) {
		// last field wasn't a float. Do nothing.
	  }
	  if (has_w) {
		if (m == 3) {
		  concept_id = syn_freq[0];
		} else {
		  concept_id = join("", syn_freq.begin(), syn_freq.end() - 2);
		}
	  } else {
		concept_id = join("", syn_freq.begin(), syn_freq.end());
	  }
	}
	return make_pair(concept_id, weight);
  }

  void read_wdict_file(const string & fname,
					   vector<string> & words,
					   WDict::wdicts_t & wdicts) {

	// optimize IO
	std::ios::sync_with_stdio(false);

	size_t N = count_lines(fname);
	vector<string>(N).swap(words);

	std::ifstream fh(fname.c_str(), ofstream::in);
	if(!fh) {
	  cerr << "Can't open " << fname << endl;
	  return;
	}

	// Parse lines of form:
	// word offset-pos:freq offset-pos:freq ...
	// Ex:
	// abandon 04135348-n:4 06081672-n:0 01663408-v:10 00451308-v:7

	string line;
	size_t line_number = 0;
	bool insertedP;
	int words_I = 0;

	try {
	  while(read_line_noblank(fh, line, line_number)) {

		vector<string> fields;

		char_separator<char> sep(" \t");
		tokenizer<char_separator<char> > tok(line, sep);
		copy(tok.begin(), tok.end(), back_inserter(fields));
		if (fields.size() == 0) continue; // blank line
		if (fields.size() < 2) {
		  cerr << "read_wdict_file error. Bad line: " << line_number << endl;
		  exit(-1);
		}
		vector<string>::const_iterator fields_it = fields.begin();
		vector<string>::const_iterator fields_end = fields.end();
		++fields_it;
		WDict::wdicts_t::iterator map_value_it;

		words[words_I] = fields[0]; // insert word
		tie(map_value_it, insertedP) = wdicts.insert(make_pair(&words[words_I], WDict_item_t()));
		if (insertedP) words_I++;

		WDict_item_t & item = map_value_it->second;

		for(; fields_it != fields_end; ++fields_it) {

		  float weight;
		  string concept_id;
		  tie(concept_id, weight) = wdict_parse_weight(*fields_it);
		  item.wsyns.push_back(concept_id);
		  if (glVars::dict::use_weight) {
			weight += glVars::dict::weight_smoothfactor;
			if (weight == 0.0)
			  throw std::runtime_error ("Error in entry " + fields[0] + ": " + *fields_it + " word has zero weight.");
			item.has_freq = 1;
			item.syns_count.push_back(weight);
		  }
		}
	  }
	} catch (std::exception & e) {
	  throw std::runtime_error("Error in read_wdict_file: " + string(e.what()) + " in line " + lexical_cast<string>(line_number));
	}
  }

  WDict::WDict() {
	read_wdict_file(glVars::dict_filename, words, m_wdicts);
  }

  WDict & WDict::instance() {
	static WDict inst;
	return inst;
  }

  WDict_entries WDict::get_entries(const std::string & word) const {
	static WDict_item_t null_entry;
	wdicts_t::const_iterator map_value_it = m_wdicts.find(&word);
	if (map_value_it == m_wdicts.end()) return WDict_entries(null_entry);
	return WDict_entries(map_value_it->second);
	//return WDict_entries(m_wdicts[&word]);
  }

  std::pair<vector<std::string>::const_iterator, vector<std::string>::const_iterator>
  WDict::get_wsyns(const std::string & word) const {
	vector<string>::const_iterator null_it;
	wdicts_t::const_iterator map_value_it = m_wdicts.find(&word);
	if (map_value_it == m_wdicts.end()) return make_pair(null_it, null_it); // null
	return make_pair(map_value_it->second.wsyns.begin(),map_value_it->second.wsyns.end());
  }


  std::pair<vector<float>::const_iterator, vector<float>::const_iterator>
  WDict::get_weights(const std::string & word) const {
	vector<float>::const_iterator null_it;
	wdicts_t::const_iterator map_value_it = m_wdicts.find(&word);
	if (map_value_it == m_wdicts.end()) return make_pair(null_it, null_it); // null
	return make_pair(map_value_it->second.syns_count.begin(),map_value_it->second.syns_count.end());
  };


  bool WDict::syn_counts(map<string, size_t> & res) const {

	map<string, size_t>().swap(res); // empty res

	wdicts_t::const_iterator m_it = m_wdicts.begin();
	wdicts_t::const_iterator m_end = m_wdicts.end();
	for(;m_it != m_end; ++m_it) {
	  const WDict_item_t & item = m_it->second;
	  if (!item.has_freq) return false;
	  assert(item.wsyns.size() == item.syns_count.size());
	  size_t m = item.wsyns.size();
	  for(size_t i = 0; i != m; ++i) {
		bool insertedP;
		map<string, size_t>::iterator map_it;
		size_t synset_n = lexical_cast<size_t>(item.syns_count[i]);
		tie(map_it, insertedP) = res.insert(make_pair(item.wsyns[i], synset_n));
		if(!insertedP) {
		  map_it->second += synset_n;
		}
	  }
	}
	return true;
  }

  //////////////////////////////////////////////////////////////
  // WDict_entries

  char WDict_entries::get_pos(size_t i) const {
	if (!glVars::input::filter_pos) return 0;
	std::string::size_type idx = _item.wsyns[i].find_last_of("-");
	if (idx == string::npos || idx == _item.wsyns[i].length() - 1)
	  throw std::runtime_error("Dictionary concept " + _item.wsyns[i] + " has no POS\n");
	return _item.wsyns[i].at(idx + 1);
  }

  float WDict_entries::get_freq(size_t i) const {
	if (_item.has_freq)
	  return _item.syns_count[i];
	return 1.0;
  }

}
