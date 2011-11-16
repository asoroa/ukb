
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
	writeV(o, item.m_wsyns);
	o << "\n";
	o << "F: ";
	writeV(o, item.m_counts);
	o << endl;
	return o;
  };


  ////////////////////////////////////////////////
  // Word2Synset stuff

  //////////////////////////////////////////////////////7
  // join


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

	float weight = 0.0f; // default weight is zero (unless glVars::dict:use_weight is false (see below))
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
		float tmp = lexical_cast<float>(syn_freq[m - 1]);
		weight = tmp;
		// reconstruct concept_id without the weight
		if (m == 3)
		  concept_id = syn_freq[0];
		else
		  concept_id = join("", syn_freq.begin(), syn_freq.end() - 2);
	  } catch (boost::bad_lexical_cast &) {
		// last field wasn't a float. Do nothing.
	  }
	}
	// if glVars::dict:use_weight is false, set weight to 1 regardless
	if (!glVars::dict::use_weight) weight = 1.0f;
	return make_pair(concept_id, weight);
  }

  char xtract_pos_cid(const string & str) {
	std::string::size_type m = str.length();
	std::string::size_type idx = str.find_last_of("-");
	if (idx == string::npos || idx == m - 1)
	  throw std::runtime_error("Dictionary concept " + str + " has no POS");
	if (m - idx > 2)
	  throw std::runtime_error("Dictionary concept " + str + " has invalid POS (more than 1 char).");
	return str.at(idx + 1);
  }


  // parse a new concep tand insert it into WDict_item_t (if previously not there)
  //
  // return true if concept (as expressed in concept_str) is new, false otherwise
  // throw if concept_str syntax is invalid

  bool WDict_item_t::parse_concept(const string & hw, const string & concept_str,
								   ccache_map_t & ccache) {
	float weight;
	string concept_id;
	char pos_c = 0;
	tie(concept_id, weight) = wdict_parse_weight(concept_str);

	// POS stuff
	if(glVars::input::filter_pos) {
	  pos_c = xtract_pos_cid(concept_id);
	}
	// Weight stuff
	if (glVars::dict::use_weight) {
	  weight += glVars::dict::weight_smoothfactor;
	  if (weight == 0.0)
		throw std::runtime_error ("Error in entry " + hw + ": " + concept_str + " word has zero weight.");
	}

	// See if concept was already there
	ccache_map_t::iterator cache_it = ccache.find(concept_id);
	if (cache_it != ccache.end()) {
	  // If not a new concept, see if previous concept had the same weight
	  if (glVars::debug::warning && weight != cache_it->second)
		cerr << "Warning in entry " + hw + ": " + concept_id + " appears twice with different weights. Skipping.\n";
	  return false;
	}

	// new concept
	m_wsyns.push_back(concept_id);
	m_thepos.push_back(pos_c);
	m_counts.push_back(weight);
	// update cache
	ccache.insert(make_pair(concept_id, weight));
	return true;
  }

  void WDict::read_wdict_file(const string & fname) {

	size_t N = count_lines(fname);
	vector<string>(N).swap(m_words);

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

	map<string, ccache_map_t> concept_cache;

	try {
	  while(read_line_noblank(fh, line, line_number)) {

		vector<string> fields;
		char_separator<char> sep(" \t");
		tokenizer<char_separator<char> > tok(line, sep);
		copy(tok.begin(), tok.end(), back_inserter(fields));
		if (fields.size() == 0) continue; // blank line
		if (fields.size() < 2) {
		  if (!glVars::dict::swallow) throw std::runtime_error("Bad line.\n");
		  cerr << "Wdict: line " << line_number <<  ": Bad line (ignoring).\n" ;
		  continue;
		}
		vector<string>::const_iterator fields_it = fields.begin();
		vector<string>::const_iterator fields_end = fields.end();
		++fields_it;
		WDict::wdicts_t::iterator map_value_it;

		m_words[words_I] = fields[0]; // insert word
		tie(map_value_it, insertedP) = m_wdicts.insert(make_pair(&m_words[words_I], WDict_item_t()));
		if (insertedP) words_I++;

		WDict_item_t & item = map_value_it->second;

		// Map to track which concepts are already stored for this headword
		map<string, ccache_map_t>::iterator cache_map_it;
		cache_map_it = concept_cache.insert(make_pair(fields[0], ccache_map_t())).first;
		ccache_map_t & ccache = cache_map_it->second;
		size_t inserted_concepts_N = ccache.size();
		for(; fields_it != fields_end; ++fields_it) {
		  try {
			if (item.parse_concept(fields[0], *fields_it, ccache)) {
			  inserted_concepts_N++;
			}
		  } catch (std::runtime_error & e) {
			if (glVars::dict::swallow) {
			  if (glVars::debug::warning) {
				cerr << "Wdict: line " << lexical_cast<string>(line_number) << ": " << e.what() << "(ignoring).\n";
			  }
			} else {
			  throw e;
			}
		  }
		}
		if (!inserted_concepts_N) {
		  // we have a headword but all the associated concepts were erroneous
		  // Erase the headword from the map
		  if (glVars::debug::warning) {
			cerr << "Wdict: line " << lexical_cast<string>(line_number) << ". Ignoring headword " << fields[0] << endl;
		  }
		  --words_I;
		  concept_cache.erase(cache_map_it);
		}
	  }
	} catch (std::exception & e) {
	  throw std::runtime_error("Error in read_wdict_file, line " + lexical_cast<string>(line_number) + "\n" + string(e.what()));
	}
  }

  WDict::WDict() {
	read_wdict_file(glVars::dict_filename);
  }

  WDict & WDict::instance() {
	static WDict inst;
	return inst;
  }

  void WDict::create_variant_map() {

	for(wdicts_t::const_iterator it = m_wdicts.begin(), end = m_wdicts.end();
		it != end; ++it) {
	  const string & hw = *(it->first);
	  const WDict_item_t & elem(it->second);
	  for(size_t i = 0, m = elem.m_wsyns.size(); i < m; ++i) {

		string hw_sense = hw + "#" + lexical_cast<string>(i + 1);

		map<string, string>::iterator celem_it = m_variants.insert(make_pair(elem.m_wsyns[i], string())).first;
		string & variant_str = celem_it->second;
		string comma = variant_str.size() ? string(", ") : string();
		variant_str.append(comma);
		variant_str.append(hw_sense);
	  }
	}
  }

  std::string WDict::variant(std::string & concept_id) const {

	static string res("Not in Dictionary");
	if (m_variants.size() == 0) {
	  // Remove const'ness
	  WDict & me = const_cast<WDict &>(*this);
	  me.create_variant_map();
	}
	map<string, string>::const_iterator it = m_variants.find(concept_id);
	if (it == m_variants.end())
	  return res;
	return it->second;
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
	return make_pair(map_value_it->second.m_wsyns.begin(),map_value_it->second.m_wsyns.end());
  }


  std::pair<vector<float>::const_iterator, vector<float>::const_iterator>
  WDict::get_weights(const std::string & word) const {
	vector<float>::const_iterator null_it;
	wdicts_t::const_iterator map_value_it = m_wdicts.find(&word);
	if (map_value_it == m_wdicts.end()) return make_pair(null_it, null_it); // null
	return make_pair(map_value_it->second.m_counts.begin(),map_value_it->second.m_counts.end());
  };


  bool WDict::syn_counts(map<string, size_t> & res) const {

	map<string, size_t>().swap(res); // empty res

	wdicts_t::const_iterator m_it = m_wdicts.begin();
	wdicts_t::const_iterator m_end = m_wdicts.end();
	for(;m_it != m_end; ++m_it) {
	  const WDict_item_t & item = m_it->second;

	  assert(item.m_wsyns.size() == item.m_counts.size());
	  size_t m = item.m_wsyns.size();
	  for(size_t i = 0; i != m; ++i) {
		bool insertedP;
		map<string, size_t>::iterator map_it;
		size_t synset_n = lexical_cast<size_t>(item.m_counts[i]);
		tie(map_it, insertedP) = res.insert(make_pair(item.m_wsyns[i], synset_n));
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
	return _item.m_thepos[i];
  }

  float WDict_entries::get_freq(size_t i) const {
	if (!glVars::dict::use_weight) return 1.0;
	return _item.m_counts[i];
  }
}
