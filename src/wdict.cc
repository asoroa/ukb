
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
	for(size_t i = 0, n = item.m_wsyns.size(); i != n; ++i) {
	  o << " " << item.m_wsyns[i];
	  if (glVars::dict::use_weight)  o << ":" << item.m_counts[i];
	}
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


  struct pw_pair_t {
	char p;
	float w;

	pw_pair_t() : p(0), w(0.0f) {}
	pw_pair_t(char pp, float ww) : p(pp), w(ww) {}

  };

  typedef std::map<Kb_vertex_t, pw_pair_t> ccache_map_t;   // conceptId -> weight

  bool parse_concept(const string & hw, const string & cstr,
					 ccache_map_t & ccache) {
	float weight;
	string concept_str;
	Kb_vertex_t concept_id;
	char pos_c = 0;
	bool aux;
	tie(concept_str, weight) = wdict_parse_weight(cstr);

	tie(concept_id, aux) = Kb::instance().get_vertex_by_name(concept_str);
	if (!aux)
	  return false;

	// POS stuff
	if(glVars::input::filter_pos) {
	  pos_c = xtract_pos_cid(concept_str);
	}
	// Weight stuff
	if (glVars::dict::use_weight) {
	  weight += glVars::dict::weight_smoothfactor;
	  if (weight == 0.0)
		throw std::runtime_error ("Error in entry " + hw + ": " + cstr + " word has zero weight.");
	}

	// See if concept was already there
	ccache_map_t::iterator cache_it = ccache.find(concept_id);
	if (cache_it != ccache.end()) {
	  // If not a new concept, see if previous concept had the same weight
	  if (glVars::debug::warning && weight != cache_it->second.w)
		cerr << "Warning in entry " + hw + ": " + concept_str + " appears twice with different weights. Skipping.\n";
	  return false;
	}
	// update cache with new concept
	ccache.insert(make_pair(concept_id, pw_pair_t(pos_c, weight)));
	return true;
  }

  void WDict::read_wdict_file(const string & fname) {

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

		map<string, ccache_map_t>::iterator cache_map_it;
		tie(cache_map_it, insertedP) = concept_cache.insert(make_pair(fields[0], ccache_map_t()));
		if (insertedP) m_words.push_back(fields[0]);
		ccache_map_t & ccache = cache_map_it->second;
		size_t inserted_concepts_N = ccache.size();
		for(; fields_it != fields_end; ++fields_it) {
		  try {
			if (parse_concept(fields[0], *fields_it, ccache)) {
			  inserted_concepts_N++;
			}
		  } catch (std::runtime_error & e) {
			if (!glVars::dict::swallow) throw e;
			if (glVars::debug::warning) {
			  cerr << "Wdict: line " << lexical_cast<string>(line_number) << ": " << e.what() << "(ignoring).\n";
			}
		  }
		}
		if (!inserted_concepts_N) {
		  // we have a headword but all the associated concepts were erroneous
		  // Erase the headword from the map
		  if (glVars::debug::warning) {
			cerr << "Wdict: line " << lexical_cast<string>(line_number) << ". Ignoring headword " << fields[0] << endl;
		  }
		  m_words.pop_back();
		  concept_cache.erase(cache_map_it);
		}
	  }
	} catch (std::exception & e) {
	  throw std::runtime_error("Error in read_wdict_file, line " + lexical_cast<string>(line_number) + "\n" + string(e.what()));
	}

	// Now fill actual dictionary
	if(m_words.size() == 0)
	  throw std::runtime_error("Error reading dict. No headwords linked to KB");
	for(vector<string>::iterator wit = m_words.begin(), wit_end = m_words.end();
		wit != wit_end; ++wit) {
	  WDict::wdicts_t::iterator map_value_it = m_wdicts.insert(make_pair(&(*wit), WDict_item_t())).first;
	  WDict_item_t & item = map_value_it->second;
	  map<string, ccache_map_t>::iterator cache_map_it = concept_cache.find(*wit);
	  for(ccache_map_t::iterator dc_it = cache_map_it->second.begin(), dc_end = cache_map_it->second.end();
	  	  dc_it != dc_end; ++dc_it) {
	  	item.m_wsyns.push_back(dc_it->first);
	  	item.m_counts.push_back(dc_it->second.w);
	  	item.m_thepos.push_back(dc_it->second.p);
	  }
	}
  }

  WDict::WDict() {
	if(glVars::dict_filename.size() == 0)
	  throw std::runtime_error("Error: no dict file\n");
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

		string synset_str = WDict_entries(elem).get_entry_str(i);
		map<string, string>::iterator celem_it = m_variants.insert(make_pair(synset_str, string())).first;
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

  //////////////////////////////////////////////////////////////
  // WDict_entries

  const std::string & WDict_entries::get_entry_str(size_t i) const {
	return Kb::instance().get_vertex_name(_item.m_wsyns[i]);
  }

  char WDict_entries::get_pos(size_t i) const {
	if (!glVars::input::filter_pos) return 0;
	return _item.m_thepos[i];
  }

  float WDict_entries::get_freq(size_t i) const {
	if (!glVars::dict::use_weight) return 1.0;
	return _item.m_counts[i];
  }


  std::ostream & operator<<(std::ostream & o, const WDict & dict) {

	for(vector<string>::const_iterator it = dict.m_words.begin(), end = dict.m_words.end();
		it != end; ++it) {
	  WDict::wdicts_t::const_iterator s_it = dict.m_wdicts.find(&(*it));
	  if (s_it == dict.m_wdicts.end()) {
		s_it = dict.m_wdicts.end();
	  }
	  assert(s_it != dict.m_wdicts.end());
	  o << *it << s_it->second << "\n";
	}
	return o;
  };


}
