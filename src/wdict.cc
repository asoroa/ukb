
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
	  throw std::runtime_error(string("Can't open ") + fname);
	}
	// First pass to count total number of words
	string line;
	size_t line_number;
	while(read_line_noblank(fh, line, line_number)) {
	  N++;
	}
	return N;
  }

  //////////////////////////////////////////////////////////////////////
  // parsing

  static size_t line_number; // global variable, used by many parsing functions

  // given a string with "concept_id:weight", extract "concept_id" and "weight"
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

  // given a concept id "concept-pos", extract pos
  static string xtract_pos_cid(const string & str) {
	std::string::size_type m = str.length();
	std::string::size_type idx = str.find_last_of("-");
	if (idx == string::npos || idx == m - 1)
	  throw ukb::wdict_error("Dictionary concept " + str + " has no POS");
	return str.substr(idx + 1);
  }

  struct concept_parse_t {
	string str;
	Kb_vertex_t u;
	string pos;
	float w;

	concept_parse_t() : str(string()), u(-1), pos(string()), w(0.0f) {}
  };

  struct ccache_map_t {
	set<Kb_vertex_t> S; // concepts already seen
	vector<concept_parse_t> V;
  };

  // given 'cstr' (concept string + weight), set a struct concept_parse_t cp:
  //
  //  cp.str -> the concept without the weight
  //  cp.u   -> vertex id correspondding to concept_str
  //  cp.pos -> the POS (empty if no pos filtering)
  //  cp.w   -> the weight (1 if no weight)

  // return code
  // 0   success
  // 1   Error: concept not in KB
  // 2   Error: parsing error (zero weight)


  static size_t parse_concept(const string & cstr,
							  concept_parse_t & cp) {
	string pos_str("");
	bool aux;
	tie(cp.str, cp.w) = wdict_parse_weight(cstr);

	//  See if concept is in KB
	tie(cp.u, aux) = Kb::instance().get_vertex_by_name(cp.str);
	if (!aux) return 1; // not in KB
	// POS stuff
	if(glVars::input::filter_pos) {
	  cp.pos = xtract_pos_cid(cp.str);
	}
	// Weight stuff
	if (glVars::dict::use_weight) {
	  cp.w += glVars::dict::weight_smoothfactor;
	  if (cp.w == 0.0f) return 2; // zero weight
	}
	return 0;
  }

  // Read a line and parse headword and concept vector
  //
  // Return codes:
  //
  // 0 -> EOF
  // 1 -> OK
  // -1 -> blank line
  // -2 -> ERROR: malformed line

  static int read_wdict_line(istream & fh,
							 vector<string> & fields) {

	string line;
	bool res = read_line_noblank(fh, line, line_number);
	if (!res) return 0; // EOF
	char_separator<char> sep(" \t");
	tokenizer<char_separator<char> > tok(line, sep);
	copy(tok.begin(), tok.end(), back_inserter(fields));
	if (fields.size() == 0) return -1; // blank line
	if (fields.size() < 2) {
	  return -2; // malformed line
	}
	return 1;
  }

  // Fill the concept vector associated with headword hw

  static size_t fill_concepts(const string & hw,
							  vector<string>::const_iterator fields_it,
							  vector<string>::const_iterator fields_end,
							  ccache_map_t & ccache) {

	static const char *concept_err_msg[] = { "(concept not in KB)",
											 "(concept with zero weight)" };
	bool aux;

	for(; fields_it != fields_end; ++fields_it) {
	  concept_parse_t cp;
	  int pc_err_status = parse_concept(*fields_it, cp);
	  if (pc_err_status != 0) {
		// deal with error
		string err_msg(string("line ") + lexical_cast<string>(line_number) + " " + cp.str
					   + " " + concept_err_msg[pc_err_status - 1]);
		if (!glVars::dict::swallow) throw ukb::wdict_error(err_msg + "\n");
		if (glVars::debug::warning) cerr << "[W] read_wdict_file: " + err_msg + " ... ignoring\n";
		continue;
	  }
	  // See if concept was already there
	  set<Kb_vertex_t>::iterator cache_it;
	  tie(cache_it, aux) = ccache.S.insert(cp.u);
	  if (aux) {
		ccache.V.push_back(cp);
	  }
	}
	return ccache.V.size();
  }

  struct pos_order {
	bool operator() (const concept_parse_t & a, const concept_parse_t & b) {
	  return a.pos < b.pos;
	}
  };

  // Given a concept_cache (which is temporary, and is used to remove possible
  // duplicated entries), fill the real dictionary.

  static void create_wdict(map<string, ccache_map_t> & concept_cache,
						   vector<std::string> & wordsV,
						   WDict::wdicts_t & m_wdicts) {

	for(vector<string>::iterator wit = wordsV.begin(), wit_end = wordsV.end();
		wit != wit_end; ++wit) {
	  WDict::wdicts_t::iterator map_value_it = m_wdicts.insert(make_pair(&(*wit), WDict_item_t())).first;
	  WDict_item_t & item = map_value_it->second;

	  map<string, ccache_map_t>::iterator cache_map_it = concept_cache.find(*wit);
	  vector<concept_parse_t> & V = cache_map_it->second.V;
	  sort(V.begin(), V.end(), pos_order());
	  size_t idx = 0;
	  size_t left = 0;
	  string old_pos("");
	  Kb_vertex_t old_concept(-1); // intialization to null vertex taken from boost .hpp
	  float old_w = 0.0;
	  for(size_t i = 0, m = V.size(); i < m; ++i) {
		if (V[i].pos != old_pos) {
		  if (left != idx) {
			item.m_pos_ranges.insert(make_pair(old_pos, wdict_range(left, idx)));
			left = idx;
		  }
		  old_pos = V[i].pos;
		}
		if (V[i].u == old_concept) {
		  if (glVars::debug::warning && V[i].w != old_w)
			cerr << "Warning in entry " + *wit + ": " + V[i].str + " appears twice with different weights. Skipping.\n";
		  continue;
		};
		// concept is different
		old_concept = V[i].u;
		old_w = V[i].w;
		item.m_wsyns.push_back(V[i].u);
		item.m_counts.push_back(V[i].w);
		idx++;
	  }
	  // insert last range
	  if (left != idx) {
		item.m_pos_ranges.insert(make_pair(old_pos, wdict_range(left, idx)));
	  }
	}
  }

  static void read_dictfile_1pass(const string & fname,
								  map<string, ccache_map_t> & concept_cache,
								  std::vector<std::string> & wordsV) {

	std::ifstream fh(fname.c_str(), ofstream::in);
	if(!fh) {
	  cerr << "Error: can't open " << fname << endl;
	  exit(-1);
	}

	// Parse lines of form:
	// word offset-pos:freq offset-pos:freq ...
	// Ex:
	// abandon 04135348-n:4 06081672-n:0 01663408-v:10 00451308-v:7

	string line;
	line_number = 0;
	bool insertedP;

	try {
	  while(true) {
		vector<string> fields;
		int rwl_res = read_wdict_line(fh, fields);
		if (rwl_res == 0) break; // EOF
		if (rwl_res == -1) continue; // blank line
		if (rwl_res == -2) {
		  // parsing error
		  if (!glVars::dict::swallow) throw ukb::wdict_error("Bad line " + lexical_cast<string>(line_number) + ".\n");
		  cerr << "[W] read_dictfile_1pass: line" << line_number <<  " is malformed (ignoring).\n" ;
		  continue;
		}
		vector<string>::const_iterator fields_it = fields.begin();
		const string & hw = *fields_it;
		++fields_it;
		map<string, ccache_map_t>::iterator cache_map_it;
		tie(cache_map_it, insertedP) = concept_cache.insert(make_pair(hw, ccache_map_t()));
		if (insertedP) wordsV.push_back(hw);
		ccache_map_t & ccache = cache_map_it->second;
		size_t concepts_N = fill_concepts(hw, fields_it, fields.end(),
										  ccache);
		if (!concepts_N) {
		  // we have a headword but all the associated concepts are erroneous
		  // Erase the headword from the map
		  if (glVars::debug::warning) {
			cerr << "[W]: line " << lexical_cast<string>(line_number) << ". Ignoring headword " << fields[0] << endl;
		  }
		  wordsV.pop_back();
		  concept_cache.erase(cache_map_it);
		}
	  }
	} catch (std::logic_error & e) {
	  throw ukb::wdict_error("[Error] read_dictfile_1pass: " + string(e.what()));
	} catch (std::exception & e) {
	  throw e; // any other exception is just thrown away
	}
  }


  void WDict::read_wdict_file(const string & fname) {

	map<string, ccache_map_t> concept_cache;

	read_dictfile_1pass(fname, concept_cache, m_words);

	if(m_words.size() == 0)
	  throw ukb::wdict_error("Error reading dict. No headwords linked to KB");

	// Now, create the actual dictionary
	create_wdict(concept_cache, m_words, m_wdicts);
  }

  // void WDict::read_wdict_alt_file(const string & fname) {

  // }

  WDict::WDict() {
	if(glVars::dict_filename.size() == 0)
	  throw std::runtime_error("Error: no dict file\n");
	read_wdict_file(glVars::dict_filename);
  }

  WDict & WDict::instance() {
	static WDict inst;
	return inst;
  }

  static void add_variant_pos(const string & hw,
							  const std::vector<Kb_vertex_t> & wsyns,
							  const string & pos,
							  size_t left,
							  size_t right,
							  std::map<std::string, std::string> & variants) {
	Kb & kb = Kb::instance();

	string hw_prefix = hw + "#";
	if (pos.size()) hw_prefix += pos + "#";
	size_t v_idx = 1;
	for(size_t i = left; i < right; ++i, ++v_idx) {

	  string hw_sense = hw_prefix + lexical_cast<string>(v_idx);
	  string synset_str = kb.get_vertex_name(wsyns[i]);
	  map<string, string>::iterator celem_it = variants.insert(make_pair(synset_str, string())).first;
	  string & variant_str = celem_it->second;
	  if(variant_str.size()) variant_str.append(", ");
	  variant_str.append(hw_sense);
	}
  }

  void WDict::create_variant_map() {

	for(wdicts_t::const_iterator it = m_wdicts.begin(), end = m_wdicts.end();
		it != end; ++it) {
	  const string & hw = *(it->first);
	  const WDict_item_t & elem(it->second);
	  std::map<std::string, wdict_range>::const_iterator rit = elem.m_pos_ranges.begin(), rend = elem.m_pos_ranges.end();
	  if(rit == rend) {
		// no pos
		add_variant_pos(hw, elem.m_wsyns,
						string(""),
						0, elem.m_wsyns.size(),
						m_variants);
	  } else
		for(;rit != rend; ++rit) {
		  add_variant_pos(hw, elem.m_wsyns,
						  rit->first,
						  rit->second.left, rit->second.right,
						  m_variants);
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

  WDict_entries WDict::get_entries(const std::string & word, const string & pos) const {
	static WDict_item_t null_entry;
	wdicts_t::const_iterator map_value_it = m_wdicts.find(&word);
	if (map_value_it == m_wdicts.end()) return WDict_entries(null_entry);
	return WDict_entries(map_value_it->second, pos);
	//return WDict_entries(m_wdicts[&word]);
  }

  //////////////////////////////////////////////////////////////
  // WDict_entries

  WDict_entries::WDict_entries(const WDict_item_t & item)
	: m_item(item), m_pos(std::string()), m_left(0), m_right(item.m_wsyns.size()) {}

  WDict_entries::WDict_entries(const WDict_item_t & item, const std::string & pos)
  : m_item(item), m_pos(pos), m_left(0), m_right(0) {
	if (!glVars::input::filter_pos || !pos.size()) {
	  string().swap(m_pos);
	  m_right = item.m_wsyns.size();
	} else {
	  std::map<std::string, wdict_range>::const_iterator it = item.m_pos_ranges.find(pos);
	  if (it != item.m_pos_ranges.end()) {
		m_left = it->second.left;
		m_right = it->second.right;
	  }
	}
  }

  size_t WDict_entries::size() const {
	return m_right - m_left;
  }

  Kb_vertex_t WDict_entries::get_entry(size_t i) const {
	return m_item.m_wsyns[i + m_left];
  }

  const std::string & WDict_entries::get_entry_str(size_t i) const {
	return Kb::instance().get_vertex_name(m_item.m_wsyns[i + m_left]);
  }

  const std::string & WDict_entries::get_pos(size_t i) const {
	return m_pos;
  }

  float WDict_entries::get_freq(size_t i) const {
	if (!glVars::dict::use_weight) return 1.0;
	return m_item.m_counts[i + m_left];
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
