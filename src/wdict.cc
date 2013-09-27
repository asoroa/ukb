
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
  //const std::string text_fname = "kb_source/enWN16";

  std::ostream & operator<<(std::ostream & o, const wdict_rhs_t & rhs) {
	Kb & kb = Kb::instance();
	for(size_t i = 0, n = rhs.m_items.size(); i != n; ++i) {
	  o << " " << kb.get_vertex_name(rhs.m_items[i].m_syn);
	  if (glVars::dict::use_weight)  o << ":" << rhs.m_items[i].m_count;
	}
	return o;
  };

  std::ostream & operator<<(std::ostream & o, const WDict_entries & entry) {
	o << entry.m_rhs;
	return o;
  }

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

	float weight = 0.0f; // default weight is zero
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

	concept_parse_t() : str(string()), u(-1), pos(""), w(0.0f) {}
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
	} else {
	  // if glVars::dict:use_weight is false, set weight to 1 regardless
	  cp.w = 1.0f;
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


  static void fill_wdict_hw(wdict_rhs_t & rhs,
							vector<concept_parse_t> & V) {
	if(glVars::input::filter_pos)
	  sort(V.begin(), V.end(), pos_order());
	vector<wdict_item_t> items;
	vector<wdict_range_t> pos_ranges;
	set<Kb_vertex_t> U;
	size_t idx = 0;
	size_t left = 0;
	string old_pos("");
	for(size_t i = 0, m = V.size(); i < m; ++i) {
	  if (!U.insert(V[i].u).second) continue; // Concept previously there
	  if (V[i].pos != old_pos) {
		if (left != idx) {
		  pos_ranges.push_back(wdict_range_t(old_pos, left, idx));
		  left = idx;
		}
		old_pos = V[i].pos;
	  }
	  items.push_back(wdict_item_t(V[i].u, V[i].w));
	  idx++;
	}
	// insert last range
	if (left != idx) {
	  pos_ranges.push_back(wdict_range_t(old_pos, left, idx));
	}
	// swap vectors from temporary, so they don't take more space than neccesary
	vector<wdict_item_t> (items).swap(rhs.m_items);
	if(glVars::input::filter_pos) {
	  vector<wdict_range_t> (pos_ranges).swap(rhs.m_pos_ranges);
	}
  }

  static size_t read_dictfile_1pass(const string & fname,
									map<string, ccache_map_t> & concept_cache) {

	std::ifstream fh(fname.c_str(), ofstream::in);
	if(!fh) {
	  cerr << "Error: can't open " << fname << endl;
	  exit(-1);
	}
	size_t N = 0;

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
		if (insertedP) N++;
		ccache_map_t & ccache = cache_map_it->second;
		size_t concepts_N = fill_concepts(hw, fields_it, fields.end(),
										  ccache);
		if (!concepts_N) {
		  // we have a headword but all the associated concepts are erroneous
		  // Erase the headword from the map
		  if (glVars::debug::warning) {
			cerr << "[W]: line " << lexical_cast<string>(line_number) << ". Ignoring headword " << fields[0] << endl;
		  }
		  N--;
		  concept_cache.erase(cache_map_it);
		}
	  }
	} catch (std::logic_error & e) {
	  throw ukb::wdict_error("[Error] read_dictfile_1pass: " + string(e.what()));
	} catch (std::exception & e) {
	  throw e; // any other exception is just thrown away
	}
	return N;
  }


  void WDict::read_wdict_file(const string & fname) {

	map<string, ccache_map_t> concept_cache;

	m_N = read_dictfile_1pass(fname, concept_cache);
	if(m_N == 0)
	  throw ukb::wdict_error("Error reading dict. No headwords linked to KB");

	// Now, create the actual dictionary given a concept_cache (which is
	// temporary, and is used to remove possible duplicated entries)

	for (map<string, ccache_map_t>::iterator it = concept_cache.begin(), end = concept_cache.end();
		 it != end; ++it) {
	  WDict::wdicts_t::iterator map_value_it = m_wdicts.insert(make_pair(it->first, wdict_rhs_t())).first;
	  fill_wdict_hw(map_value_it->second, it->second.V);
	}
  }

  void WDict::read_alternate_file(const string & fname) {

	map<string, ccache_map_t> concept_cache;
	bool insertedP;

	read_dictfile_1pass(fname, concept_cache);

	for (map<string, ccache_map_t>::iterator it = concept_cache.begin(), end = concept_cache.end();
		it != end; ++it) {
	  wdicts_t::iterator mit;
	  wdict_rhs_t empty_rhs;
	  tie (mit, insertedP) = m_wdicts.insert(make_pair(it->first, empty_rhs));
	  if (!insertedP) {
		// already in map, so empty existing concepts (will be replaced by alternate dict)
		mit->second.swap(empty_rhs);
	  } else {
		m_N++; // New headword
	  }
	  fill_wdict_hw(mit->second, it->second.V);
	}
  }

  WDict::WDict() : m_N(0) {
	if(glVars::dict::text_fname.size() == 0 and glVars::dict::bin_fname.size() == 0)
	  throw std::runtime_error("[E] WDict: no dict file\n");
	if (glVars::dict::text_fname.size()) read_wdict_file(glVars::dict::text_fname);
	else {
	  read_wdict_binfile(glVars::dict::bin_fname);
	}
  }

  WDict & WDict::instance() {
	static WDict inst;
	return inst;
  }

  size_t WDict::size() {
	return m_N;
  }

  void WDict::size_bytes() {
	long D = 0;
	long C = 0;
	long O_str = 0; // overhead
	long O_rhs = 0; // overhead
	long O_items = 0; // overhead
	long O_ranges = 0; // overhead
	long O_hash = 0; // overhead
	long V = 0;
	O_hash += sizeof(m_wdicts);
	for (WDict::wdicts_t::const_iterator it = m_wdicts.begin(), end = m_wdicts.end();
		 it != end; ++it) {
	  O_hash += sizeof(*it);
	  O_str += sizeof(it->first);
	  D += it->first.size();
	  C += it->first.capacity();
	  const wdict_rhs_t & rhs = it->second;
	  O_rhs += sizeof(rhs);
	  D += rhs.m_items.size() * sizeof(Kb_vertex_t);
	  C += rhs.m_items.capacity() * sizeof(Kb_vertex_t);
	  O_items += sizeof(rhs.m_items);
	  O_ranges += sizeof(rhs.m_pos_ranges);
	  for(std::vector<wdict_range_t>::const_iterator pit = rhs.m_pos_ranges.begin();
		  pit != rhs.m_pos_ranges.end(); ++pit) {
		O_ranges += sizeof(*pit);
		D += pit->pos.size();
		C += pit->pos.capacity();
		D += sizeof(pit->left);
		D += sizeof(pit->right);
	  }
	  C += (rhs.m_pos_ranges.capacity() - rhs.m_pos_ranges.size()) * sizeof(wdict_range_t);
	}
	for(map<std::string, std::string>::const_iterator it = m_variants.begin(), end = m_variants.end();
		 it != end; ++it) {
	  V += it->first.size();
	  V += it->second.size();
	}
	cout << "Dict: " << D << " " << "Capacity: " << C << " Overhead: " << O_hash + O_str + O_rhs + O_items + O_ranges << "\n";
	cout << "O_hash " << O_hash << " ";
	cout << "O_rhs " << O_rhs << " ";
	cout << "O_str " << O_rhs << " ";
	cout << "O_items " << O_items << " ";
	cout << "O_ranges " << O_ranges << "\n";

  }


  static void add_variant_pos(const string & hw,
							  const std::vector<wdict_item_t> & items,
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
	  string synset_str = kb.get_vertex_name(items[i].m_syn);
	  map<string, string>::iterator celem_it = variants.insert(make_pair(synset_str, string())).first;
	  string & variant_str = celem_it->second;
	  if(variant_str.size()) variant_str.append(", ");
	  variant_str.append(hw_sense);
	}
  }

  void WDict::create_variant_map() {

	for(wdicts_t::const_iterator it = m_wdicts.begin(), end = m_wdicts.end();
		it != end; ++it) {
	  const string & hw = it->first;
	  const wdict_rhs_t & rhs(it->second);
	  std::vector<wdict_range_t>::const_iterator rit = rhs.m_pos_ranges.begin(), rend = rhs.m_pos_ranges.end();
	  if(rit == rend) {
		// no pos
		add_variant_pos(hw, rhs.m_items,
						string(""),
						0, rhs.m_items.size(),
						m_variants);
	  } else
		for(;rit != rend; ++rit) {
		  add_variant_pos(hw, rhs.m_items,
						  rit->pos,
						  rit->left, rit->right,
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
	static wdict_rhs_t null_entry;
	wdicts_t::const_iterator map_value_it = m_wdicts.find(word);
	if (map_value_it == m_wdicts.end()) return WDict_entries(null_entry);
	return WDict_entries(map_value_it->second, pos);
	//return WDict_entries(m_wdicts[&word]);
  }

  //////////////////////////////////////////////////////////////
  // WDict_entries

  WDict_entries::WDict_entries(const wdict_rhs_t & rhs)
	: m_rhs(rhs), m_pos(std::string()), m_left(0), m_right(rhs.m_items.size()) {}


  struct wdict_range_pos_P {
	bool operator() (const wdict_range_t & a) {
	  return a.pos == m_p;
	}
	wdict_range_pos_P(const string & p) : m_p(p) {}
	const string & m_p;
  };

  WDict_entries::WDict_entries(const wdict_rhs_t & rhs, const std::string & pos)
  : m_rhs(rhs), m_pos(pos), m_left(0), m_right(0) {
	if (!glVars::input::filter_pos || !pos.size()) {
	  string().swap(m_pos);
	  m_right = rhs.m_items.size();
	} else {
	  std::vector<wdict_range_t>::const_iterator end = rhs.m_pos_ranges.end();
	  std::vector<wdict_range_t>::const_iterator it = std::find_if(rhs.m_pos_ranges.begin(), end, wdict_range_pos_P(pos));
	  if (it != end) {
		m_left = it->left;
		m_right = it->right;
	  }
	}
  }

  size_t WDict_entries::size() const {
	return m_right - m_left;
  }

  Kb_vertex_t WDict_entries::get_entry(size_t i) const {
	return m_rhs.m_items[i + m_left].m_syn;
  }

  const std::string & WDict_entries::get_entry_str(size_t i) const {
	return Kb::instance().get_vertex_name(this->get_entry(i));
  }

  const std::string & WDict_entries::get_pos(size_t i) const {
	return m_pos;
  }

  float WDict_entries::get_freq(size_t i) const {
	if (!glVars::dict::use_weight) return 1.0;
	return m_rhs.m_items[i + m_left].m_count;
  }

  std::ostream & operator<<(std::ostream & o, const WDict & dict) {

	for (WDict::wdicts_t::const_iterator it = dict.m_wdicts.begin(), end = dict.m_wdicts.end();
		 it != end; ++it) {
	  o << it->first << it->second << "\n";
	}
	return o;
  };

  ////////////////////////////////////////////////////////////////////////////////
  // Streaming

  static const size_t magic_id_dict_v0 = 0x130926;


  void WDict::read_wdict_binfile(const string & fname) {

	ifstream fi(fname.c_str(), ifstream::binary|ifstream::in);
	if (!fi) {
	  throw runtime_error("[E] reading serialized dictionary: Can not open " + fname);
	}
	read_dict_from_stream(fi);
  }

  void WDict::write_wdict_binfile(const string & fname) {

	ofstream fo(fname.c_str(),  ofstream::binary|ofstream::out);
	if (!fo) {
	  throw runtime_error("[E] writing serialized dict: Can not create " + fname);
	}
	write_dict_to_stream(fo);
  }

  ostream & write_posRangeM_to_stream (std::ostream & os, const std::vector<wdict_range_t> & pr) {

	size_t m = pr.size();
	os.write(reinterpret_cast<const char *>(&m), sizeof(m));
	if(m) {
	  for(std::vector<wdict_range_t>::const_iterator it = pr.begin(), end = pr.end();
		  it != end; ++it) {
		write_atom_to_stream(os, it->pos);
		write_atom_to_stream(os, it->left);
		write_atom_to_stream(os, it->right);
	  }
	}
	return os;
  }

  void read_posRangeM_from_stream (std::istream & is, std::vector<wdict_range_t> & pr) {
	string hw;
	size_t m;
	vector<wdict_range_t> auxV;

	read_atom_from_stream(is, m);
	if (!is) return;
	for (size_t i = 0; i < m; ++i) {
	  wdict_range_t r;
	  read_atom_from_stream(is, r.pos);
	  read_atom_from_stream(is, r.left);
	  read_atom_from_stream(is, r.right);
	  auxV.push_back(r);
	}
	vector<wdict_range_t> (auxV).swap(pr);
  }

  ostream & write_itemV_to_stream (std::ostream & os, const vector<wdict_item_t> & vi) {
	size_t m = vi.size();
	os.write(reinterpret_cast<const char *>(&m), sizeof(m));
	if(m) {
	  for(std::vector<wdict_item_t>::const_iterator it = vi.begin(), end = vi.end();
		  it != end; ++it) {
		write_atom_to_stream(os, it->m_syn);
		write_atom_to_stream(os, it->m_count);
	  }
	}
	return os;
  }

  void read_itemV_from_stream (std::istream & is, vector<wdict_item_t> & vi) {

	size_t m;
	vector<wdict_item_t> auxV;

	read_atom_from_stream(is, m);
	if (!is) return;
	for (size_t i = 0; i < m; ++i) {
	  wdict_item_t i;
	  read_atom_from_stream(is, i.m_syn);
	  read_atom_from_stream(is, i.m_count);
	  auxV.push_back(i);
	}
	vector<wdict_item_t> (auxV).swap(vi);
  }

  ostream & write_rhs_to_stream (std::ostream & os, const wdict_rhs_t & rhs) {

	write_itemV_to_stream(os, rhs.m_items);
	write_posRangeM_to_stream(os, rhs.m_pos_ranges);

	return os;
  }

  void read_rhs_from_stream (std::istream & is, wdict_rhs_t & rhs) {

	read_itemV_from_stream(is, rhs.m_items);
	read_posRangeM_from_stream(is, rhs.m_pos_ranges);
  }

  ostream & WDict::write_dict_to_stream (std::ostream & os) const {

	write_atom_to_stream(os, magic_id_dict_v0);
	// Write map
	os.write(reinterpret_cast<const char *>(&m_N), sizeof(m_N));
	if(m_N) {
	  for(wdicts_t::const_iterator it = m_wdicts.begin(), end = m_wdicts.end();
		  it != end; ++it) {
		write_atom_to_stream(os, it->first);
		write_rhs_to_stream(os, it->second);
	  }
	}
	return os;
  }

  void WDict::read_dict_from_stream (std::istream & is) {

	size_t id;
	read_atom_from_stream(is, id);
	if (id != magic_id_dict_v0) {
	  throw runtime_error("[E] reading serialized dictionary: invalid id (same platform used to compile the KB?)");
	}
	string hw;
	read_atom_from_stream(is, m_N);
	if (!is) return;
	for (size_t i = 0; i < m_N; ++i) {
	  read_atom_from_stream(is, hw);
	  wdict_rhs_t e;
	  read_rhs_from_stream(is, e);
	  m_wdicts.insert(make_pair(hw, e));
	}
  }
}
