// -*-C++-*-

#ifndef WDICT_H
#define WDICT_H

#include <string>
#include <iterator>
#include <map>
#include <set>
#include <vector>
#include <iosfwd>
#include <stdexcept>
#include <tr1/unordered_map>

////////////////////////////////////////

#include "kbGraph.h"

namespace ukb {

  class wdict_error : public std::logic_error {
  public:
	wdict_error(const std::string& msg = "") : std::logic_error(msg) {}
  };


  // specicies a range [left, right)
  struct wdict_range_t {
	std::string pos;
	size_t left;
	size_t right;
	wdict_range_t() : pos(), left(0), right(0) {}
	wdict_range_t(std::string p, size_t a, size_t b) : pos(p), left(a), right(b) {}
  };

  // type of dict items
  struct WDict_item_t {
	std::vector<Kb_vertex_t> m_wsyns;
	std::vector<float> m_counts;
	std::vector<wdict_range_t> m_pos_ranges;

	WDict_item_t() {}

	void swap(WDict_item_t & o) {
	  m_wsyns.swap(o.m_wsyns);
	  m_counts.swap(o.m_counts);
	  m_pos_ranges.swap(o.m_pos_ranges);
	}

  };

  // Accessor class for WDict entries associated to a word

  class WDict_entries {


  public:
	WDict_entries(const WDict_item_t & item);
	WDict_entries(const WDict_item_t & item, const std::string & pos);
	~WDict_entries() {}

	size_t size() const;
	Kb_vertex_t get_entry(size_t i) const;
	const std::string & get_entry_str(size_t i) const;
	float get_freq(size_t i) const;
	const std::string & get_pos(size_t i) const;
	size_t dist_pos() const;

	friend std::ostream & operator<<(std::ostream & o, const WDict_entries & item);

  private:
	const WDict_item_t & m_item;
	std::string m_pos;
	size_t m_left;
	size_t m_right;

  };

  std::ostream & operator<<(std::ostream & o, const WDict_item_t & item);

  class WDict {
  public:
	// Singleton
	static WDict & instance();

	size_t size();
	WDict_entries get_entries(const std::string & word, const std::string & pos = std::string()) const;

	//const std::vector<std::string> & headwords() const { return m_words; }

	std::string variant(std::string &concept_id) const;
	std::string variant(Kb_vertex_t v) const;

	void read_alternate_file(const std::string & fname);
	friend std::ostream& operator<<(std::ostream & o, const WDict & dict);

	void write_wdict_binfile(const std::string & fname);

	// Debug
	void  size_bytes();

  private:

	WDict();
	WDict(const WDict &);
	WDict & operator=(const WDict &);

	void read_wdict_file(const std::string & fname);

	void create_variant_map();

	// Streaming
	void read_dict_from_stream (std::istream & is);
	void read_wdict_binfile(const std::string & fname);
	std::ostream & write_dict_to_stream (std::ostream & os) const;

  public:
	typedef std::tr1::unordered_map<std::string, WDict_item_t > wdicts_t;

  private:
	wdicts_t m_wdicts;
	size_t m_N; // number of headwords
	std::map<std::string, std::string> m_variants;
  };
}
#endif
