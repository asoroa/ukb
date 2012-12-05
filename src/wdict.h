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

////////////////////////////////////////

#include "kbGraph.h"

namespace ukb {

  class wdict_error : public std::logic_error {
  public:
	wdict_error(const std::string& msg = "") : std::logic_error(msg) {}
  };


  // specicies a range [left, right)
  struct wdict_range {
	size_t left;
	size_t right;
	wdict_range(size_t a, size_t b) : left(a), right(b) {}
  };

  // type of dict items
  struct WDict_item_t {
	std::vector<Kb_vertex_t> m_wsyns;
	std::vector<float> m_counts;
	std::map<std::string, wdict_range> m_pos_ranges;

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

	WDict_entries get_entries(const std::string & word, const std::string & pos = std::string()) const;

	const std::vector<std::string> & headwords() const { return m_words; }

	std::string variant(std::string &concept_id) const;
	std::string variant(Kb_vertex_t v) const;

	void read_alternate_file(const std::string & fname);
	friend std::ostream& operator<<(std::ostream & o, const WDict & dict);

  private:

	WDict();
	WDict(const WDict &);
	WDict & operator=(const WDict &);

	void read_wdict_file(const std::string & fname);

	void create_variant_map();

  public:
	typedef std::map<const std::string, WDict_item_t > wdicts_t;

  private:
	wdicts_t m_wdicts;
	std::vector<std::string> m_words;
	std::map<std::string, std::string> m_variants;
  };
}
#endif
