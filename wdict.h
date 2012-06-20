// -*-C++-*-

#ifndef WDICT_H
#define WDICT_H

#include <string>
#include <iterator>
#include <map>
#include <set>
#include <vector>
#include <iosfwd>

////////////////////////////////////////

#include "kbGraph.h"

namespace ukb {

  // type of dict items


  struct WDict_item_t {
	std::vector<Kb_vertex_t> m_wsyns;
	std::vector<float> m_counts;
	std::vector<char> m_thepos;

	WDict_item_t() {}
  };

  // Accessor class for WDict entries associated to a word

  class WDict_entries {

	const WDict_item_t & _item;

  public:
	WDict_entries(const WDict_item_t & item) : _item(item) {}
	~WDict_entries() {}

	size_t size() const { return _item.m_wsyns.size(); }
	Kb_vertex_t get_entry(size_t i) const { return _item.m_wsyns[i]; }
	const std::string & get_entry_str(size_t i) const;
	float get_freq(size_t i) const;
	char get_pos(size_t i) const;
	size_t dist_pos() const;
  };


  std::ostream & operator<<(std::ostream & o, const WDict_item_t & item);

  class WDict {
  public:
	// Singleton
	static WDict & instance();

	WDict_entries get_entries(const std::string & word) const;

	const std::vector<std::string> & headwords() const { return m_words; }

	std::string variant(std::string &concept_id) const;
	std::string variant(Kb_vertex_t v) const;

	friend std::ostream& operator<<(std::ostream & o, const WDict & dict);

  private:

	WDict();
	WDict(const WDict &);
	WDict & operator=(const WDict &);

	void read_wdict_file(const std::string & fname);

	void create_variant_map();

  public:
	// functor for comparing keys

	struct Mycomp {
	  bool operator() (const std::string * a, const std::string * b) const {
		return *a < *b;
	  }
	};

	typedef std::map<const std::string *, WDict_item_t, Mycomp > wdicts_t;

  private:
	wdicts_t m_wdicts;
	std::vector<std::string> m_words;
	std::map<std::string, std::string> m_variants;
  };
}
#endif
