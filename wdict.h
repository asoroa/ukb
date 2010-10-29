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

namespace ukb {

  // type of dict items

  struct WDict_item_t {
	std::vector<std::string> wsyns;
	std::vector<float> syns_count;
	std::vector<char> m_thepos;
	std::set<char> m_distpos;

	WDict_item_t() {}

  };

  // Accessor class for WDict entries associated to a word

  class WDict_entries {

	const WDict_item_t & _item;

  public:
	WDict_entries(const WDict_item_t & item) : _item(item) {}
	~WDict_entries() {}

	size_t size() const { return _item.wsyns.size(); }
	const std::string & get_entry(size_t i) const { return _item.wsyns[i]; }
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


	// old stuff
	std::pair<std::vector<std::string>::const_iterator, std::vector<std::string>::const_iterator>
	get_wsyns(const std::string & word) const;

	std::pair<std::vector<float>::const_iterator, std::vector<float>::const_iterator>
	get_weights(const std::string & word) const;

	bool syn_counts(std::map<std::string, size_t> & res) const;

	const std::vector<std::string> & get_wordlist() const { return m_words; }

  private:

	WDict();
	WDict(const WDict &);
	WDict & operator=(const WDict &);

	void read_wdict_file(const std::string & fname);

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
  };
}
#endif
