// -*-C++-*-

#ifndef CSENTENCE_H
#define CSENTENCE_H

#include "kbGraph.h"
#include <string>
#include <vector>
#include <iosfwd>
#include <boost/graph/graph_traits.hpp>

namespace ukb {

  class CSentence; // forward declaration

  class CWord {

  public:

	// Type of CWords:
	//
	// cw_ctxword: context word. Affects initial PV calculation but is not disambiguated.
	// cw_tgtword: target word: Affects initial PV calculation and is disambiguated.
	// cw_concept: concept. Affects initial PV calculation and is not diambiguated
	// cw_tgtword_nopv: target 'nopv' word. Does not affect initial PV calculation and will be disambiguated.

	enum cwtype {
	  cw_ctxword = 0,
	  cw_tgtword = 1,
	  cw_concept = 2,
	  cw_tgtword_nopv = 3,
	  cw_error
	};

	typedef std::vector<std::string>::const_iterator const_iterator;
	typedef std::vector<std::string>::iterator iterator;
	typedef std::vector<std::string>::reference reference;
	typedef std::vector<std::string>::const_reference const_reference;
	typedef std::vector<std::string>::value_type value_type;
	typedef std::vector<std::string>::size_type size_type;

	explicit CWord() : m_pos(0), m_weight(1.0), m_type(cw_error), m_disamb(false) {};
	CWord(const std::string & w_, const std::string & id, char pos, cwtype type, float wght_ = 1.0);
	CWord & operator=(const CWord & cw_);
	~CWord() {};

	// Attach synsets of a new lemma to the CWord.
	//
	// Note, the CWord does not track the lemma

	void attach_lemma(const std::string & lemma, char pos = 0);

	iterator begin() {return m_syns.begin();}
	iterator end() {return m_syns.end();}
	const_iterator begin() const {return m_syns.begin();}
	const_iterator end() const {return m_syns.end();}
	size_type size() const {return m_syns.size(); }

	const std::string & syn(size_t i) const { return m_syns[i];}
	float rank(size_t i) const { return m_ranks[i];}

	std::string word() const { return w; }

	std::string wpos() const;

	std::string id() const {return m_id;}
	char get_pos() const {return m_pos;}

	float get_weight() const { return m_weight;}
	float get_linkw_factor() const { return m_linkw_factor; }
	void set_weight(float w) { m_weight = w;}

	bool is_tgtword() const { return (m_type == cw_tgtword); }
	bool is_disambiguated() const { return m_disamb; }
	bool is_monosemous() const { return (1 == m_syns.size()); }
	bool is_synset() const { return m_type == cw_concept; }

	cwtype type() const { return m_type; }

	void empty_synsets() {
	  std::vector<std::string>().swap(m_syns);
	  std::vector<std::pair<Kb_vertex_t, float> >().swap(m_V);
	  std::vector<float>().swap(m_ranks);
	  m_disamb = false;
	}
	std::vector<std::string> & get_syns_vector() { return m_syns; }
	const std::vector<std::pair<Kb_vertex_t, float> > & V_vector() const { return m_V; }

	template <typename Map>
	void rank_synsets(Map rankMap, bool use_prior) {
	  size_t n = m_syns.size();
	  size_t i;
	  if (!n) return; // No synsets
	  for(i = 0; i != n; ++i) {
		m_ranks[i] = rankMap[m_V[i].first];
		if (use_prior) m_ranks[i] *= m_V[i].second * m_linkw_factor;
	  }
	}

	// Used in disambGraph
	template <typename G, typename Map>
	void rank_synsets(G & g, Map rankMap) {
	  size_t n = m_syns.size();
	  size_t i;
	  if (!n) return; // No synsets
	  for(i = 0; i != n; ++i)
		m_ranks[i] = rankMap[g.get_vertex_by_name(m_syns[i]).first];
	}

	void disamb_cword();

	friend std::ostream& operator<<(std::ostream & o, const CWord & cw_);
	std::ostream & print_cword_simple(std::ostream & o) const;
	std::ostream & print_cword_aw(std::ostream & o) const;
	std::ostream & print_cword_semcor_aw(std::ostream & o) const;
	friend class CSentence;

	// Debug

	std::ostream & debug(std::ostream & o) const;

  private:

	size_t link_dict_concepts(const std::string & lemma, char pos);
	void read_from_stream (std::ifstream & is);
	std::ofstream & write_to_stream(std::ofstream & o) const;
	void shuffle_synsets();

	std::string w;
	std::string m_id;
	char m_pos; // 'n', 'v', 'a', 'r' or 0 (no pos)
	float m_weight;     // Initial weight for PPV
	std::vector<std::string> m_syns;
	std::vector<std::pair<Kb_vertex_t, float> > m_V;
	std::vector<float> m_ranks;
	float m_linkw_factor; // 1 / (sum of all link weights)
	cwtype m_type;
	bool m_disamb;      // If word is disambiguated, that is, if the synset
	// are ordered according to their ranks
  };

  class CSentence {

  public:

	typedef std::vector<CWord>::const_iterator const_iterator;
	typedef std::vector<CWord>::iterator iterator;
	typedef std::vector<CWord>::reference reference;
	typedef std::vector<CWord>::const_reference const_reference;
	typedef std::vector<CWord>::value_type value_type;
	typedef std::vector<CWord>::size_type size_type;

	CSentence() {};
	CSentence(const std::vector<std::string> & sent_);

	CSentence(const CSentence & cs_) : v(cs_.v) , cs_id(cs_.cs_id) {};
	CSentence & operator=(const CSentence & cs_);

	void append(const CSentence & cs_);

	iterator begin() {return v.begin();}
	const_iterator begin() const {return v.begin();}
	iterator end() {return v.end();}
	const_iterator end() const {return v.end();}
	void push_back(const CWord & cw_) { v.push_back(cw_); }
	reference back() {return v.back();}
	const_reference back() const {return v.back();}

	const CWord & operator[]( int i ) const {
	  return v[i];
	}

	CWord & operator[]( int i ) {
	  return v[i];
	}

	size_type size() const {return v.size();}
	std::string id() const {return cs_id;}

	void distinguished_synsets(std::vector<std::string> & res) const;

	std::istream & read_aw(std::istream & is, size_t & l_n);

	void write_to_binfile (const std::string & fName) const;
	void read_from_binfile (const std::string & fName);
	friend std::ostream& operator<<(std::ostream & o, const CSentence & cs_);
	std::ostream & print_csent_aw(std::ostream & o) const;
	std::ostream & print_csent_semcor_aw(std::ostream & o) const;

	std::ostream & print_csent_simple(std::ostream & o) const;

  private:
	void read_from_stream (std::ifstream & is);
	std::ofstream & write_to_stream(std::ofstream & o) const;
	std::vector<CWord> v;
	std::string cs_id;
  };

  bool calculate_kb_ppr(const CSentence & cs,
						std::vector<float> & res);

  bool calculate_kb_ppr_by_word(const CSentence & cs,
								CSentence::const_iterator tgtw_it,
								std::vector<float> & ranks);

  int calculate_kb_ppr_by_word_and_disamb(CSentence & cs);

  bool calculate_kb_ppv_csentence(CSentence & cs, std::vector<float> & res);

  void disamb_csentence_kb(CSentence & cs,
							const std::vector<float> & ranks);


  // Functions for calculating initial PV given a CSentence

  int pv_from_cs_onlyC(const CSentence & cs,
					   std::vector<float> & pv,
					   CSentence::const_iterator exclude_word_it);

}
#endif
