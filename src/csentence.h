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

		// cw_ctxword: context word. Affects initial PV calculation but is not
		//             disambiguated.
		//
		// cw_tgtword: target word: Affects initial PV calculation and is
		//             disambiguated.
		//
		// cw_concept: concept. Affects initial PV calculation and is not
		//             diambiguated
		//
		// cw_tgtword_nopv: target 'nopv' word. Does not affect initial PV
		//                  calculation and will be disambiguated.
		//
		// cw_ctxword_nopv: context 'nopv' word. Does not affect initial PV
		//                  calculation and will not be disambiguated. Useful for
		//                  'grouping' concepts.

		enum cwtype {
			cw_ctxword = 0,
			cw_tgtword = 1,
			cw_concept = 2,
			cw_tgtword_nopv = 3,
			cw_ctxword_nopv = 4,
			cw_error
		};

		typedef std::vector<std::string>::const_iterator const_iterator;
		typedef std::vector<std::string>::iterator iterator;
		typedef std::vector<std::string>::reference reference;
		typedef std::vector<std::string>::const_reference const_reference;
		typedef std::vector<std::string>::value_type value_type;
		typedef std::vector<std::string>::size_type size_type;

		explicit CWord() : m_pos(0), m_weight(1.0), m_linkw_factor(1.0), m_type(cw_error), m_disamb(false) {};
		CWord(const std::string & w_, const std::string & id, const std::string & pos, cwtype type, float wght_ = 1.0);
		CWord & operator=(const CWord & cw_);
		~CWord() {};

		// Attach synsets of a new lemma to the CWord.
		//
		// Note, the lemma of the CWord remains untouched

		void attach_lemma(const std::string & lemma, const std::string & pos = std::string());

		// set concepts attached to the cword (used in tgt_noppv type cwords)
		size_t set_concepts(std::map<std::string, float> & C);

		iterator begin() {return m_syns.begin();}
		iterator end() {return m_syns.end();}
		const_iterator begin() const {return m_syns.begin();}
		const_iterator end() const {return m_syns.end();}

		size_type size() const {return m_syns.size(); }

		const std::string & syn(size_t i) const { return m_syns[i];}
		float rank(size_t i) const { return m_ranks[i];}

		std::string word() const { return m_w; }

		std::string wpos() const;

		std::string id() const {return m_id;}

		float get_weight() const { return m_weight;}
		float get_linkw_factor() const { return m_linkw_factor; }
		void set_weight(float w) { m_weight = w;}

		bool is_tgtword() const { return (m_type == cw_tgtword || m_type == cw_tgtword_nopv); }
		bool is_disambiguated() const { return m_disamb; }
		bool is_monosemous() const { return (1 == m_syns.size()); }
		bool is_synset() const { return m_type == cw_concept; }

		bool has_concept(const std::string & str);

		cwtype type() const { return m_type; }

		void empty_synsets();
		const std::vector<std::pair<Kb_vertex_t, float> > & V_vector() const { return m_V; }

		template <typename Map>
		void rank_synsets(Map rankMap, bool use_prior) {
			size_t n = m_syns.size();
			size_t i;
			if (!n) return; // No synsets
			if (rankMap.size()) return; // No ranks
			for(i = 0; i != n; ++i) {
				m_ranks[i] = rankMap[m_V[i].first];
				if (use_prior) m_ranks[i] *= m_V[i].second * m_linkw_factor;
			}
		}

		// Used when only one tw. Apply MFS.
		void rank_synsets_one_tw(bool use_prior) {
			size_t n = m_syns.size();
			size_t i;
			if (!n) return; // No synsets
			float factor = 1.0 / (float) n;
			for(i = 0; i != n; ++i) {
				m_ranks[i] = factor;
				if(use_prior) m_ranks[i] *= m_V[i].second * m_linkw_factor;
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
		std::ostream & print_cword(std::ostream & o) const;
		friend class CSentence;

		// Debug

		std::ostream & debug(std::ostream & o) const;

	private:

		size_t link_dict_concepts(const std::string & lemma, const std::string & pos);
		void shuffle_synsets();

		std::string m_w;
		std::string m_id;
		std::string m_pos; // 'n', 'v', 'a', 'r' or 0 (no pos)
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

		CSentence() : m_tgtN(0) {};
		CSentence(const std::vector<std::string> & sent_);

		CSentence(const CSentence & cs_) : m_tgtN(0), m_v(cs_.m_v) , m_id(cs_.m_id) {};
		CSentence & operator=(const CSentence & cs_);

		void append(const CSentence & cs_);

		iterator begin() {return m_v.begin();}
		const_iterator begin() const {return m_v.begin();}
		iterator end() {return m_v.end();}
		const_iterator end() const {return m_v.end();}
		void push_back(const CWord & cw_) { m_v.push_back(cw_); }
		reference back() {return m_v.back();}
		const_reference back() const {return m_v.back();}

		const CWord & operator[]( int i ) const {
			return m_v[i];
		}

		CWord & operator[]( int i ) {
			return m_v[i];
		}

		size_type size() const {return m_v.size();}
		size_type has_tgtwords() const { return m_tgtN; }
		std::string id() const {return m_id;}

		void distinguished_synsets(std::vector<std::string> & res) const;

		std::istream & read_aw(std::istream & is, size_t & l_n);

		void write_to_binfile (const std::string & fName) const;
		void read_from_binfile (const std::string & fName);
		friend std::ostream& operator<<(std::ostream & o, const CSentence & cs_);
		std::ostream & print_csent(std::ostream & o) const;
		std::ostream & debug(std::ostream & o) const;
	private:
		size_t m_tgtN; // number of target words in Csentence
		std::vector<CWord> m_v;
		std::string m_id;
	};

	bool calculate_kb_ppr(const CSentence & cs,
						  std::vector<float> & res);

	bool calculate_kb_ppr_by_word(const CSentence & cs,
								  CSentence::const_iterator tgtw_it,
								  std::vector<float> & ranks);

	int calculate_kb_ppr_by_word_and_disamb(CSentence & cs);

	bool calculate_kb_ppv_csentence(CSentence & cs, std::vector<float> & res);

	bool disamb_csentence_kb(CSentence & cs,
							 const std::vector<float> & ranks);


	// Functions for calculating initial PV given a CSentence

	int pv_from_cs_onlyC(const CSentence & cs,
						 std::vector<float> & pv,
						 CSentence::const_iterator exclude_word_it);

}
#endif
