// -*-C++-*-

#ifndef WDICT_H
#define WDICT_H

#include "globalVars.h"
#include "wdict_vector.h"

#include <string>
#include <iterator>
#include <map>
#include <set>
#include <vector>
#include <iosfwd>
#include <stdexcept>
#include <boost/unordered_map.hpp>
#include <boost/iterator/iterator_facade.hpp>

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

	struct wdict_item_t {
		Kb::vertex_descriptor m_syn;
		float m_count;

		wdict_item_t() {}
		wdict_item_t(Kb::vertex_descriptor syn, float count) : m_syn(syn), m_count(count) {}

		void swap(wdict_item_t & o) {
			std::swap(m_syn, o.m_syn);
			std::swap(m_count, o.m_count);
		}
	};

	// type of dict items
	struct wdict_rhs_t {
		wdict_vector<wdict_item_t> m_items;
		wdict_vector<wdict_range_t> m_pos_ranges;

		wdict_rhs_t() {}

		void swap(wdict_rhs_t & o) {
			m_items.swap(o.m_items);
			m_pos_ranges.swap(o.m_pos_ranges);
		}

	};

	// Accessor class for WDict entries associated to a word

	class WDict_entries {

	public:
		WDict_entries(const wdict_rhs_t & item);
		WDict_entries(const wdict_rhs_t & item, const std::string & pos);
		~WDict_entries() {}

		const wdict_item_t *begin() const;
		const wdict_item_t *end() const;
		size_t size() const;
		Kb::vertex_descriptor get_entry(size_t i) const;
		const std::string & get_entry_str(size_t i) const;
		float get_freq(size_t i) const;
		const std::string & get_pos(size_t i) const;
		size_t dist_pos() const;

		friend std::ostream & operator<<(std::ostream & o, const WDict_entries & item);

		class freq_const_iterator:
			public boost::iterator_facade<
			freq_const_iterator,
			float,
			boost::bidirectional_traversal_tag,
			const float &> {

		public:

			typedef const wdict_item_t * value_type;

			freq_const_iterator() : m_current() {}
			explicit freq_const_iterator(value_type p) : m_current(p) {}

		private:
			friend class boost::iterator_core_access;

			const float & dereference() const {
				static float one = 1.0f;
				if (!glVars::dict::use_weight) return one;
				return m_current->m_count;
			}
			void increment() { ++m_current; }
			void decrement() { --m_current; }
			void advance(int n) { m_current += n; }
			std::ptrdiff_t distance_to(const freq_const_iterator & o) {
				return o.m_current - m_current;
			}
			bool equal(const freq_const_iterator & o) const {
				return m_current == o.m_current;
			}

			value_type m_current;
		};

	private:
		const wdict_rhs_t & m_rhs;
		std::string m_pos;
		size_t m_left;
		size_t m_right;

	};

	// inverse dictionary

	struct winvdict_item_t {
		std::string m_word;
		float m_count;

		winvdict_item_t() {}
		winvdict_item_t(std::string word, float count) : m_word(word), m_count(count) {}

		void swap(winvdict_item_t & o) {
			std::swap(m_word, o.m_word);
			std::swap(m_count, o.m_count);
		}
	};

	// type of inverse dict items
	typedef std::vector<winvdict_item_t> winvdict_rhs_t;

	class WInvdict_entries {

	public:
		WInvdict_entries(const winvdict_rhs_t & item) : m_rhs(item) {};
		~WInvdict_entries() {}

		size_t size() const { return m_rhs.size(); }
		const std::string & get_word(size_t i) const { return m_rhs[i].m_word; }
		float get_prob(size_t i) const { return m_rhs[i].m_count; };

		const winvdict_rhs_t::const_iterator begin() const;
		const winvdict_rhs_t::const_iterator end() const;

		//friend std::ostream & operator<<(std::ostream & o, const WInvdict_entries & item);
		class freq_const_iterator:
			public boost::iterator_facade<
			freq_const_iterator,
			float,
			boost::bidirectional_traversal_tag,
			const float &> {

		public:
			typedef winvdict_rhs_t::const_iterator value_type;

			freq_const_iterator() : m_current() {}
			explicit freq_const_iterator(value_type p) : m_current(p) {}

		private:
			friend class boost::iterator_core_access;

			const float & dereference() const {
				static float one = 1.0f;
				if (!glVars::dict::use_weight) return one;
				return m_current->m_count;
			}
			void increment() { ++m_current; }
			void decrement() { --m_current; }
			void advance(int n) { m_current += n; }
			std::ptrdiff_t distance_to(const freq_const_iterator & o) {
				return o.m_current - m_current;
			}
			bool equal(const freq_const_iterator & o) const {
				return m_current == o.m_current;
			}

			value_type m_current;
		};

	private:
		const winvdict_rhs_t & m_rhs;
	};

	std::ostream & operator<<(std::ostream & o, const wdict_rhs_t & item);

	class WDict {

	public:
		typedef boost::unordered_map<std::string, wdict_rhs_t > wdict_t;
		typedef boost::unordered_map<Kb::vertex_descriptor, winvdict_rhs_t > winvdict_t;

		// Singleton
		static WDict & instance();

		size_t size() const;
		size_t size_inv() const;

		WDict_entries get_entries(const std::string & word, const std::string & pos = std::string()) const;
		const wdict_t & wdict() const  { return m_wdict; };

		//const std::vector<std::string> & headwords() const { return m_words; }

		std::string variant(std::string &concept_id) const;
		std::string variant(Kb::vertex_descriptor v) const;

		WInvdict_entries words(Kb::vertex_descriptor u) const;
		WInvdict_entries words(const std::string & ustr) const;

		void read_alternate_file(const std::string & fname);
		friend std::ostream& operator<<(std::ostream & o, const WDict & dict);

		void write_wdict_binfile(const std::string & fname);

		// Debug
		void  size_bytes();

		friend class WDictHeadwords;

	private:

		WDict();
		WDict(const WDict &);
		WDict & operator=(const WDict &);

		void read_wdict_file(const std::string & fname);

		void create_variant_map();
		void create_inverse_dict() const;

		// Streaming
		void read_dict_from_stream (std::istream & is);
		void read_wdict_binfile(const std::string & fname);
		std::ostream & write_dict_to_stream (std::ostream & os) const;

	private:
		wdict_t m_wdict;
		mutable winvdict_t m_wdict_inv;
		size_t m_N; // number of headwords
		std::map<std::string, std::string> m_variants;
	};

	class WDictHeadwords {

	public:

		WDictHeadwords(const WDict & wdict);

		const std::string & hw(size_t i) const;
		WDict_entries rhs(size_t i) const;
		size_t size() const;

	private:

		WDictHeadwords(const WDictHeadwords &);
		WDictHeadwords & operator=(const WDictHeadwords &);

		typedef const std::pair<const std::string, wdict_rhs_t> * pair_ptr_t;
		std::vector<pair_ptr_t> m_V;
	};

	// return concept priors:
	//
	// P[concept] = freq  ==> how many times 'concept' has been used in the dictionary
	// N => total number of concept counts

	float concept_priors(boost::unordered_map<Kb::vertex_descriptor, float> & P);

}



#endif
