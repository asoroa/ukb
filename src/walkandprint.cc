#include "walkandprint.h"
#include "wdict.h"
#include "common.h"
#include "globalVars.h"
#include "kbGraph.h"
#include <string>
#include <iostream>
#include <fstream>

// Random

#include <boost/random/random_device.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/discrete_distribution.hpp>

namespace ukb {

	using namespace std;
	using namespace boost;

	template<class It>
	int rnumber_distribution(It beg, It end) {
		boost::random::discrete_distribution<> dist(beg, end);
		return dist(glVars::rnd::urng);
	}

	// returns a random floating-point value uniformly distributed in the range [0..1)
	static float rnumber_01() {
		boost::random::uniform_01<> dist;
		return dist(glVars::rnd::urng);
	}

	// generate between 0 and b, both inclusive
	static int rnumber(int b) {
		boost::random::uniform_int_distribution<> dist(0, b);
		return dist(glVars::rnd::urng);
	}

	// generate between 0 and b, both inclusive
	static size_t rnumber(size_t b) {
		boost::random::uniform_int_distribution<> dist(0, b);
		return dist(glVars::rnd::urng);
	}

	static float rnumber(float b) {
		boost::random::uniform_real_distribution<> dist(0.0f, b);
		return dist(glVars::rnd::urng);
	}

	///////////////////////////////////////////////
	// vsampling stuff

	struct vsampling_t::sort_t {

		const vector<float> & m_v ;
		sort_t(const vector<float> & score) : m_v(score) {}
		int operator()(const int & i, const int & j) {
			// Descending order
			return m_v[i] > m_v[j];
		}
	};

	void vsampling_t::debug() {
		const vector<float> & sprank = Kb::instance().static_prank();
		cerr << std::setprecision(4) << 1.0f / static_cast<float>(m_bucket_N) << "\n";
		for(size_t i = 0; i < m_intervals.size(); i++) {
			size_t isize = m_intervals[i].second - m_intervals[i].first ;
			cerr << "Interval (" << isize << ")" << i << ": ";
			for (int j = m_intervals[i].first; j != m_intervals[i].second; j++) {
				string v = Kb::instance().get_vertex_name(m_idx[j]);
				float rank = sprank[m_idx[j]];
				cerr << v << ":" << std::setprecision(4) << rank << " ";
			}
			cerr << "\n";
		}
	}

	vsampling_t::vsampling_t(size_t buckets) :
		m_N( Kb::instance().size() ), m_bucket_N(buckets) {
		if (!m_bucket_N) return ;
		init(Kb::instance().static_prank());
	}

	// Create buckets according to ranks vector, which has to be a probability
	// vector

	vsampling_t::vsampling_t(size_t buckets, const vector<float> & ranks) :
		m_N( Kb::instance().size() ), m_bucket_N(buckets) {
		if (!m_bucket_N) return ;
		init(ranks);
	}

	// Init buckets according to ranks vector, which has to be a probability
	// vector

	void vsampling_t::init(const vector<float> & ranks)  {
		vector<int>(m_N).swap(m_idx);
		for(size_t i = 0; i < m_N; i++) m_idx[i] = i;
		// sort idx according to rank vector
		sort(m_idx.begin(), m_idx.end(), sort_t(ranks));
		float factor = 1.0f / static_cast<float>(m_bucket_N) ;
		float accum = 0.0;
		int left = 0;
		for(size_t i = 0; i < m_N; i++) {
			accum += ranks[ m_idx[i] ];
			if (accum > factor) {
				m_intervals.push_back(make_pair(left, i));
				if (m_intervals.size() == m_bucket_N) break;
				left = i;
				accum = 0.0;
			}
		}
		m_intervals.back().second = m_N; // put remaining nodes in last bucket
		m_bucket_N = m_intervals.size(); // safety measure
	}

	// Method to sample a vertex
	// First, decide which bucket, then select one at random from the bucket
	int vsampling_t::sample() {
		if (!m_bucket_N) {
			return rnumber(m_N - 1);
		}
		size_t bucket_i = rnumber(m_bucket_N - 1);
		int left = m_intervals[bucket_i].first;
		int m = m_intervals[bucket_i].second - left;
		assert(m > 0);
		size_t off = rnumber(m - 1);
		return m_idx[left + off];
	}

	///////////////////////////////////////////////
	// sampling according to strong components

	void vsampling_components_t::debug() {
		for(size_t i = 0; i < m_intervals.size(); i++) {
			cerr << "Interval [ " << m_intervals[i].first << ", " << m_intervals[i].second << ") W: "<< std::setprecision(4) << m_compW[i] << " : ";
			// for (int j = m_intervals[i].first; j != m_intervals[i].second; j++) {
			//	string v = Kb::instance().get_vertex_name(m_idx[j]);
			//	cerr << v << " ";
			// }
			cerr << "\n";
		}
	}

	struct vsampling_components_t::sort_t {

		const vector<size_t> & m_components ;
		const vector<int> & m_comp_sizes ;
		sort_t(const vector<size_t> & components, const vector<int> & comp_sizes) : m_components(components), m_comp_sizes(comp_sizes) {}
		int operator()(const int & i, const int & j) {
			// Descending order
			size_t comp_i = m_components[i];
			size_t comp_j = m_components[j];
			int size_delta = m_comp_sizes[ comp_i ] - m_comp_sizes[ comp_j ];
			if (size_delta > 0) return true;  // comp_i is bigger that comp_j
			if (size_delta < 0) return false; // comp_i is smaller that comp_j
			return comp_i < comp_j; // if equal size, sort descending order
		}
	};

	vsampling_components_t::vsampling_components_t() : m_N(0), m_component_N(0) {

		Kb & kb = Kb::instance();
		m_N = kb.size();
		if (!m_N) return ;
		vector<size_t> components;
		m_component_N = kb.components(components);
		if (m_component_N == 1) return;
		vector<int> comp_sizes(m_component_N, 0);
		for(vector<size_t>::const_iterator it = components.begin(), end = components.end();
			it != end; ++it) {
			comp_sizes[*it]++;
		}
		vector<int>(m_N).swap(m_idx);
		for(size_t i = 0; i < m_N; i++) m_idx[i] = i;
		// sort idx according to component
		sort(m_idx.begin(), m_idx.end(), sort_t(components, comp_sizes));

		// calculate component intervals
		// and fill probabilities of components
		float prob_factor = 1.0f / static_cast<float>(m_N);
		size_t current_comp = components[ m_idx[0] ];
		int left = 0;
		for(size_t i = 1; i < m_N; i++) {
			size_t new_comp = components[ m_idx[i] ];
			if (new_comp != current_comp) {
				m_intervals.push_back(make_pair(left, i));
				m_compW.push_back(comp_sizes[current_comp] * prob_factor);
				left = i;
				current_comp = new_comp;
			}
		}
		m_intervals.push_back(make_pair(left, m_N));
		m_compW.push_back(comp_sizes[current_comp] * prob_factor);
		// accumulate component probabilities
		float accum = 0.0;
		for(size_t i = 0; i < m_component_N - 1; i++) {
			accum += m_compW[i];
			m_compW[i] = accum;
		}
		m_compW[m_component_N - 1] = 1.0;
	}

	// Method to sample a vertex
	// First, decide which component (according to size), then select one at random from the component
	int vsampling_components_t::sample() {
		if (m_component_N < 2) {
			return rnumber(m_N - 1);
		}
		float toss = rnumber_01();
		size_t comp_i;
		for(comp_i = 0; comp_i < m_component_N; ++comp_i)
			if (m_compW[comp_i] > toss) break;
		assert(comp_i < m_component_N);

		int left = m_intervals[comp_i].first;
		int m = m_intervals[comp_i].second - left;
		assert(m > 0);
		if (m == 1) return m_idx[left];
		size_t off = rnumber(m - 1);
		return m_idx[left + off];
	}

	static bool emit_word_vertex_mono(Kb::vertex_descriptor current, string & emit_word,
									  vector<float> & vertex2word_tweight) {

		WInvdict_entries words = WDict::instance().words(current);
		if (!words.size()) return false;
		float & total_weight = vertex2word_tweight[current];
		if (total_weight == 0.0f) {
			float T = 0.0f;
			WInvdict_entries::freq_const_iterator freq_it( words.begin() );
			WInvdict_entries::freq_const_iterator freq_end( words.end() );
			for(;freq_it != freq_end; ++freq_it) T += *freq_it;
			total_weight = T;
		}
		if (total_weight == 0.0f) return false; // Should never happen
		float rand_value = rnumber(total_weight);
		float w_accum = 0;
		for(size_t i = 0; i < words.size(); ++i) {
			w_accum += words.get_prob(i);
			emit_word = words.get_word(i);
			if (rand_value < w_accum) break;
		}
		return true;
	}

	static pair<string, bool> xtract_lang(const string & str) {
		size_t m = str.size();
		if (m < 3) return make_pair(string(), false);
		if (str[m-3] != '#')
			return make_pair(string(), false);
		return make_pair(string(str, m - 2), true);
	}

	struct wap_langselect_sort_t {
		wap_langselect_sort_t() {}
		int operator()(const pair<string, float> & a, const pair<string, float> & b) {
			string langA, langB;
			bool ok;
			tie(langB, ok) = xtract_lang(b.first);
			if (!ok) return 1;
			tie(langA, ok) = xtract_lang(a.first);
			if (!ok) return 0;
			return langA < langB;
		}
	};

	// select words of one language at random (languages are equiprobable)
	static pair<vector<pair<string, float> >::iterator,
				vector<pair<string, float> >::iterator>
	subvec_language(WInvdict_entries & words,
					vector<pair<string, float> > & words_mono)
	{
		size_t N = words.size();
		for(size_t i = 0; i < N; ++i) {
			words_mono.push_back(make_pair(words.get_word(i), words.get_prob(i)));
		}
		sort(words_mono.begin(), words_mono.end(), wap_langselect_sort_t());
		vector<pair<string, float> >::iterator it = words_mono.begin(), end = words_mono.end();
		bool ok;
		string stored_lang, current_lang;
		tie(stored_lang, ok) = xtract_lang(it->first);
		if (!ok) return make_pair(words_mono.begin(), end);
		vector<vector<pair<string, float> >::iterator> lang_idx;
		lang_idx.push_back(it);
		++it;
		while(it != end) {
			tie(current_lang, ok) = xtract_lang(it->first);
			if (!ok) return make_pair(words_mono.begin(), end);
			if(stored_lang != current_lang) {
				lang_idx.push_back(it);
				stored_lang = current_lang;
			}
			++it;
		}
		int langN = lang_idx.size();
		if (langN == 1)
			return make_pair(words_mono.begin(), end);
		int idx = rnumber(langN - 1);
		it = lang_idx[idx];
		if (idx == langN - 1)
			end = words_mono.end();
		else
			end = lang_idx[idx + 1];
		return make_pair(it, end);
	}

	static bool emit_word_vertex_multi(Kb::vertex_descriptor current, string & emit_word) {

		WInvdict_entries words = WDict::instance().words(current);
		if (!words.size()) return false;
		vector<pair<string, float> > words_mono;
		vector<pair<string, float> >::iterator beg, end;
		tie(beg, end) = subvec_language(words, words_mono);
		float total_weight = 0.0f;
		for(vector<pair<string, float> >::iterator it = beg;
			it != end; ++it) {
			total_weight += it->second;
		}
		if (total_weight == 0.0f) return false; // Should never happen

		float rand_value = rnumber(total_weight);
		float w_accum = 0;
		for(; beg != end; ++beg) {
			w_accum += beg->second;
			emit_word = beg->first;
			if (rand_value < w_accum) break;
		}
		return true;
	}

	static bool emit_word_vertex(Kb::vertex_descriptor current, string & emit_word,
								 vector<float> & vertex2word_tweight) {
		Kb & kb = Kb::instance();

		float toss = rnumber_01();
		if (toss >= glVars::wap::wemit_prob) {
			emit_word = kb.get_vertex_name(current);
			return true;
		}
		if(glVars::wap::multilang)
			return emit_word_vertex_multi(current, emit_word);
		else
			return emit_word_vertex_mono(current, emit_word, vertex2word_tweight);
	}

	static bool select_next_vertex_eweight(Kb::vertex_descriptor & current,
										   vector<float> & vertex_out_tweight) {

		Kb & kb = Kb::instance();
		Kb::vertex_descriptor previous = current;

		Kb::out_edge_iterator out_it, out_end;

		tie(out_it, out_end) = kb.out_neighbors(previous);
		float & total_weight = vertex_out_tweight[previous];
		if (total_weight == 0.0f) {
			float T = 0.0f;
			for(Kb::out_edge_iterator auxit = out_it; auxit < out_end; ++auxit) {
				T += kb.get_edge_weight(*out_it);
			}
			total_weight = T;
		}
		if (total_weight == 0.0f) return false; // danglink link. The RW is over.
		float rand_value = rnumber(total_weight);
		float w_accum = 0.0f;
		for(; out_it < out_end; ++out_it) {
			w_accum += kb.get_edge_weight(*out_it);
			current = kb.edge_target(*out_it);
			if (rand_value < w_accum) break;
		}
		return true;
	}

	static bool select_next_vertex_degree(Kb::vertex_descriptor & current,
										  vector<float> & vertex_out_tweight) {

		Kb & kb = Kb::instance();
		Kb::boost_graph_t & G = kb.graph();
		Kb::vertex_descriptor previous = current;

		Kb::out_edge_iterator out_it, out_end;

		tie(out_it, out_end) = kb.out_neighbors(previous);
		float & total_weight = vertex_out_tweight[previous];
		if (total_weight == 0.0f) {
			float T = 0.0f;
			for(Kb::out_edge_iterator auxit = out_it; auxit < out_end; ++auxit) {
				Kb::vertex_descriptor uu = kb.edge_target(*auxit);
				T += in_degree(uu, G);
			}
			total_weight = T;
		}
		if (total_weight == 0.0f) return false; // danglink link. The RW is over.
		float rand_value = rnumber(total_weight);
		float w_accum = 0.0f;
		for(; out_it < out_end; ++out_it) {
			Kb::vertex_descriptor uu = kb.edge_target(*out_it);
			w_accum += in_degree(uu, G);
			current = uu;
			if (rand_value < w_accum) break;
		}
		return true;
	}

	static bool select_next_vertex(Kb::vertex_descriptor & current,
								   vector<float> & vertex_out_tweight) {
		if (glVars::wap::prefer_indegree) return select_next_vertex_degree(current,
																		   vertex_out_tweight);
		return select_next_vertex_eweight(current,
										  vertex_out_tweight);
	}

	// perform one complete rw starting from v

	static void do_complete_mc(Kb::vertex_descriptor v,
							   vector<string> & emited_words,
							   vector<float> & vertex2word_tweight,
							   vector<float> & vertex_out_tweight) {

		Kb::vertex_descriptor current = v;  //Start the iteration in the V vertex

		for (float r = rnumber_01(); r <= glVars::prank::damping; r = rnumber_01() ) {
			// emit word from current vertex
			string emit_word;
			if (emit_word_vertex(current, emit_word, vertex2word_tweight))
				emited_words.push_back(emit_word);
			// Select next vertex to jump
			if (!select_next_vertex(current, vertex_out_tweight)) break;
		}
	}

	// perform one complete rw starting from v (fixed path length)

	static void do_complete_length(Kb::vertex_descriptor current,
								   size_t t,
								   vector<string> & emited_words,
								   vector<float> & vertex2word_tweight,
								   vector<float> & vertex_out_tweight) {

		if (!t) return;

		vector<string>().swap(emited_words);

		for (size_t i = 0; i < t; ++i ) {
			// emit word from current vertex
			string emit_word;
			if (emit_word_vertex(current, emit_word, vertex2word_tweight))
				emited_words.push_back(emit_word);
			// Select next vertex to jump
			if (!select_next_vertex(current, vertex_out_tweight)) break;
		}
	}

	// Main algorithm for walk&print

	bool Wap::next(vector<string> & emited_words) {

		vector<string>().swap(emited_words);

		if (m_n && m_i >= m_n) return false;

		int idx = m_vsampler.sample();
		Kb::vertex_descriptor u(idx);
		if (!m_cache_init) {
			if (!glVars::wap::multilang)
				vector<float>(Kb::instance().size(), 0.0f).swap(m_vertex2word_tweight);
			vector<float>(Kb::instance().size(), 0.0f).swap(m_vertex_out_tweight);
			m_cache_init = true;
		}
		do_complete_mc(u, emited_words, m_vertex2word_tweight, m_vertex_out_tweight);
		m_i++;
		return true;
	}

	// Main algorithm for walk&print

	bool WapComponents::next(vector<string> & emited_words) {

		vector<string>().swap(emited_words);

		if (m_n && m_i >= m_n) return false;

		int idx = m_vsampler.sample();
		Kb::vertex_descriptor u(idx);
		if (!m_cache_init) {
			if (!glVars::wap::multilang)
				vector<float>(Kb::instance().size(), 0.0f).swap(m_vertex2word_tweight);
			vector<float>(Kb::instance().size(), 0.0f).swap(m_vertex_out_tweight);
			m_cache_init = true;
		}
		do_complete_mc(u, emited_words, m_vertex2word_tweight, m_vertex_out_tweight);
		m_i++;
		return true;
	}

	// Main algorithm for Walk&Print starting from a word

	bool WapWord::next(vector<string> & emited_words) {

		vector<string>().swap(emited_words);

		if (m_n && m_i >= m_n) return false;

		static boost::unordered_map<const std::string *, float > dweight_cache;
		boost::unordered_map<const std::string *, float >::iterator dweight_it;
		bool P;

		// select synset at random
		tie(dweight_it, P) = dweight_cache.insert(std::make_pair(&m_seed, 0.0f));
		float & total_weight = dweight_it->second;
		if (P) {
			for(size_t i = 0; i < m_synsets.size(); ++i) {
				total_weight += m_synsets.get_freq(i);
			}
		}
		if (total_weight == 0.0f) return false;; // word has no attached vertices
		float rand_value = rnumber(total_weight);
		float w_accum = 0;
		Kb::vertex_descriptor synset = 0;
		for(size_t i = 0; i < m_synsets.size(); ++i) {
			w_accum += m_synsets.get_freq(i);
			synset = m_synsets.get_entry(i);
			if (rand_value < w_accum) break;
		}
		if (!m_cache_init) {
			if (!glVars::wap::multilang)
				vector<float>(Kb::instance().size(), 0.0f).swap(m_vertex2word_tweight);
			vector<float>(Kb::instance().size(), 0.0f).swap(m_vertex_out_tweight);
			m_cache_init = true;
		}
		do_complete_mc(synset, emited_words, m_vertex2word_tweight, m_vertex_out_tweight);
		m_i++;
		return true;
	}

	bool DeepWalk::next(vector<string> & emited_words) {

		vector<string>().swap(emited_words);

		if (m_g >= m_gamma) return false;

		Kb::vertex_descriptor u(m_i);
		if (!m_cache_init) {
			if (!glVars::wap::multilang)
				vector<float>(Kb::instance().size(), 0.0f).swap(m_vertex2word_tweight);
			vector<float>(Kb::instance().size(), 0.0f).swap(m_vertex_out_tweight);
			m_cache_init = true;
		}
		do_complete_length(u, m_t, emited_words, m_vertex2word_tweight, m_vertex_out_tweight);
		if (++m_i >= m_N) {
			m_i = 0;
			m_g++;
		}
		return true;
	}

	// Not used

	// void wap_do_mc_words(size_t n) {

	//	WDictHeadwords dicthws(WDict::instance());

	//	static boost::unordered_map<const std::string &, float > dweight_cache;
	//	boost::unordered_map<const std::string &, float >::iterator dweight_it;
	//	bool P;
	//	size_t m = dicthws.size();

	//	for(size_t i = 0; i < n; ++i) {
	//		int idx = rnumber((int) m - 1);
	//		const string & hw = dicthws.hw(idx);

	//		// select synset at random
	//		WDict_entries synsets(dicthws.rhs(idx));
	//		tie(dweight_it, P) = dweight_cache.insert(make_pair(hw, 0.0f));
	//		float & total_weight = dweight_it->second;
	//		if (P) {
	//			for(size_t i = 0; i < synsets.size(); ++i) {
	//				total_weight += synsets.get_freq(i);
	//			}
	//		}
	//		float rand_value = rnumber(total_weight);
	//		float w_accum = 0;
	//		Kb::vertex_descriptor synset = 0;
	//		for(size_t i = 0; i < synsets.size(); ++i) {
	//			w_accum += synsets.get_freq(i);
	//			synset = synsets.get_entry(i);
	//			if (rand_value < w_accum) break;
	//		}
	//		// print RW
	//		vector<string> emited_words;
	//		do_mc_complete(synset, emited_words);
	//		// print RW
	//		if(emited_words.size()) {
	//			cout << hw << "\t";
	//			vector<string>::iterator it = emited_words.begin();
	//			vector<string>::iterator end = emited_words.end();
	//			if (it != end) {
	//				--end;
	//				for(;it != end; ++it) {
	//					cout << *it << " ";
	//				}
	//				cout << *end << "\n";
	//			}
	//		}
	//	}
	// }

};
