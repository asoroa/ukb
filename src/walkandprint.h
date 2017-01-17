// -*-C++-*-

#ifndef WALKPRINT_H
#define WALKPRINT_H

#include <string>
#include <vector>
#include <utility> // pair
#include "kbGraph.h"
#include "wdict.h"


namespace ukb {

	// bucket sampling
	struct vsampling_t {

		struct sort_t; // predicate for sorting

		vsampling_t(size_t buckets);
		// Create buckets according to ranks vector, which has to be a probability
		// vector
		vsampling_t(size_t buckets, const std::vector<float> & ranks);

		int sample();
		void debug();

	private:
		void init(const std::vector<float> & ranks);
		std::vector<int> m_idx;
		size_t m_N; // size of graph
		size_t m_bucket_N;
		std::vector<std::pair<int, int> > m_intervals;

	};

	struct vsampling_components_t {

		struct sort_t; // predicate for sorting

		vsampling_components_t();
		int sample();
		void debug();

		size_t m_N; // size of graph
		size_t m_component_N; // number of components in graph
		std::vector<float> m_compW; // accumulated probabilities of components
		std::vector<int> m_idx; // vertices ordered by component (from larger to smaller)
		std::vector<std::pair<int, int> > m_intervals;

	};


	// Walk&Print class
	class Wap {

	public:

		Wap(size_t n, size_t bucket_size, std::vector<float> & P) : m_n(n), m_bsize(bucket_size), m_vsampler(bucket_size, P), m_i(0), m_cache_init(false) {}
		Wap(size_t n, size_t bucket_size) : m_n(n), m_bsize(bucket_size), m_vsampler(bucket_size), m_i(0), m_cache_init(false) {}
		//Wap(size_t n, size_t bucket_size, const vector<float> & vpriors) : m_n(n), m_bsize(bucket_size), m_vsampler(bucket_size, priors), m_i(0) {}
		~Wap() {};

		// perform an iteration leaving the result context in C
		// return false if iteration is over

		bool next(std::vector<std::string> & C);

	private:

		size_t m_n;     // number of contexts to produce
		size_t m_bsize; // bucket size
		vsampling_t m_vsampler; // sampling from buckets
		size_t m_i;     // number of context produced so far

		// cache vertex2word
		bool m_cache_init;
		std::vector<float> m_vertex2word_tweight;
		// cache vertex2vertex
		std::vector<float> m_vertex_out_tweight;

	};

	// Walk&Print class
	class WapComponents {

	public:

		WapComponents(size_t n) : m_n(n), m_vsampler(), m_i(0), m_cache_init(false) {}
		~WapComponents() {};

		// perform an iteration leaving the result context in C
		// return false if iteration is over

		bool next(std::vector<std::string> & C);

	private:

		size_t m_n;     // number of contexts to produce
		vsampling_components_t m_vsampler; // sampling from buckets
		size_t m_i;     // number of context produced so far

		// cache vertex2word
		bool m_cache_init;
		std::vector<float> m_vertex2word_tweight;
		// cache vertex2vertex
		std::vector<float> m_vertex_out_tweight;
	};


	// Deepwalk algorithm
	class DeepWalk {
	public:
		DeepWalk(size_t gamma, size_t t)
			: m_N(Kb::instance().size()), m_i(0), m_gamma(gamma), m_g(0), m_t(t), m_cache_init(false) {}
		~DeepWalk() {}

		// perform an iteration leaving the result context in C
		// return false if iteration is over

		bool next(std::vector<std::string> & C);

	private:
		size_t m_N;     // number of graph vertices
		size_t m_i;     // number of context produced so far
		float m_gamma;  // number of rw per vertex
		float m_g;      // number of gamma iterations so far
		size_t m_t;     // context size

		// cache vertex2word
		bool m_cache_init;
		std::vector<float> m_vertex2word_tweight;
		// cache vertex2vertex
		std::vector<float> m_vertex_out_tweight;

	};



	// Walk&Print class starting from a word
	class WapWord {

	public:

		WapWord(std::string & hw, size_t n) : m_seed(hw), m_n(n), m_i(0), m_synsets(WDict::instance().get_entries(hw)), m_cache_init(false) {}
		~WapWord() {};

		// perform an iteration leaving the result context in C
		// return false if iteration is over

		bool next(std::vector<std::string> & C);

	private:

		std::string & m_seed;     // seed word
		size_t m_n;     // number of contexts to produce
		size_t m_i;     // number of context produced so far
		WDict_entries m_synsets;

		// cache vertex2word
		bool m_cache_init;
		std::vector<float> m_vertex2word_tweight;
		// cache vertex2vertex
		std::vector<float> m_vertex_out_tweight;
	};

	// Walk&Print starting from a word
	void wap_do_mc_word(const std::string & hw, size_t n);

}


#endif
