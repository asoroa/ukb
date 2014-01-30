// -*-C++-*-

#ifndef GLVARS_H
#define GLVARS_H

#include<iosfwd>
#include<cstddef>
#include<string>
#include<vector>

// For srand & friends

#include <boost/random/mersenne_twister.hpp>

namespace ukb {
  namespace glVars {

	extern char ukb_version[];

	// debug
	extern int verbose;

	namespace debug {
	  extern bool warning;
	}

	extern std::vector<std::string> rel_source;

	namespace csentence {
	  extern bool concepts_in;
	  extern bool disamb_minus_static;
	  extern bool mult_priors; // multiply priors to tw synsets if --ppr_w2w and --dict_weight are selected
	}

	namespace dict {
	  extern std::string text_fname; // The name of the dictionary file (text file)
	  extern std::string bin_fname; // The name of the dictionary file (compiled version)
	  extern bool use_weight; // Use weights when linking words to concepts
	  extern float weight_smoothfactor;
	  extern bool use_shuffle; // Use random shuffle at reading words from dictionary
	  extern bool swallow; // When reading the dictionary, go ahead even if malformed input
	}

	namespace chsq {
	  extern size_t cooc_min;
	  extern float threshold;
	}

	namespace prank {
	  extern bool use_weight;   // Use weights in pagerank calculations
	  extern int num_iterations;
	  extern float threshold;
	  extern float damping;
	}

	// Input
	namespace input {
	  // Wether input words must be filtered by pos when attaching them
	  // to the KB. It also has effects on dictionary.

	  extern bool filter_pos;
	  extern bool swallow;        // tries to swallow malformed input
	  extern bool weight;
	}


	// Output stuff
	namespace output {
	  extern bool allranks;       // print all ranks (with weights or just best ranks)
	  extern bool monosemous;     // print monosemous words
	  extern bool ties;           // print also if ties
	  extern bool norm_ranks;     // normalize relative ranks of a CWord
	}

	namespace kb {
	  extern bool keep_reltypes; // Wether edges locally keep the relation types
	  extern bool v1_kb; // Wether input has v1 format
	  extern bool filter_src; // Wether input relations should be filtered by relation source
	  extern bool keep_directed; // Wether we will allow directed edges (default true)
	}

	namespace dGraph {
	  extern int max_depth;
	  extern bool stopCosenses;
	}

	extern boost::mt19937 rand_generator;

	// Maiztasunak

	// Ezaugarriak

	enum RankAlg {
	  no_alg, // error
	  pageRank,
	  degree
	};


	extern RankAlg rAlg;

	RankAlg get_algEnum(const std::string & alg);

	//   enum VertexRelate { // How to compute relateness between vertices in cooc. graph
	//     freq,
	//     chsq,         // Chi square
	//     rel_error
	//   };

	//  extern VertexRelate vRelate;

	//   enum HubsAlg {
	//     veronis,
	//     pageRank,
	//     hits,
	//     degree,
	//     mcl,
	//     random
	//   };

	//   extern HubsAlg hubsAlg;


	//   struct Param {
	//     // coocurrences with a frequency of 5 or more are retained
	//     static size_t edgeFreqMin;

	//     // Words with less than 10 occurrences are discarded

	//     static size_t vertexFreqMin;

	//     // Edges with a weight below 0.1 are arbitrarily (totally) elminiated
	//     static float ertzPisuMin;

	//     // Context containing fewer than 4 words are discarded
	//     static size_t ctxSizeMin;
	//   };

	//   struct VeronisParam {
	//     // Good candidate for Hub:
	//     // The mean of the weights btw. the candidate node and its 6 most frequent
	//     // neighbors must be below 0.8
	//     static size_t hubNeighborSize;
	//     static float hubNeighborWeightMin;
	//     static float hubFreqThreshold;
	//   };

	//   struct PrankParam {
	//     static float   freqThreshold;
	//     static int     siblingsCannotBeHubs;
	//   };

	//   struct MclParam {
	//     static float Inflate;
	//   };

  }

  void show_global_variables (std::ostream & o);
}
#endif
