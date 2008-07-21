// -*-C++-*-

#ifndef GLVARS_H
#define GLVARS_H

#include<iosfwd>
#include<cstddef>
#include<string>
#include<vector>

// For srand & friends 

#include <boost/random/mersenne_twister.hpp>

namespace glVars {
  // debug
  extern int verbose;

  extern std::vector<std::string> rel_source;
  extern std::string w2s_filename;

  extern bool mcr_with_freqs;

  namespace csentence {
    extern bool word_weight;
  }

  namespace chsq {
    extern size_t cooc_min;
    extern float threshold;
  }

  namespace prank {
    extern size_t num_iterations;
  }

  // Output stuff
  namespace output {
    extern bool allranks; // print all ranks (with weights or just best ranks)
    extern bool monosemous; // print monosemous words
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

//     // Edges with a weigth below 0.1 are arbitrarily (totally) elminiated
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

};

void show_global_variables (std::ostream & o);

#endif
