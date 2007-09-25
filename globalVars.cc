#include "globalVars.h"
#include "common.h"
#include <iostream>
#include <ctime>

using namespace glVars;

namespace glVars {
  int verbose = 0;

  std::vector<std::string> rel_source;
  std::string w2s_filename = "../Data/Preproc/wn1.6_index.sense_freq";

  bool mcr_with_freqs = false;
  bool output_monosemous = true;

  boost::mt19937 rand_generator(static_cast<unsigned int>(std::time(0)));

//   bool word_norep = 0;
//   int hub_neighborhood = 0;
//   float context_threshold = 0;

//   bool FeaturesCoocurrence = true;
//   size_t numHubsFixed = 0;

//   VertexRelate vRelate = freq;
//   HubsAlg hubsAlg = veronis;

//   size_t Param::edgeFreqMin = 5;
//   size_t Param::vertexFreqMin = 10;
//   float  Param::ertzPisuMin = 0.1;
//   size_t Param::ctxSizeMin = 4;

//   size_t VeronisParam::hubNeighborSize = 6;
//   float  VeronisParam::hubNeighborWeightMin = 0.2;
//   float  VeronisParam::hubFreqThreshold = 0.001;

//   float  PrankParam::freqThreshold = 0.002; // 
//   int    PrankParam::siblingsCannotBeHubs = 0;

//   float MclParam::Inflate = 0.0f;
}

void show_global_variables(std::ostream & o) {

  o << "General options" << std::endl << "******" << std::endl;
  o << "verbose:" << verbose << '\t' ;
  o << "PPV:" << mcr_with_freqs << std::endl;
  o << "Rels: ";
  writeV(o, rel_source);
  o << std::endl;
//   o <<  "Hub_neighborhood: " << hub_neighborhood << '\t';
//   //o << "Context_threshold: " << context_threshold << '\t';
//   o << "FeaturesCoocurrence: " << FeaturesCoocurrence << '\t';
//   o << "word_norep: " << word_norep << std::endl;

//   o << "Hub algorithm: ";

//   switch (hubsAlg) {
//   case akOrok::veronis:
//     o << "veronis";
//     break;
//   case akOrok::pageRank:
//     o << "pageRank";
//     break;
//   case akOrok::degree:
//     o << "degree";
//     break;
//   case akOrok::hits:
//     o << "hits";
//     break;
//   case akOrok::mcl:
//     o << "mcl";
//     break;
//   case akOrok::random:
//     o << "random";
//     break;
//   default:
//     o << "unknow (ERROR)";
//   }
//   if ((hubsAlg == akOrok::pageRank) || (hubsAlg == akOrok::hits)) {
//     o << '\t' << "Stop Criterion: ";
//     if (akOrok::numHubsFixed) o << "Fixed num hubs";
//     else o << "Freq. threshold";
//   }
//   o << std::endl;
//   o << "Relatedness: ";
//   switch (vRelate) {
//   case akOrok::freq:
//     o << "freq";
//     break;
//   case akOrok::chsq:
//     o << "chsq";
//     break;
//   default:
//     o << "unknow (ERROR)";
//   }
//   o << std::endl;
//   o << "Param::edgeFreqMin: " << Param::edgeFreqMin << '\t';
//   o << "Param::vertexFreqMin: " << Param::vertexFreqMin << std::endl;
//   o << "Param::ertzPisuMin: " << Param::ertzPisuMin << '\t';
//   o << "Param::ctxSizeMin: " << Param::ctxSizeMin << std::endl;
//   if (akOrok::numHubsFixed) {
//     o << "akOrok::numHubsFixed: " << akOrok::numHubsFixed << std::endl;
//   } else {
//     if (hubsAlg == veronis) {
//       o << "VeronisParam::hubNeighborSize: " << VeronisParam::hubNeighborSize << '\t';
//       o << "VeronisParam::hubNeighborWeightMin: " << VeronisParam::hubNeighborWeightMin << std::endl;
//       o << "VeronisParam::hubFreqThreshold: " << VeronisParam::hubFreqThreshold << std::endl;
//     } else {
//       o << "PrankParam::freqThreshold: " << PrankParam::freqThreshold << '\t';
//       o << "PrankParam::siblingsCannotBeHubs: " << PrankParam::siblingsCannotBeHubs << std::endl;
//     }
//   }
//   if (hubsAlg == mcl) {
//     o << "MclParam::Inflate: " << MclParam::Inflate << std::endl;
//   }
}
