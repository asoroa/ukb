#include "globalVars.h"
#include "common.h"
#include "version.h"
#include <iostream>
#include <ctime>


#ifdef UKB_UKB_VERSION
#define UKB_VERSION UKB_UKB_VERSION
#else
#define UKB_VERSION "unknown"
#endif

namespace ukb {
	namespace glVars {

		char ukb_version[] = UKB_VERSION;

		int verbose = 0;

		namespace debug {
			bool warning = false;
		}

		std::vector<std::string> rel_source;

		boost::mt19937 rand_generator(static_cast<unsigned int>(std::time(0)));

		namespace csentence {
			bool concepts_in = true;
			bool disamb_minus_static = false;
			bool mult_priors = true;
		}

		namespace dict {
			std::string text_fname;
			std::string bin_fname;
			bool use_weight = false; // Use weights when linking words to concepts
			float weight_smoothfactor = 1.0;
			bool use_shuffle = false; // Use random shuffle at reading words from dictionary
			bool swallow = true;
		}

		namespace chsq {
			size_t cooc_min = 5;
			float threshold = 3.84146; // 95.0% confidence;
		}

		namespace prank {
			bool use_weight = false;
			int num_iterations = 30; // Conservative, but stop if threshold is reached. If zero, just use threshold.
			float threshold = 0.0001; // If zero just use num_iterations
			float damping = 0.85; // damping factor
			PrankImpl impl = pm; // default is power method
			float nibble_epsilon = 0.0000005;
		}

		namespace input {
			bool filter_pos = true;
			bool swallow = false;
			bool weight = true;
		}

		namespace output {
			bool allranks = false;
			bool monosemous = true;
			bool ties = false;
			bool norm_ranks = true;
		}

		RankAlg rAlg = pageRank;

		namespace kb {
			bool keep_reltypes = false;
			bool keep_directed = true;
			bool v1_kb = true;
			bool filter_src = true;
		}

		namespace dGraph {
			int max_depth = 6;
			bool stopCosenses = false;
		}

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

		RankAlg get_algEnum(const std::string & alg) {

			if ("pageRank" == alg) return pageRank;
			if ("degree" == alg) return degree;
			return no_alg;
		}
	}

	using namespace glVars;
	void show_global_variables(std::ostream & o) {

		o << "General options" << std::endl << "******" << std::endl;
		o << "verbose:" << verbose << '\t' ;
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
}
