#include "common.h"
#include "globalVars.h"

#include <boost/numeric/conversion/cast.hpp>

namespace ukb {

	/////////////////////////////////////////////////////////////////////
	// Random

	boost::minstd_rand global_mersenne;

	int g_randTarget(int Target) {
		// Return rand up to, but not including, Target
		boost::uniform_int<> local_int_dist(0, Target-1);
		//boost::variate_generator<boost::mt19937&, boost::uniform_int<> > local_uni_int(global_mersenne, local_int_dist);
		boost::variate_generator<boost::minstd_rand, boost::uniform_int<> > local_uni_int(global_mersenne, local_int_dist);
		return local_uni_int();
	}

	//////////////////////////////////////////////////////7
	// split

	std::vector<std::string> split(const std::string & str, const std::string & delims) {
		std::string::size_type start_index, end_index;
		std::vector<std::string> ret;

		// Skip leading delimiters, to get to the first token
		start_index = str.find_first_not_of(delims);

		// While found a beginning of a new token
		//
		while (start_index != std::string::npos)
			{
				// Find the end of this token
				end_index = str.find_first_of(delims, start_index);

				// If this is the end of the std::string
				// If this is the end of the string
				if (end_index == std::string::npos) {
					ret.push_back(str.substr(start_index, str.size() - start_index));
					break;
				}

				ret.push_back(str.substr(start_index, end_index - start_index));

				// Find beginning of the next token
				start_index = str.find_first_not_of(delims, end_index);
			}

		return ret;
	}

	std::string join(const std::string &delim,
					 std::vector<std::string>::const_iterator it,
					 std::vector<std::string>::const_iterator end) {

		std::string res("");

		if (it == end) return res;
		--end;
		for(;it != end; ++it) {
			res+=*it;
			res+=delim;
		}
		res+=*end;
		return res;
	}

	std::string join(const std::string & delim, const std::vector<std::string> & V) {
		return join(delim, V.begin(), V.end());
	}


	void trim_spaces(std::string &l) {

		std::string res("");
		std::string::size_type start = l.find_first_not_of(" \t\r");
		std::string::size_type end = l.find_last_not_of(" \t\r");
		if (start != std::string::npos) {
			res = l.substr(start, end - start + 1);
		}
		l.swap(res);
	}

	std::istream & read_line_noblank(std::istream & is, std::string & line, size_t & l_n) {
		do {
			std::getline(is, line, '\n');
			trim_spaces(line);
			l_n++;
		} while(is && !line.size());
		return is;
	}

	// Manage iterations/threshold

	void set_pr_convergence(size_t iterations, float thresh) {

		if (thresh < 0.0 || thresh > 1.0) {
			std::cerr << "Warning: invalid prank_thresh value " << thresh << " set to zero.\n";
			thresh = 0.0f;
		}
		if (iterations == 0 && thresh == 0.0f)
			throw std::runtime_error(std::string("Error: both iterations and threshold are set to zero."));
		if (iterations && thresh != 0.0) {
			// user specified both
			glVars::prank::num_iterations = boost::numeric_cast<int>(iterations);
			glVars::prank::threshold = thresh;
		} else {
			if (iterations) {
				glVars::prank::num_iterations = boost::numeric_cast<int>(iterations);
				glVars::prank::threshold = 0.0;
			}
			if (thresh != 0.0) {
				glVars::prank::num_iterations = 0;
				glVars::prank::threshold = thresh;
			}
		}
	}
}
