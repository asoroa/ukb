#include "common.h"

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
      if (end_index == std::string::npos)
	end_index = str.length();

      ret.push_back(str.substr(start_index, end_index - start_index));

      // Find beginning of the next token                                                                 
      start_index = str.find_first_not_of(delims, end_index);
    }

  return ret;
}

