// -*-C++-*-

#ifndef FILEELEM_H
#define FILEELEM_H

#include <string>
#include <vector>

namespace ukb {

	using std::string;

	std::vector<std::string> extract_input_files(const std::string & fullname,
												 const std::string & extension = std::string());

	bool exists_file(const std::string & fname);

	std::string basename(const std::string & fname);
	std::string get_fname_absolute(const std::string & fname);

	struct File_elem {

		File_elem(const string & fname);
		File_elem(const string & fullname_in,
				  const string & out_dir,
				  const string & new_ext = string());

		void fill(const string & str);

		void set_path(const string & out_dir);

		string get_fname() const; // gett full filename (name + extension)

		string path;  // path to file
		string fname; // filename without extension
		string ext;   // extension
	};
}
#endif
