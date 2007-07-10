// -*-C++-*-

#ifndef FILEELEM_H
#define FILEELEM_H

#include <string>

// Basename & friends
#include <boost/filesystem/operations.hpp>
#include "boost/filesystem/path.hpp"

using std::string;

struct File_elem {

  File_elem(const string & fname);
  File_elem(const string & fullname_in,
	    const string & out_dir,
	    const string & new_ext = string());

  void fill(const string & str);

  void set_path(const string & out_dir);

  string get_fname() const;

  string path;
  string fname;
  string ext;
};

#endif
