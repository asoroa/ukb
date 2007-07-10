#include "fileElem.h"
#include <iostream>

using std::cerr;
using std::cout;
using std::endl;

/////////////////////////////////////////////////////////////
// Filesystem stuff (paths, extensions, etc)

string::const_iterator find_last(string::const_iterator sit,
				 string::const_iterator sit_end, 
				 char delim) {
  string::const_iterator sit_found = sit_end;
  for(;sit != sit_end;++sit) 
    if (*sit == delim) sit_found = sit;
  return sit_found;
}



File_elem::File_elem(const string & fname) {
  fill(fname);
}

File_elem::File_elem(const string & fullname_in,
		     const string & out_dir,
		     const string & new_ext) {

  fill(fullname_in);

  set_path(out_dir);
    
  if (new_ext.size())
    ext = new_ext;
}

void File_elem::fill(const string & str) {

  boost::filesystem::path p(str);
  p.normalize();

  path = p.branch_path().string();
  if (path == "") path = ".";

  string file_fname = p.leaf(); // name + extension

  string::const_iterator beg = file_fname.begin();
  string::const_iterator end = file_fname.end();
  string::const_iterator it = find_last(file_fname.begin(), file_fname.end(), '.');

  ext.assign(it, end);
  fname.assign(beg, it); // without extension
}

void File_elem::set_path(const string & out_dir) {
  size_t m = out_dir.size();
  if (m) {
    if(!boost::filesystem::exists(out_dir)) {
      // Better die
      //boost::filesystem::create_directory( out_dir_dir );
      cerr << "Error: " << out_dir << " doesn't exist.\n" << endl;
      exit(-1);
    }
    if (out_dir[m-1] == '/') --m;
    path.assign(out_dir, 0, m);
  }
}

string File_elem::get_fname() const {
  string res(path);
  res.append("/");
  res.append(fname);
  res.append(ext);
  return res;
}
