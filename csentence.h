// -*-C++-*-

#ifndef CSENTENCE_H
#define CSENTENCE_H

//#include "mcrGraph.h"
#include <string>
#include <vector>
#include <iosfwd>

//typedef std::vector<Mcr_vertex_t> CWord;

typedef unsigned int Vertex_t; // Achtung!

class CSentence; // forward declaration

class CWord {

public:

  typedef std::vector<Vertex_t>::const_iterator const_iterator;
  typedef std::vector<Vertex_t>::iterator iterator;
  typedef std::vector<Vertex_t>::reference reference;
  typedef std::vector<Vertex_t>::const_reference const_reference;
  typedef std::vector<Vertex_t>::value_type value_type;
  typedef std::vector<Vertex_t>::size_type size_type;

  explicit CWord() {};
  CWord(const std::string & w_);
  CWord(const std::string & w_, const std::string & id, char pos, bool is_dist);
  CWord & operator=(const CWord & cw_);
  ~CWord() {};

  iterator begin() {return syns.begin();}
  iterator end() {return syns.end();}
  const_iterator begin() const {return syns.begin();}
  const_iterator end() const {return syns.end();}
  std::vector<Vertex_t> & get_syns_vector() { return syns; }
  size_type size() const {return syns.size(); }
  std::string word() const { return w; }
  bool is_distinguished() const { return distinguished; }
  std::string id() const {return cw_id;}

  friend std::ostream& operator<<(std::ostream & o, const CWord & cw_);
  friend class CSentence;
  
private:

  void read_from_stream (std::ifstream & is);
  std::ofstream & write_to_stream(std::ofstream & o) const;

  std::string w;
  std::string cw_id;
  char pos; // 'n', 'v', 'a', 'r' or 0 (no pos)
  std::vector<Vertex_t> syns;
  bool distinguished;
};

class CSentence {

public:

  typedef std::vector<CWord>::const_iterator const_iterator;
  typedef std::vector<CWord>::iterator iterator;
  typedef std::vector<CWord>::reference reference;
  typedef std::vector<CWord>::const_reference const_reference;
  typedef std::vector<CWord>::value_type value_type;
  typedef std::vector<CWord>::size_type size_type;

  CSentence() {};
  CSentence(const std::vector<std::string> & sent_);

  CSentence(const CSentence & cs_) : v(cs_.v) , cs_id(cs_.cs_id) {};
  CSentence & operator=(const CSentence & cs_);

  void append(const CSentence & cs_);

  iterator begin() {return v.begin();}
  const_iterator begin() const {return v.begin();}
  iterator end() {return v.end();}
  const_iterator end() const {return v.end();}
  void push_back(const CWord & cw_) { v.push_back(cw_); }
  reference back() {return v.back();}
  const_reference back() const {return v.back();}

  size_type size() const {return v.size();}
  std::string id() const {return cs_id;}

  void distinguished_vertices(std::vector<Vertex_t> & res) const;
  
  std::istream & read_aw(std::istream & is);

  void write_to_binfile (const std::string & fName) const;
  void read_from_binfile (const std::string & fName);
  friend std::ostream& operator<<(std::ostream & o, const CSentence & cs_);
  
private:  

  void read_from_stream (std::ifstream & is);
  std::ofstream & write_to_stream(std::ofstream & o) const;


  std::vector<CWord> v;
  std::string cs_id;
};

#endif
