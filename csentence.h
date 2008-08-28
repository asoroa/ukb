// -*-C++-*-

#ifndef CSENTENCE_H
#define CSENTENCE_H

//#include "mcrGraph.h"
#include <string>
#include <vector>
#include <iosfwd>
#include <boost/graph/graph_traits.hpp>

//typedef std::vector<Mcr_vertex_t> CWord;

//typedef unsigned int Vertex_t; // Achtung!

class CSentence; // forward declaration

class CWord {

public:

  typedef std::vector<std::string>::const_iterator const_iterator;
  typedef std::vector<std::string>::iterator iterator;
  typedef std::vector<std::string>::reference reference;
  typedef std::vector<std::string>::const_reference const_reference;
  typedef std::vector<std::string>::value_type value_type;
  typedef std::vector<std::string>::size_type size_type;

  explicit CWord() : m_weight(1.0), distinguished(false), ranks_equal(true), disamb(false) {};
  CWord(const std::string & w_);
  CWord(const std::string & w_, const std::string & id, char pos, bool is_dist);
  CWord & operator=(const CWord & cw_);
  ~CWord() {};

  iterator begin() {return syns.begin();}
  iterator end() {return syns.end();}
  const_iterator begin() const {return syns.begin();}
  const_iterator end() const {return syns.end();}
  size_type size() const {return syns.size(); }

  const std::string & syn(size_t i) const { return syns[i];}
  float rank(size_t i) const { return ranks[i];}

  std::string word() const { return w; }

  std::string wpos() const;

  std::string id() const {return cw_id;}
  char get_pos() const {return pos;}

  float get_weight() const { return m_weight;}
  void set_weight(float w) { m_weight = w;}

  bool is_distinguished() const { return distinguished; }
  bool is_disambiguated() const { return disamb; }
  bool is_monosemous() const { return (1 == syns.size()); }

  void empty_synsets() { 
    std::vector<std::string>().swap(syns); 
    std::vector<float>().swap(ranks); 
  }
  std::vector<std::string> & get_syns_vector() { return syns; }

  template <typename G, typename Map> 
  void rank_synsets(G & g, Map rankMap) {
    size_t n = syns.size();
    size_t i;
    if (!n) return; // No synsets
    for(i = 0; i != n; ++i) 
      ranks[i] = rankMap[g.get_vertex_by_name(syns[i]).first];
    
    for(i = 1; i != n; ++i) {
      if(ranks[i] != ranks[i-1]) {
	ranks_equal = false;
	break;
      }
    }
  }

  void disamb_cword();

  friend std::ostream& operator<<(std::ostream & o, const CWord & cw_);
  std::ostream & print_cword_simple(std::ostream & o) const;
  std::ostream & print_cword_aw(std::ostream & o) const;
  std::ostream & print_cword_semcor_aw(std::ostream & o) const;
  friend class CSentence;
  
private:

  void read_from_stream (std::ifstream & is);
  std::ofstream & write_to_stream(std::ofstream & o) const;
  void shuffle_synsets();

  std::string w;
  std::string cw_id;
  char pos; // 'n', 'v', 'a', 'r' or 0 (no pos)
  float m_weight;     // Initial weight for PPV
  std::vector<std::string> syns;
  std::vector<float> ranks;
  bool distinguished;
  bool ranks_equal; // If all ranks have the same value (for disambiguating)
  bool disamb;      // If word is disambiguated, that is, if the synset
		    // are ordered according to their ranks
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

  void distinguished_synsets(std::vector<std::string> & res) const;

  std::istream & read_aw(std::istream & is);

  void write_to_binfile (const std::string & fName) const;
  void read_from_binfile (const std::string & fName);
  friend std::ostream& operator<<(std::ostream & o, const CSentence & cs_);
  std::ostream & print_csent_aw(std::ostream & o) const;
  std::ostream & print_csent_semcor_aw(std::ostream & o) const;

  std::ostream & print_csent_simple(std::ostream & o) const;
  
private:  
  void read_from_stream (std::ifstream & is);
  std::ofstream & write_to_stream(std::ofstream & o) const;
  std::vector<CWord> v;
  std::string cs_id;
};

bool calculate_mcr_hr(const CSentence & cs,
		      std::vector<float> & res,
		      bool with_weight);

void calculate_mcr_hr_by_word_and_disamb(CSentence & cs,
					 bool with_weight);

bool calculate_mcr_ppv_csentence(CSentence & cs, std::vector<float> & res);

void disamb_csentence_mcr(CSentence & cs,
                          std::vector<float> & ranks);

#endif
