// -*-C++-*-

#ifndef W2SYN_H
#define W2SYN_H

#include <string>
#include <iterator>
#include <map>
#include <vector>
#include <iosfwd>

// Global variables

extern const std::string w2s_filename;

// Word2Syns stuff (singleton)


struct W2Syn_item {
  std::vector<std::string> wsyns;
  std::vector<size_t> syns_count;
  bool has_freq;
};

std::ostream & operator<<(std::ostream & o, const W2Syn_item & item);

class W2Syn {
public:
  // Singleton
  static W2Syn & instance();

  std::pair<std::vector<std::string>::const_iterator, std::vector<std::string>::const_iterator> 
  get_wsyns(const std::string & word) const;
  bool syn_counts(std::map<std::string, size_t> & res) const;
private:
  W2Syn();
  W2Syn(const W2Syn &);
  W2Syn & operator=(const W2Syn &);
  std::map<std::string, W2Syn_item > m_w2syns;
};

#endif
