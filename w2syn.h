// -*-C++-*-

#ifndef W2SYN_H
#define W2SYN_H

#include <string>
#include <iterator>
#include <map>
#include <vector>
#include <fstream>

// Global variables

extern const std::string w2s_filename;

// Word2Syns stuff (singleton)

class W2Syn {
public:
  // Singleton
  static W2Syn & instance();

  const std::vector<std::string> & get_syns(const std::string & word) const;
  void w2syn_reverse(std::map<std::string, std::string> & rev) const;
private:
  W2Syn();
  W2Syn(const W2Syn &);
  W2Syn & operator=(const W2Syn &);
  std::map<std::string, std::vector<std::string> > m_w2syns;
};

#endif
