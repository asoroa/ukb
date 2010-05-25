#include "csentence.h"
#include "common.h"
#include "globalVars.h"
#include "kbGraph.h"
#include "wdict.h"

// Tokenizer
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

#include<boost/tuple/tuple.hpp> // for "tie"

namespace ukb {
  using namespace std;
  using namespace boost;

  static CWord::cwtype cast_int_cwtype(int i) {
	CWord::cwtype res;

	switch(i) {
	case CWord::cwtoken:
	  res = CWord::cwtoken;
	  break;
	case CWord::cwdist:
	  res = CWord::cwdist;
	  break;
	case CWord::cwsynset:
	  res = CWord::cwsynset;
	  break;
	default:
	  res = CWord::cwerror;
	  break;
	}
	return res;
  }

  bool CWord::tie_to_kb() {

	Kb & kb = ukb::Kb::instance();
	bool existP;
	map<string, pair<Kb_vertex_t, float> > str2kb;
	vector<string>::const_iterator str_it;
	vector<string>::const_iterator str_end;

	WDict_entries entries = WDict::instance().get_entries(w);

	for(size_t i= 0; i < entries.size(); ++i) {

	  const string & syn_str = entries.get_entry(i);
	  float syn_freq = glVars::dict::use_weight ? entries.get_freq(i) : 1.0;

	  if(m_pos && glVars::input::filter_pos) {
		// filter synsets by pos
		char synpos = entries.get_pos(i);
		if(!synpos) {
		  throw std::runtime_error("CWord: dictionary concept " + syn_str + " has no POS\n");
		}
		if (m_pos != synpos) continue;
	  }

	  Kb_vertex_t syn_v;
	  tie(syn_v, existP) = kb.get_vertex_by_name(syn_str, Kb::is_concept);
	  if (existP) {
		m_syns.push_back(syn_str);
		str2kb[syn_str] = make_pair(syn_v, syn_freq);
	  } else {
		if (glVars::debug::warning) {
		  cerr << "W:CWord: synset " << syn_str << " of word " << w << " is not in KB" << endl;
		}
	  }
	}

	if (m_syns.size() == 0) return false;

	// Shuffle synsets string vector
	boost::random_number_generator<boost::mt19937, long int> rand_dist(glVars::rand_generator);
	std::random_shuffle(m_syns.begin(), m_syns.end(), rand_dist);

	// Update ranks
	vector<float>(m_syns.size(), 0.0).swap(m_ranks);

	// Update KB with CWord
	Kb_vertex_t word_v;
	Kb_vertex_t wpos_v;
	Kb_vertex_t w_v;

	// Insert word in KB

	word_v = kb.find_or_insert_word(word());
	w_v = word_v;
	// If pos then insert wpos and link to word
	if(m_pos && glVars::input::filter_pos) {
	  wpos_v = kb.find_or_insert_word(wpos());
	  kb.find_or_insert_edge(word_v, wpos_v, 1.0);
	  w_v = wpos_v;
	}

	// Insert related concepts/synsets

	for(vector<string>::iterator it = m_syns.begin(), end = m_syns.end();
		it != end; ++it) {
	  pair<Kb_vertex_t, float> uf = str2kb[*it];
	  m_V.push_back(uf.first);
	  // tie word to synsets
	  kb.find_or_insert_edge(w_v, uf.first, uf.second);
	}
	return true;
  }

  CWord::CWord(const string & w_, const string & id_, char pos_, cwtype type_, float wght_)
	: w(w_), m_id(id_), m_pos(pos_), m_weight(wght_), m_type(type_) {
	switch(m_type) {
	case cwsynset:

	  Kb_vertex_t u;
	  bool P;
	  tie(u, P) = ukb::Kb::instance().get_vertex_by_name(w, Kb::is_concept);
	  if (!P) {
		throw std::runtime_error("CWord concept " + w + " not in KB");
	  }
	  m_syns.push_back(w);
	  m_V.push_back(u);
	  m_ranks.push_back(0.0f);

	  break;
	case cwtoken:
	case cwdist:
	  tie_to_kb();
	  m_disamb = (1 == m_syns.size()); // monosemous words are disambiguated
	  break;
	default:
	  break;
	}
  }

  CWord & CWord::operator=(const CWord & cw_) {
	if (&cw_ != this) {
	  w = cw_.w;
	  m_id = cw_.m_id;
	  m_weight = cw_.m_weight;
	  m_pos = cw_.m_pos;
	  m_syns = cw_.m_syns;
	  m_V = cw_.m_V;
	  m_ranks = cw_.m_ranks;
	  m_disamb = cw_.m_disamb;
	}
	return *this;
  }

  string CWord::wpos() const {

	if (is_synset()) {
	  std::runtime_error("CWoord::wpos: can't get wpos of Cword synset " + w);
	}
	string wpos(w);
	char pos = get_pos();
	if(!glVars::input::filter_pos || pos == 0) return wpos;
	wpos.append("#");
	wpos.append(1,pos);
	return wpos;
  };

  ostream & CWord::debug(ostream & o) const  {

	o << "w: " << w << " \n";
	o <<  "m_id: " << m_id << string(" \n");
	o << "m_pos: "  << m_pos << string(" \n");
	o << "m_weight: "  << lexical_cast<string>(m_weight) << string(" \n");
	o << "m_distinguished: "  << lexical_cast<int>(m_type) << string(" \n");
	o << "m_disamb: "  << lexical_cast<bool>(m_disamb) << string(" \n");
	o << "m_syns: ";
	writeV(o, m_syns);
	o << string(" \n");
	o << "m_V: ";
	writeV(o, m_V);
	o << string(" \n");
	o << "m_ranks: ";
	writeV(o, m_ranks);
	o << string(" \n");
	return o;
  }


  struct CWSort {

	CWSort(const vector<float> & _v) : v(_v) {}
	int operator () (const int & i, const int & j) {
	  // Descending order
	  return v[i] > v[j];
	}
	const vector<float> & v;
  };

  void CWord::disamb_cword() {

	size_t n = m_syns.size();
	if (!n) return;
	if (n == 1) {
	  m_disamb = true; // Monosemous words are disambiguated
	  return;
	}

	vector<string> syns(n);
	vector<float> ranks(n);
	vector<int> idx(n);

	for(size_t i=0; i < n; ++i)
	  idx[i] = i;
	sort(idx.begin(), idx.end(), CWSort(m_ranks));

	bool ranks_equal = true;
	for(size_t i=0; i < n; ++i) {
	  syns[i]  = m_syns[idx[i]];
	  ranks[i] = m_ranks[idx[i]];
	  if (ranks[i] != ranks[0])
		ranks_equal = false;
	}
	if (ranks_equal) return; // If all ranks have same value the word is not disambiguated
	syns.swap(m_syns);
	ranks.swap(m_ranks);
	m_disamb = true;
  }


  std::ostream& operator<<(std::ostream & o, const CWord & cw_) {

	//KbGraph & g = ukb::Kb::instance().graph();

	o << cw_.w;
	if (cw_.m_pos)
	  o << "-" << cw_.m_pos;
	if(glVars::csentence::concepts_in)
	  o << "#" << cw_.m_weight;
	o << "#" << cw_.m_id << "#" << cw_.m_type << " " << cw_.m_disamb;
	o << '\n';
	for(size_t i = 0; i < cw_.m_syns.size(); ++i) {
	  assert(i < cw_.m_ranks.size());
	  o << cw_.m_syns[i] << ":" << cw_.m_ranks[i] << '\n';
	}
	return o;
  }

  ostream & cw_aw_print_best(ostream & o,
							 const vector<string> & syns,
							 const vector<float> & ranks) {
	size_t n = syns.size();
	o << " " << syns[0];
	for(size_t i = 1; i != n; ++i) {
	  if (ranks[i] != ranks[0]) break;
	  o << " " << syns[i];
	}
	return o;
  }

  ostream & cw_aw_print_all(ostream & o,
							const vector<string> & syns,
							const vector<float> & ranks) {
	float rsum = 0.0;
	for(size_t i = 0; i != syns.size(); ++i) {
	  rsum += ranks[i];
	}
	float norm_factor = 1.0 / rsum;
	for(size_t i = 0; i != syns.size(); ++i) {
	  o << " " << syns[i] << "/" << ranks[i]*norm_factor;
	}
	return o;
  }

  ostream & CWord::print_cword_aw(ostream & o) const {

	vector<string> id_fields(split(m_id, "."));
	assert(id_fields.size() > 0);
	o << id_fields[0] << " " << m_id;
	if(!glVars::output::allranks) cw_aw_print_best(o, m_syns, m_ranks);
	else cw_aw_print_all(o, m_syns, m_ranks);
	o << " !! " << w << "\n";
	return o;
  }

  ostream & CWord::print_cword_semcor_aw(ostream & o) const {

	if ((1 == m_syns.size()) && !glVars::output::monosemous) return o; // Don't output monosemous words

	vector<string> id_fields(split(m_id, "."));
	assert(id_fields.size() > 0);
	o << id_fields[0] << "." << id_fields[1] << " " << m_id;
	if(!glVars::output::allranks) cw_aw_print_best(o, m_syns, m_ranks);
	else cw_aw_print_all(o, m_syns, m_ranks);
	o << " !! " << w << "\n";
	return o;
  }

  ostream & CWord::print_cword_simple(ostream & o) const {

	o << m_id << " ";
	if(!glVars::output::allranks) cw_aw_print_best(o, m_syns, m_ranks);
	else cw_aw_print_all(o, m_syns, m_ranks);
	o << " !! " << w << "\n";
	return o;
  }

  ////////////////////////////////////////////////////////
  // Streaming

  std::ofstream & CWord::write_to_stream(std::ofstream & o) const {

	write_atom_to_stream(o,w);
	write_atom_to_stream(o,m_id);
	write_atom_to_stream(o,m_pos);
	write_vector_to_stream(o,m_syns);
	write_atom_to_stream(o,m_type);
	return o;

  };

  void CWord::read_from_stream(std::ifstream & i) {

	read_atom_from_stream(i,w);
	read_atom_from_stream(i,m_id);
	read_atom_from_stream(i,m_pos);
	read_vector_from_stream(i,m_syns);
	vector<float>(m_syns.size()).swap(m_ranks); // Init ranks vector
	read_atom_from_stream(i,m_type);

	m_disamb = (1 == m_syns.size());
  };

  ////////////////////////////////////////////////////////////////
  // CSentence

  struct ctw_parse_t {
	string lemma;
	string pos;
	string id;
	int dist; // 0 -> no dist; 1 -> dist; 2 -> concept (no dist)
	float w;

	ctw_parse_t() : lemma(), pos(), id(), dist(0), w(1.0) {}

  };

  ctw_parse_t parse_ctw(const string & word) {

	ctw_parse_t res;

	char_separator<char> wsep("#", "", keep_empty_tokens);
	tokenizer<char_separator<char> > wtok(word, wsep);
	tokenizer<char_separator<char> >::iterator it = wtok.begin();
	tokenizer<char_separator<char> >::iterator end = wtok.end();
	if (it == end) return res; // empty line
	vector<string> fields;
	copy(it, end, back_inserter(fields));
	size_t m = fields.size();
	if (m != 4 && m != 5) {
	  throw std::runtime_error(word + " : too few fields.");
	}
	res.lemma = fields[0];
	res.pos = fields[1];
	res.id = fields[2];
	try {
	  res.dist = lexical_cast<int>(fields[3]);
	  if (m == 5)
		res.w = lexical_cast<float>(fields[4]);
	} catch (boost::bad_lexical_cast &) {
	  throw std::runtime_error(word + " : Parsing error.");
	}

	if (res.w < 0.0) {
	  throw std::runtime_error(word + " : Negative weight.");
	}
	return res;
  }

  // AW file read (create a csentence from a context)

  istream & CSentence::read_aw(istream & is, size_t & l_n) {

  string line;

	if(read_line_noblank(is, line, l_n)) {

	  // first line is id
	  char_separator<char> sep(" \t");
	  vector<string> ctx;

	  tokenizer<char_separator<char> > tok_id(line, sep);
	  copy(tok_id.begin(), tok_id.end(), back_inserter(ctx));
	  if (ctx.size() == 0) return is; // blank line or EOF
	  cs_id = ctx[0];
	  vector<string>().swap(ctx);
	  // next comes the context
	  if(!read_line_noblank(is, line, l_n)) return is;

	  tokenizer<char_separator<char> > tok_ctx(line, sep);
	  copy(tok_ctx.begin(), tok_ctx.end(), back_inserter(ctx));
	  if (ctx.size() == 0) return is; // blank line or EOF
	  int i = 0;
	  try {
		for(vector<string>::const_iterator it = ctx.begin(); it != ctx.end(); ++i, ++it) {

		  ctw_parse_t ctwp = parse_ctw(*it);
		  if (ctwp.lemma.size() == 0) continue;
		  char pos(0);
		  if (ctwp.pos.size() && glVars::input::filter_pos) pos = ctwp.pos[0];
		  CWord new_cw(ctwp.lemma, ctwp.id, pos, cast_int_cwtype(ctwp.dist), ctwp.w);

		  if (new_cw.size()) {
			v.push_back(new_cw);
		  } else {
			// No synset for that word.
			if (glVars::debug::warning) {
			  cerr << "W: " << ctwp.lemma;
			  if (glVars::input::filter_pos && pos)
				cerr << "-" << pos;
			  cerr << " can't be mapped to KB." << endl;
			}
		  }
		}
	  } catch (std::exception & e) {
		throw std::runtime_error("Context error in line " + lexical_cast<string>(l_n) + "\n" + e.what() );
	  }
	}
	return is;
  }


  // CSentence::CSentence(const vector<string> & sent) : cs_id(string()) {

  //   set<string> wordS;

  //   vector<string>::const_iterator s_it = sent.begin();
  //   vector<string>::const_iterator s_end = sent.end();
  //   for(;s_it != s_end; ++s_it) {
  //     bool insertedP;
  //     set<string>::iterator it;
  //     tie(it, insertedP) = wordS.insert(*s_it);
  //     if (insertedP) {
  //       v.push_back(CWord(*s_it));
  //     }
  //   }
  // }

  CSentence & CSentence::operator=(const CSentence & cs_) {
	if (&cs_ != this) {
	  v = cs_.v;
	  cs_id = cs_.cs_id;
	}
	return *this;
  }

  void CSentence::append(const CSentence & cs_) {
	vector<CWord> tenp(v.size() + cs_.v.size());
	vector<CWord>::iterator aux_it;
	aux_it = copy(v.begin(), v.end(), tenp.begin());
	copy(cs_.v.begin(), cs_.v.end(), aux_it);
	v.swap(tenp);
  }

  void CSentence::distinguished_synsets(vector<string> & res) const {

	vector<CWord>::const_iterator cw_it, cw_end;
	cw_it = v.begin();
	cw_end = v.end();
	for(; cw_it != cw_end; ++cw_it) {
	  if (!cw_it->is_distinguished()) continue;
	  copy(cw_it->begin(), cw_it->end(), back_inserter(res));
	}
  }

  std::ostream& operator<<(std::ostream & o, const CSentence & cs_) {
	o << cs_.cs_id << endl;
	copy(cs_.begin(), cs_.end(), ostream_iterator<CWord>(o, "\n"));
	return o;
  }

  std::ostream & CSentence::print_csent_aw(std::ostream & o) const {

	vector<CWord>::const_iterator cw_it = v.begin();
	vector<CWord>::const_iterator cw_end = v.end();

	for(; cw_it != cw_end; ++cw_it) {
	  if (cw_it->size() == 0) continue;
	  if (!cw_it->is_distinguished()) continue;
	  if (!cw_it->is_disambiguated() && !glVars::output::ties) continue;
	  if (cw_it->is_monosemous() && !glVars::output::monosemous) return o; // Don't output monosemous words

	  cw_it->print_cword_aw(o);
	}
	return o;
  }

  std::ostream & CSentence::print_csent_semcor_aw(std::ostream & o) const {

	vector<CWord>::const_iterator cw_it = v.begin();
	vector<CWord>::const_iterator cw_end = v.end();

	for(; cw_it != cw_end; ++cw_it) {
	  if (cw_it->size() == 0) continue;
	  if (!cw_it->is_distinguished()) continue;
	  if (!cw_it->is_disambiguated() && !glVars::output::ties) continue;
	  if (cw_it->is_monosemous() && !glVars::output::monosemous) return o; // Don't output monosemous words

	  cw_it->print_cword_semcor_aw(o);
	}
	return o;
  }

  std::ostream & CSentence::print_csent_simple(std::ostream & o) const {

	vector<CWord>::const_iterator cw_it = v.begin();
	vector<CWord>::const_iterator cw_end = v.end();

	for(; cw_it != cw_end; ++cw_it) {
	  if (cw_it->size() == 0) continue;
	  if (!cw_it->is_distinguished()) continue;
	  if (!cw_it->is_disambiguated() && !glVars::output::ties) continue;
	  if (cw_it->is_monosemous() && !glVars::output::monosemous) return o; // Don't output monosemous words
	  o << cs_id << " ";
	  cw_it->print_cword_simple(o);
	}
	return o;
  }

  ////////////////////////////////////////////////////////////////////////////////
  //
  // PageRank in Kb


  // Get personalization vector giving an csentence
  // 'light' wpos

  int cs_pv_vector_w(const CSentence & cs,
					 vector<float> & pv,
					 CSentence::const_iterator exclude_word_it) {

	Kb & kb = ukb::Kb::instance();
	vector<float> (kb.size(), 0.0).swap(pv);
	bool aux;
	set<string> S;
	vector<float> ranks;
	int inserted_i = 0;

	float K = 0.0;
	// put pv to the synsets of words except exclude_word
	for(CSentence::const_iterator it = cs.begin(), end = cs.end();
		it != end; ++it) {
	  if (it == exclude_word_it) continue;
	  unsigned char sflags = it->is_synset() ? Kb::is_concept : Kb::is_word;
	  string wpos = it->wpos();

	  set<string>::iterator aux_set;
	  tie(aux_set, aux) = S.insert(wpos);
	  if (aux) {
		Kb_vertex_t u;
		tie(u, aux) = kb.get_vertex_by_name(wpos, sflags);
		float w = glVars::csentence::pv_no_weight ? 1.0 : it->get_weight();
		if (aux && w != 0.0) {
		  pv[u] += w;
		  K +=w;
		  inserted_i++;
		}
	  }
	}
	if (!inserted_i) return 0;
	// Normalize PPV vector
	float div = 1.0 / static_cast<float>(K);
	for(vector<float>::iterator it = pv.begin(), end = pv.end();
		it != end; ++it) *it *= div;
	return inserted_i;
  }


  // Given a CSentence apply Personalized PageRank and obtain obtain it's
  // Personalized PageRank Vector (PPV)
  //
  // Initial V is computed by activating the words of the context

  bool calculate_kb_ppr(const CSentence & cs,
						vector<float> & res) {

	Kb & kb = ukb::Kb::instance();
	vector<float> pv;
	int aux = cs_pv_vector_w(cs, pv, cs.end());
	if (!aux) return false;
	// Execute PageRank
	kb.pageRank_ppv(pv, res);
	return true;
  }


  // Given 2 vectors (va, vb) return the vector going from va to vb
  // res[i] = vb[i] - va[1]

  struct va2vb {
	va2vb(const vector<float> & va, const vector<float> & vb) :
	  m_va(va), m_vb(vb) {}
	float operator[](size_t i) const {
	  return m_vb[i] - m_va[i];
	}
  private:
	const vector<float> & m_va;
	const vector<float> & m_vb;
  };

  // given a word,
  // 1. put a ppv in the synsets of the rest of words.
  // 2. Pagerank
  // 3. use rank for disambiguating word

  void calculate_kb_ppr_by_word_and_disamb(CSentence & cs) {

	Kb & kb = ukb::Kb::instance();
	vector<float> pv;
	vector<float> ranks;

	vector<CWord>::iterator cw_it = cs.begin();
	vector<CWord>::iterator cw_end = cs.end();
	for(; cw_it != cw_end; ++cw_it) {

	  // Target word must be distinguished.
	  if(!cw_it->is_distinguished()) continue;

	  int aux = cs_pv_vector_w(cs, pv, cw_it);
	  // Execute PageRank
	  if (aux) {
		kb.pageRank_ppv(pv, ranks);
		// disambiguate cw_it
		if (glVars::csentence::disamb_minus_static) {
		  struct va2vb newrank(ranks, kb.static_prank());
		  cw_it->rank_synsets(newrank);
		} else {
		  cw_it->rank_synsets(ranks);
		}
	  }
	  cw_it->disamb_cword();
	}
  }

  //
  // Given a previously disambiguated CSentence (all synsets of words
  // have a rank), calculate a kb prgaRank where PPV is formed by
  // 'activating' just the synsets of CWords with appropiate
  // (normalized) rank
  //

  bool calculate_kb_ppv_csentence(CSentence & cs, vector<float> & res) {

	Kb & kb = ukb::Kb::instance();
	bool aux;

	// Initialize result vector
	vector<float> (kb.size(), 0.0).swap(res);

	// Initialize PPV vector
	vector<float> ppv(kb.size(), 0.0);

	vector<CWord>::iterator cw_it = cs.begin();
	vector<CWord>::iterator cw_end = cs.end();
	float K = 0.0;
	for(; cw_it != cw_end; ++cw_it) {
	  for(size_t i = 0; i != cw_it->size(); ++i) {
		Kb_vertex_t u;
		tie(u, aux) = kb.get_vertex_by_name(cw_it->syn(i), Kb::is_concept);
		if (aux) {
		  ppv[u] = cw_it->rank(i);
		  K += cw_it->rank(i);
		}
	  }
	}

	if (K == 0.0) return false;
	// Normalize PPV vector
	float div = 1.0 / K;
	for(vector<float>::iterator rit = ppv.begin(); rit != ppv.end(); ++rit)
	  *rit *= div;

	// Execute PageRank
	kb.pageRank_ppv(ppv, res);
	return true;
  }


  //
  // Disambiguate a CSentence given a vector of ranks
  //

  void disamb_csentence_kb(CSentence & cs,
							const vector<float> & ranks) {

	Kb & kb = ukb::Kb::instance();

	vector<CWord>::iterator cw_it = cs.begin();
	vector<CWord>::iterator cw_end = cs.end();
	for(; cw_it != cw_end; ++cw_it) {
	  if (!cw_it->is_distinguished()) continue;
	  if (glVars::csentence::disamb_minus_static) {
		struct va2vb newrank(ranks, kb.static_prank());
		cw_it->rank_synsets(newrank);
	  } else {
		cw_it->rank_synsets(ranks);
	  }
	  cw_it->disamb_cword();
	}
  }


  ////////////////////////////////////////////////////////
  // Streaming

  const size_t magic_id = 0x070704;

  ofstream & CSentence::write_to_stream (std::ofstream & o) const {
	size_t vSize = v.size();

	write_atom_to_stream(o, magic_id);
	write_atom_to_stream(o,vSize);
	for(size_t i=0; i< vSize; ++i) {
	  v[i].write_to_stream(o);
	}
	write_atom_to_stream(o,cs_id);
	return o;
  }

  void CSentence::read_from_stream (std::ifstream & is) {

	size_t vSize;
	size_t id;

	// Read a vector of CWords
	read_atom_from_stream(is, id);
	if (id != magic_id) {
	  cerr << "Invalid file. Not a csentence." << endl;
	  exit(-1);
	}
	read_atom_from_stream(is, vSize);
	v.resize(vSize);
	for(size_t i=0; i<vSize; ++i) {
	  v[i].read_from_stream(is);
	}
	read_atom_from_stream(is,cs_id);
  }

  void CSentence::write_to_binfile (const std::string & fName) const {
	ofstream fo(fName.c_str(),  ofstream::binary|ofstream::out);
	if (!fo) {
	  cerr << "Error: can't create" << fName << endl;
	  exit(-1);
	}
	write_to_stream(fo);
  }


  void CSentence::read_from_binfile (const std::string & fName) {
	ifstream fi(fName.c_str(), ifstream::binary|ifstream::in);
	if (!fi) {
	  cerr << "Error: can't open " << fName << endl;
	  exit(-1);
	}
	read_from_stream(fi);
  }
}
