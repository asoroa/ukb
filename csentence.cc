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
	case 0:
	  res = CWord::cw_ctxword;
	  break;
	case 1:
	  res = CWord::cw_tgtword;
	  break;
	case 2:
	  res = CWord::cw_concept;
	  break;
	case 3:
	  res = CWord::cw_tgtword_nopv;
	  break;
	default:
	  res = CWord::cw_error;
	  break;
	}
	return res;
  }

  size_t CWord::link_dict_concepts(const string & lemma, char pos) {

	size_t new_c = 0;

	set<Kb_vertex_t> syns_S;
	for(size_t i = 0, m = m_V.size(); i < m; ++i) {
	  syns_S.insert(m_V[i].first);
	}

	WDict_entries entries = WDict::instance().get_entries(lemma);
	if (!entries.size()) return 0;

	vector<size_t> sidxV; // synset index vector
	for(size_t i = 0, m = entries.size(); i < m; ++i)
	  sidxV.push_back(i);

	if(glVars::dict::use_shuffle) {
	  // Shuffle index vector
	  boost::random_number_generator<boost::mt19937, long int> rand_dist(glVars::rand_generator);
	  std::random_shuffle(sidxV.begin(), sidxV.end(), rand_dist);
	}

	for(size_t i= 0, m = sidxV.size(); i < m; ++i) {
	  size_t idx = sidxV[i];
	  if(pos) {
		// filter synsets by pos
		char synpos = entries.get_pos(idx);
		if (pos != synpos) continue;
	  }
	  Kb_vertex_t syn_v = entries.get_entry(idx);
	  if (!syns_S.insert(syn_v).second) continue; // Synset previously there
	  const string & syn_str = entries.get_entry_str(idx);
	  float syn_freq = entries.get_freq(idx);
	  m_syns.push_back(syn_str);
	  m_V.push_back(make_pair(syn_v, syn_freq));
	  new_c++;
	}

	if (new_c == 0) return 0;

	// (Re)calculate m_linkw_factor
	float wlink = 0.0;
	for(vector<pair<Kb_vertex_t, float> >::iterator it = m_V.begin(), end = m_V.end();
		it != end; ++it) wlink += it->second;
	if (wlink == 0) return 0;
	m_linkw_factor = 1.0 / wlink;

	// Update ranks
	vector<float>(m_syns.size(), 0.0).swap(m_ranks);
	return new_c;
  }

  CWord::CWord(const string & w_, const string & id_, char pos_, cwtype type_, float wght_)
	: w(w_), m_id(id_), m_pos(pos_), m_weight(wght_), m_type(type_) {

	// cw_tgtword_nopv are internally cw_tgtword with zero weight
	if (m_type == cw_tgtword_nopv) {
	  m_type = cw_tgtword;
	  m_weight = 0;
	}

	switch(m_type) {
	case cw_concept:

	  Kb_vertex_t u;
	  bool P;
	  tie(u, P) = ukb::Kb::instance().get_vertex_by_name(w);
	  if (!P) {
		throw std::runtime_error("CWord concept " + w + " not in KB");
	  }
	  m_syns.push_back(w);
	  m_V.push_back(make_pair(u, 1.0f));
	  m_ranks.push_back(0.0f);
	  m_linkw_factor = 1.0;
	  break;
	case cw_ctxword:
	case cw_tgtword:
	  if (!link_dict_concepts(w, m_pos)) {
		// empty CWord
		std::vector<std::string>().swap(m_syns);
	  }
	  m_disamb = (1 == m_syns.size()); // monosemous words are disambiguated
	  break;
	default:
	  break;
	}
  }

  void CWord::attach_lemma(const string & lemma, char pos) {
	if (!w.size())
	  throw std::runtime_error("CWord::attach_lemma error: can't attach lemma to an empty CWord.");
	if (m_type == cw_concept)
	  throw std::runtime_error("CWord::attach_lemma error: can't attach lemma to a CWord of type cw_concept.");
	if (glVars::input::filter_pos && !pos)
	  throw std::runtime_error("CWord::attach_lemma error: no POS.");
	link_dict_concepts(lemma, pos);
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

	if (is_synset()) return word();
	string wpos(word());
	char pos = get_pos();
	if(pos == 0) return wpos;
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
	for(vector<pair<Kb_vertex_t, float> >::const_iterator it = m_V.begin(), end = m_V.end();
		it != end; ++it) {
	  o << "[" << it->first << ", " << it->second << "], ";
	}
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

	o << cw_.w << "#";
	if (cw_.m_pos)
	  o << cw_.m_pos;
	o << "#" << cw_.m_id << "#" << cw_.m_type;
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
	float norm_factor = 1.0;
	if (glVars::output::norm_ranks) {
	  float rsum = 0.0;
	  for(size_t i = 0; i != syns.size(); ++i) {
		rsum += ranks[i];
	  }
	  norm_factor *= 1.0 / rsum;
	}
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

		  try {
			ctw_parse_t ctwp = parse_ctw(*it);
			if (ctwp.lemma.size() == 0) continue;
			char pos(0);
			CWord::cwtype cw_type = cast_int_cwtype(ctwp.dist);
			if (cw_type == CWord::cw_error) {
			  throw std::runtime_error(*it + " fourth field is invalid.");
			}
			if (cw_type != CWord::cw_concept && glVars::input::filter_pos) {
			  if (!ctwp.pos.size()) throw std::runtime_error(*it + " has no POS.");
			  if (ctwp.pos.size() != 1) throw std::runtime_error(*it + " has invalid POS (more than 1 char).");
			  pos = ctwp.pos[0];
			}
			if(!glVars::input::weight)
			  ctwp.w = 1.0;
			CWord new_cw(ctwp.lemma, ctwp.id, pos, cw_type, ctwp.w);

			if (new_cw.size()) {
			  v.push_back(new_cw);
			} else {
			  // No synset for that word.
			  if (glVars::debug::warning)
				cerr << "W:" << *it << " can't be mapped to KB.";
			}
		  } catch (std::exception & e) {
			string msg(e.what());
			if (!glVars::input::swallow) throw std::runtime_error(msg);
			if (glVars::debug::warning) {
			  cerr << msg << "\n";
			}
		  }
		}
	  } catch (std::exception & e) {
		throw std::runtime_error("Context error in line " + lexical_cast<string>(l_n) + "\n" + e.what());
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
	  if (!cw_it->is_tgtword()) continue;
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
	  if (!cw_it->is_tgtword()) continue;
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
	  if (!cw_it->is_tgtword()) continue;
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
	  if (!cw_it->is_tgtword()) continue;
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


  // Update pv of synsets pointed by CWord of token type

	// cw is the cword node
	// For each v pointed in V, (i.e., cw->v should be in the graph), init the PV with:
	//
	// PV[v] += normalized_cw_w * e[cw->v] / Sum_{cw->u}(e[cw->u]) =
    //        = factor * e[cw->v]
	//
	// Note:
	//
	// Can not be used when glVars::prank::poslightw is set


  size_t update_pv_cw(const vector<pair<Kb_vertex_t, float> > & m_V,
					  float factor,
					  vector<float> & pv) {

	// Sum edge weights

	size_t inserted = 0;
	// Uppdate PV
	for(vector<pair<Kb_vertex_t, float> >::const_iterator it = m_V.begin(), end = m_V.end();
		it != end; ++it) {
	  inserted++;
	  pv[it->first] += it->second * factor;
	}
	return inserted;
  }

  // Get personalization vector giving an csentence. onlyC variant.

  int pv_from_cs_onlyC(const CSentence & cs,
					   vector<float> & pv,
					   CSentence::const_iterator exclude_word_it) {

	Kb & kb = ukb::Kb::instance();

	if (kb.size() == pv.size()) {
	  std::fill(pv.begin(), pv.end(), 0.0);
	} else {
	  vector<float> (kb.size(), 0.0).swap(pv);
	}

	Kb_vertex_t u;
	int inserted_i = 0;

	float CW_w = 0.0;

	vector<const CWord *> cs_uniq;
	set<string> S;
	bool aux;

	if(exclude_word_it != cs.end()) {
	  S.insert(exclude_word_it->wpos());
	}

	// Remove duplicated tokens and get words weight normalization factor
	for(CSentence::const_iterator it = cs.begin(), end = cs.end();
		it != end; ++it) {
	  if (it == exclude_word_it) continue;
	  float w = it->get_weight();
	  if (w == 0.0) continue;
	  set<string>::iterator aux_set;
	  tie(aux_set, aux) = S.insert(it->wpos());
	  if (!aux) continue; // already inserted
	  cs_uniq.push_back(&(*it));
	  CW_w += w;
	}
	float CW_w_factor = 1.0f / CW_w;

	// put pv to the synsets of words
	for(vector<const CWord *>::const_iterator it = cs_uniq.begin(), end = cs_uniq.end();
		it != end; ++it) {
	  const CWord & cw = **it;
	  float cw_w = cw.get_weight() * CW_w_factor;
	  if (cw.type() == CWord::cw_concept) {
		tie(u, aux) = kb.get_vertex_by_name(cw.word());
		assert(aux);
		pv[u] += cw_w;
		inserted_i++;
	  } else {
		inserted_i += update_pv_cw(cw.V_vector(),
								   cw_w * cw.get_linkw_factor(),
								   pv);
	  }
	}
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
	int aux = pv_from_cs_onlyC(cs, pv, cs.end());
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


  // given a word (pointed by tgtw_it),
  // 1. put a ppv in the synsets of the rest of words.
  // 2. Pagerank
  // 3. use rank for disambiguating word

  bool calculate_kb_ppr_by_word(const CSentence & cs,
								CSentence::const_iterator tgtw_it,
								vector<float> & ranks) {

	Kb & kb = ukb::Kb::instance();
	vector<float> pv;

	int aux = pv_from_cs_onlyC(cs, pv, tgtw_it);
	// Execute PageRank
	if (aux) {
	  kb.pageRank_ppv(pv, ranks);
	}
	return aux;
  }

  // given a CSentence
  // for each target word
  //   1. put a ppv in the synsets of the rest of words.
  //   2. Pagerank
  //   3. use rank for disambiguating word

  int calculate_kb_ppr_by_word_and_disamb(CSentence & cs) {

	Kb & kb = ukb::Kb::instance();
	vector<float> ranks;
	int success_n = 0;

	bool use_prior = glVars::dict::use_weight && glVars::csentence::mult_priors; // If --dict-weight and ppr_w2w, use priors when ranking synsets

	vector<CWord>::iterator cw_it = cs.begin();
	vector<CWord>::iterator cw_end = cs.end();
	for(; cw_it != cw_end; ++cw_it) {

	  // Target word must be distinguished.
	  if(!cw_it->is_tgtword()) continue;

	  if (calculate_kb_ppr_by_word(cs, cw_it, ranks)) {
		success_n++;
		if (glVars::csentence::disamb_minus_static) {
		  struct va2vb newrank(ranks, kb.static_prank());
		  cw_it->rank_synsets(newrank, use_prior);
		} else {
		  cw_it->rank_synsets(ranks, use_prior);
		}
	  }
	  cw_it->disamb_cword();
	}
	return success_n;
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
		tie(u, aux) = kb.get_vertex_by_name(cw_it->syn(i));
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
	  if (!cw_it->is_tgtword()) continue;
	  if (glVars::csentence::disamb_minus_static) {
		struct va2vb newrank(ranks, kb.static_prank());
		cw_it->rank_synsets(newrank, false);
	  } else {
		cw_it->rank_synsets(ranks, false);
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
