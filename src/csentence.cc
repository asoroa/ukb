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
	case 4:
	  res = CWord::cw_ctxword_nopv;
	  break;
	default:
	  res = CWord::cw_error;
	  break;
	}
	return res;
  }

  void CWord::empty_synsets() {
	std::vector<std::string>().swap(m_syns);
	std::vector<std::pair<Kb_vertex_t, float> >().swap(m_V);
	std::vector<float>().swap(m_ranks);
	m_linkw_factor = 1;
	m_disamb = false;
  }

  size_t CWord::link_dict_concepts(const string & lemma, const string & pos) {

	size_t new_c = 0;

	set<Kb_vertex_t> syns_S;
	for(size_t i = 0, m = m_V.size(); i < m; ++i) {
	  syns_S.insert(m_V[i].first);
	}

	WDict_entries entries = WDict::instance().get_entries(lemma, pos);
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

  CWord::CWord(const string & w_, const string & id_, const string & pos_, cwtype type_, float wght_)
	: m_w(w_), m_id(id_), m_pos(pos_), m_weight(wght_), m_type(type_) {

	string empty_str;

	if (m_type == cw_concept) {
	  Kb_vertex_t u;
	  bool P;
	  tie(u, P) = ukb::Kb::instance().get_vertex_by_name(m_w);
	  if (!P) {
		throw std::logic_error("CWord concept " + m_w + " not in KB");
	  }
	  m_syns.push_back(m_w);
	  m_V.push_back(make_pair(u, 1.0f));
	  m_ranks.push_back(0.0f);
	  m_linkw_factor = 1.0;
	  m_pos.swap(empty_str); // concepts have no POS
	  return;
	}
	// empty POS string when no pos filtering
	if(!glVars::input::filter_pos)
	  m_pos.swap(empty_str);
	m_linkw_factor = 0;

	if (m_type == cw_tgtword_nopv or m_type == cw_ctxword_nopv) return; // if nopv we are done

	// m_type == cw_ctxword or cw_tgtword:
	if (!link_dict_concepts(m_w, m_pos)) {
	  // empty CWord
	  empty_synsets();
	}
	m_disamb = (1 == m_syns.size()); // monosemous words are disambiguated
  }

  void CWord::attach_lemma(const string & lemma, const string & pos) {
	if (!m_w.size())
	  throw std::logic_error("CWord::attach_lemma error: can't attach lemma to an empty CWord.");
	if (m_type == cw_concept)
	  throw std::logic_error("CWord::attach_lemma error: can't attach lemma to a CWord of type cw_concept.");
	if (glVars::input::filter_pos && !pos.size())
	  throw std::logic_error("CWord::attach_lemma error: no POS.");
	link_dict_concepts(lemma, pos);
  }

  // (Re)Set cword and link it to new concepts.
  // (used in tgt_noppv type cwords)
  // return: number of new concepts

  size_t CWord::set_concepts(map<string, float> & C) {

	this->empty_synsets();

	Kb_vertex_t u;
	bool P;
	float total_w = 0.0f;
	size_t N = 0;
	for(map<string, float>::iterator it = C.begin(), end = C.end();
		it != end; ++it) {
	  N++;
	  tie(u, P) = ukb::Kb::instance().get_vertex_by_name(it->first);
	  if (!P) {
		throw std::logic_error("reset_concepts: " + it->first + " not in KB");
	  }
	  m_syns.push_back(it->first);
	  m_V.push_back(make_pair(u, it->second));
	  total_w += it->second;
	}
	if (total_w == 0.0f) return 0;
	m_linkw_factor = 1.0 / total_w;
	// Update ranks
	vector<float>(m_syns.size(), 0.0).swap(m_ranks);
	m_disamb = (N == 1); // monosemous words are disambiguated
	return N;
  }

  CWord & CWord::operator=(const CWord & cw_) {
	if (&cw_ != this) {
	  m_w = cw_.m_w;
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
	if(m_pos.size() == 0) return wpos;
	wpos.append("#");
	wpos.append(m_pos);
	return wpos;
  };

  ostream & CWord::debug(ostream & o) const  {

	o << "w: " << m_w << " \n";
	o <<  "m_id: " << m_id << string(" \n");
	o << "m_pos: "  << m_pos << string(" \n");
	o << "m_weight: "  << lexical_cast<string>(m_weight) << string(" \n");
	o << "m_distinguished: "  << lexical_cast<int>(m_type) << string(" \n");
	o << "m_disamb: "  << lexical_cast<bool>(m_disamb) << string(" \n");
	o << "m_syns: ";
	writeV(o, m_syns);
	o << string(" \n");
	o << "m_linkw_factor: "  << lexical_cast<string>(m_linkw_factor) << string(" \n");
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

	for(size_t i=0; i < n; ++i) {
	  syns[i]  = m_syns[idx[i]];
	  ranks[i] = m_ranks[idx[i]];
	}
	if(ranks[0] == ranks[n-1]) return; // If all ranks have same value the word is not disambiguated
	syns.swap(m_syns);
	ranks.swap(m_ranks);
	m_disamb = true;
  }


  std::ostream& operator<<(std::ostream & o, const CWord & cw_) {

	//KbGraph & g = ukb::Kb::instance().graph();

	o << cw_.m_w << "#";
	if (cw_.m_pos.size())
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
	  if (rsum) norm_factor *= 1.0 / rsum;
	}
	for(size_t i = 0; i != syns.size(); ++i) {
	  o << " " << syns[i] << "/" << ranks[i]*norm_factor;
	}
	return o;
  }

  ostream & CWord::print_cword(ostream & o) const {

	o << m_id << " ";
	if(!glVars::output::allranks) cw_aw_print_best(o, m_syns, m_ranks);
	else cw_aw_print_all(o, m_syns, m_ranks);
	o << " !! " << m_w << "\n";
	return o;
  }

  ////////////////////////////////////////////////////////////////
  // CSentence

  struct ctw_parse_t {
	string lemma;
	string pos;
	string id;
	int dist; // See cwtype enum.
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
	  throw std::logic_error(word + " : too few fields.");
	}
	res.lemma = fields[0];
	res.pos = fields[1];
	res.id = fields[2];
	try {
	  res.dist = lexical_cast<int>(fields[3]);
	  if (m == 5)
		res.w = lexical_cast<float>(fields[4]);
	} catch (boost::bad_lexical_cast &) {
	  throw std::logic_error(word + " : Parsing error.");
	}

	if (res.w < 0.0) {
	  throw std::logic_error(word + " : Negative weight.");
	}
	return res;
  }

  static void push_ctx(const vector<string> & ctx,
					   vector<CWord> & cws,
					   size_t & tgtN) {

	CWord *last_nopv = 0;
	map<string, float> nopv_concepts;
	vector<string>::const_iterator end = ctx.end();
	for(vector<string>::const_iterator it = ctx.begin();
		it != end or last_nopv; ++it) {
	  try {
		if (it == end) {
		  // last check to fill last nopv
		  if (!last_nopv->set_concepts(nopv_concepts)) { // false means no concepts attached
			if (last_nopv->is_tgtword()) {
			  tgtN--;
			}
			cws.pop_back();
		  }
		  break;
		}
		ctw_parse_t ctwp = parse_ctw(*it);
		if (ctwp.lemma.size() == 0) {
		  throw std::logic_error(*it + " has no lemma.");
		}
		string pos("");
		CWord::cwtype cw_type = cast_int_cwtype(ctwp.dist);
		if (cw_type == CWord::cw_error) {
		  throw std::logic_error(*it + " fourth field is invalid.");
		}
		if ((cw_type == CWord::cw_ctxword || cw_type == CWord::cw_tgtword) && glVars::input::filter_pos) {
		  if (!ctwp.pos.size()) throw std::logic_error(*it + " has no POS.");
		  pos = ctwp.pos;
		}
		if(!glVars::input::weight)
		  ctwp.w = 1.0;

		CWord new_cw(ctwp.lemma, ctwp.id, pos, cw_type, ctwp.w);

		if (last_nopv) {
		  // there is a nopv element 'active'
		  if (cw_type == CWord::cw_concept) {
			// attach concept to nopv
			nopv_concepts.insert(make_pair(ctwp.lemma, ctwp.w));
			continue;
		  }
		  // new elem is not concept, so set last nopv with new concepts
		  if (!last_nopv->set_concepts(nopv_concepts)) {// false means no concepts attached
			if (last_nopv->is_tgtword()) {
			  tgtN--;
			}
			cws.pop_back();
		  }
		  last_nopv = 0; // and reset last_nopv
		}

		if (cw_type == CWord::cw_tgtword_nopv or cw_type == CWord::cw_ctxword_nopv) {
		  // New nopv cword
		  map<string, float>().swap(nopv_concepts);
		  cws.push_back(new_cw);
		  last_nopv = &cws.back();
		  if (new_cw.is_tgtword()) {
			tgtN++;
		  }
		  continue;
		}
		if (!new_cw.size()) {
		  if (glVars::debug::warning)
			// No synset for that word.
			cerr << "W:" << *it << " can't be mapped to KB.";
		  continue;
		}
		cws.push_back(new_cw);
		if (new_cw.is_tgtword()) {
		  tgtN++;
		}
	  } catch (ukb::wdict_error & e) {
		throw e;
	  } catch (std::logic_error & e) {
		string msg(e.what());
		if (!glVars::input::swallow) throw std::runtime_error(msg);
		if (glVars::debug::warning) {
		  cerr << msg << "\n";
		}
	  }
	}
  }

  // AW file read (create a csentence from a context)

  istream & CSentence::read_aw(istream & is, size_t & l_n) {

	string line;
	map<int, map<string, float> > Ties;

	if(read_line_noblank(is, line, l_n)) {

	  // first line is id
	  char_separator<char> sep(" \t");
	  vector<string> ctx;

	  tokenizer<char_separator<char> > tok_id(line, sep);
	  copy(tok_id.begin(), tok_id.end(), back_inserter(ctx));
	  if (ctx.size() == 0) return is; // blank line or EOF
	  m_id = ctx[0];
	  vector<string>().swap(ctx);
	  // next comes the context
	  if(!read_line_noblank(is, line, l_n)) return is;
	  tokenizer<char_separator<char> > tok_ctx(line, sep);
	  copy(tok_ctx.begin(), tok_ctx.end(), back_inserter(ctx));
	  if (ctx.size() == 0) return is; // blank line or EOF
	  try {
		push_ctx(ctx, m_v, m_tgtN);
	  } catch (ukb::wdict_error & e) {
		throw e;
	  } catch (std::exception & e) {
		throw std::runtime_error("Context error in line " + lexical_cast<string>(l_n) + "\n" + e.what());
	  }
	}
	return is;
  }


  // CSentence::CSentence(const vector<string> & sent) : m_id(string()) {

  //   set<string> wordS;

  //   vector<string>::const_iterator s_it = sent.begin();
  //   vector<string>::const_iterator s_end = sent.end();
  //   for(;s_it != s_end; ++s_it) {
  //     bool insertedP;
  //     set<string>::iterator it;
  //     tie(it, insertedP) = wordS.insert(*s_it);
  //     if (insertedP) {
  //       m_v.push_back(CWord(*s_it));
  //     }
  //   }
  // }

  CSentence & CSentence::operator=(const CSentence & cs_) {
	if (&cs_ != this) {
	  m_tgtN = cs_.m_tgtN;
	  m_v = cs_.m_v;
	  m_id = cs_.m_id;
	}
	return *this;
  }

  void CSentence::append(const CSentence & cs_) {
	vector<CWord> tenp(m_v.size() + cs_.m_v.size());
	vector<CWord>::iterator aux_it;
	aux_it = copy(m_v.begin(), m_v.end(), tenp.begin());
	copy(cs_.m_v.begin(), cs_.m_v.end(), aux_it);
	m_v.swap(tenp);
	m_tgtN += cs_.m_tgtN;
  }

  void CSentence::distinguished_synsets(vector<string> & res) const {

	vector<CWord>::const_iterator cw_it, cw_end;
	cw_it = m_v.begin();
	cw_end = m_v.end();
	for(; cw_it != cw_end; ++cw_it) {
	  if (!cw_it->is_tgtword()) continue;
	  for(size_t i = 0, m = cw_it->size();
		  i < m; ++i)
		res.push_back(cw_it->syn(i));
	}
  }

  std::ostream& operator<<(std::ostream & o, const CSentence & cs_) {
	o << cs_.m_id << endl;
	copy(cs_.begin(), cs_.end(), ostream_iterator<CWord>(o, " "));
	o << "\n";
	return o;
  }

  std::ostream & CSentence::debug(std::ostream & o) const {
	o << m_id << endl;
	for(vector<CWord>::const_iterator it = m_v.begin(), end=m_v.end();
		it != end; ++it) {
	  o << "**CWord" << endl;
	  it->debug(o);
	}
	return o;
  }

  std::ostream & CSentence::print_csent(std::ostream & o) const {

	if (!m_tgtN) return o;

	vector<CWord>::const_iterator cw_it = m_v.begin();
	vector<CWord>::const_iterator cw_end = m_v.end();

	for(; cw_it != cw_end; ++cw_it) {
	  if (cw_it->size() == 0) continue;
	  if (!cw_it->is_tgtword()) continue;
	  if (!cw_it->is_disambiguated() && !glVars::output::ties) continue;
	  if (cw_it->is_monosemous() && !glVars::output::monosemous) continue; // Don't output monosemous words
	  o << m_id << " ";
	  cw_it->print_cword(o);
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
		u = cw.V_vector().at(0).first;
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
  // * Initial V is computed by activating the words of the context
  //
  // * Dispatch to proper implementation using glVars::pprImpl global variable

  bool calculate_kb_ppr(const CSentence & cs,
						vector<float> & res,
						CSentence::const_iterator tgtw_it) {

	if (tgtw_it == CSentence::const_iterator()) tgtw_it = cs.end();
	Kb & kb = ukb::Kb::instance();
	vector<float> pv;
	int aux = pv_from_cs_onlyC(cs, pv, tgtw_it);
	const std::vector<float> pers = pv;
	if (!aux) return false;
	switch(glVars::prank::impl) {
	case glVars::pm:
	  kb.pageRank_ppv(pv, res); // power method
	  break;
	case glVars::mc_complete:
	  kb.monte_carlo_complete(glVars::prank::damping,
							  pers,
							  glVars::prank::mc_m,
							  res); // monte carlo complete
	  break;
	case glVars::mc_end:
	  kb.monte_carlo_end_point_cyclic(glVars::prank::damping,
									  pers,
									  glVars::prank::mc_m,
									  res); // monte carlo endpoint
	  break;
	}
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

	if (!cs.has_tgtwords()) return false; // no target words
	return calculate_kb_ppr(cs, ranks, tgtw_it);
  }

  // given a CSentence
  // for each target word
  //   1. put a ppv in the synsets of the rest of words.
  //   2. Pagerank
  //   3. use rank for disambiguating word

  int calculate_kb_ppr_by_word_and_disamb(CSentence & cs) {

	if (!cs.has_tgtwords()) return 0; // no target words

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

	if (!cs.has_tgtwords()) return false; // no target words

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

  bool disamb_csentence_kb(CSentence & cs,
						   const vector<float> & ranks) {

	if (!cs.has_tgtwords()) return false; // no target words

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
	return true;
  }
}
