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

	static bool cwtype_is_tgtword(CWord::cwtype t) {
		return t == CWord::cw_tgtword || t == CWord::cw_tgtword_nopv;
	}

	bool CWord::is_tgtword() const {
		return cwtype_is_tgtword(m_type);
	}

	void CWord::empty_synsets() {
		std::vector<std::string>().swap(m_syns);
		std::vector<std::pair<Kb::vertex_descriptor, float> >().swap(m_V);
		std::vector<float>().swap(m_ranks);
		m_linkw_factor = 1;
		m_disamb = false;
	}

	void CWord::swap(CWord & o) {
		m_w.swap(o.m_w);
		m_id.swap(o.m_id);
		m_pos.swap(o.m_pos);
		std::swap(m_weight, o.m_weight);
		m_syns.swap(o.m_syns);
		m_V.swap(o.m_V);
		m_ranks.swap(o.m_ranks);
		std::swap(m_linkw_factor, o.m_linkw_factor);
		std::swap(m_type, o.m_type);
		std::swap(m_disamb, o.m_disamb);
	}

	size_t CWord::link_dict_concepts(const string & lemma, const string & pos) {

		size_t new_c = 0;

		set<Kb::vertex_descriptor> syns_S;
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
			boost::random_number_generator<boost::mt19937, long int> rand_dist(glVars::rnd::urng);
			std::random_shuffle(sidxV.begin(), sidxV.end(), rand_dist);
		}

		for(size_t i= 0, m = sidxV.size(); i < m; ++i) {
			size_t idx = sidxV[i];
			Kb::vertex_descriptor syn_v = entries.get_entry(idx);
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
		for(vector<pair<Kb::vertex_descriptor, float> >::iterator it = m_V.begin(), end = m_V.end();
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
			Kb::vertex_descriptor u;
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

		std::vector<std::string> syns;
		std::vector<std::pair<Kb::vertex_descriptor, float> > V;
		Kb::vertex_descriptor u;
		bool P;
		float total_w = 0.0f;
		size_t N = 0;
		for(map<string, float>::iterator it = C.begin(), end = C.end();
			it != end; ++it) {
			tie(u, P) = ukb::Kb::instance().get_vertex_by_name(it->first);
			if (!P) {
				if (glVars::debug::warning)
					// No synset for that word.
					cerr << "reset_concepts: " + it->first + " not in KB\n";
				continue;
			}
			N++;
			syns.push_back(it->first);
			V.push_back(make_pair(u, it->second));
			total_w += it->second;
		}
		if (total_w == 0.0f) {
			this->empty_synsets();
			return 0;
		}

		// Update concepts
		m_linkw_factor = 1.0 / total_w;
		m_syns.swap(syns);
		m_V.swap(V);
		vector<float>(N, 0.0).swap(m_ranks);
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

		o << m_w << "#";
		if (m_pos.size())
			o << m_pos;
		o << "#" << m_id << "#" << m_type;
		o << "#" << lexical_cast<string>(m_weight) << "#" << lexical_cast<string>(m_linkw_factor) << "\t";
		for(size_t i = 0; i < m_syns.size(); i++) {
			o << "[" << m_syns[i] << ", " << m_V[i].first << ", " << m_V[i].second << "] ";
		}
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

		//Kb::boost_graph_t & g = ukb::Kb::instance().graph();

		return cw_.write(o, cw_.m_id, cw_.m_type);
	}

	std::ostream & CWord::write(std::ostream & o, const string & id_, CWord::cwtype t_) const {

		o << m_w << "#";
		if (m_pos.size())
			o << m_pos;
		o << "#" << id_ << "#" << t_;
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

	ostream & CWord::print_cword(ostream & o, const std::string &id) const {

		o << id << " ";
		if(!glVars::output::allranks) cw_aw_print_best(o, m_syns, m_ranks);
		else cw_aw_print_all(o, m_syns, m_ranks);
		o << " !! " << m_w << "\n";
		return o;
	}

	////////////////////////////////////////////////////////////////
	// CSentence

	CSentence::CSentence(const std::string & id, const std::string & ctx_str) :
		m_tgtN(0), m_weight(0.0), m_w_factor(0.0f), m_id(id) {
		if (m_id.empty()) throw std::runtime_error(string("empty id"));
		vector<string> ctx;
		vector<CWord> V;
		char_separator<char> sep(" \t");
		tokenizer<char_separator<char> > tok_ctx(ctx_str, sep);
		copy(tok_ctx.begin(), tok_ctx.end(), back_inserter(ctx));
		if (ctx.size() == 0) return;
		try {
			push_ctx(ctx);
		} catch (ukb::wdict_error & e) {
			throw e;
		} catch (std::exception & e) {
			throw std::runtime_error(string("context error: ") + e.what());
		}
	}

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


	// NOTE: destroys new_cw
	void CSentence::push_cw(CWord & new_cw,
							 map<string, size_t > & CW,
							 bool is_nopv) {
		map<string, size_t >::iterator cwit;
		bool P;
		string id = new_cw.id();
		CWord::cwtype type = new_cw.type();
		bool is_tgtword = new_cw.is_tgtword();
		float w = new_cw.m_weight;

		tie(cwit, P) = CW.insert(make_pair(new_cw.wpos(), m_vuniq.size()));
		if (P) {
			m_vuniq.push_back(CWord());
			m_vuniq.back().swap(new_cw);
		} else {
			CWord & old = m_vuniq[cwit->second];
			// check type compatibility
			bool old_is_nopv = old.type() == CWord::cw_tgtword_nopv || old.type() == CWord::cw_ctxword_nopv;
			if ((is_nopv && !old_is_nopv) ||
				(!is_nopv && old_is_nopv))
				throw std::logic_error(string("push error: word " + old.wpos() + " has incompatible types"));
			// if current is target word, set original to target word
			if (is_tgtword)
				old.m_type = type;
			old.m_weight += w; // aggregate weights
		}
		m_tokens.push_back(cwtoken_t(id, type, cwit->second));
		if (is_tgtword) m_tgtN++;
		m_weight += w;
	}

	void CSentence::push_ctx(const vector<string> & ctx) {

		map<string, size_t > CW; // words inserted so far
		bool last_is_nopv = false;
		CWord last_nopv;
		map<string, float> nopv_concepts;
		vector<string>::const_iterator end = ctx.end();
		for(vector<string>::const_iterator it = ctx.begin();
			it != end or last_is_nopv; ++it) {
			try {
				if (it == end) {
					// last check to fill last nopv
					last_is_nopv = false;
					if (last_nopv.set_concepts(nopv_concepts)) { // false means no concepts attached
						map<string, float>().swap(nopv_concepts);
						push_cw(last_nopv, CW, true);
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
				if (glVars::input::filter_pos & (cw_type == CWord::cw_ctxword || cw_type == CWord::cw_tgtword)) {
					if (!ctwp.pos.size()) throw std::logic_error(*it + " has no POS.");
					pos = ctwp.pos;
				}
				if(!glVars::input::weight)
					ctwp.w = 1.0;

				if (last_is_nopv) {
					// there is a nopv element 'active'
					if (cw_type == CWord::cw_concept) {
						// attach concept to nopv
						nopv_concepts.insert(make_pair(ctwp.lemma, ctwp.w));
						continue;
					}
					last_is_nopv = false;
					// new elem is not concept, so push last nopv
					if (last_nopv.set_concepts(nopv_concepts)) {
						map<string, float>().swap(nopv_concepts);
						--it; // if push_cw throws, make sure current word is processed again
						push_cw(last_nopv, CW, true);
						++it; // did not throw anyway
					}
				}

				// last_is_nopv == false
				CWord new_cw(ctwp.lemma, ctwp.id, pos, cw_type, ctwp.w);
				if (cw_type == CWord::cw_tgtword_nopv or cw_type == CWord::cw_ctxword_nopv) {
					// New nopv cword
					new_cw.swap(last_nopv);
					last_is_nopv = true;
					continue;
				}
				if (!new_cw.size()) {
					if (glVars::debug::warning)
						// No synset for that word.
						cerr << "W:" << *it << " can't be mapped to KB.\n";
					continue;
				}
				push_cw(new_cw, CW, false);
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
		if (m_vuniq.size()) {
			if (m_weight == 0.0)
				throw std::logic_error(string("Context with zero weight"));
			m_w_factor = 1.0 / m_weight;
		}
	}

	std::ostream& operator<<(std::ostream & o, const CSentence & cs_) {
		o << cs_.m_id << endl;
		vector<CSentence::cwtoken_t>::const_iterator cw_it = cs_.m_tokens.begin();
		vector<CSentence::cwtoken_t>::const_iterator cw_end = cs_.m_tokens.end();

		cw_end--;
		for(; cw_it != cw_end; ++cw_it) {
			cs_.m_vuniq[cw_it->idx].write(o, cw_it->id, cw_it->t);
			o << " ";
		}
		cs_.m_vuniq[cw_end->idx].write(o, cw_end->id, cw_end->t);
		o << " ";
		return o;
	}

	std::ostream & CSentence::debug(std::ostream & o) const {
		o << m_id << " tgtN:" << lexical_cast<string>(m_tgtN) <<
			" W:" << lexical_cast<string>(m_weight) << "(" << lexical_cast<string>(m_w_factor) << ")" << endl;
		o << "Unique:\n";
		for(vector<CWord>::const_iterator it = m_vuniq.begin(), end=m_vuniq.end();
			it != end; ++it) {
			o << "  ";
			it->debug(o);
			o << endl;
		}
		o << "Tokens:" << endl;
		vector<cwtoken_t>::const_iterator cw_it = m_tokens.begin();
		vector<cwtoken_t>::const_iterator cw_end = m_tokens.end();
		for(; cw_it != cw_end; ++cw_it) {
			o << "[" << cw_it->idx << ", " << cw_it->id << ", " << cw_it->t << "] ";
		}
		o << endl;
		return o;
	}

	std::ostream & CSentence::print_csent(std::ostream & o) const {

		if (!m_tgtN) return o;

		vector<cwtoken_t>::const_iterator cw_it = m_tokens.begin();
		vector<cwtoken_t>::const_iterator cw_end = m_tokens.end();

		for(; cw_it != cw_end; ++cw_it) {
			const CWord & cw = m_vuniq[cw_it->idx];
			if (cw.size() == 0) continue;
			if (!cwtype_is_tgtword(cw_it->t)) continue;
			if (!cw.is_disambiguated() && !glVars::output::ties) continue;
			if (cw.is_monosemous() && !glVars::output::monosemous) continue; // Don't output monosemous words
			o << m_id << " ";
			cw.print_cword(o, cw_it->id);
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


	size_t update_pv_cw(const vector<pair<Kb::vertex_descriptor, float> > & m_V,
						float factor,
						vector<float> & pv) {

		// Sum edge weights

		size_t inserted = 0;
		// Uppdate PV
		for(vector<pair<Kb::vertex_descriptor, float> >::const_iterator it = m_V.begin(), end = m_V.end();
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

		Kb::vertex_descriptor u;
		int inserted_i = 0;

		// put pv to the synsets of words
		for(vector<CWord>::const_iterator it = cs.ubegin(), end = cs.uend();
			it != end; ++it) {
			if (it == exclude_word_it) continue;
			const CWord & cw = *it;
			float cw_w = cw.get_weight() * cs.weigth_factor();
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
	// Initial V is computed by activating the words of the context

	bool calculate_kb_ppr(const CSentence & cs,
						  vector<float> & ranks) {

		return calculate_kb_ppr_by_word(cs, cs.uend(), ranks);
	}

	// given a word (pointed by tgtw_it),
	// 1. put a ppv in the synsets of the rest of words.
	// 2. Pagerank
	//
	// Note: tgtw_it can point to end(), effectively turning this function into
	// a 'standard' PPR

	bool calculate_kb_ppr_by_word(const CSentence & cs,
								  CSentence::const_iterator tgtw_it,
								  vector<float> & ranks) {

		Kb & kb = ukb::Kb::instance();
		vector<float> pv;
		int aux = pv_from_cs_onlyC(cs, pv, tgtw_it);
		// Execute PageRank
		if (aux) {
			kb.pageRank_ppv(pv, ranks);
			if (glVars::csentence::disamb_minus_static) {
				const vector<float> & staticV = kb.static_prank();
				for(size_t i = 0, n = staticV.size();
					i != n; ++i) {
					ranks[i] = staticV[i] - ranks[i];
				}
			}
		}
		return aux;
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

		vector<CWord>::iterator cw_it = cs.ubegin();
		vector<CWord>::iterator cw_end = cs.uend();
		float K = 0.0;
		for(; cw_it != cw_end; ++cw_it) {
			for(size_t i = 0; i != cw_it->size(); ++i) {
				Kb::vertex_descriptor u;
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

		vector<CWord>::iterator cw_it = cs.ubegin();
		vector<CWord>::iterator cw_end = cs.uend();
		for(; cw_it != cw_end; ++cw_it) {
			if (!cw_it->is_tgtword()) continue;
			cw_it->rank_synsets(ranks, false);
			cw_it->disamb_cword();
		}
		return true;
	}
}
