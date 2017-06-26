#include "kbGraph.h"
#include "common.h"
#include "globalVars.h"
#include "prank.h"

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <string>
#include <map>
#include <iterator>
#include <algorithm>
#include <ostream>

// Tokenizer
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

// Stuff for generating random numbers

#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>

// bfs

#include <boost/graph/visitors.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/pending/indirect_cmp.hpp>

#if BOOST_VERSION > 104400
#include <boost/range/irange.hpp>
#else
#include <boost/pending/integer_range.hpp>
#endif

#include <boost/graph/graph_utility.hpp> // for boost::make_list

// dijkstra

#include <boost/graph/dijkstra_shortest_paths.hpp>

// strong components

#include <boost/graph/strong_components.hpp>


namespace ukb {

	using namespace std;
	using namespace boost;

	////////////////////////////////////////////////////////////////////////////////
	// Class Kb


	////////////////////////////////////////////////////////////////////////////////
	// Singleton stuff

	Kb* Kb::p_instance = 0;

	Kb *Kb::create() {

		static Kb theKb;
		return &theKb;
	}

	Kb & Kb::instance() {
		if (!p_instance) {
			throw runtime_error("KB not initialized");
		}
		return *p_instance;
	}

	void Kb::create_from_txt(const string & synsFileName,
							 const std::set<std::string> & src_allowed) {
		if (p_instance) return;
		Kb *tenp = create();
		tenp->read_from_txt(synsFileName, src_allowed);
		p_instance = tenp;
	}

	void Kb::create_from_txt(std::istream & is,
							 const std::set<std::string> & src_allowed) {
		if (p_instance) return;
		Kb *tenp = create();
		tenp->read_from_txt(is, src_allowed);
		p_instance = tenp;
	}

	void Kb::create_from_binfile(const std::string & fname) {

		if (p_instance) return;
		Kb *tenp = create();

		if (!fname.size())
			throw std::runtime_error(string("[E] loading KB: no KB name"));

		ifstream fi(fname.c_str(), ifstream::binary|ifstream::in);
		if (!fi)
			throw std::runtime_error(string("[E] loading KB: can not open ") + fname);

		tenp->read_from_stream(fi);
		p_instance = tenp;
	}


	void Kb::create_from_kbgraph16(Kb16 & kbg) {
		if (p_instance) return;
		Kb *tenp = create();
		precsr_t precsr16;

		Kb16::boost_graph_t oldg = kbg.g;
		graph_traits<Kb16::boost_graph_t>::edge_iterator eit, eend;
		tie(eit, eend) = edges(oldg);
		for(; eit != eend; ++eit) {
			string ustr(get(vertex_name, oldg, source(*eit, oldg)));
			string vstr(get(vertex_name, oldg, target(*eit, oldg)));
			precsr16.insert_edge(ustr, vstr, get(edge_weight, oldg, *eit), get(edge_rtype, oldg, *eit));
		}

		boost_graph_t *new_g = new boost_graph_t(boost::edges_are_unsorted_multi_pass,
												 precsr16.E.begin(), precsr16.E.end(),
												 precsr16.eProp.begin(),
												 precsr16.m_vsize);

		tenp->m_g.reset(new_g);

		BGL_FORALL_VERTICES(v, *(tenp->m_g), Kb::boost_graph_t) {
			(*(tenp->m_g))[v].name = precsr16.vProp[v].name;
		}

		tenp->m_vertexN = num_vertices(*(tenp->m_g));
		tenp->m_edgeN = num_edges(*(tenp->m_g));
		// relation sources
		std::set<std::string>(kbg.relsSource).swap(tenp->m_relsSource);
		// vertex map
		tenp->m_synsetMap.swap(precsr16.m_vMap);
		// relation types
		tenp->m_rtypes.m_strtypes.swap(kbg.rtypes);
		// Notes
		tenp->m_notes = kbg.notes;
		tenp->m_notes.push_back("--");
		tenp->m_notes.push_back("converted_to_2.0");

		p_instance = tenp;
	}


	////////////////////////////////////////////////////////////////////////////////


	void Kb::add_comment(const string & str) {
		m_notes.push_back(str);
	}

	const vector<string> & Kb::get_comments() const {
		return m_notes;
	}

	////////////////////////////////////////////////////////////////////////////////
	// bfs

	struct kb_bfs_init:public base_visitor<kb_bfs_init> {
	public:
		kb_bfs_init(Kb::vertex_descriptor *v):m_v(v) { }
		typedef on_initialize_vertex event_filter;
		inline void operator()(Kb::vertex_descriptor u, const Kb::boost_graph_t & g)
		{
			m_v[u] = u;
		}
		Kb::vertex_descriptor *m_v;
	};

	struct kb_bfs_pred:public base_visitor<kb_bfs_pred> {
	public:
		kb_bfs_pred(Kb::vertex_descriptor *v):m_v(v) { }
		typedef on_tree_edge event_filter;
		inline void operator()(Kb::edge_descriptor e, const Kb::boost_graph_t & g) {
			m_v[target(e, g)] = source(e, g);
		}
		Kb::vertex_descriptor *m_v;
	};


	// Note:
	//
	// after bfs, if (parents[v] == v) and (v != u), then u and v are not
	// connected in the graph.

	bool Kb::bfs(Kb::vertex_descriptor src,
				 std::vector<Kb::vertex_descriptor> & parents) const {

		size_t m = num_vertices(*m_g);
		if(parents.size() == m) {
			std::fill(parents.begin(), parents.end(), Kb::vertex_descriptor());
		} else {
			vector<Kb::vertex_descriptor>(m).swap(parents);  // reset parents
		}

		breadth_first_search(*m_g,
							 src,
							 boost::visitor(boost::make_bfs_visitor
											(boost::make_list(kb_bfs_init(&parents[0]),
															  kb_bfs_pred(&parents[0])))));
		return true;
	}


	bool Kb::dijkstra (Kb::vertex_descriptor src,
					   std::vector<Kb::vertex_descriptor> & parents) const {

		size_t m = num_vertices(*m_g);
		if(parents.size() == m) {
			std::fill(parents.begin(), parents.end(), Kb::vertex_descriptor());
		} else {
			vector<Kb::vertex_descriptor>(m).swap(parents);  // reset parents
		}

		// Hack to remove const-ness
		Kb & me = const_cast<Kb &>(*this);
		vector<float> w;
		vector<float> dist(m);
		property_map<Kb::boost_graph_t, boost::vertex_index_t>::type indexmap = get(vertex_index, *m_g);
		property_map<Kb::boost_graph_t, float edge_prop_t::*>::type wmap = get(&edge_prop_t::weight, *(me.m_g));

		dijkstra_shortest_paths(*m_g,
								src,
								predecessor_map(make_iterator_property_map(parents.begin(),
																		   get(vertex_index, *m_g))).
								distance_map(make_iterator_property_map(dist.begin(),
																		get(vertex_index, *m_g))).
								weight_map(wmap).
								vertex_index_map(indexmap));

		return true;
	}

	////////////////////////////////////////////////////////////////////////////////
	// Get shortest subgraphs

	class bfs_subg_terminate : public std::exception {};

	struct subg {
		vector<Kb::vertex_descriptor> V;
		vector<vector<Kb::vertex_descriptor> > E;
	};


	class bfs_subg_visitor : public default_bfs_visitor {

	public:
		bfs_subg_visitor(subg & s_, Kb::vertex_descriptor u, int limit)
			: m_sg(s_), m_idx(), m_i(0), m_t(0), m_max(limit) {
			add_v(u);
		}

		void tree_edge(Kb::edge_descriptor e, const Kb::boost_graph_t & g)
		{

			Kb::vertex_descriptor u = source(e,g);
			Kb::vertex_descriptor v = target(e,g);

			// vertex v is new, but yet undiscovered
			int v_i = add_v(v);
			if (v_i == -1) return; // max limit reached.
			int u_i = get_v(u);
			add_e(u_i, v_i);

			Kb::edge_descriptor aux;
			bool existsP;
			tie(aux, existsP) = edge(v, u, g);
			if (existsP) add_e(v_i, u_i); // as this edge is no more traversed.
		}

		void non_tree_edge(Kb::edge_descriptor e, const Kb::boost_graph_t & g)
		{
			// cross edge. source is previously stored for sure. target probably
			// is, unless max limit was reached
			int v_i = get_v(target(e,g));
			if (v_i == -1) return; // target vertex not stored because max limit.
			int u_i = get_v(source(e,g));

			add_e(u_i, v_i);
		}

		void discover_vertex(Kb::vertex_descriptor u, const Kb::boost_graph_t & g)
		{
			if (m_t == m_max) throw bfs_subg_terminate();
			++m_t;
		}

		int add_v(Kb::vertex_descriptor v)
		{
			if(m_i == m_max) return -1;
			m_sg.V.push_back(v);
			m_idx[v] = m_i;
			m_sg.E.push_back(vector<Kb::vertex_descriptor>());
			int res = m_i;
			++m_i;
			return res;
		}

		int get_v(Kb::vertex_descriptor v) {
			map<Kb::vertex_descriptor, int>::iterator it=m_idx.find(v);
			if(it == m_idx.end()) return -1;
			return it->second;
		}

		void add_e(int u_i, int v_i) {
			m_sg.E[u_i].push_back(m_sg.V[v_i]);
		}

	private:
		subg & m_sg;
		map<Kb::vertex_descriptor, int> m_idx;
		int m_i; // num of inserted vertices
		int m_t; // time
		int m_max;
	};


	void Kb::get_subgraph(const string & src,
						  vector<string> & V,
						  vector<vector<string> > & E,
						  size_t limit) {

		Kb::vertex_descriptor u;
		bool aux;
		tie(u,aux) = get_vertex_by_name(src);
		if(!aux) return;

		subg sg;
		bfs_subg_visitor vis(sg, u, limit);

		try {
			breadth_first_search(*m_g, u, boost::visitor(vis));
		} catch (bfs_subg_terminate & ) {}

		size_t N = sg.V.size();
		vector<string>(N).swap(V);
		vector<vector<string> >(N).swap(E);

		for(size_t i=0; i < N; ++i) {
			V[i] = (*m_g)[sg.V[i]].name;
			size_t m = sg.E[i].size();
			vector<string> l(m);
			for(size_t j=0; j < m; ++j) {
				l[j] =  (*m_g)[sg.E[i].at(j)].name;
			}
			E[i].swap(l);
		}
	}

	bool Kb::get_shortest_paths(const std::string & src,
								const std::vector<std::string> & targets,
								std::vector<std::vector<std::string> > & paths) {
		vector<Kb::vertex_descriptor> parents;
		Kb::vertex_descriptor u;
		bool aux;
		tie(u,aux) = get_vertex_by_name(src);
		if(!aux) return false;
		std::vector<std::vector<std::string> >().swap(paths);
		this->bfs(u, parents);
		for(std::vector<std::string>::const_iterator it = targets.begin(), end = targets.end();
			it != end; ++it) {
			Kb::vertex_descriptor v;
			tie(v,aux) = get_vertex_by_name(*it);
			if (!aux) continue;
			if (parents[v] == v) continue; // either (u == v) or v is not connected to u.
			paths.push_back(vector<string>());
			vector<string> & P = paths.back();
			// iterate until source is met
			P.push_back(get_vertex_name(v));
			while(1) {
				v = parents[v];
				P.push_back(get_vertex_name(v));
				if (v == u) break;
			}
			std::reverse(P.begin(), P.end());
		}
		return paths.size();
	}


	////////////////////////////////////////////////////////////////////////////////
	// strings <-> vertex_id

	pair<Kb::vertex_descriptor, bool> Kb::get_vertex_by_name(const std::string & str) const {
		map<string, Kb::vertex_descriptor>::const_iterator it;

		it = m_synsetMap.find(str);
		if (it != m_synsetMap.end()) return make_pair(it->second, true);
		return make_pair(Kb::vertex_descriptor(), false);
	}

	void Kb::edge_add_reltype(Kb::edge_descriptor e, const string & rel) {
		m_rtypes.add_type(rel, (*m_g)[e].etype);
	}

	std::vector<std::string> Kb::edge_reltypes(Kb::edge_descriptor e) const {
		return m_rtypes.tvector((*m_g)[e].etype);
	}

	float Kb::get_edge_weight(Kb::edge_descriptor e) const {
		if (!glVars::prank::use_weight) return 1.0f;
		return (*m_g)[e].weight;
	}

	void Kb::set_edge_weight(Kb::edge_descriptor e, float w) {
		(*m_g)[e].weight = w;
	}

	std::pair<Kb::out_edge_iterator, Kb::out_edge_iterator> Kb::out_neighbors(Kb::vertex_descriptor u) {
		return out_edges(u, *m_g);
	}

	std::pair<Kb::in_edge_iterator, Kb::in_edge_iterator> Kb::in_neighbors(Kb::vertex_descriptor u) {
		return in_edges(u, *m_g);
	}


	////////////////////////////////////////////////////////////////////////////////
	// Query and retrieval

	// filter_mode
	//   0 -> no filter
	//   1 -> only words
	//   2 -> only concepts


	void vname_filter(const map<string, Kb::vertex_descriptor> & theMap,
					  const vector<float> & ranks,
					  const Kb::boost_graph_t & g,
					  vector<float> & outranks,
					  vector<string> & vnames) {

		size_t v_m = theMap.size();

		vector<Kb::vertex_descriptor> V(v_m);
		// empty output vectors
		vector<float>(v_m).swap(outranks);
		vector<string>(v_m).swap(vnames);

		map<string, Kb::vertex_descriptor>::const_iterator m_it = theMap.begin();
		map<string, Kb::vertex_descriptor>::const_iterator m_end = theMap.end();
		size_t v_i = 0;

		// Fill vertices index vector
		for(; m_it != m_end; ++m_it, ++v_i) {
			V[v_i] = m_it->second;
		}
		// Sort to guarantee uniqueness
		sort(V);

		for(v_i = 0; v_i < v_m; ++v_i) {
			outranks[v_i] = ranks[V[v_i]];
			vnames[v_i] = g[V[v_i]].name;
		}
	}


	void Kb::filter_ranks_vnames(const vector<float> & ranks,
								 vector<float> & outranks,
								 vector<string> & vnames,
								 int filter_mode) const {

		size_t v_i, v_m;

		// No filtering
		v_m = ranks.size();
		outranks.resize(v_m);
		vnames.resize(v_m);
		for(v_i = 0; v_i < v_m; ++v_i) {
			outranks[v_i] = ranks[v_i];
			vnames[v_i] = (*m_g)[v_i].name;
		}
	}

	////////////////////////////////////////////////////////////////////////////////
	// Get static pageRank vector

	const std::vector<float> & Kb::static_prank() const {
		if (m_static_ppv.size()) return m_static_ppv;

		// Hack to remove const-ness
		Kb & me = const_cast<Kb &>(*this);

		if (m_vertexN == 0) return m_static_ppv; // empty graph
		vector<float> pv(m_vertexN, 1.0/static_cast<float>(m_vertexN));
		me.pageRank_ppv(pv, me.m_static_ppv);
		return m_static_ppv;
	}

	////////////////////////////////////////////////////////////////////////////////
	// Random

	Kb::vertex_descriptor Kb::get_random_vertex() const {

		int r = g_randTarget(num_vertices(*m_g));

		return r;
	}

	////////////////////////////////////////////////////////////////////////////////
	// read from textfile and create graph


	// Line format:
	//
	// u:synset v:synset t:rel i:rel s:source d:directed w:weight
	//
	// u: source vertex. Mandatory.
	// v: target vertex. Mandatory.
	// t: relation type (hyperonym, meronym, etc) of edge u->v. Optional.
	// i: (inverse) relation type of edge v->u (hyponym, etc). Optional. Useless on undirected graphs.
	// s: source of relation (wn30, kb17, etc). Optional.
	// d: wether the relation is directed. Optional, default is undirected.
	// w: relation weight. Must be positive. Optional.


	struct rel_parse {
		string u;
		string v;
		string rtype;
		string irtype;
		string src;
		float w;
		bool directed;

		rel_parse() : u(), v(), rtype(), irtype(), src(), w(0.0), directed(false) {}

	};

	bool parse_line(const string & line, rel_parse & out) {

		rel_parse res;

		char_separator<char> sep(" \t");
		tokenizer<char_separator<char> > tok(line, sep);
		tokenizer<char_separator<char> >::iterator it = tok.begin();
		tokenizer<char_separator<char> >::iterator end = tok.end();
		if (it == end) return false; // empty line
		for(;it != end; ++it) {

			string str = *it;
			if (str.length() < 3 || str[1] != ':') {
				throw runtime_error("parse_line error. Malformed line: " + line);
			}
			char f = str[0];
			string val = str.substr(2);
			if(!val.size()) continue;
			switch (f) {
			case 'u':
				res.u = val;
				break;
			case 'v':
				res.v = val;
				break;
			case 't':
				res.rtype = val;
				break;
			case 'i':
				res.irtype = val;
				break;
			case 's':
				res.src = val;
				break;
			case 'w':
				res.w = lexical_cast<float>(val);
				break;
			case 'd':
				res.directed = glVars::kb::keep_directed && lexical_cast<bool>(val);
				break;
			default:
				throw runtime_error("parse_line error. Unknown value " + str);
				break;
			}
		}
		if (!res.u.size()) throw runtime_error("parse_line error. No source vertex.");
		if (!res.v.size()) throw runtime_error("parse_line error. No target vertex.");
		out = res;
		return true;
	}

	void Kb::read_from_txt(istream & kbFile,
						   const set<string> & src_allowed) {
		string line;
		size_t line_number = 0;
		precsr_t csr_pre;

		set<string>::const_iterator srel_end = src_allowed.end();
		while(kbFile) {
			vector<string> fields;
			read_line_noblank(kbFile, line, line_number);
			if(!kbFile) continue;
			if (line[0] == '#') continue;
			rel_parse f;
			try {
				if (!parse_line(line, f)) continue;

				if (glVars::kb::filter_src) {
					if (src_allowed.find(f.src) == srel_end) continue; // Skip this relation
				}

				if (f.u == f.v) continue; // no self-loops

				if (f.src.size()) {
					this->add_relSource(f.src);
				}

				float w = f.w ? f.w : 1.0;
				// add edge

				// relation type
				// empty f.rtype unless glVars::kb::keep_reltypes

				if (!glVars::kb::keep_reltypes)
					string().swap(f.rtype);

				csr_pre.insert_edge(f.u, f.v, w, f.rtype);

				// Insert v->u if undirected relation

				if (!f.directed || !glVars::kb::keep_directed) {
					csr_pre.insert_edge(f.v, f.u, w, f.rtype);
				}
			} catch (std::exception & e) {
				string msg(string(e.what()) + " in line " + lexical_cast<string>(line_number));
				if(!glVars::input::swallow) throw std::runtime_error(msg);
				if (glVars::debug::warning) {
					cerr << msg << " (Skipping)\n";
				}
			}
		}

		Kb::boost_graph_t *new_g = new Kb::boost_graph_t(boost::edges_are_unsorted_multi_pass,
														 csr_pre.E.begin(), csr_pre.E.end(),
														 csr_pre.eProp.begin(),
														 csr_pre.m_vsize);

		m_g.reset(new_g);
		// add_edges(csr_pre.E.begin(), csr_pre.E.end(),
		//		  // csr_pre.eProp.begin(),
		//		  // csr_pre.eProp.end(),
		//		  g);

		BGL_FORALL_VERTICES(v, *m_g, Kb::boost_graph_t) {
			(*m_g)[v].name = csr_pre.vProp[v].name;
		}

		m_vertexN = num_vertices(*m_g);
		m_edgeN = num_edges(*m_g);
		// Init vertex map
		m_synsetMap.swap(csr_pre.m_vMap);
	}

	void Kb::read_from_txt(const std::string & synsFileName,
						   const set<string> & src_allowed) {

		std::ifstream input_ifs(synsFileName.c_str(), ofstream::in);
		if (!input_ifs) {
			throw runtime_error("Kb::read_from_txt error: Can't open " + synsFileName);
		}
		if(glVars::kb::v1_kb) {
			throw runtime_error(synsFileName + " :sorry, the binary representation has an old format.");
		} else {
			std::istream input_is(input_ifs.rdbuf());
			this->read_from_txt(input_is, src_allowed);
		}
	}

	void Kb::display_info(std::ostream & o) const {

		o << "Relation sources: ";
		writeS(o, m_relsSource);
		if (m_notes.size()) {
			o << "\nM_Notes: ";
			writeV(o, m_notes);
		}
		size_t edge_n = num_edges(*m_g);

		o << "\n" << num_vertices(*m_g) << " vertices and " << edge_n << " edges.\n(Note that if graph is undirected you should divide the edge number by 2)" << endl;
		if (m_rtypes.size()) {
			o << "Relations:";
			writeV(o, m_rtypes.m_strtypes);
			o << endl;
		}
	}


	std::pair<size_t, size_t> Kb::indeg_maxmin() const {

		size_t m = std::numeric_limits<size_t>::max();
		size_t M = std::numeric_limits<size_t>::min();

		size_t d = 0;

		graph_traits<Kb::boost_graph_t>::vertex_iterator it, end;
		tie(it, end) = vertices(*m_g);
		for(; it != end; ++it) {
			d = in_degree(*it, *m_g);
			if (d > M) M = d;
			if (d < m) m = d;
		}
		return make_pair(M, m);
	}

	std::pair<size_t, size_t> Kb::outdeg_maxmin() const {

		size_t m = std::numeric_limits<size_t>::max();
		size_t M = std::numeric_limits<size_t>::min();

		size_t d;

		graph_traits<Kb::boost_graph_t>::vertex_iterator it, end;
		tie(it, end) = vertices(*m_g);
		for(; it != end; ++it) {
			d = out_degree(*it, *m_g);
			if (d > M) M = d;
			if (d < m) m = d;
		}
		return make_pair(M, m);
	}


	int Kb::components() const {

		vector<size_t> v(num_vertices(*m_g));
		boost::iterator_property_map<
			std::vector<size_t>::iterator,
			boost::property_map<Kb::boost_graph_t, boost::vertex_index_t>::type> pm(v.begin(), get(boost::vertex_index, *m_g));

		int i = boost::strong_components(*m_g, pm);

		return i;

	}

	int Kb::components(vector<size_t> & v) const {

		vector<size_t> (num_vertices(*m_g)).swap(v);
		boost::iterator_property_map<
			std::vector<size_t>::iterator,
			boost::property_map<Kb::boost_graph_t, boost::vertex_index_t>::type> pm(v.begin(), get(boost::vertex_index, *m_g));

		int i = boost::strong_components(*m_g, pm);

		return i;

	}

	void Kb::ppv_weights(const vector<float> & ppv) {

		graph_traits<Kb::boost_graph_t>::edge_iterator it, end;

		tie(it, end) = edges(*m_g);
		for(; it != end; ++it) {
			(*m_g)[*it].weight = ppv[target(*it, *m_g)];
		}
	}

	////////////////////////////////////////////////////////////////////////////////
	// PageRank in KB


	// PPV version

	void Kb::pageRank_ppv(const vector<float> & ppv_map,
						  vector<float> & ranks) {

		typedef graph_traits<Kb::boost_graph_t>::edge_descriptor edge_descriptor;
		property_map<Kb::boost_graph_t, float edge_prop_t::*>::type weight_map = get(&edge_prop_t::weight, *m_g);
		prank::constant_property_map <edge_descriptor, float> cte_weight(1.0); // always return 1

		if (0 == m_out_coefs.size()) {
			vector<float>(m_vertexN, 0.0f).swap(m_out_coefs);
			if (glVars::prank::use_weight) {
				prank::init_out_coefs(*m_g,  &m_out_coefs[0], weight_map);
			} else {
				prank::init_out_coefs(*m_g,  &m_out_coefs[0], cte_weight);
			}
		}
		if (m_vertexN == ranks.size()) {
			std::fill(ranks.begin(), ranks.end(), 0.0);
		} else {
			vector<float>(m_vertexN, 0.0).swap(ranks); // Initialize rank vector
		}
		vector<float> rank_tmp(m_vertexN, 0.0);    // auxiliary rank vector

		switch(glVars::prank::impl) {
		  case glVars::pm:
			  if (glVars::prank::use_weight) {
				  prank::do_pageRank(*m_g, m_vertexN, &ppv_map[0],
									 weight_map, &ranks[0], &rank_tmp[0],
									 glVars::prank::num_iterations,
									 glVars::prank::threshold,
									 glVars::prank::damping,
									 m_out_coefs);
			  } else {
				  prank::do_pageRank(*m_g, m_vertexN, &ppv_map[0],
									 cte_weight, &ranks[0], &rank_tmp[0],
									 glVars::prank::num_iterations,
									 glVars::prank::threshold,
									 glVars::prank::damping,
									 m_out_coefs);
			  }
			  break;
		  case glVars::nibble:
			  prank::pageRank_nibble_lazy(*m_g, ppv_map, m_out_coefs, glVars::prank::damping, glVars::prank::nibble_epsilon, ranks);
			  break;
		default:
			cerr << "Error! undefined method for PageRank calculation.\n";
			exit(1);
			break;
		}
	}


	////////////////////////////////////////////////////////////////////////////////
	// Debug

	ostream & Kb::dump_graph(std::ostream & o) const {
		o << "Sources: ";
		writeS(o, m_relsSource);
		o << endl;
		graph_traits<Kb::boost_graph_t>::vertex_iterator it, end;
		tie(it, end) = vertices(*m_g);
		for(;it != end; ++it) {
			o << (*m_g)[*it].name;
			graph_traits<Kb::boost_graph_t>::out_edge_iterator e, e_end;
			tie(e, e_end) = out_edges(*it, *m_g);
			if (e != e_end)
				o << "\n";
			for(; e != e_end; ++e) {
				o << "  ";
				vector<string> r = edge_reltypes(*e);
				writeV(o, r);
				o << " " << (*m_g)[target(*e, *m_g)].name;
				o << " (" << (*m_g)[*e].weight << ")\n";
			}
		}
		return o;
	}

	////////////////////////////////////////////////////////////////////////////////
	// Streaming

	static const size_t magic_id_v1 = 0x070201;
	static const size_t magic_id = 0x080826;
	static const size_t magic_id_csr = 0x110501;

	// CSR read

	vertex_prop_t read_vertex_prop_from_stream(istream & is) {

		string name;

		read_atom_from_stream(is, name);
		return vertex_prop_t(name);
	}

	edge_prop_t read_edge_prop_from_stream(istream & is) {

		float w;
		etype_t::value_type etype;

		read_atom_from_stream(is, w);
		read_atom_from_stream(is, etype);

		return edge_prop_t(w, etype);

	}

	void  Kb::read_from_stream (std::istream & is) {

		size_t vertex_n;
		size_t edge_n;
		size_t id;
		Kb::boost_graph_t *new_g;

		try {
			read_atom_from_stream(is, id);
			if (id != magic_id_csr) {
				if (id == magic_id_v1 || id == magic_id)
					throw runtime_error("Old (pre 2.0) binary serialization format. Convert the graph to new format using the \"convert2.0\" utility.");
				else
					throw runtime_error("Invalid id (same platform used to compile the KB?)");
			}
			read_set_from_stream(is, m_relsSource);
			m_rtypes.read_from_stream(is);
			read_map_from_stream(is, m_synsetMap);

			read_atom_from_stream(is, id);
			if (id != magic_id_csr) {
				throw runtime_error("Invalid id after reading maps");
			}

			read_atom_from_stream(is, edge_n);
			read_atom_from_stream(is, vertex_n);
			read_atom_from_stream(is, id);

			if (id != magic_id_csr) {
				throw runtime_error("Invalid id after reading graph sizes");
			}
			new_g = new Kb::boost_graph_t();

			new_g->m_forward.m_rowstart.resize(0);
			read_vector_from_stream(is, new_g->m_forward.m_rowstart);
			read_vector_from_stream(is, new_g->m_forward.m_column);
			new_g->m_backward.m_rowstart.resize(0);
			read_vector_from_stream(is, new_g->m_backward.m_rowstart);
			read_vector_from_stream(is, new_g->m_backward.m_column);
			read_vector_from_stream(is, new_g->m_backward.m_edge_properties);

			for(size_t i = 0; i != vertex_n; ++i) {
				new_g->vertex_properties().m_vertex_properties.push_back(read_vertex_prop_from_stream(is));
			}

			for(size_t i = 0; i != edge_n; ++i) {
				new_g->m_forward.m_edge_properties.push_back(read_edge_prop_from_stream(is));
			}

			read_atom_from_stream(is, id);
			if (id != magic_id_csr) {
				throw runtime_error("Invalid id after reading graph");
			}
			read_vector_from_stream(is, m_notes);
		} catch (std::exception & e) {
			throw runtime_error(string("Error when reading serialized graph: ") + e.what());
		}

		m_g.reset(new_g);
		vector<float>().swap(m_static_ppv); // empty static rank vector

		m_vertexN = vertex_n;
		m_edgeN = edge_n;
		assert(num_vertices(*m_g) == m_vertexN);
		assert(num_edges(*m_g) == m_edgeN);
	}

	// write

	ostream & write_vertex_prop_to_stream(ostream & o,
										  const vertex_prop_t & p) {
		write_atom_to_stream(o, p.name);
		return o;
	}

	ostream & write_edge_prop_to_stream(ostream & o,
										const edge_prop_t & ep) {
		write_atom_to_stream(o, ep.weight);
		write_atom_to_stream(o, ep.etype);
		return o;
	}

	ostream & Kb::write_to_stream(ostream & o) const {

		// First write maps

		assert(m_vertexN == num_vertices(*m_g));
		assert(m_edgeN == num_edges(*m_g));

		write_atom_to_stream(o, magic_id_csr);

		write_vector_to_stream(o, m_relsSource);
		m_rtypes.write_to_stream(o);
		write_map_to_stream(o, m_synsetMap);

		write_atom_to_stream(o, magic_id_csr);

		write_atom_to_stream(o, m_edgeN);
		write_atom_to_stream(o, m_vertexN);

		write_atom_to_stream(o, magic_id_csr);

		write_vector_to_stream(o, m_g->m_forward.m_rowstart);
		write_vector_to_stream(o, m_g->m_forward.m_column);
		write_vector_to_stream(o, m_g->m_backward.m_rowstart);
		write_vector_to_stream(o, m_g->m_backward.m_column);
		write_vector_to_stream(o, m_g->m_backward.m_edge_properties);

		//	write_vector_to_stream(o, m_g->vertex_properties().m_vertex_properties);
		//	write_vector_to_stream(o, m_g->m_forward.m_edge_properties);

		size_t vProp_n = m_g->vertex_properties().m_vertex_properties.size();
		assert(vProp_n == m_vertexN);
		for(size_t i = 0; i != vProp_n; ++i) {
			write_vertex_prop_to_stream(o, m_g->vertex_properties().m_vertex_properties[i]);
		}

		size_t eProp_n = m_g->m_forward.m_edge_properties.size();
		assert(eProp_n == m_edgeN);
		for(size_t i = 0; i != eProp_n; ++i) {
			write_edge_prop_to_stream(o, m_g->m_forward.m_edge_properties[i]);
		}

		write_atom_to_stream(o, magic_id_csr);

		write_vector_to_stream(o, m_notes);
		return o;
	}


	void Kb::write_to_binfile (const string & fName) {

		ofstream fo(fName.c_str(),  ofstream::binary|ofstream::out);
		if (!fo) {
			cerr << "Error: can't create" << fName << endl;
			exit(-1);
		}
		write_to_stream(fo);
	}

	// text write

	ostream & Kb::write_to_textstream(ostream & o) const {

		graph_traits<Kb::boost_graph_t>::edge_iterator it, end;
		tie(it, end) = edges(*m_g);
		for(;it != end; ++it) {
			const string & u_str = (*m_g)[source(*it, *m_g)].name;
			const string & v_str = (*m_g)[target(*it, *m_g)].name;
			vector<string> r = edge_reltypes(*it);
			if (r.size()) {
				for(vector<string>::const_iterator rit = r.begin(), rend = r.end();
					rit != rend; ++rit) {
					o << "u:" << u_str << " v:" << v_str << " s:" << *rit << " d:1\n";
				}
			} else {
				o << "u:" << u_str << " v:" << v_str << " d:1\n";
			}
		}
		return o;
	}

	void Kb::write_to_textfile (const string & fName) const {

		ofstream fo(fName.c_str(),  ofstream::out);
		if (!fo) {
			cerr << "Error: can't create" << fName << endl;
			exit(-1);
		}
		write_to_textstream(fo);
	}
}
