#include "kbGraph_v16.h"
#include "common.h"
#include "globalVars.h"
#include "wdict.h"
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


	static vector<string>::size_type get_reltype_idx(const string & rel,
													 vector<string> & rtypes) {

		vector<string>::iterator it = rtypes.begin();
		vector<string>::iterator end = rtypes.end();
		vector<string>::size_type idx = 0;

		for(;it != end; ++it) {
			if (*it == rel) break;
			++idx;
		}
		if (it == end) {
			// new relation type
			rtypes.push_back(rel);
		}
		if (idx > 31) {
			throw runtime_error("get_rtype_idx error: too many relation types !");
		}
		return idx;
	}

	////////////////////////////////////////////////////////////////////////////////
	// Class Kb


	////////////////////////////////////////////////////////////////////////////////
	// Singleton stuff

	Kb16* Kb16::p_instance = 0;

	Kb16 *Kb16::create() {

		static Kb16 theKb;
		return &theKb;
	}

	Kb16 & Kb16::instance() {
		if (!p_instance) {
			throw runtime_error("KB not initialized");
		}
		return *p_instance;
	}

	void Kb16::create_from_binfile(const std::string & fname) {

		if (p_instance) return;
		Kb16 *tenp = create();

		ifstream fi(fname.c_str(), ifstream::binary|ifstream::in);
		if (!fi) {
			cerr << "Error: can't open " << fname << endl;
			exit(-1);
		}
		try {
			tenp->read_from_stream(fi);
		} catch(std::exception& e) {
			cerr << e.what() << "\n";
			exit(-1);
		}

		p_instance = tenp;
	}


	////////////////////////////////////////////////////////////////////////////////
	// strings <-> vertex_id

	pair<Kb16_vertex_t, bool> Kb16::get_vertex_by_name(const std::string & str,
													   unsigned char flags) const {
		map<string, Kb16_vertex_t>::const_iterator it;

		if(flags & Kb16::is_concept) {
			it = synsetMap.find(str);
			if (it != synsetMap.end()) return make_pair(it->second, true);
		}
		// is it a word ?
		if(flags & Kb16::is_word) {
			it = wordMap.find(str);
			if (it != wordMap.end()) return make_pair(it->second, true);
		}
		return make_pair(Kb16_vertex_t(), false);
	}

	Kb16_vertex_t Kb16::InsertNode(const string & name, unsigned char flags) {
		coef_status = 0; // reset out degree coefficients
		if (static_ppv.size()) vector<float>().swap(static_ppv); // empty static rank vector
		Kb16_vertex_t u = add_vertex(g);
		put(vertex_name, g, u, name);
		put(vertex_flags, g, u, flags);
		return u;
	}

	Kb16_vertex_t Kb16::find_or_insert_synset(const string & str) {
		map<string, Kb16_vertex_t>::iterator it;
		bool insertedP;
		tie(it, insertedP) = synsetMap.insert(make_pair(str, Kb16_vertex_t()));
		if(insertedP) {
			// new vertex
			it->second = InsertNode(str, Kb16::is_concept);
		}
		return it->second;
	}

	Kb16_vertex_t Kb16::find_or_insert_word(const string & str) {
		map<string, Kb16_vertex_t>::iterator it;
		bool insertedP;
		tie(it, insertedP) = wordMap.insert(make_pair(str, Kb16_vertex_t()));
		if(insertedP) {
			// new vertex
			it->second = InsertNode(str, Kb16::is_word);
		}
		return it->second;
	}

	Kb16_edge_t Kb16::find_or_insert_edge(Kb16_vertex_t u, Kb16_vertex_t v,
										  float w) {

		Kb16_edge_t e;
		bool existsP;

		if (u == v)
			throw runtime_error("Can't insert self loop !");
		//if (w != 1.0) ++w; // minimum weight is 1
		tie(e, existsP) = edge(u, v, g);
		if(!existsP) {
			coef_status = 0; // reset out degree coefficients
			if (static_ppv.size()) vector<float>().swap(static_ppv); // empty static rank vector
			e = add_edge(u, v, g).first;
			put(edge_weight, g, e, w);
			put(edge_rtype, g, e, static_cast<boost::uint32_t>(0));
		}
		return e;
	}

	void Kb16::edge_add_reltype(Kb16_edge_t e, const string & rel) {
		boost::uint32_t m = get(edge_rtype, g, e);
		vector<string>::size_type idx = get_reltype_idx(rel, rtypes);
		m |= (1UL << idx);
		put(edge_rtype, g, e, m);
	}

	std::vector<std::string> Kb16::edge_reltypes(Kb16_edge_t e) const {
		vector<string> res;
		boost::uint32_t m = get(edge_rtype, g, e);
		vector<string>::size_type idx = 0;
		boost::uint32_t i = 1;
		while(idx < 32) {
			if (m & i) {
				res.push_back(rtypes[idx]);
			}
			idx++;
			i <<= 1;
		}
		return res;
	}

	////////////////////////////////////////////////////////////////////////////////
	// Query and retrieval

	bool Kb16::vertex_is_synset(Kb16_vertex_t u) const {
		return !vertex_is_word(u);
	}

	bool Kb16::vertex_is_word(Kb16_vertex_t u) const {
		return (get(vertex_flags, g, u) & Kb16::is_word);
	}

	////////////////////////////////////////////////////////////////////////////////
	// Streaming

	const size_t magic_id_v1 = 0x070201;
	const size_t magic_id = 0x080826;

	// read

	Kb16_vertex_t read_vertex_from_stream_v1(ifstream & is,
											 Kb16Graph & g) {

		string name;

		read_atom_from_stream(is, name);
		Kb16_vertex_t v = add_vertex(g);
		put(vertex_name, g, v, name);
		put(vertex_flags, g, v, 0);
		return v;
	}

	Kb16_edge_t read_edge_from_stream_v1(ifstream & is,
										 Kb16Graph & g) {

		size_t sIdx;
		size_t tIdx;
		float w = 0.0;
		//size_t source;
		bool insertedP;
		Kb16_edge_t e;

		read_atom_from_stream(is, sIdx);
		read_atom_from_stream(is, tIdx);
		read_atom_from_stream(is, w);
		//read_atom_from_stream(is, id);
		//read_atom_from_stream(is, source);
		tie(e, insertedP) = add_edge(sIdx, tIdx, g);
		assert(insertedP);
		put(edge_weight, g, e, w);
		//put(edge_source, g, e, source);

		return e;
	}

	Kb16_vertex_t read_vertex_from_stream(ifstream & is,
										  Kb16Graph & g) {

		string name;
		string gloss;

		read_atom_from_stream(is, name);
		read_atom_from_stream(is, gloss);
		Kb16_vertex_t v = add_vertex(g);
		put(vertex_name, g, v, name);
		put(vertex_flags, g, v, static_cast<unsigned char>(Kb16::is_concept));
		return v;
	}

	Kb16_edge_t read_edge_from_stream(ifstream & is,
									  Kb16Graph & g) {

		size_t sIdx;
		size_t tIdx;
		float w = 0.0;
		boost::uint32_t rtype;
		bool insertedP;
		Kb16_edge_t e;

		read_atom_from_stream(is, sIdx);
		read_atom_from_stream(is, tIdx);
		read_atom_from_stream(is, w);
		read_atom_from_stream(is, rtype);
		//read_atom_from_stream(is, source);
		tie(e, insertedP) = add_edge(sIdx, tIdx, g);
		assert(insertedP);
		put(edge_weight, g, e, w);
		put(edge_rtype, g, e, rtype);

		return e;
	}

	void  Kb16::read_from_stream (std::ifstream & is) {

		size_t vertex_n;
		size_t edge_n;
		size_t i;
		size_t id;

		std::map<std::string, int> relMap_aux;     // Obsolete map from relation name to relation id

		try {
			coef_status = 0;
			vector<float>().swap(static_ppv); // empty static rank vector
			read_atom_from_stream(is, id);
			if (id == magic_id_v1) {

				// Backward compatibility with binary v1 format

				read_set_from_stream(is, relsSource);
				read_map_from_stream(is, relMap_aux);

				read_map_from_stream(is, synsetMap);
				read_map_from_stream(is, wordMap);
				//read_map_from_stream(is, sourceMap);

				read_atom_from_stream(is, id);
				if(id != magic_id_v1) {
					cerr << "Error: invalid id after reading maps" << endl;
					exit(-1);
				}

				read_atom_from_stream(is, vertex_n);
				for(i=0; i<vertex_n; ++i) {
					read_vertex_from_stream_v1(is, g);
				}

				read_atom_from_stream(is, id);
				if(id != magic_id_v1) {
					cerr << "Error: invalid id after reading vertices" << endl;
					exit(-1);
				}

				read_atom_from_stream(is, edge_n);
				for(i=0; i<edge_n; ++i) {
					read_edge_from_stream_v1(is, g);
				}

				read_atom_from_stream(is, id);
				if(id != magic_id_v1) {
					cerr << "Error: invalid id after reading edges" << endl;
					exit(-1);
				}
				read_vector_from_stream(is, notes);
				if(id != magic_id_v1) {
					cerr << "Error: invalid id (filename is a kbGraph?)" << endl;
					exit(-1);
				}
			} else {
				// Normal case
				read_set_from_stream(is, relsSource);
				read_vector_from_stream(is, rtypes);

				read_map_from_stream(is, synsetMap);
				read_map_from_stream(is, wordMap);

				read_atom_from_stream(is, id);
				if(id != magic_id) {
					cerr << "Error: invalid id after reading maps" << endl;
					exit(-1);
				}

				read_atom_from_stream(is, vertex_n);
				for(i=0; i<vertex_n; ++i) {
					read_vertex_from_stream(is, g);
				}

				read_atom_from_stream(is, id);
				if(id != magic_id) {
					cerr << "Error: invalid id after reading vertices" << endl;
					exit(-1);
				}

				read_atom_from_stream(is, edge_n);
				for(i=0; i<edge_n; ++i) {
					read_edge_from_stream(is, g);
				}

				read_atom_from_stream(is, id);
				if(id != magic_id) {
					cerr << "Error: invalid id after reading edges" << endl;
					exit(-1);
				}
				read_vector_from_stream(is, notes);
				if(id != magic_id) {
					cerr << "Error: invalid id (filename is a kbGraph?)" << endl;
					exit(-1);
				}
			}
		} catch (...) {
			throw runtime_error("Error when reading serialized graph (same platform used to compile the KB?)\n");
		}

		map<string, Kb16_vertex_t>::iterator m_it(wordMap.begin());
		map<string, Kb16_vertex_t>::iterator m_end(wordMap.end());
		for(; m_it != m_end; ++m_it) {
			put(vertex_flags, g, m_it->second,
				get(vertex_flags, g, m_it->second) || Kb16::is_word);
		}
	}

	// write

	//
	// Auxiliary functions for removing isolated vertices
	//


	static size_t vdelta_isolated = numeric_limits<size_t>::max();

	static size_t get_vdeltas(const Kb16Graph & g,
							  vector<size_t> & vdeltas) {

		size_t d = 0;

		graph_traits<Kb16Graph>::vertex_iterator vit, vend;
		tie(vit, vend) = vertices(g);
		for(;vit != vend; ++vit) {
			if (out_degree(*vit, g) + in_degree(*vit, g) == 0) {
				// isolated vertex
				vdeltas[*vit] = vdelta_isolated;
				++d;
			} else {
				vdeltas[*vit] = d;
			}
		}
		return d;
	}

	static void map_update(const vector<size_t> & vdelta,
						   map<string, Kb16_vertex_t> & theMap) {

		map<string, Kb16_vertex_t>::iterator it = theMap.begin();
		map<string, Kb16_vertex_t>::iterator end = theMap.end();

		while(it != end) {
			if (vdelta[it->second] == vdelta_isolated) {
				// erase element
				theMap.erase(it++);
			} else {
				// update vertex id
				it->second -= vdelta[it->second];
				++it;
			}
		}
	}


	// write functions

	ofstream & write_vertex_to_stream(ofstream & o,
									  const Kb16Graph & g,
									  const vector<size_t> & vdelta,
									  const Kb16_vertex_t & v) {
		string name;

		if (vdelta[v] != vdelta_isolated) {
			write_atom_to_stream(o, get(vertex_name, g, v));
			write_atom_to_stream(o, get(vertex_gloss, g, v));
		}
		return o;
	}


	ofstream & write_edge_to_stream(ofstream & o,
									const Kb16Graph & g,
									const vector<size_t> & vdelta,
									const Kb16_edge_t & e) {

		size_t uIdx = get(vertex_index, g, source(e,g));
		uIdx -= vdelta[uIdx];
		size_t vIdx = get(vertex_index, g, target(e,g));
		vIdx -= vdelta[vIdx];

		float w = get(edge_weight, g, e);
		boost::uint32_t rtype = get(edge_rtype, g, e);

		o.write(reinterpret_cast<const char *>(&uIdx), sizeof(uIdx));
		o.write(reinterpret_cast<const char *>(&vIdx), sizeof(vIdx));
		o.write(reinterpret_cast<const char *>(&w), sizeof(w));
		o.write(reinterpret_cast<const char *>(&rtype), sizeof(rtype));
		return o;
	}

	ofstream & Kb16::write_to_stream(ofstream & o) {

		// First remove isolated vertices and
		// - get delta vector
		// - remove from map

		// - get deltas
		vector<size_t> vdelta(num_vertices(g), 0);
		size_t visol_size = get_vdeltas(g, vdelta);

		// - update the maps

		if (visol_size) {
			map_update(vdelta, synsetMap);
			map_update(vdelta, wordMap);
		}

		// Write maps

		write_atom_to_stream(o, magic_id);

		write_vector_to_stream(o, relsSource);
		write_vector_to_stream(o, rtypes);

		write_map_to_stream(o, synsetMap);
		write_map_to_stream(o, wordMap);
		//write_map_to_stream(o, sourceMap);

		write_atom_to_stream(o, magic_id);

		size_t vertex_n = num_vertices(g) - visol_size;

		write_atom_to_stream(o, vertex_n);
		graph_traits<Kb16Graph>::vertex_iterator v_it, v_end;

		tie(v_it, v_end) = vertices(g);
		for(; v_it != v_end; ++v_it) {
			write_vertex_to_stream(o, g, vdelta, *v_it);
		}

		write_atom_to_stream(o, magic_id);

		size_t edge_n = num_edges(g);

		write_atom_to_stream(o, edge_n);
		graph_traits<Kb16Graph>::edge_iterator e_it, e_end;

		tie(e_it, e_end) = edges(g);
		for(; e_it != e_end; ++e_it) {
			write_edge_to_stream(o, g, vdelta, *e_it);
		}
		write_atom_to_stream(o, magic_id);
		if(notes.size()) write_vector_to_stream(o, notes);
		return o;
	}

	void Kb16::write_to_binfile (const string & fName) {

		ofstream fo(fName.c_str(),  ofstream::binary|ofstream::out);
		if (!fo) {
			cerr << "Error: can't create" << fName << endl;
			exit(-1);
		}
		write_to_stream(fo);
	}

	// text write

	ofstream & write_to_textstream(const Kb16Graph & g, ofstream & o) {

		graph_traits<Kb16Graph>::edge_iterator e_it, e_end;

		tie(e_it, e_end) = edges(g);
		for(; e_it != e_end; ++e_it) {
			o << "u:" << get(vertex_name, g, source(*e_it, g)) << " ";
			o << "v:" << get(vertex_name, g, target(*e_it, g)) << " d:1\n";
		}
		return o;
	}

}
