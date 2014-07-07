#include "kbGraph_common.h"

namespace ukb {

	using namespace std;
	using namespace boost;

	////////////////////////////////////////////////////////////////////////////////
	// etype_t

	etype_t & etype_t::operator=(const etype_t & o) {
		if (&o != this) {
			m_strtypes = o.m_strtypes;
		}
		return *this;
	}

	void etype_t::swap(etype_t & o) {
		m_strtypes.swap(o.m_strtypes);
	}

	size_t etype_t::size() const { return m_strtypes.size(); }

	void etype_t::add_type(const std::string & tstr, value_type & val) {
		int pos = strpos(tstr);
		if (pos == -1)
			pos = stradd(tstr);
		val |= (value_type(1) << pos);
	}

	bool etype_t::has_type(const std::string & tstr, value_type val) const {
		int pos = strpos(tstr);
		if (pos == -1)
			return false;
		return (val & (value_type(1) << pos));
	}

	std::vector<std::string> etype_t::tvector(value_type val) const {
		std::vector<std::string> res;
		if (m_strtypes.size() == 0) return res;

		std::vector<std::string>::size_type idx = 0;
		value_type i(1);
		while(idx < 32) {
			if (val & i) {
				res.push_back(m_strtypes[idx]);
			}
			idx++;
			i <<= 1;
		}
		return res;
	}

	void etype_t::read_from_stream (std::istream & o) {
		read_vector_from_stream(o, m_strtypes);
	}
	std::ostream & etype_t::write_to_stream(std::ostream & o) const {
		write_vector_to_stream(o, m_strtypes);
		return o;
	}

	int etype_t::strpos(const std::string & tstr) const {
		int res = 1;
		std::vector<std::string>::const_iterator it = m_strtypes.begin();
		std::vector<std::string>::const_iterator end = m_strtypes.end();
		for(;it != end; ++it) {
			if (*it == tstr) break;
			res++;
		}
		if(it == end) return -1;
		return res;
	}

	size_t etype_t::stradd(const std::string & tstr) {
		m_strtypes.push_back(tstr);
		size_t pos = m_strtypes.size();
		if (pos > 32)
			throw std::runtime_error("etype_t:::add_type error: too many relation types !");
		return pos;
	}

	// precsr_t


	size_t precsr_t::insert_vertex(const std::string & ustr) {

		bool insertedP;
		vertex_map_t::iterator vit;

		boost::tie(vit, insertedP) = m_vMap.insert(std::make_pair(ustr, m_vsize));
		if(insertedP) {
			vProp.push_back(vertex_prop_t(ustr));
			++m_vsize;
		}
		return vit->second;
	}

	size_t precsr_t::insert_edge(const std::string & ustr,
								 const std::string & vstr,
								 float w,
								 etype_t::value_type etype) {

		size_t u = insert_vertex(ustr);
		size_t v = insert_vertex(vstr);

		edge_map_t::iterator eit;
		edge_map_t::key_type k = std::make_pair(u, v);
		bool insertedP;

		boost::tie(eit, insertedP) = m_eMap.insert(std::make_pair(k, m_esize));
		if(insertedP) {
			E.push_back(k);
			eProp.push_back(edge_prop_t());
			++m_esize;
		}
		eProp[eit->second].etype |= etype;
		eProp[eit->second].weight = w;

		return eit->second;
	}

	size_t precsr_t::insert_edge(const std::string & ustr,
								 const std::string & vstr,
								 float w,
								 const std::string & rtype) {
		size_t eidx = insert_edge(ustr, vstr, w, etype_t::value_type(0));
		// add edge type
		if (rtype.size())
			m_rtypes.add_type(rtype,eProp[eidx].etype);
		return eidx;
	}
}
