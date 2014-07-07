// -*-C++-*-

#ifndef WDICT_VECTOR_H
#define WDICT_VECTOR_H

#include <algorithm>
#include <boost/cstdint.hpp> // for uint32_t

template<class T>
class wdict_vector {

public:
	typedef  T *      iterator;
	typedef const T * const_iterator;


	wdict_vector() : m_ptr(0), m_size(0) {}

	~wdict_vector() {
		if (m_ptr) {
			try {
				delete [] m_ptr;
			} catch (...) {}
		}
	}

	// Copy ctor
	wdict_vector(const wdict_vector & o) : m_ptr(0), m_size(o.m_size) {
		if (m_size) {
			m_ptr = new T[ m_size ];
			std::copy(o.m_ptr, o.m_ptr + m_size, m_ptr);
		}
	}

	// Assignment operator
	wdict_vector & operator=(const wdict_vector & o) {
		wdict_vector aux(o);
		this->swap(aux);
		return *this;
	}

	// Assign from vector
	template<class Vector>
	void assign(Vector & v) {
		if (m_ptr) {
			delete [] m_ptr;
			m_ptr = 0;
		}
		m_size = v.size();
		if (m_size) {
			m_ptr = new T[ m_size ];
			std::copy(v.begin(), v.end(), m_ptr);
		}
	}

	T & operator[](size_t i) {
		return *(m_ptr + i);
	}

	const T & operator[](size_t i) const {
		return *(m_ptr + i);
	}

	iterator begin() {
		return m_ptr;
	}

	const_iterator begin() const {
		return m_ptr;
	}

	iterator end() {
		return m_ptr + m_size;
	}

	const_iterator end() const {
		return m_ptr + m_size;
	}

	void swap(wdict_vector & o) {
		std::swap(m_size, o.m_size);
		std::swap(m_ptr, o.m_ptr);
	}

	size_t size() const {
		return m_size;
	}

	size_t capacity() const {
		return m_size;
	}

private:

	T *m_ptr;
	boost::uint32_t m_size;

};

#endif
