// -*-C++-*-

#ifndef PRANK_H
#define PRANK_H

#include <boost/graph/graph_concepts.hpp>
#include<boost/tuple/tuple.hpp> // for "tie"
#include <iosfwd>

/////////////////////////////////////////////////////////////////////
// pageRank
//
// it solves pageRank with the so called power method see
// @article{langville04,
// title = {{Deeper Inside PageRank}},
// author = {A.N. Langville and C.D. Meyer},
// journal = {Internet Mathematics},
// number = {3},
// pages = {335--380},
// volume = {1},
// year = {2004},
// }
//
// Note: it correctly handles dangling nodes

namespace ukb {

  namespace prank {

	//
	// a constant property map that always returns the same value
	//


	template<typename K, typename V>
	class constant_property_map
	  : public boost::put_get_helper<V, constant_property_map<K, V> > {
	public:
	  typedef K key_type;
	  typedef V value_type;
	  typedef V reference;
	  typedef boost::readable_property_map_tag category;

	  constant_property_map(V value) : store(value) {}

	  inline value_type operator[](const key_type& v) const { return store; }
	private:
	  V store;
	};

	////////////////////////////////
	//
	// Init out_coefs so that out_coefs[v] has the sum of weights of
	// out-edges.
	//

	template<typename G, typename coefmap_t, typename wmap_t>
	void init_out_coefs(const G & g,
						const std::vector<typename graph_traits<G>::vertex_descriptor> & V,
						coefmap_t W,
						wmap_t wmap) {
	  typedef typename graph_traits<G>::vertex_descriptor vertex_t;
	  typename std::vector<vertex_t>::const_iterator v = V.begin();
	  typename std::vector<vertex_t>::const_iterator end = V.end();
	  for (; v != end; ++v) {
		typename graph_traits<G>::out_edge_iterator e, e_end;
		tie(e, e_end) = out_edges(*v, g);
		float total_w = 0.0;
		for(; e != e_end; ++e) {
		  if(target(*e, g) == *v) continue; // Don't count self loops
		  total_w += wmap[*e];
		}
		W[*v] = 1.0f / total_w;
	  }
	}

	//
	// Apply one step of pageRank algorithm
	//

	template<typename G, typename ppvMap_t, typename wMap_t, typename map1_t, typename map2_t>
	float update_pRank(G & g,
					   std::vector<typename graph_traits<G>::vertex_descriptor> & V,
					   float damping,
					   ppvMap_t ppv_V,
					   const std::vector<float> & out_coef,
					   wMap_t & wmap,
					   const map1_t rank_map1,
					   map2_t rank_map2) {

	  typedef typename graph_traits<G>::vertex_descriptor vertex_descriptor;

	  typename std::vector<vertex_descriptor>::iterator v = V.begin();
	  typename std::vector<vertex_descriptor>::iterator end = V.end();
	  float norm = 0.0;
	  for (; v != end; ++v) {
		double rank = 0.0;
		typename graph_traits<G>::in_edge_iterator e, e_end;
		tie(e, e_end) = in_edges(*v, g);
		for(; e != e_end; ++e) {
		  vertex_descriptor u = source(*e, g);
		  if (*v == u) continue; // No self-loops
		  rank += rank_map1[u] * wmap[*e] * static_cast<double>(out_coef[u]);
		}
		double dangling_factor = 0.0;
		if (0.0 == out_coef[*v]) {
		  // dangling link
		  dangling_factor = static_cast<double>(damping)*rank_map1[*v];
		}
		rank_map2[*v] = damping*rank + (dangling_factor + static_cast<double>(1.0 - damping) )*ppv_V[*v];
		norm += fabs(rank_map2[*v] - rank_map1[*v]);
	  }
	  return norm;
	}

	//
	// Initialize rank and iterate
	//

	template<typename G, typename ppvMap_t, typename wMap_t, typename map1_t, typename map2_t>
	void do_pageRank(G & g,
					 std::vector<typename graph_traits<G>::vertex_descriptor> & V,
					 ppvMap_t ppv_V,
					 wMap_t & wmap,
					 map1_t rank_map1,
					 map2_t rank_map2,
					 int iterations,
					 float threshold,
					 const std::vector<float> & out_coef) {

	  if (iterations == 0 && threshold == 0.0)
		throw std::runtime_error("prank error: iterations and threshold are set to zero!\n");
	  if (!iterations) iterations = std::numeric_limits<int>::max();

	  float damping = 0.85;
	  // Initialize rank_map1 appropriately
	  {
		const double init_value = 1.0f/static_cast<double>(V.size());
		typename graph_traits<G>::vertex_iterator v, end;
		for (tie(v, end) = vertices(g); v != end; ++v) {
		  rank_map1[*v]=init_value;
		}
	  }

	  // Continue iterating until the termination condition is met

	  bool to_map_2 = true;
	  float residual = 0.0;
	  while(iterations--) {
		// Update to the appropriate rank map
		if (to_map_2)
		  residual = update_pRank(g, V, damping, ppv_V, out_coef, wmap, rank_map1, rank_map2);
		else
		  residual = update_pRank(g, V, damping, ppv_V, out_coef, wmap, rank_map2, rank_map1);
		// The next iteration will reverse the update mapping
		to_map_2 = !to_map_2;
		if (residual < threshold) break;
	  }

	  // If we stopped after writing the latest results to rank_map2,
	  // copy the results back to rank_map1 for the caller
	  if (!to_map_2) {
		typename graph_traits<G>::vertex_iterator v, end;
		for (tie(v, end) = vertices(g); v != end; ++v) {
		  rank_map1[*v] = rank_map2[*v];
		}
	  }
	}


	/////////////////////////////////////////////////////////////////
	// PageRank iteration
	//

	//
	// main entry point
	//

	template<typename G, typename ppvMap_t, typename wMap_t, typename map1_t, typename map2_t>
	void pageRank_iterate(G & g,
						  std::vector<typename graph_traits<G>::vertex_descriptor> & V,
						  ppvMap_t ppv_V,
						  wMap_t & wmap,
						  map1_t rank_map1,
						  map2_t rank_map2,
						  int iterations,
						  float threshold) {

	  std::vector<float> out_coef(num_vertices(g), 0.0f);
	  // Initialize out_coef
	  init_out_coefs(g, V, &out_coef[0], wmap);
	  do_pageRank(g, V, ppv_V, wmap, rank_map1, rank_map2, iterations, threshold, out_coef);
	}

	//
	// PageRank algorithm without weights
	//
	// main entry point
	//

	template<typename G, typename ppvMap_t, typename map1_t, typename map2_t>
	void pageRank_iterate_now(G & g,
							  std::vector<typename graph_traits<G>::vertex_descriptor> & V,
							  ppvMap_t ppv_V,
							  map1_t rank_map1,
							  map2_t rank_map2,
							  int iterations,
							  float threshold) {

	  typedef typename graph_traits<G>::edge_descriptor edge_descriptor;
	  constant_property_map <edge_descriptor, float> cte_weight(1); // always return 1
	  pageRank_iterate(g, V, ppv_V, cte_weight, rank_map1, rank_map2, iterations, threshold);
	}

	/////////////////////////////////////////////////////////////////////////

	template<typename G, typename coefmap_t, typename wmap_t>
	void init_degree(const G & g,
					 coefmap_t In,
					 wmap_t wmap) {

	  typename graph_traits<G>::vertex_iterator v, v_end;

	  tie(v, v_end) = vertices(g);
	  for(; v != v_end; ++v) {
		In[*v] = out_degree(*v, g);
	  }
	}
  }
}
#endif
