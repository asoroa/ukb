// -*-C++-*-

#ifndef PRANK_H
#define PRANK_H

/////////////////////////////////////////////////////////////////////
// pageRank
//


// PageRank iteration
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
// Note: 
//    * as G is an undirected graph, there is no dangling node


template<typename G, typename ppvMap_t, typename wMap_t, typename map1_t, typename map2_t>
void update_pRank(G & g,
		  std::vector<typename graph_traits<G>::vertex_descriptor> & V,
		  float dfactor,
		  ppvMap_t ppv_V,
		  const std::vector<float> & out_coef, // 1/N
		  wMap_t & wmap,
		  const map1_t rank_map1,
		  map2_t rank_map2) {
  
  typedef typename graph_traits<G>::vertex_descriptor vertex_descriptor;
  
  typename std::vector<vertex_descriptor>::iterator v = V.begin();
  typename std::vector<vertex_descriptor>::iterator end = V.end();
  for (; v != end; ++v) {
    float rank=0.0;
    typename graph_traits<G>::in_edge_iterator e, e_end;
    tie(e, e_end) = in_edges(*v, g);
    for(; e != e_end; ++e) {
      vertex_descriptor u = source(*e, g);
      rank += rank_map1[u] * wmap[*e] * out_coef[u];
    }
    rank_map2[*v] = dfactor*rank + (1-dfactor)*ppv_V[*v];
  }
}

template<typename G, typename ppvMap_t, typename wMap_t, typename map1_t, typename map2_t>
void pageRank_iterate(G & g,
		      std::vector<typename graph_traits<G>::vertex_descriptor> & V,
		      ppvMap_t ppv_V,
		      const std::vector<float> & out_coef,
		      wMap_t & wmap,
		      map1_t rank_map1,
		      map2_t rank_map2,
		      int iterations) {
  int erpin_kop = V.size();
  float damping = 0.85;
  // Initialize rank_map1 appropriately
  {
    const float init_value = 1.0/(float) erpin_kop;
    typename graph_traits<G>::vertex_iterator v, end;
    for (tie(v, end) = vertices(g); v != end; ++v) {
      rank_map1[*v]=init_value;
    }
  }

  // Continue iterating until the termination condition is met

  bool to_map_2 = true;
  while(iterations--) {
    // Update to the appropriate rank map
    if (to_map_2)
      update_pRank(g, V, damping, ppv_V, out_coef, wmap, rank_map1, rank_map2);
    else
      update_pRank(g, V, damping, ppv_V, out_coef, wmap, rank_map2, rank_map1);
    // The next iteration will reverse the update mapping
    to_map_2 = !to_map_2;
  }

  // We stopped after writing the latest results to rank_map2, so copy the
  // results back to rank_map1 for the caller
  if (!to_map_2) {
    typename graph_traits<G>::vertex_iterator v, end;
    for (tie(v, end) = vertices(g); v != end; ++v) {
      rank_map1[*v] = rank_map2[*v];
    }
  }
}

template<typename G, typename wmap_t>
void init_out_coefs(const G & g,
		    const std::vector<typename graph_traits<G>::vertex_descriptor> & V,
		    std::vector<float> & W,
		    wmap_t wmap) {
  typedef typename graph_traits<G>::vertex_descriptor vertex_t;
  typename std::vector<vertex_t>::const_iterator v = V.begin();
  typename std::vector<vertex_t>::const_iterator end = V.end();
  for (; v != end; ++v) {
    typename graph_traits<G>::out_edge_iterator e, e_end;
    tie(e, e_end) = out_edges(*v, g);
    float total_w = 0.0;
    for(; e != e_end; ++e) {
      total_w += wmap[*e];
    }
    W[*v] = 1.0f / total_w;
  }
}

#endif
