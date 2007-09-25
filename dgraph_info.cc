#include "disambGraph.h"
#include <string>
#include <iostream>
#include <fstream>

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

using namespace std;
using namespace boost;

int main(int argc, char *argv[]) {

  if(argc < 2) {
    cerr << "Usage: dgraph_info dgraph_file" <<endl;
    return -1;
  }

  DisambGraph dg;
  dg.read_from_binfile(argv[1]);

  DisambGraph::boost_graph_type & g = dg.graph();

  size_t N = num_vertices(g);

  cout << "V: " << N << " E:" << num_edges(g) << endl;

  std::vector<int> component(N);
  int num = connected_components(g, &component[0]); 

  //int ccomp = strong_components(g, make_iterator_property_map(ccompV.begin(), get(vertex_index, g), ccompV[0]));

  cout << "Strong components: " << num << endl;

  return 0;
}

/*
 * Local Variables:
 * mode: c++
 * compile-command: "make -k dgraph_info"
 * End:
 */
