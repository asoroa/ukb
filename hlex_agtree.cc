#include "hlex_agtree.h"

// IO

namespace ukb {

  using namespace std;
  using namespace boost;


  long AgZuhaitza::magic_id =  0x231270L;
  long AgZuhaitza::magic_idG = 0x010999L;

  AgZuhaitza::vertex_descriptor read_vertex_from_stream(ifstream & is,
														AgTree & g) {

	size_t len;
	string izena;
	size_t hDist;
	size_t hI;
	size_t vO;

	is.read(reinterpret_cast<char *> (&len), sizeof(len));
	if (len) {
	  char * iz=new char [len+1]; // tenporala
	  is.read(iz, len);
	  iz[len] = 0;
	  izena = string(iz);
	  delete[] iz;
	}
	is.read(reinterpret_cast<char *> (&hDist), sizeof(hDist));
	is.read(reinterpret_cast<char *> (&hI), sizeof(hI));
	is.read(reinterpret_cast<char *> (&vO), sizeof(vO));
	AgZuhaitza::vertex_descriptor u = add_vertex(g);
	put(vertex_name, g, u, izena);
	put(vertex_hubDist, g, u, hDist);
	put(vertex_hubN, g, u, hI);
	put(vertex_vOrig, g, u, vO);
	return u;
  }

  AgZuhaitza::edge_descriptor read_edge_from_stream(ifstream & is,
													AgTree & g) {
	size_t sIdx;
	size_t tIdx;
	is.read(reinterpret_cast<char *> (&sIdx), sizeof(sIdx));
	is.read(reinterpret_cast<char *> (&tIdx), sizeof(tIdx));
	//is.read(reinterpret_cast<char *> (&pisua), sizeof(pisua));
	//is.read(reinterpret_cast<char *> (&maizt), sizeof(maizt));
	//vertex_descriptor s(sIdx);
	//vertex_descriptor t(tIdx);
	AgZuhaitza::edge_descriptor e;
	bool sartu_da;
	//tie(e,sartu_da) = add_edge(s,t,g);
	tie(e,sartu_da) = add_edge(sIdx,tIdx,g);
	//put(edge_weight, g, e, pisua);
	//put(edge_freq, g, e, maizt);
	return e;
  }

  ifstream & AgZuhaitza::read_from_stream(ifstream & is)  {

	size_t erpin_kop;
	size_t ertz_kop;
	size_t len;
	long id;

	is.read(reinterpret_cast<char *> (&id), sizeof(id));
	if (id != magic_id) {
	  cerr << "Errorea zuhaitza irakurtzerakoan (fitxategia zuhaitza da ?)" << endl;
	  exit(-1);
	}
	is.read(reinterpret_cast<char *> (&len), sizeof(len));
	if (len) {
	  char * iz=new char [len+1]; // tenporala
	  is.read(iz, len);
	  iz[len] = 0;
	  targetWord = string(iz);
	  delete[] iz;
	}
	is.read(reinterpret_cast<char *> (&erpin_kop), sizeof(erpin_kop));
	is.read(reinterpret_cast<char *> (&ertz_kop), sizeof(ertz_kop));
	is.read(reinterpret_cast<char *> (&targetVertex), sizeof(vertex_descriptor));

	//   if (erpin_kop) {
	//     if (is) {
	//       targetVertex = read_vertex_from_stream(is, g);
	//       --erpin_kop;
	//     } else {
	//       cerr << "Errorea grafoa irakurtzerakoan (hasierako erpina)" << endl;
	//       exit(-1);
	//     }
	//   }
	if (is) {
	  for(size_t i=0;i != erpin_kop; ++i)
		read_vertex_from_stream(is, g);
	} else {
	  cerr << "Errorea zuhaitza irakurtzerakoan (tarteko erpin bat)" << endl;
	  exit(-1);
	}
	if (is) {
	  for(size_t i=0;i != ertz_kop; ++i)
		read_edge_from_stream(is, g);
	} else {
	  cerr << "Errorea zuhaitza irakurtzerakoan (ertzak)" << endl;
	  exit(-1);
	}
	// hitzMapa hash-a sortu
	graph_traits<AgTree>::vertex_iterator vIt, vItEnd;
	tie(vIt, vItEnd) = vertices(g);
	for(;vIt != vItEnd; ++vIt)
	  if (*vIt == targetVertex) {
		targetWord = get(vertex_name, g, *vIt);
	  } else
		hitzMapa.insert(make_pair<string, vertex_descriptor> (get(vertex_name, g, *vIt),
															  *vIt));
	//gDisplay(g);
	// root zenbakia lortu
	hubsN=0;
	graph_traits<AgTree>::out_edge_iterator eIt, eItEnd;
	tie(eIt, eItEnd) = out_edges(targetVertex, g);
	for(;eIt != eItEnd; ++eIt) {
	  hubs.push_back(get(vertex_name, g, target(*eIt, g)));
	  ++hubsN;
	}
	return is;
  }

  void AgZuhaitza::read_from_file(const string & fizena) {

	//  ofstream fo(fozena.c_str(),  ofstream::binary|ofstream::out);

	ifstream fi(fizena.c_str(), ifstream::binary|ifstream::in);

	if (!fi) {
	  cerr << "Errorea:" << fizena << " fitxategia ezin dut ireki !" << endl;
	  exit(-1);
	}
	read_from_stream(fi);
  }
}
