#include "common.h"
#include "wdict.h"
#include "globalVars.h"
#include "configFile.h"
#include "kbGraph.h"
#include <string>
#include <iostream>
#include <fstream>

// Boost libraries

#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>


// timer

#include <boost/timer.hpp>


using namespace std;
using namespace boost;
using namespace ukb;

static bool opt_variants = false;

void query (const string & str) {

	Kb & kb = Kb::instance();
	KbGraph & g = kb.graph();

	bool aux;
	Kb_vertex_t u;

	tie(u, aux) = kb.get_vertex_by_name(str);
	if (aux) {
		cout << kb.get_vertex_name(u);
		cout << "\n";
		graph_traits<Kb::boost_graph_t>::out_edge_iterator it , end;
		tie(it, end) = out_edges(u, g);
		for(;it != end; ++it) {
			cout << "  ";
			cout << kb.get_vertex_name(kb.edge_target(*it));
			cout << ":" << g[*it].weight << "\n";
		}
	}
}

bool parse_csv(const string & str,
			   vector<string> & V) {

	typedef tokenizer<char_separator<char> > tokenizer_t;

	char_separator<char> sep(" \t,");
	tokenizer_t tok(str.begin(), str.end(), sep);
	copy(tok.begin(), tok.end(), back_inserter(V));
	return (V.size() > 0);
}

void set_source_rels(const string & str,
					 set<string> & S) {

	vector<string> V;
	set<string>().swap(S); // empty set

	parse_csv(str, V);

	copy(V.begin(), V.end(),
		 inserter(S, S.end())); // Remove duplicates

}

void print_iquery_v(KbGraph & g, Kb_vertex_t u, float w = 0, int sp = 0) {

	string hw(g[u].name);

	for (int i = 0; i < sp; ++i)
		cout << "  ";
	cout << hw;
	if (w)
		cout << ":" << w;
	if (opt_variants)
		cout << " " << WDict::instance().variant(hw);
	cout << "\n";
}

void show_bfs_path(string & str) {

	Kb & kb = Kb::instance();
	KbGraph & g = kb.graph();
	bool aux;
	Kb_vertex_t u, v;
	vector<Kb_vertex_t> p;
	vector<string> path_str;

	path_str = split(str, " ");
	if(path_str.size() < 2) {
		cout << "\n"<< " p needs two concepts!";
		return;
	}
	tie(u, aux) = kb.get_vertex_by_name(path_str[0]);
	if(!aux) {
		cout << path_str[0] << " not present!";
		return;
	}

	tie(v, aux) = kb.get_vertex_by_name(path_str[1]);
	if(!aux) {
		cout << path_str[1] << " not present!";
		return;
	}
	kb.bfs(u, p);
	Kb_vertex_t w = v;
	while(w != u) {
		print_iquery_v(g, w, 0, 1);
		w = p[w];
	}
	print_iquery_v(g, u, 0, 1);
}

void show_neighbours(string & str) {

	Kb & kb = Kb::instance();
	KbGraph & g = kb.graph();

	bool inv = false;
	bool aux;
	Kb_vertex_t u;

	if (str.substr(0, 2) == "i:") {
		inv=true;
		str = str.substr(2);
	}
	tie(u, aux) = kb.get_vertex_by_name(str);
	if (aux) {
		print_iquery_v(g, u);
		if (!inv) {
			graph_traits<Kb::boost_graph_t>::out_edge_iterator it , end;
			tie(it, end) = out_edges(u, g);
			for(;it != end; ++it) {
				print_iquery_v(g, target(*it, g), g[*it].weight, 2);
			}
		} else {
			graph_traits<Kb::boost_graph_t>::in_edge_iterator iit , iend;
			tie(iit, iend) = in_edges(u, g);
			for(;iit != iend; ++iit) {
				print_iquery_v(g, source(*iit, g), g[*iit].weight, 2);
			}
		}
	} else {
		cout << "\n"<< str << " not present!\n";
	}
}

void show_dict_entries(const string & str) {

	WDict_entries entries = WDict::instance().get_entries(str);

	if (!entries.size()) {
		cout << str << " not in dictionary.\n";
		return;
	}
	cout << str << " ->";
	for(size_t i = 0; i < entries.size(); i++) {
		cout << " " << entries.get_entry_str(i) << ":" << entries.get_freq(i);
	}
	cout << "\n";
}

void iquery() {

	while	(cin) {
		string str;
		size_t l;
		cout << "\nQ:";
		read_line_noblank(cin, str, l);;
		if (str == "q") break;
		if (str.substr(0, 2) == "p:") {
			str = str.substr(2);
			show_bfs_path(str);
			continue;
		}
		if (str.substr(0, 2) == "d:") {
			str = str.substr(2);
			show_dict_entries(str);
			continue;
		}
		show_neighbours(str);
	}
}

void subg(const string & initV, size_t N) {

	std::vector<std::string> V;
	std::vector<std::vector<std::string> > E;
	Kb::instance().get_subgraph(initV, V, E, N);
	cout << "graph UMLS_subg {\n";
	for (size_t i = 0, im = V.size(); i != im; ++i) {
		for(size_t j = 0, jm = E[i].size(); j != jm; ++j) {
			cout << V[i] << " -- " << E[i][j] << ";\n";
		}
	}
	cout << "}\n";
}

void sPath(const string & sPathV) {

	string source;
	vector<string> aux,targets;
	std::vector<std::vector<std::string> > paths;
	aux = split(sPathV, "#");
	if (aux.size() < 2) {
		cerr << "Spath error: you must at least specify two nodes\n.";
		exit(-1);
	}
	source = aux[0];
	set<string> S(aux.begin() + 1, aux.end());
	copy(S.begin(), S.end(), back_inserter(targets));
	Kb::instance().get_shortest_paths(source, targets, paths);
	for(std::vector<std::vector<std::string> >::iterator it = paths.begin(), end = paths.end();
		it != end; ++it) {
		writeV(cout, *it);
		cout << "\n";
	}
}

int main(int argc, char *argv[]) {

	srand(3);

	timer load;

	bool opt_dbfile = false;
	bool opt_info = false;
	bool opt_Info = false;
	bool opt_query = false;
	bool opt_iquery = false;
	bool opt_dump = false;

	// subgraph options
	string subg_init;
	size_t subgN = 100;

	string fullname_out("kb_wnet.bin");
	string kb_file;
	string query_vertex;
	string sPathV;

	glVars::kb::v1_kb = false; // Use v2 format
	glVars::kb::filter_src = false; // by default, don't filter relations by src

	set<string> src_allowed;

	string cmdline("cmd: !! -v ");
	cmdline += glVars::ukb_version;
	for (int i=0; i < argc; ++i) {
		cmdline += " ";
		cmdline += argv[i];
	}

	const char desc_header[] = "compile_kb: create a serialized image of the KB\n"
		"Usage:\n"
		"compile_kb -o output.bin [-f \"src1, src2\"] kb_file.txt kb_file.txt ... -> Create a KB image reading relations textfiles.\n"
		"compile_kb -o dict.bin -D dict_textfile --serialize_dict kb_file.bin -> Create a dictionary image reading text dictionary.\n"
		"compile_kb -i kb_file.bin -> Get info of a previously compiled KB.\n"
		"compile_kb -q concept-id kb_file.bin -> Query a node on a previously compiled KB.\n"
		"Options:";

	using namespace boost::program_options;

	options_description po_desc("General options");
	po_desc.add_options()
		("help,h", "This help page.")
		("version", "Show version.")
		("verbose,v", "Be verbose.")
		("serialize_dict", "Serialize dictionary.")
		;

	options_description po_desc_create("Options for creating binary graphs");
	po_desc_create.add_options()
		("output,o", value<string>(), "Output file name.")
		("filter_src,f", value<string>(), "Filter relations according to their sources.")
		("undirected,U", "Force undirected graph.")
		("rtypes,r", "Keep relation types on edges.")
		("minput", "Do not die when dealing with malformed input.")
		("nopos", "Don't filter words by Part of Speech when reading dict.")
		("note", value<string>(), "Add a comment to the graph.")
		;

	options_description po_desc_query("Options for querying over binary graphs");
	po_desc_query.add_options()
		("info,i", "Give info about some Kb binfile.")
		("Info,I", "Give more info about Kb binfile. This option can be computationally expensive.")
		("dump", "Dump a serialized graph. Warning: very verbose!.")
		("query,q", value<string>(), "Given a vertex name, display its relations.")
		("iquery,Q", "Interactively query graph.")
		("subG,S", value<string>(), "Get a subgraph starting at this vertex. See subG_depth.")
		("subG_N", value<size_t>(), "Max. number of nodes in subgraph (see --subG). Default is 100.")
		("sPaths", value<string>(), "Get shortest paths. Value is a vector of nodes, separated by character '#'. First node is source.")
		("dict_file,D", value<string>(), "Dictionary text file. Use only when querying (--quey or --iquery) or when creating serialized dict (--serialize_dict).")
		;

	options_description po_hidden("Hidden");
	po_hidden.add_options()
		("input-file",value<string>(), "Input files.")
		;
	options_description po_visible(desc_header);
	po_visible.add(po_desc).add(po_desc_create).add(po_desc_query);

	options_description po_desc_all("All options");
	po_desc_all.add(po_visible).add(po_hidden);

	positional_options_description po_optdesc;
	po_optdesc.add("input-file", -1);

	variables_map vm;

	try {
		store(command_line_parser(argc, argv).
			  options(po_desc_all).
			  positional(po_optdesc).
			  run(), vm);
		notify(vm);

		// If asked for help, don't do anything more

		if (vm.count("help")) {
			cout << po_visible << endl;
			exit(0);
		}

		if (vm.count("version")) {
			cout << glVars::ukb_version << endl;
			exit(0);
		}

		if (vm.count("verbose")) {
			glVars::verbose = 1;
		}

		if (vm.count("nopos")) {
			glVars::input::filter_pos = false;
		}

		if (vm.count("serialize_dict")) {
			opt_dbfile = true;
		}

		if (vm.count("info")) {
			opt_info = true;
		}

		if (vm.count("Info")) {
			opt_Info = true;
		}

		if (vm.count("iquery")) {
			opt_iquery = true;
		}

		if (vm.count("subG")) {
			subg_init = vm["subG"].as<string>();
		}

		if (vm.count("subG_depth")) {
			subgN = vm["subG_depth"].as<size_t>();
		}

		if (vm.count("sPaths")) {
			sPathV = vm["sPaths"].as<string>();
		}

		if (vm.count("dict_file")) {
			glVars::dict::text_fname = vm["dict_file"].as<string>();
			opt_variants = true;
		}

		if (vm.count("note")) {
		}

		if (vm.count("query")) {
			opt_query = true;
			query_vertex = vm["query"].as<string>();
		}

		if (vm.count("filter_src")) {
			glVars::kb::filter_src = true;
			set_source_rels(vm["filter_src"].as<string>(), src_allowed);
		}

		if (vm.count("dump")) {
			opt_dump = true;
		}

		if (vm.count("undirected")) {
			glVars::kb::keep_directed = false;
		}

		if (vm.count("rtypes")) {
			glVars::kb::keep_reltypes = true;
		}

		if (vm.count("minput")) {
			glVars::input::swallow = true;
		}

		if (vm.count("input-file")) {
			kb_file = vm["input-file"].as<string>();
		}

		if (vm.count("output")) {
			fullname_out = vm["output"].as<string>();
		}
	}
	catch(std::exception& e) {
		cerr << e.what() << "\n";
		exit(-1);
	}

	if (!kb_file.size()) {
		cerr << po_visible << "\n";
		cerr << "Error: no input files\n" << endl;
		exit(1);
	}

	if (opt_info) {
		Kb::create_from_binfile(kb_file);
		Kb::instance().display_info(cout);
		return 0;
	}

	if (opt_Info) {
		Kb::create_from_binfile(kb_file);

		int id_m, id_M;
		int od_m, od_M;
		int comp;

		Kb::instance().display_info(cout);

		tie(id_m, id_M) = Kb::instance().indeg_maxmin();
		tie(od_m, od_M) = Kb::instance().outdeg_maxmin();
		comp = Kb::instance().components();

		cout << "In degree (max, min):  (" << id_m << ", " << id_M << ")\n";
		cout << "Out degree (max, min): (" << od_m << ", " << od_M << ")\n";
		cout << "Number of (strong) components: " << comp << endl;

		return 0;
	}

	if (opt_iquery) {
		Kb::create_from_binfile(kb_file);
		iquery();
		return 0;
	}

	if (opt_dump) {
		Kb::create_from_binfile(kb_file);
		Kb::instance().dump_graph(cout);
		return 0;
	}

	if (opt_query) {
		Kb::create_from_binfile(kb_file);
		query(query_vertex);
		return 0;
	}

	if(subg_init.size()) {
		Kb::create_from_binfile(kb_file);
		subg(subg_init, subgN);
		return 0;
	}

	if(sPathV.size()) {
		Kb::create_from_binfile(kb_file);
		sPath(sPathV);
		return 0;
	}

	if (glVars::verbose) {
		show_global_variables(cerr);
	}

	if (glVars::verbose)
		cerr << "Reading relations"<< endl;

	if (opt_dbfile) {
		// Serialize dict
		if (glVars::dict::text_fname.size() == 0) {
			cerr << "--serialize_dict error: -D option missing.\n";
		}
		if (kb_file.size() == 0) {
			cerr << "--serialize_dict error: graph missing.\n";
			exit(-1);
		}
		Kb::create_from_binfile(kb_file);
		WDict::instance().write_wdict_binfile(fullname_out);
		exit(0);
	}

	try {
		// If first input file is "-", open std::cin
		if (kb_file == "-") {
			cmdline += " <STDIN>";
			Kb::create_from_txt(std::cin, src_allowed );
		} else {
			Kb::create_from_txt(kb_file, src_allowed );
		}
	}  catch(std::exception& e) {
		cerr << e.what() << "\n";
		exit(-1);
	}


	if (glVars::verbose)
		cerr << "Writing binary file: "<< fullname_out<< endl;
	Kb::instance().add_comment(cmdline);
	Kb::instance().write_to_binfile(fullname_out);
	if (glVars::verbose)
		cerr << "Wrote " << num_vertices(Kb::instance().graph()) << " vertices and " << num_edges(Kb::instance().graph()) << " edges" << endl;

	return 0;
}

/*
 * Local Variables:
 * mode: c++
 * End:
 */
