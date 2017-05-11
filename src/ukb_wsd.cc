#include "wdict.h"
#include "common.h"
#include "globalVars.h"
#include "kbGraph.h"
#include "disambGraph.h"
#include "fileElem.h"
#include <string>
#include <iostream>
#include <fstream>
#include <syslog.h>

#include "ukbServer.h"

// Basename & friends
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

// Program options

#include <boost/program_options.hpp>

using namespace ukb;
using namespace std;
using namespace boost;

/////////////////////////////////////////////////////////////
// Global variables

enum main_method {
	m_ppr,       // ppr, whole graph
	m_ppr_w2w,   // ppr_w2w, whole graph
	m_static,    // static, whole graph
	//	m_degree,    // degree, whole graph
	m_bfs,       // bfs subgraph
	m_dfs
};

enum dgraph_rank {
	r_ppr,
	r_ppr_w2w,
	r_static,
	r_degree
};


dgraph_rank dgraph_rank_method = r_ppr;

main_method opt_dmethod = m_ppr_w2w;
string cmdline;
bool opt_daemon = false;
bool opt_dump_dgraph = false;

// Program options stuff

/* Auxiliary functions for checking input for validity. */

/* Function used to check that 'opt1' and 'opt2' are not specified
   at the same time. */
void conflicting_options(const boost::program_options::variables_map& vm,
						 const char* opt1, const char* opt2)
{
	if (vm.count(opt1) && !vm[opt1].defaulted()
		&& vm.count(opt2) && !vm[opt2].defaulted())
		throw logic_error(string("Conflicting options '")
						  + opt1 + "' and '" + opt2 + "'.");
}

/* Function used to check that of 'for_what' is specified, then
   'required_option' is specified too. */
void option_dependency(const boost::program_options::variables_map& vm,
					   const char* for_what, const char* required_option)
{
	if (vm.count(for_what) && !vm[for_what].defaulted())
		if (vm.count(required_option) == 0 || vm[required_option].defaulted())
			throw logic_error(string("Option '") + for_what
							  + "' requires option '" + required_option + "'.");
}

///////////////////////////////////////

// Main functions

// Disambiguate using disambiguation graph (dgraph) method

static void build_dgraph(CSentence & cs, DisambGraph & dgraph) {

	switch (opt_dmethod) {
	case m_dfs:
		if (glVars::dGraph::stopCosenses)
			build_dgraph_dfs_nocosenses(cs, dgraph);
		else build_dgraph_dfs(cs, dgraph);
		break;
	case m_bfs:
		build_dgraph_bfs(cs, dgraph);
		break;
	default:
		cerr << "build_dgraph: [E] logic error\n";
		exit(1);
	}
}

static bool rank_dgraph(const CSentence & cs,
						DisambGraph & dg,
						vector<float> & ranks) {

	bool ok = false;
	switch(dgraph_rank_method) {
	case r_ppr:
		//ok = (opt_dmethod == m_mention) ? dgraph_mention_ppr(cs, dg, ranks) : dgraph_ppr(cs, dg, ranks);
		ok = dgraph_ppr(cs, dg, ranks);
	break;
	case r_degree:
		ok = dgraph_degree(dg, ranks);
		break;
	case r_static:
		ok = dgraph_static(dg, ranks);
		break;
	case r_ppr_w2w:
		cerr << "rank_dgraph: [E] can't use ppr_w2w\n";
		exit(1);
		break;
	}
	return ok;
}

static void disamb_dgraph(CSentence & cs,
						  DisambGraph & dg,
						  vector<float> & ranks) {
		disamb_csentence_dgraph(cs, dg, ranks);
}


void dgraph_w2w_csent(CSentence & cs) {

	// fall back to static if csentence has only one word
	if (cs.size() == 1) {
		if (glVars::debug::warning)
			cerr << "dis_csent: using static for context " << cs.id() << endl;
		const vector<float> ranks = Kb::instance().static_prank();
		disamb_csentence_kb(cs, ranks);
	} else {
		DisambGraph dgraph;
		build_dgraph(cs, dgraph);
		vector<float> ranks;
		for(CSentence::iterator cw_it = cs.ubegin(), cw_end = cs.uend();
			cw_it != cw_end; ++cw_it) {
			if(!cw_it->is_tgtword()) continue;
			bool ok = dgraph_ppr(cs, dgraph, ranks, cw_it);
			if (!ok && glVars::debug::warning) {
				cerr << "dgraph_w2w_csent: [W] No ranks for sentence " << cs.id() << "\n";
				return;
			}
			disamb_cword_dgraph(cw_it, dgraph, ranks);
		}
	}
}

void dgraph_csent(CSentence & cs) {

	if (dgraph_rank_method == r_ppr_w2w) {
		dgraph_w2w_csent(cs);
		return;
	}

	// fall back to static if csentence has only one word
	if (cs.size() == 1) {
		if (glVars::debug::warning)
			cerr << "dis_csent: using static for context " << cs.id() << endl;
		const vector<float> & ranks = Kb::instance().static_prank();
		disamb_csentence_kb(cs, ranks);
	} else {
		DisambGraph dgraph;
		build_dgraph(cs, dgraph);
		vector<float> ranks;
		bool ok = rank_dgraph(cs, dgraph, ranks);
		if (!ok && glVars::debug::warning) {
			cerr << "dgraph_csent: [W] No ranks for sentence " << cs.id() << "\n";
			return;
		}
		disamb_dgraph(cs, dgraph, ranks);
	}
}

void ppr_csent(CSentence & cs) {

	vector<float> ranks;
	bool ok = calculate_kb_ppr(cs,ranks);
	if (!ok && glVars::debug::warning) {
		std::cerr << "ppr_csent: [W] Error in sentence " << cs.id() << "\n";
		return;
	}
	disamb_csentence_kb(cs, ranks);
}

// w2w approach
// for each target word
//   1. init pv with the synsets of the rest of words
//   2. run Personalized Pagerank
//   3. use rank for disambiguating word

void ppr_w2w_csent(CSentence & cs) {

	vector<float> ranks;
	int success_n = 0;

	vector<CWord>::iterator cw_it = cs.ubegin();
	vector<CWord>::iterator cw_end = cs.uend();
	for(; cw_it != cw_end; ++cw_it) {
		// Target word must be distinguished.
		if(!cw_it->is_tgtword()) continue;
		if (!cw_it->is_monosemous() &&
			calculate_kb_ppr_by_word(cs, cw_it, ranks)) {
			success_n++;
			cw_it->rank_synsets(ranks, glVars::csentence::mult_priors);
		}
		cw_it->disamb_cword();
	}
	if (!success_n && glVars::debug::warning) {
		std::cerr << "ppr_w2w_csent: [W] Error in sentence " << cs.id() << "\n";
		return;
	}
}

void static_csent(CSentence &cs) {

	static vector<float> ranks;
	if (!ranks.size()) {
		ranks = Kb::instance().static_prank();
	}
	disamb_csentence_kb(cs, ranks);
}


void dispatch_run_cs(CSentence & cs) {

	switch(opt_dmethod) {
	case m_bfs:
	case m_dfs:
		dgraph_csent(cs);
		break;
	case m_ppr:
		ppr_csent(cs);
		break;
	case m_ppr_w2w:
		ppr_w2w_csent(cs);
		break;
	case m_static:
		static_csent(cs);
		break;
	};
}

void dispatch_run(istream & is, ostream & os) {

	size_t l_n = 0;
	string cid, ctx;
	while (read_ukb_ctx(is, l_n, cid, ctx)) {
		try {
			CSentence cs(cid, ctx);
			if(ctx.size()) {
				dispatch_run_cs(cs);
				cs.print_csent(os);
			} else {
				if (glVars::debug::warning) {
					cerr << "[W] empty context " << cs.id() + " in line " + lexical_cast<string>(l_n) + "\n";
				}
			}
		} catch (ukb::wdict_error & e) {
			throw e;
		} catch (std::logic_error & e) {
			string msg = "[E] Bad context in line " + lexical_cast<string>(l_n) + "\n" + e.what();
			if (!glVars::input::swallow) throw std::runtime_error(msg);
			if (glVars::debug::warning) {
				cerr << msg << "\n";
			}
		}
	}
}

///////////////////////////////////////////////
// Server/clien functions

// Return FALSE means kill server

#ifdef UKB_SERVER
bool handle_server_read(sSession & session) {
	string ctx_id;
	string ctx;
	try {
		session.receive(ctx);
		if (ctx == "stop") return false;
		session.send(cmdline);
		while(1) {
			if (!session.receive(ctx_id)) break;
			if (!session.receive(ctx)) break;
			CSentence cs(ctx_id, ctx);
			dispatch_run_cs(cs);
			ostringstream oss;
			cs.print_csent(oss);
			string oss_str(oss.str());
			if (!oss_str.length()) oss_str = "#"; // special line if not output
			session.send(oss_str);
		}
	} catch (std::exception& e)	{
		// send error and close the session.
		// Note: the server is still alive for new connections.
		session.send(e.what());
	}
	return true;
}

bool client(istream & is, ostream & os, unsigned int port) {
	// connect to ukb port and send data to it.
	sClient client("localhost", port);
	string server_cmd;
	string go("go");
	if (client.error()) {
		std::cerr << "Error when connecting: " << client.error_str() << std::endl;
		return false;
	}
	string id, ctx, out;
	size_t l_n = 0;
	try {
		client.send(go);
		client.receive(server_cmd);
		os << server_cmd << std::endl;
		while(read_line_noblank(is, id, l_n)) {
			if(!read_line_noblank(is, ctx, l_n)) return false;
			client.send(id);
			client.send(ctx);
			client.receive(out);
			if (out == "#") continue; // empty output for that context
			os << out;
			os.flush();
		}
	} catch (std::exception& e)	{
		std::cerr << e.what() << std::endl;
		return false;
	}
	return true;
}


bool client_stop_server(unsigned int port) {
	// connect to ukb port and tell it to stop
	sClient client("localhost", port);
	string stop("stop");
	if (client.error()) {
		std::cerr << "client_stop_server: [E] Error when connecting: " << client.error_str() << std::endl;
		return false;
	}
	try {
		client.send(stop);
	} catch (std::exception& e)	{
		std::cerr << e.what() << std::endl;
		return false;
	}
	return true;
}
#endif

void load_kb_and_dict(bool from_daemon) {

	if (from_daemon) {
		string aux("Loading KB ");
		aux += glVars::kb::fname;
		syslog(LOG_INFO | LOG_USER, "%s", aux.c_str());
	} else if (glVars::verbose) {
		cout << "Loading KB " + glVars::kb::fname + "\n";
	}
	Kb::create_from_binfile(glVars::kb::fname);
	// Explicitly load dictionary only if:
	// - there is a dictionary name (textual or binary)
	// - from_daemon is set
	if (!from_daemon) return;
	if (!(glVars::dict::text_fname.size() + glVars::dict::bin_fname.size())) return;
	string aux("Loading Dict ");
	aux += glVars::dict::text_fname.size() ? glVars::dict::text_fname : glVars::dict::bin_fname;
	syslog(LOG_INFO | LOG_USER, "%s", aux.c_str());
	// looking for "fake_entry" causes dictionary to be loaded in memory
	WDict_entries fake_entry = WDict::instance().get_entries("kaka", "");
	if (glVars::dict::altdict_fname.size()) {
		WDict::instance().read_alternate_file(glVars::dict::altdict_fname);
	}
}

void test() {

	size_t l_n = 0;

	glVars::dict::use_shuffle = false;
	string cid, ctx;
	read_ukb_ctx(cin, l_n, cid, ctx);
	CSentence cs(cid, ctx);
	cs.debug(cerr);
}

int main(int argc, char *argv[]) {


	map<string, dgraph_rank> map_dgraph_ranks;

	map_dgraph_ranks["ppr"] = r_ppr;
	map_dgraph_ranks["ppr_w2w"] = r_ppr_w2w;
	map_dgraph_ranks["degree"] = r_degree;
	map_dgraph_ranks["static"] = r_static;

	bool opt_do_test = false;
	bool opt_client = false;
	bool opt_shutdown = false;

	cmdline = string("!! -v ");
	cmdline += glVars::ukb_version;
	for (int i=0; i < argc; ++i) {
		cmdline += " ";
		cmdline += argv[i];
	}

	vector<string> input_files;
	string fullname_in;
	ifstream input_ifs;

#ifdef UKB_SERVER
	unsigned int port = 10000;
#endif
	size_t iterations = 0;
	float thresh = 0.0;
	bool check_convergence = false;

	using namespace boost::program_options;

	const char desc_header[] = "ukb_wsd: perform WSD with KB based algorithm\n"
		"Usage examples:\n"
		"ukb_wsd -D dict.txt -K kb.bin --ppr input.txt    -> Disambiguate input.txt using PPR technique according to graph kb.bin and dictionary dict.txt\n"
		"ukb_wsd -D dict.txt -K kb.bin --dgraph_dfs input.txt -> Disambiguate input.txt using disambiguation graph technique, according to graph kb.bin and dictionary dict.txt\n"
		"Options";

	//options_description po_desc(desc_header);

	options_description po_desc("General");

	po_desc.add_options()
		("help,h", "This page")
		("version", "Show version.")
		("kb_binfile,K", value<string>(), "Binary file of KB (see compile_kb).")
		("dict_file,D", value<string>(), "Dictionary text file.")
		("dict_binfile", value<string>(), "Dictionary binary file.")
		;

	options_description po_desc_wsd("WSD methods");
	po_desc_wsd.add_options()
		("ppr", "Given a text input file, disambiguate context using Personalized PageRank method.")
		("ppr_w2w", "Given a text input file, disambiguate context using Personalized PageRank method word by word (see README).")
		("static", "Given a text input file, disambiguate context using static pageRank over kb.")
		("dgraph_bfs", "Given a text input file, disambiguate context using disambiguation graph mehod (bfs).")
		("dgraph_dfs", "Given a text input file, disambiguate context using disambiguation graph mehod (dfs).")
		("nostatic", "Substract static ppv to final ranks.")
		("noprior", "Don't multiply priors to target word synsets if --ppr_w2w and --dict_weight are selected.")
		;

	options_description po_desc_input("Input options");
	po_desc_input.add_options()
		("nopos", "Don't filter words by Part of Speech.")
		("minput", "Do not die when dealing with malformed input.")
		("ctx_noweight", "Do not use weights of input words (defaut is use context weights).")
		;

	options_description po_desc_prank("pageRank general options");
	po_desc_prank.add_options()
		("prank_nibble", "Use PageRank approximation (PageRank-nibble). See also nibble_epsilon.")
		("prank_weight,w", "Use weights in pageRank calculation. Serialized graph edges must have some weight.")
		("prank_iter", value<size_t>(), "Number of iterations in pageRank. Default is 30.")
		("prank_threshold", value<float>(), "Threshold for stopping PageRank. Default is zero. Good value is 0.0001.")
		("prank_damping", value<float>(), "Set damping factor in PageRank equation. Default is 0.85.")
		("dgraph_rank", value<string>(), "Set disambiguation method for dgraphs. Options are: ppr(default), ppr_w2w, coherence, static, degree.")
		("dgraph_maxdepth", value<size_t>(), "If --dgraph_dfs is set, specify the maximum depth (default is 6).")
		("dgraph_nocosenses", "If --dgraph_dfs, stop DFS when finding one co-sense of target word in path.")
		("nibble_epsilon", value<float>(), "Error for approximate pageRank as computed by the nibble algorithm.")
		;

	options_description po_desc_dict("Dictionary options");
	po_desc_dict.add_options()
		("altdict", value<string>(), "Provide an alternative dictionary overriding the values of default dictionary.")
		("dict_weight", "Use weights when linking words to concepts (dict file has to have weights). This is the default setting.")
		("dict_noweight", "Do not use weights when linking words to concepts.")
		("smooth_dict_weight", value<float>(), "Smoothing factor to be added to every weight in dictionary concepts. Default is 1.")
		("dict_strict", "Be strict when reading the dictionary and stop when any error is found.")
		;

	options_description po_desc_output("Output options");
	po_desc_output.add_options()
		("allranks", "Write key file with all synsets associated with ranks.")
		("verbose,v", "Be verbose.")
		("no-monosemous", "Don't output anything for monosemous words.")
		("ties", "Output also in case of ties.")
		("rank_nonorm", "Do not normalize the ranks of target words.")
		("dump_dgraph", "Dump subgraphs using gaphviz format.")
		;

	options_description po_desc_server("Client/Server options");
	po_desc_server.add_options()
		("daemon", "Start a daemon listening to port. Assumes --port")
		("port", value<unsigned int>(), "Port to listen/send information.")
		("client", "Use client mode to send contexts to the ukb daemon. Bare in mind that the configuration is that of the server.")
		("shutdown", "Shutdown ukb daemon.")
		;

	options_description po_visible(desc_header);
	po_visible.add(po_desc).add(po_desc_wsd).add(po_desc_prank).add(po_desc_input).add(po_desc_dict).add(po_desc_output).add(po_desc_server);

	options_description po_hidden("Hidden");
	po_hidden.add_options()
		("bcomp_kb_binfile,M", value<string>(), "Backward compatibility with -K.")
		("bcomp_dictfile,W", value<string>(), "Backward compatibility with -D.")
		("only_ctx_words,C", "Backward compatibility with -C.")
		("concept_graph,G", "Backward compatibility with -G.")
		("dgraph", "Backward compatibility with --dgraph.")
		("nodict_weight", "alias of --dict_noweight")
		("test,t", "(Internal) Do a test.")
		("input-file",value<string>(), "Input file.")
		;
	options_description po_all("All options");
	po_all.add(po_visible).add(po_hidden);

	positional_options_description po_optdesc;
	po_optdesc.add("input-file", 1);
	//    po_optdesc.add("output-file", 1);

	try {

		variables_map vm;
		store(command_line_parser(argc, argv).
			  options(po_all).
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

		if (vm.count("nopos")) {
			glVars::input::filter_pos = false;
		}

		if (vm.count("minput")) {
			glVars::input::swallow = true;
		}

		if (vm.count("ctx_noweight")) {
			glVars::input::weight = false;
		}

		if (vm.count("ppr")) {
			opt_dmethod = m_ppr;
		}

		if (vm.count("ppr_w2w")) {
			opt_dmethod = m_ppr_w2w;
		}

		if (vm.count("static")) {
			opt_dmethod = m_static;
		}

		if (vm.count("dgraph_bfs")) {
			opt_dmethod = m_bfs;
		}

		if (vm.count("dgraph")) {
			opt_dmethod = m_bfs;
		}

		if (vm.count("dgraph_dfs")) {
			opt_dmethod = m_dfs;
		}

		if (vm.count("prank_iter")) {
			iterations = vm["prank_iter"].as<size_t>();
			check_convergence = true;
		}

		if (vm.count("prank_threshold")) {
			thresh = vm["prank_threshold"].as<float>();
			check_convergence = true;
		}

		if (vm.count("prank_damping")) {
			float dp = vm["prank_damping"].as<float>();
			if (dp <= 0.0 || dp > 1.0) {
				cerr << "Error: invalid prank_damping value " << dp << "\n";
				exit(-1);
			}
			glVars::prank::damping = dp;
		}

		if (vm.count("prank_nibble")) {
			glVars::prank::impl = glVars::nibble;
		}

		if (vm.count("nibble_epsilon")) {
			float dp = vm["nibble_epsilon"].as<float>();
			if (dp <= 0.0 || dp > 1.0) {
				cerr << "Error: invalid nibble_epsilon value " << dp << "\n";
				exit(-1);
			}
			glVars::prank::nibble_epsilon = dp;
		}

		if (vm.count("dgraph_maxdepth")) {
			size_t md = vm["dgraph_maxdepth"].as<size_t>();
			if (md == 0) {
				cerr << "Error: invalid dgraph_maxdepth of zero\n";
				exit(-1);
			}
			glVars::dGraph::max_depth = md;
		}

		if (vm.count("dgraph_nocosenses")) {
			glVars::dGraph::stopCosenses = true;
		}

		if (vm.count("dgraph_rank")) {
			string str = vm["dgraph_rank"].as<string>();
			map<string, dgraph_rank>::iterator it = map_dgraph_ranks.find(str);
			if (it == map_dgraph_ranks.end()) {
				cerr << "Error: invalid dgraph_rank method. Should be one of: ";
				for(map<string, dgraph_rank>::iterator iit = map_dgraph_ranks.begin();
					iit != map_dgraph_ranks.end(); ++iit)
					cerr << " " << iit->first;
				cerr << "\n";
				exit(-1);
			}
			dgraph_rank_method = it->second;
		}

		if (vm.count("nostatic")) {
			glVars::csentence::disamb_minus_static = true;
		}

		if (vm.count("noprior")) {
			glVars::csentence::mult_priors = false;
		}

		if (vm.count("bcomp_dictfile")) {
			glVars::dict::text_fname = vm["bcomp_dictfile"].as<string>();
		}

		if (vm.count("dict_file")) {
			glVars::dict::text_fname = vm["dict_file"].as<string>();
		}

		if (vm.count("dict_binfile")) {
			glVars::dict::bin_fname = vm["dict_binfile"].as<string>();
		}

		if (vm.count("dict_strict")) {
			glVars::dict::swallow = false;
		}

		if (vm.count("dict_weight")) {
			glVars::dict::use_weight = true;
		}

		if (vm.count("dict_noweight")) {
			glVars::dict::use_weight = false;
		}

		if (vm.count("nodict_weight")) {
			glVars::dict::use_weight = false;
		}

		if (vm.count("smooth_dict_weight")) {
			glVars::dict::weight_smoothfactor = vm["smooth_dict_weight"].as<float>();
		}

		if (vm.count("bcomp_kb_binfile")) {
			glVars::kb::fname = vm["bcomp_kb_binfile"].as<string>();
		}

		if (vm.count("kb_binfile")) {
			glVars::kb::fname = vm["kb_binfile"].as<string>();
		}

		if (vm.count("rank_alg")) {
			glVars::RankAlg alg = glVars::get_algEnum(vm["rank_alg"].as<string>());
			if (alg == glVars::no_alg) {
				cerr << "Error: Undefined rank algorithm " << vm["rank_alg"].as<string>() << endl;
				exit(-1);
			}
			glVars::rAlg = alg;
		}

		if (vm.count("prank_weight")) {
			glVars::prank::use_weight = true;
		}

		if (vm.count("input-file")) {
			fullname_in = vm["input-file"].as<string>();
		}

		if (vm.count("verbose")) {
			glVars::verbose = 1;
			glVars::debug::warning = 1;
		}

		if (vm.count("allranks")) {
			glVars::output::allranks = true;
		}

		if (vm.count("dump_dgraph")) {
			opt_dump_dgraph = true;
		}

		if (vm.count("test")) {
			opt_do_test = true;
		}

		if (vm.count("no-monosemous")) {
			glVars::output::monosemous = false;
		}

		if (vm.count("ties")) {
			glVars::output::ties = true;
		}

		if (vm.count("rank_nonorm")) {
			glVars::output::norm_ranks = false;
		}

		if (vm.count("altdict")) {
			glVars::dict::altdict_fname = vm["altdict"].as<string>();
		}

		if (vm.count("daemon")) {
#ifdef UKB_SERVER
			opt_daemon = true;
#else
		cerr << "[E] server not available (compile ukb without -DUKB_SERVER switch)\n";
		exit(1);
#endif
		}

		if (vm.count("port")) {
#ifdef UKB_SERVER
			port = vm["port"].as<unsigned int>();;
#else
		cerr << "[E] server not available (compile ukb without -DUKB_SERVER switch)\n";
		exit(1);
#endif
		}

		if (vm.count("client")) {
#ifdef UKB_SERVER
			opt_client = true;
#else
			cerr << "[E] server not available (compile ukb without -DUKB_SERVER switch)\n";
			exit(1);
#endif
		}

		if (vm.count("shutdown")) {
#ifdef UKB_SERVER
			opt_shutdown = true;
#else
			cerr << "[E] server not available (compile ukb without -DUKB_SERVER switch)\n";
			exit(1);
#endif
		}

	} catch(std::exception& e) {
		cerr << e.what() << "\n";
		exit(-1);
	}

	if(opt_shutdown) {
#ifdef UKB_SERVER
		if (client_stop_server(port)) {
			cerr << "Stopped UKB daemon on port " << lexical_cast<string>(port) << "\n";
			return 0;
		} else {
			cerr << "Can not stop UKB daemon on port " << lexical_cast<string>(port) << "\n";
			return 1;
		}
#endif
	}

	// if not daemon, check input files (do it early before loading KB and dictionary)
	if (!fullname_in.size() and !opt_daemon) {
		cout << po_visible << endl;
		cout << "Error: No input" << endl;
		exit(-1);
	}

	if (check_convergence) set_pr_convergence(iterations, thresh);

	// if --daemon, fork server process (has to be done before loading KB and dictionary)
	if (opt_daemon) {
#ifdef UKB_SERVER
		try {
			// Get absolute names of KB and dict
			glVars::kb::fname =  get_fname_absolute(glVars::kb::fname);
			glVars::dict::text_fname = get_fname_absolute(glVars::dict::text_fname);
			glVars::dict::altdict_fname = get_fname_absolute(glVars::dict::altdict_fname);
		} catch(std::exception& e) {
			cerr << e.what() << "\n";
			return 1;
		}
		if (!glVars::kb::fname.size()) {
			cerr << "Error: no KB file\n";
			return 1;
		}
		// accept malformed contexts, as we don't want the daemon to die.
		glVars::input::swallow = true;
		cout << "Starting UKB WSD daemon on port " << lexical_cast<string>(port) << " ... \n";
		return start_daemon(port, &load_kb_and_dict, &handle_server_read);
#endif
	}

	// if not --client, load KB
	if (!opt_client) {
		if (!glVars::kb::fname.size()) {
			cerr << "Error: no KB file\n";
			exit(1);
		}
		try {
			load_kb_and_dict(false);
		} catch (std::exception & e) {
			cerr << e.what() << "\n";
			return 1;
		}
	}

	// create stream from input file
	if (fullname_in == "-" ) {
		// read from <STDIN>
		cmdline += " <STDIN>";
		fullname_in = "<STDIN>";
	} else {
		input_ifs.open(fullname_in.c_str(), ofstream::in);
		if (!input_ifs) {
			cerr << "[E] Can't open " << fullname_in << endl;
			exit(-1);
		}
		// redirect std::cin to read from file
		std::cin.rdbuf(input_ifs.rdbuf());
	}

	if (opt_client) {
#ifdef UKB_SERVER
		// TODO :
		// - check parameters
		// - cmdline
		return !client(std::cin, std::cout, port);
#endif
	}

	if (opt_do_test) {
		test();
		return 0;
	}

	cout << cmdline << "\n";

	try {
		dispatch_run(std::cin, std::cout);
	} catch (std::exception & e) {
		cerr << "[E] Error reading " << fullname_in << "\n" << e.what() << "\n";
		return 0;
	}

	return 0;
}
