#include "globalVars.h"
#include "walkandprint.h"
#include "fileElem.h"
#include "ukbServer.h"
#include <string>
#include <iostream>
#include <fstream>
#include <syslog.h>

/*
  #include "wdict.h"
  #include "common.h"
  #include "kbGraph.h"
  #include "disambGraph.h"
*/

// Basename & friends
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

// Program options

#include <boost/program_options.hpp>

using namespace ukb;
using namespace std;
using namespace boost;

// Lots of global variables
// (I should fix this)

bool opt_deepwalk = false;
size_t opt_deepwalk_gamma = 80; // as for (Perozzi et al., 2014)
size_t opt_deepwalk_length = 10; // as for (Perozzi et al., 2014)
string seed_word;
bool opt_vcomponents = false;
bool opt_dictbucket = false;
size_t opt_bucket_size = 0;

bool opt_daemon = false;
IWap * walker = 0;
string cmdline;

static void print_ctx(const vector<string> & emited_words, ostream & o) {

	if(emited_words.size() > 1) {
		vector<string>::const_iterator it = emited_words.begin();
		vector<string>::const_iterator end = emited_words.end();
		if (it != end) {
			--end;
			for(;it != end; ++it) {
				o << *it << " ";
			}
			o << *end << "\n";
		}
	}
}

IWap *create_walker(size_t N) {
	if(opt_deepwalk) return new DeepWalk(opt_deepwalk_gamma, opt_deepwalk_length);
	if (seed_word.size()) return new WapWord(seed_word, N);
	if (opt_vcomponents) return new WapComponents(N);
	if (opt_dictbucket) {
		vector<float> Priors(Kb::instance().size());
		{
			boost::unordered_map<Kb::vertex_descriptor, float> P;
			float N = concept_priors(P);
			for(boost::unordered_map<Kb::vertex_descriptor, float>::iterator it = P.begin(), end = P.end();
				it != end; ++it) {
				Priors[ it->first ] = it->second / N;
			}
		}
		return new Wap(N, opt_bucket_size, Priors);
	}
	return new Wap(N, opt_bucket_size);
}

///////////////////////////////////////////////
// Server/client functions

// Return FALSE means kill server

#ifdef UKB_SERVER
bool handle_server_read(sSession & session) {
	string ctx;
	vector<string> rwctx;
	try {
		if (!session.receive(ctx)) return true;
		if (ctx == "stop") {
			// release walker
			if (walker) delete walker;
			walker = 0;
			return false;
		}
		if (!walker) {
			syslog(LOG_DEBUG | LOG_USER, "%s", "UKB daemon: creating walker");
			walker = create_walker(0);
		}
		session.send(cmdline);
		while(1) {
			if (!session.receive(ctx)) break;
			vector<string>().swap(rwctx);
			walker->next(rwctx);
			ostringstream oss;
			print_ctx(rwctx, oss);
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

bool client(istream & is, ostream & os, unsigned int port, size_t N) {
	// connect to ukb port and send data to it.
	sClient client("localhost", port);
	string server_cmd;
	string go("go");
	if (client.error()) {
		std::cerr << "Error when connecting: " << client.error_str() << std::endl;
		return false;
	}
	string id, ctx, out;
	try {
		client.send(go);
		client.receive(server_cmd);
		os << server_cmd << std::endl;
		while(N--) {
			client.send(go);
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
	// - there is a dictionary name
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

int main(int argc, char *argv[]) {

	cmdline = ("!! -v ");
	cmdline += glVars::ukb_version;
	for (int i=0; i < argc; ++i) {
		cmdline += " ";
		cmdline += argv[i];
	}

	vector<string> input_files;
	string N_str;
	ifstream input_ifs;

	glVars::input::filter_pos = false;
	int opt_srand = 0;

	bool opt_client = false;
	bool opt_shutdown = false;
#ifdef UKB_SERVER
	unsigned int port = 10000;
#endif

	using namespace boost::program_options;

	const char desc_header[] = "ukb_walkandprint: generate contexts for creating neural network embeddings\n"
		"Usage examples:\n"
		"ukb_embedding -D dict.txt -K graph.bin --dict_weigth 10000 -> Create context for 10000 random walks\n"

		"Options";

	//options_description po_desc(desc_header);

	options_description po_desc("General");

	po_desc.add_options()
		("help,h", "This page")
		("version", "Show version.")
		("verbose,v", "Be verbose.")
		("kb_binfile,K", value<string>(), "Binary file of KB (see compile_kb).")
		("dict_file,D", value<string>(), "Dictionary text file.")
		;

	options_description po_desc_waprint("walk and print options");
	po_desc_waprint.add_options()
		("deepwalk", "Use deepwalk algorithm.")
		("deepwalk_walks", value<size_t>(), "Deekwalk walks per vertex (gamma).")
		("deepwalk_length", value<size_t>(), "Deekwalk walk length (t).")
		("indeg", "Prefer vertices with higher indegree when walking.")
		("srand", value<int>(), "Seed number for random number generator.")
		("vsample", "Sample vertices according to static prank.")
		("dsample", "Sample vertices according to prior calculated from dictionary.")
		("csample", "Sample vertices according to graph components.")
		("multilang", "Whether dictionary contains words in many languages.")
		("buckets", value<size_t>(), "Number of buckets used in vertex sampling (default is 10).")
		("wemit_prob", value<float>(), "Probability to emit a word when on an vertex (emit vertex name instead). Default is 1.0 (always emit word).")
		("seed_word", value<string>(), "Select concepts associate to the word to start the random walk. Default is start at any vertex at random.")
		;

	options_description po_desc_prank("pageRank general options");
	po_desc_prank.add_options()
		("prank_damping", value<float>(), "Set damping factor in PageRank equation. Default is 0.85.")
		;

	options_description po_desc_dict("Dictionary options");
	po_desc_dict.add_options()
		("dict_weight", "Use weights when linking words to concepts (dict file has to have weights). This is the default setting.")
		("dict_noweight", "Do not use weights when linking words to concepts.")
		("smooth_dict_weight", value<float>(), "Smoothing factor to be added to every weight in dictionary concepts. Default is 1.")
		("dict_strict", "Be strict when reading the dictionary and stop when any error is found.")
		;

	options_description po_desc_server("Client/Server options");
	po_desc_server.add_options()
		("daemon", "Start a daemon listening to port. Assumes --port")
		("port", value<unsigned int>(), "Port to listen/send information.")
		("client", "Use client mode to send contexts to the ukb daemon. Bare in mind that the configuration is that of the server.")
		("shutdown", "Shutdown ukb daemon.")
		;

	options_description po_visible(desc_header);
	po_visible.add(po_desc).add(po_desc_waprint).add(po_desc_prank).add(po_desc_dict).add(po_desc_server);

	options_description po_hidden("Hidden");
	po_hidden.add_options()
		("test,t", "(Internal) Do a test.")
		("nodict_weight", "alias of --dict_noweight")
		("N",value<string>(), "Number of RW.")
		;
	options_description po_all("All options");
	po_all.add(po_visible).add(po_hidden);

	positional_options_description po_optdesc;
	po_optdesc.add("N", 1);
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

		if (vm.count("wemit_prob")) {
			glVars::wap::wemit_prob = vm["wemit_prob"].as<float>();
			if (glVars::wap::wemit_prob < 0.0f || glVars::wap::wemit_prob > 1.0f) {
				cerr << "Error: invalid wemit_prob value:" << glVars::wap::wemit_prob << "\n";
				exit(1);
			}
		}

		if (vm.count("prank_damping")) {
			float dp = vm["prank_damping"].as<float>();
			if (dp <= 0.0 || dp > 1.0) {
				cerr << "Error: invalid prank_damping value " << dp << "\n";
				exit(1);
			}
			glVars::prank::damping = dp;
		}

		if (vm.count("dict_file")) {
			glVars::dict::text_fname = vm["dict_file"].as<string>();
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

		if (vm.count("deepwalk")) {
			opt_deepwalk = true;
		}

		if (vm.count("deepwalk_walks")) {
			opt_deepwalk_gamma = vm["deepwalk_walks"].as<size_t>();
			if (!opt_deepwalk_gamma) {
				cerr << "Error: --deepwalk_walks has to be possitive.\n";
				exit(1);
			}
		}

		if (vm.count("deepwalk_length")) {
			opt_deepwalk_length = vm["deepwalk_length"].as<size_t>();
			if (!opt_deepwalk_length) {
				cerr << "Error: --deepwalk_length has to be possitive.\n";
				exit(1);
			}
		}

		if (vm.count("srand")) {
			opt_srand = vm["srand"].as<int>();
			if (opt_srand == 0) {
				cerr << "Error: --srand parameter can not be zero.\n";
				exit(1);
			}
		}

		if (vm.count("vsample")) {
			opt_bucket_size = 10 ; // default value
		}

		if (vm.count("dsample")) {
			opt_bucket_size = 10 ; // default value
			opt_dictbucket = true;
		}

		if (vm.count("buckets")) {
			opt_bucket_size = vm["buckets"].as<size_t>();
		}

		if (vm.count("csample")) {
			opt_vcomponents = true;
		}

		if (vm.count("multilang")) {
			glVars::wap::multilang = true;
		}

		if (vm.count("kb_binfile")) {
			glVars::kb::fname = vm["kb_binfile"].as<string>();

		}

		if (vm.count("seed_word")) {
			seed_word = vm["seed_word"].as<string>();
		}

		if (vm.count("N")) {
			N_str = vm["N"].as<string>();
		}

		if (vm.count("indeg")) {
			glVars::wap::prefer_indegree = true;
		}

		if (vm.count("verbose")) {
			glVars::verbose = 1;
			glVars::debug::warning = 1;
		}

		if (vm.count("test")) {
			//opt_do_test = true;
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

	if (!opt_srand) {
		glVars::rnd::init_random_device();
	} else {
		glVars::rnd::init_random_device(static_cast<int>(opt_srand));
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

	// if --daemon, fork server process (has to be done before loading KB and dictionary)
	if (opt_daemon) {
#ifdef UKB_SERVER
		if (opt_deepwalk) {
			cerr << "[E] deepwalk not supported as client/server (sorry)\n";
			exit(1);
		}
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
		cout << "Starting UKB walk&print daemon on port " << lexical_cast<string>(port) << " ... \n";
		return start_daemon(port, &load_kb_and_dict, &handle_server_read);
#endif
	}

	// if not --client, load KB
	if (!opt_client) {
		if (!glVars::kb::fname.size()) {
			cerr << po_visible << endl;
			cout << "Error: no KB file\n";
			exit(1);
		}
		try {
			load_kb_and_dict(false);
		} catch (std::exception & e) {
			cerr << e.what() << "\n";
			return 1;
		}
	}
	size_t N = 0;
	if (opt_client) {
#ifdef UKB_SERVER
		if (!N_str.size()) {
			cout << po_visible << endl;
			cout << "Please specify the number of random walks" << endl;
			exit(-1);
		}
		return !client(std::cin, std::cout, port, lexical_cast<size_t>(N_str));
#endif
	}

	try {
		cout << cmdline << "\n";
		if (!opt_deepwalk) {
			if (!N_str.size()) {
				cout << po_visible << endl;
				cout << "Please specify the number of random walks" << endl;
				exit(-1);
			}
			N = lexical_cast<size_t>(N_str);
		}
		walker = create_walker(N);
		vector<string> ctx;
		while(walker->next(ctx))
			print_ctx(ctx, cout);
	} catch(std::exception& e) {
		cerr << e.what() << "\n";
		exit(-1);
	}
}
