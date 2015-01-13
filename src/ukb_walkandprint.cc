#include "wdict.h"
#include "common.h"
#include "globalVars.h"
#include "kbGraph.h"
#include "disambGraph.h"
#include <string>
#include <iostream>
#include <fstream>

// Basename & friends
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

// Program options

#include <boost/program_options.hpp>
#include <boost/random/random_device.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/discrete_distribution.hpp>

using namespace ukb;
using namespace std;
using namespace boost;

boost::random::mt19937 mt;

float opt_wemit_prob = 1.0; // whether to emit words

template<class It>
int rnumber_distribution(It beg, It end) {
	boost::random::discrete_distribution<> dist(beg, end);
	return dist(mt);
}

// generate between 0 and b, both inclusive
int rnumber(int b) {
    boost::random::uniform_int_distribution<> dist(0, b);
    return dist(mt);
}

// generate between 0 and b, both inclusive
size_t rnumber(size_t b) {
    boost::random::uniform_int_distribution<> dist(0, b);
    return dist(mt);
}

float rnumber(float b) {
    boost::random::uniform_real_distribution<> dist(0.0f, b);
    return dist(mt);
}

///////////////////////////////////////////////


struct vpart_sort_t {

	const vector<float> & m_v ;
	vpart_sort_t(const vector<float> & score) : m_v(score) {}
	int operator()(const int & i, const int & j) {
		// Descending order
		return m_v[i] > m_v[j];
	}
};

struct vsampling_t {

	vector<int> m_idx;
	size_t m_N; // size of graph
	size_t m_bucket_N;
	vector<pair<int, int> > m_intervals;

	void debug() {
		const vector<float> & sprank = Kb::instance().static_prank();
		cout << std::setprecision(4) << 1.0f / static_cast<float>(m_bucket_N) << "\n";
		for(size_t i = 0; i < m_intervals.size(); i++) {
			cout << "Interval " << i << ": ";
			for (int j = m_intervals[i].first; j != m_intervals[i].second; j++) {
				string v = Kb::instance().get_vertex_name(m_idx[j]);
				float rank = sprank[m_idx[j]];
				cout << v << ":" << std::setprecision(4) << rank << " ";
			}
			cout << "\n";
		}
	}

	vsampling_t(size_t buckets = 10) : m_bucket_N(buckets) {

		Kb & kb = Kb::instance();
		m_N = kb.size();
		if (!buckets) return ;
		const vector<float> & sprank = kb.static_prank();
		vector<int>(m_N).swap(m_idx);
		for(size_t i = 0; i < m_N; i++) m_idx[i] = i;
		// sort idx according to static prank
		sort(m_idx.begin(), m_idx.end(), vpart_sort_t(sprank));
		float factor = 1.0f / static_cast<float>(m_bucket_N) ;
		float accum = 0.0;
		pair<int, int> interval = make_pair(0,0);
		for(size_t i = 0; i < m_N; i++) {
			interval.second++;
			accum += sprank[ m_idx[i] ];
			if (accum > factor) {
				m_intervals.push_back(interval);
				if (m_intervals.size() == m_bucket_N) break;
				interval = make_pair(interval.second, interval.second);
				accum = 0.0;
			}
		}
		m_intervals.back().second = m_N;
		m_bucket_N = m_intervals.size();
	}

	// Method to sample a vertex
	// First, decide which bucket, then select one at random from the bucket
	int sample() {
		if (!m_bucket_N) {
			return rnumber(m_N);
		}
		size_t bucket_i = rnumber(m_bucket_N - 1);
		int left = m_intervals[bucket_i].first;
		int m = m_intervals[bucket_i].second - left;
		assert(m > 0);
		size_t off = rnumber(m);
		return m_idx[left + off];
	}
};

///////////////////////////////////////////////


bool emit_word_vertex(Kb_vertex_t current, string & emit_word) {

	static vector<float> vertex2word_tweight;

	Kb & kb = Kb::instance();

	float toss = rnumber(1.0f);
	if (toss >= opt_wemit_prob) {
		emit_word = kb.get_vertex_name(current);
		return true;
	}

	if (!vertex2word_tweight.size()) vector<float>(kb.size(), 0.0f).swap(vertex2word_tweight);

	WInvdict_entries words = WDict::instance().words(current);
	if (!words.size()) return false;
	float & total_weight = vertex2word_tweight[current];
	if (total_weight == 0.0f) {
		WInvdict_entries::freq_const_iterator freq_it( words.begin() );
		WInvdict_entries::freq_const_iterator freq_end( words.end() );
		for(;freq_it != freq_end; ++freq_it) total_weight += *freq_it;
	}
	if (total_weight == 0.0f) return false; // Should never happen
	float rand_value = rnumber(total_weight);
	float w_accum = 0;
	for(size_t i = 0; i < words.size(); ++i) {
		w_accum += words.get_prob(i);
		emit_word = words.get_word(i);
		if (rand_value < w_accum) break;
	}
	return true;
}

bool select_next_vertex(Kb_vertex_t & current) {

	static vector<float> vertex_out_tweight;

	Kb & kb = Kb::instance();
	Kb_vertex_t previous = current;

	if (!vertex_out_tweight.size()) vector<float>(kb.size(), 0.0f).swap(vertex_out_tweight);

	Kb_out_edge_iter_t out_it, out_end;

	tie(out_it, out_end) = kb.out_neighbors(previous);
	float & total_weight = vertex_out_tweight[previous];
	if (total_weight == 0.0f) {
		for(; out_it < out_end; ++out_it) {
			total_weight += kb.get_edge_weight(*out_it);
		}
		tie(out_it, out_end) = kb.out_neighbors(previous);
	}
	if (total_weight == 0.0f) return false; // danglink link. The RW is over.
	float rand_value = rnumber(total_weight);
	float w_accum = 0.0f;
	for(; out_it < out_end; ++out_it) {
		w_accum += kb.get_edge_weight(*out_it);
		current = kb.edge_target(*out_it);
		if (rand_value < w_accum) break;
	}
	return true;
}

void do_mc_complete(Kb_vertex_t v, vector<string> & emited_words) {

	Kb_vertex_t current = v;  //Start the iteration in the V vertex
	vector<string>().swap(emited_words);

	for (float r = rnumber(1.0f); r <= glVars::prank::damping; r = rnumber(1.0f) ) {
		// emit word from current vertex
		string emit_word;
		if (emit_word_vertex(current, emit_word))
			emited_words.push_back(emit_word);
		// Select next vertex to jump
		if (!select_next_vertex(current)) break;
	}
}

void do_mc_synsets(size_t n, size_t bucket_size) {

	Kb & kb = Kb::instance();
	size_t m = kb.size();
	vector<string> emited_words;
	vsampling_t vsampler(bucket_size);
	if (!m) return;
	m--;
	for(size_t i = 0; i < n; ++i) {
		int idx = vsampler.sample();
		Kb_vertex_t u(idx);
		do_mc_complete(u, emited_words);
		if(emited_words.size() > 1) {
			vector<string>::iterator it = emited_words.begin();
			vector<string>::iterator end = emited_words.end();
			if (it != end) {
				--end;
				for(;it != end; ++it) {
					cout << *it << " ";
				}
				cout << *end << "\n";
			}
		}
	}
}

void do_mc_word(const string & hw, size_t n) {

	WDictHeadwords dicthws(WDict::instance());

	static std::tr1::unordered_map<const std::string &, float > dweight_cache;
	std::tr1::unordered_map<const std::string &, float >::iterator dweight_it;
	bool P;

	WDict_entries synsets(WDict::instance().get_entries(hw));

	for(size_t i = 0; i < n; ++i) {
		// select synset at random
		tie(dweight_it, P) = dweight_cache.insert(make_pair(hw, 0.0f));
		float & total_weight = dweight_it->second;
		if (P) {
			for(size_t i = 0; i < synsets.size(); ++i) {
				total_weight += synsets.get_freq(i);
			}
		}
		if (total_weight == 0.0f) break; // word has no attached vertices
		float rand_value = rnumber(total_weight);
		float w_accum = 0;
		Kb_vertex_t synset;
		for(size_t i = 0; i < synsets.size(); ++i) {
			w_accum += synsets.get_freq(i);
			synset = synsets.get_entry(i);
			if (rand_value < w_accum) break;
		}
		// print RW
		vector<string> emited_words;
		do_mc_complete(synset, emited_words);
		// print RW
		if(emited_words.size() > 1) {
			cout << hw << "\t";
			vector<string>::iterator it = emited_words.begin();
			vector<string>::iterator end = emited_words.end();
			if (it != end) {
				--end;
				for(;it != end; ++it) {
					cout << *it << " ";
				}
				cout << *end << "\n";
			}
		}
	}
}

void do_mc_words(size_t n) {

	WDictHeadwords dicthws(WDict::instance());

	static std::tr1::unordered_map<const std::string &, float > dweight_cache;
	std::tr1::unordered_map<const std::string &, float >::iterator dweight_it;
	bool P;
	size_t m = dicthws.size();

	for(size_t i = 0; i < n; ++i) {
		int idx = rnumber((int) m - 1);
		const string & hw = dicthws.hw(idx);

		// select synset at random
		WDict_entries synsets(dicthws.rhs(idx));
		tie(dweight_it, P) = dweight_cache.insert(make_pair(hw, 0.0f));
		float & total_weight = dweight_it->second;
		if (P) {
			for(size_t i = 0; i < synsets.size(); ++i) {
				total_weight += synsets.get_freq(i);
			}
		}
		float rand_value = rnumber(total_weight);
		float w_accum = 0;
		Kb_vertex_t synset;
		for(size_t i = 0; i < synsets.size(); ++i) {
			w_accum += synsets.get_freq(i);
			synset = synsets.get_entry(i);
			if (rand_value < w_accum) break;
		}
		// print RW
		vector<string> emited_words;
		do_mc_complete(synset, emited_words);
		// print RW
		if(emited_words.size()) {
			cout << hw << "\t";
			vector<string>::iterator it = emited_words.begin();
			vector<string>::iterator end = emited_words.end();
			if (it != end) {
				--end;
				for(;it != end; ++it) {
					cout << *it << " ";
				}
				cout << *end << "\n";
			}
		}
	}
}


int main(int argc, char *argv[]) {

	string kb_binfile("");

	string cmdline("!! -v ");
	cmdline += glVars::ukb_version;
	for (int i=0; i < argc; ++i) {
		cmdline += " ";
		cmdline += argv[i];
	}

	vector<string> input_files;
	string N_str;
	size_t N;
	ifstream input_ifs;

	bool opt_do_test = false;
	glVars::input::filter_pos = false;
	size_t opt_bucket_size = 0;

	using namespace boost::program_options;

	const char desc_header[] = "ukb_walkandprint: generate contexts for creating neural network embeddings\n"
		"Usage examples:\n"
		"ukb_embedding -D dict.txt -K graph.bin --dict_weigth 10000 -> Create context for 10000 random walks\n"

		"Options";

	//options_description po_desc(desc_header);

	options_description po_desc("General");
	string seed_word;

	po_desc.add_options()
		("help,h", "This page")
		("version", "Show version.")
		("verbose,v", "Be verbose.")
		("kb_binfile,K", value<string>(), "Binary file of KB (see compile_kb).")
		("dict_file,D", value<string>(), "Dictionary text file.")
		;

	options_description po_desc_waprint("walk and print options");
	po_desc_waprint.add_options()
		("vsample", "Sample vertices according to static prank.")
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
		("dict_weight", "Use weights when linking words to concepts (dict file has to have weights).")
		("smooth_dict_weight", value<float>(), "Smoothing factor to be added to every weight in dictionary concepts. Default is 1.")
		("dict_strict", "Be strict when reading the dictionary and stop when any error is found.")
		;

	options_description po_visible(desc_header);
	po_visible.add(po_desc).add(po_desc_waprint).add(po_desc_prank).add(po_desc_dict);

	options_description po_hidden("Hidden");
	po_hidden.add_options()
		("test,t", "(Internal) Do a test.")
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
			opt_wemit_prob = vm["wemit_prob"].as<float>();
			if (opt_wemit_prob < 0.0f || opt_wemit_prob > 1.0f) {
				cerr << "Error: invalid wemit_prob value:" << opt_wemit_prob << "\n";
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

		if (vm.count("smooth_dict_weight")) {
			glVars::dict::weight_smoothfactor = vm["smooth_dict_weight"].as<float>();
		}

		if (vm.count("vsample")) {
			opt_bucket_size = 10 ; // default value
		}

		if (vm.count("buckets")) {
			opt_bucket_size = vm["buckets"].as<size_t>();
		}

		if (vm.count("kb_binfile")) {
			kb_binfile = vm["kb_binfile"].as<string>();
		}

		if (vm.count("seed_word")) {
			seed_word = vm["word"].as<string>();
		}

		if (vm.count("N")) {
			N_str = vm["N"].as<string>();
		}

		if (vm.count("verbose")) {
			glVars::verbose = 1;
			glVars::debug::warning = 1;
		}

		if (vm.count("test")) {
			opt_do_test = true;
		}
	} catch(std::exception& e) {
		cerr << e.what() << "\n";
		exit(-1);
	}

	boost::random_device rd;
	mt.seed(rd());

	if (!N_str.size()) {
		cout << po_visible << endl;
		cout << "Please specify the number of random walks" << endl;
		exit(-1);
	}

	N = lexical_cast<size_t>(N_str);
	Kb::create_from_binfile(kb_binfile);

	cout << cmdline << "\n";

	if (seed_word.size()) {
		do_mc_word(seed_word, N);
	} else {
		do_mc_synsets(N, opt_bucket_size);
	}

}
