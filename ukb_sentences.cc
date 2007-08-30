#include "common.h"
#include "w2syn.h"
#include "globalVars.h"
#include "configFile.h"
#include "fileElem.h"
#include "mcrGraph.h"
#include "disambGraph.h"

#include <string>
#include <iostream>
#include <fstream>

// Program options

#include <boost/program_options.hpp>

// timer

#include <boost/timer.hpp>

// bfs 

#include <boost/graph/breadth_first_search.hpp>
#include <boost/pending/indirect_cmp.hpp>
#include <boost/pending/integer_range.hpp>

using namespace std;
using namespace boost;

const char *mcr_default_binfile = "mcr_wnet.bin";

void read_skip_relations (const string & file,
			  set<string> & skip_rels) {

  const string delims(" \t");
  ifstream fi(file.c_str(), ifstream::binary|ifstream::in);
  if (!fi) {
    cerr << "Error: can't open " << file << endl;
    exit(-1);
  }
  string line;
  while(fi) {
    getline(fi, line);
    vector<string> rels = split(line, delims);
    for(vector<string>::iterator it = rels.begin(); it != rels.end(); ++it)
      skip_rels.insert(*it);
  }
}

void add_words_to_mcr() {

}

void compute_sentence_vectors(string & fullname_in, string & out_dir) {

  Mcr & mcr = Mcr::instance();
  if (glVars::verbose) 
    cerr << "Reading words ...\n";
  W2Syn & w2syn = W2Syn::instance();

  ifstream fh_in(fullname_in.c_str());

  if (!fh_in) {
    cerr << "Can't open " << fullname_in << endl;
    exit(-1);
  }

  vector<CSentence> vcs;
  CSentence cs;

  // add words and word#pos into Mcr

  vector<string>::const_iterator word_it = w2syn.get_wordlist().begin();
  vector<string>::const_iterator word_end = w2syn.get_wordlist().end();

  if (glVars::verbose) 
    cerr << "Adding words to Mcr ...";

  size_t word_i = 0;
  for(; word_it != word_end; ++word_it) {
    vector<string>::const_iterator syn_it, syn_end;
    tie(syn_it, syn_end) = w2syn.get_wsyns(*word_it);
    mcr.add_tokens(*word_it, syn_it, syn_end);
    ++word_i;
    if (glVars::verbose && (word_i % 1000 == 0))
      cerr << word_i << " ... ";
  }

  if (glVars::verbose) 
    Mcr::instance().display_info(cerr);
  if (glVars::verbose) 
    cerr << "\nReady for pageRank!\n";

  // Read sentences and compute rank vectors

  try {
    while (cs.read_aw(fh_in)) {

      // Initialize rank vector
      vector<float> ranks;

      // Initialize PPV vector
      vector<float> ppv(mcr.size(), 0.0);
      CSentence::iterator it = cs.begin();
      CSentence::iterator end = cs.end();
      size_t K = 0;
      for(;it != end; ++it) {
	if(!it->is_distinguished()) continue;
	
	string wpos(it->word());
	wpos.append("#");
	char pos = it->get_pos();
	wpos.append(1,pos);
	bool aux;
	Mcr_vertex_t u;
	tie(u, aux) = mcr.getVertexByName(wpos);
	ppv[u] = 1;
	++K;
      }
      if (!K) continue;
      // Normalize PPV vector
      float div = 1.0 / static_cast<float>(K);
      for(vector<float>::iterator rit = ppv.begin(); rit != ppv.end(); ++rit) 
	*rit *= div;      

      // Execute PageRank
      mcr.pageRank_ppv(ppv, ranks);

      File_elem fout(fullname_in, out_dir, ".ppv");
      fout.fname = cs.id();

      ofstream fo(fout.get_fname().c_str(),  ofstream::out);
      if (!fo) {
	cerr << "Error: can't create" << fout.get_fname() << endl;
	exit(-1);
      }
      copy(ranks.begin(), ranks.end(), ostream_iterator<float>(fo, "\n"));
      cs = CSentence();
    }
  } 
  catch (string & e) {
    cerr << "Errore reading " << fullname_in << ":" << e << "\n";
    throw(e);    
  }
}


int main(int argc, char *argv[]) {

  srand(3);

  timer load;

  bool opt_param = false;

  string mcr_binfile(mcr_default_binfile);

  string out_dir;
  string fullname_in;

  const char desc_header[] = "Usage:\n"
    "ukb_sentences context_file.txt\n"
    "  Creates one file per sentence (.ppv extension) with the vector of the sentence"
    "Options:";
  
  using namespace boost::program_options;

  options_description po_desc(desc_header);

  po_desc.add_options()
    ("help,h", "This help page.")
    ("mcr_binfile,M", value<string>(), "Binary file of MCR (see create_mcrbin). Default is mcr_wnet.bin")
    ("out_dir,O", value<string>(), "Directory for leaving output files.")
    ("w2syn_file,W", value<string>(), "Word to synset map file (default is ../Data/Preproc/wn1.6_index.sense_freq).")
    ("param,p", value<string>(), "Specify parameter file.")
    ("verbose,v", "Be verbose.")
    ;
  options_description po_desc_hide("Hidden");
  po_desc_hide.add_options()
    ("input-file",value<string>(), "Input file.")
    ("output-file",value<string>(), "Output file.")
    ;
  options_description po_desc_all("All options");
  po_desc_all.add(po_desc).add(po_desc_hide);
  
  positional_options_description po_optdesc;
  po_optdesc.add("input-file", 1);
  po_optdesc.add("output-file", 1);  

  try {
    variables_map vm;
    store(command_line_parser(argc, argv).
	  options(po_desc_all).
	  positional(po_optdesc).
	  run(), vm);
    notify(vm);

    // If asked for help, don't do anything more

    if (vm.count("help")) {
      cout << po_desc << endl;
      exit(0);
    }

    // verbosity

    if (vm.count("verbose")) {
      glVars::verbose = 1;
    }


    // Config params first, so that they can be overriden by switch options

    if (vm.count("mcr_binfile")) {
      mcr_binfile = vm["mcr_binfile"].as<string>();
    }

    if (vm.count("param")) {
      parse_config(vm["param"].as<string>());
      opt_param = true;
    }

    if (vm.count("out_dir")) {
      out_dir = vm["out_dir"].as<string>();
    }

    if (vm.count("w2syn_file")) {
      glVars::w2s_filename = vm["w2syn_file"].as<string>();
    }

    if (vm.count("input-file")) {
      fullname_in = vm["input-file"].as<string>();
    }

  }
  catch(exception& e) {
    cerr << e.what() << "\n";
    throw(e);
  }

  if(fullname_in.size() == 0) {
    cout << po_desc << endl;
    return 1;
  }

  if (glVars::verbose) 
    cerr << "Reading binary mcr file " << mcr_binfile;
  Mcr::create_from_binfile(mcr_binfile);
  if (glVars::verbose) 
    Mcr::instance().display_info(cerr);

  compute_sentence_vectors(fullname_in, out_dir);

  return 0;
}

/*
 * Local Variables:
 * mode: c++
 * compile-command: "make -k ukb_sentences"
 * End:
 */
