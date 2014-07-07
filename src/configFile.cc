#include "configFile.h"
#include "globalVars.h"

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <algorithm>
#include <cctype> // tolower


#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

namespace ukb {

	using namespace std;
	using namespace boost;

	typedef tokenizer<char_separator<char> > tokenizer_t;

	static bool init = false;

	static map<string, int> varMap;

	//static map<string, glVars::HubsAlg> hubsAlgMap;

	//static map<string, glVars::VertexRelate> relateMap;

	static int l_number;

	void init_varMap () {
		varMap["verbose"] = 0;
		varMap["rel_sources"] = 1;
		varMap["text_fname"] = 2;
	}


	void c_error(const string & str) {

		cerr << "Error in line " << l_number << ": " << str << endl;
		exit(-1);

	}

	void c_warn(const string & str) {

		cerr << "Warning in line " << l_number << ": " << str << endl;

	}

	bool parse_csv(const string & str, 
				   vector<string> & V) {

		char_separator<char> sep(" \t,");
		tokenizer_t tok(str.begin(), str.end(), sep);
		copy(tok.begin(), tok.end(), back_inserter(V));
		return (V.size() > 0);
	}

	void parseVarValue(string & variable, const string & value, 
					   int isValue) {

		bool is_negated = false;
		string::const_iterator sit = variable.begin();
		string::const_iterator sit_end = variable.end();
  
		if (sit == sit_end) return;
		if (*sit == '!') {
			is_negated = true;
			++sit;
		}

		string var(sit, sit_end);

		map<string, int>::const_iterator it = varMap.find(var);
		if (it == varMap.end()) {
			if (glVars::verbose) c_warn("[" + var + "] config variable undefined.");
			return;
		}

		//map<string, glVars::HubsAlg>::const_iterator mit_alg;
		//map<string, glVars::VertexRelate>::const_iterator mit_rel;

		switch (it->second) {
		case 0: // verbose
			glVars::verbose = (is_negated == false);
			break;
		case 1: // rel_sources
			if (!parse_csv(value, glVars::rel_source)) {
				c_warn("zero-sized vector");
			}
			break;
		case 2: // text_fname
			glVars::dict::text_fname = value;
			break;
		}
	}


	static bool skip_line(const string &l) {

		if (l.size() == 0) return true;
		string::const_iterator sIt = l.begin();
		string::const_iterator sItEnd = l.end();
		while(sIt != sItEnd && isspace(*sIt)) sIt++;
		if (sIt == sItEnd) return true;
		if (*sIt == '#') return true;
		return false;
	}

	int count_tokens(tokenizer_t & tok) {

		int ema = 0;
		tokenizer_t::const_iterator tit = tok.begin();
		tokenizer_t::const_iterator titEnd = tok.end();
		while(tit != titEnd) {
			++ema;
			++tit;
		}
		return ema;
	}

	// See if the line corresponds to a config parameter
	bool is_config_param(tokenizer_t & tok) {
		string str = *(tok.begin());

		// If the first character is not a number, it is a config parameter
		return !isdigit(str[0]);
  
	}


	void parse_config_fh(istream & fi) {

		vector<string> parValues;

		l_number = -1;
		while (fi) {
			string l;
			string var;
			string value;

			l_number++;
			getline(fi,l,'\n');

			if(skip_line(l)) continue;

			char_separator<char> sep(" \t=");
			tokenizer_t tok(l.begin(), l.end(), sep);
			tokenizer_t::const_iterator tit = tok.begin();
			tokenizer_t::const_iterator titEnd = tok.end();
			assert(tit != titEnd);

			int n = count_tokens(tok);

			if ( n == 0) continue;

			if (is_config_param(tok)) {
				// var = value
				var = *tit;
				transform(var.begin(), var.end(), var.begin(), (int(*)(int))tolower); // lowercase var string
				--n;
				++tit;
				if (n) {
					value = *tit;
					transform(value.begin(), value.end(), value.begin(), (int(*)(int))tolower); // lowercase value string
				} else {
					value = string("");
				}
				parseVarValue(var, value, n);
			} else {
				parValues.clear(); // empty (probably previous used) parameter values vector
				copy(tit, titEnd, back_inserter(parValues));
			}
		}
		//parseParamValues(parValues);
	}

	void parse_config(const std::string & fizena) {
		ifstream fi(fizena.c_str());

		if (!init) {
			init_varMap();
			init = true;
		}
		l_number = 0;
		if (!fi) {
			cerr << "cannot open " << fizena << " parameter file." << endl;
			exit(-1);
		}
		parse_config_fh(fi);
	}
}
