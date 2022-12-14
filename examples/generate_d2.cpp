//--------------------------------------------------------------------
//
// This file is part of PEACE.
// 
// PEACE is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// PEACE is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with PEACE.  If not, see <http://www.gnu.org/licenses/>.
// 
// Miami University makes no representations or warranties about the
// suitability of the software, either express or implied, including
// but not limited to the implied warranties of merchantability,
// fitness for a particular purpose, or non-infringement.  Miami
// University shall not be liable for any damages suffered by licensee
// as a result of using, result of using, modifying or distributing
// this software or its derivatives.
//
// By using or copying this Software, Licensee agrees to abide by the
// intellectual property laws, and all other applicable laws of the
// U.S., and the terms of GNU General Public License (version 3).
//
// Authors:   John Karro               karroje@muohio.edu
//
//---------------------------------------------------------------------

#include "ESTAnalyzerFactory.h"
#include "HeuristicChain.h"
#include "ParameterSetManager.h"

#include "ESTAnalyzer.h"
#include "Heuristic.h"
#include "ESTList.h"

#include "generate_d2.h"

#include <fstream>
#include <cstdlib>
#include <ctime>
#include <unistd.h>

const int max_id_length = 1000;
const int max_string_length = 100000;

void parse_args(int argc, char** argv);

// Parameters for d2 and heurstics
double window_fraction = 0; 
int window_length = 100;
int word_length = 6;

int uv_u = 8;
int uv_skip = 8;

int tv_t = 40;

bool use_two_pass = true;

// Parameters for simulaton
int segmentLength = 500;
int min_overlap = 1;
int max_overlap = -1;
int num_trials = 10;
double error_rate = 0.03;
double N_rate = 0.00;

// Parameters for input
long seed = -1;
string file = "all_zf_cdnas.reduced.fa";

// Parameters for output
string output_file = "";
bool unixCRLF = false;
bool header = true;
string header_string = "overlap\td2\t2pd2\tuvtv";
bool header_only = false;

bool help = false;
bool find_error = false;

string add_errors(double error_rate, string& s) {
    for (string::iterator i=s.begin(); i != s.end(); i++) {
        if (drand48() < N_rate) {
            cout << "HERE";
            *i = 'N';
        }
        else if (drand48() < error_rate) {
            int shift = (int)(((int)3*drand48()) + 1);
            for (int j=0; j < shift; j++) {
                switch (*i) {
                case 'A' : *i = 'C'; break;
                case 'C' : *i = 'G'; break;
                case 'G' : *i = 'T'; break;
                case 'T' : *i = 'A'; break;
                case 'a' : *i = 'c'; break;
                case 'c' : *i = 'g'; break;
                case 'g' : *i = 't'; break;
                case 't' : *i = 'a'; break;
                }
            }
        }
    }
    return s;
}

string make_id(int i) {
    char b[100];
    sprintf(b, "%d", i);
    return (string)b;
}

int main(int argc, char** argv) {

    parse_args(argc, argv);

    //*****************************
    // Setup
    ifstream fin(file.c_str());
    assert(fin);

    ostream* out = output_file=="" ? &cout : new ofstream(output_file.c_str());

    if (seed < 0) {
        seed = (long) time(NULL) + getpid();
    }
    srand48(seed);

    string eol = unixCRLF ? "\n" : "\r\n";

    //*****************************
    // Print header in file (if necessary0
    if (header) 
        *out << header_string << eol;
    if (header_only) 
        return 0;


    //*****************************
    // Read in sequence files
    char id_str[max_id_length];
    char s[max_string_length];
    vector<string> seqs;


    // NOTE: Assumes that the sequence is contained on one line!
    while (fin.getline(id_str, max_id_length)) {
        assert(fin.getline(s, max_string_length));
        string seq = (string) s;
        if ((int)seq.length() > 2*segmentLength)
            seqs.push_back(s);
    }

    int id = 0;
    vector<int> start_coords;
    ESTList estList;
  
    //*************************************
    // Generate the overlapping sequences
    for (int overlap=min_overlap; overlap <= max_overlap; overlap++) {
        for (int i=0; i < num_trials; i++) {
            int index = (int)(seqs.size()*drand48());
            string& seq = seqs[index];
      
            int start1 = (int)((seq.length() - 2*segmentLength + overlap)*drand48());
            int start2 = start1 + segmentLength - overlap;

            start_coords.push_back(start1);
            start_coords.push_back(start2);

            string s1 = seq.substr(start1, segmentLength);
            s1 = add_errors(error_rate, s1);
            estList.add(id, make_id(id), s1);

            string s2 = seq.substr(start2, segmentLength);
            s2 = add_errors(error_rate, s2);
            estList.add(id+1, make_id(id+1), s2);

            id += 2;
        }
    }
    int id_overlap = id;

    //******************************************
    // Generate the non-overlapping sequences
    for (int i=0; i < num_trials; i++) {
        int index1, index2;
    
        do {
            index1 = (int)(seqs.size()*drand48());
            index2 = (int)(seqs.size()*drand48());
        } while (index1 == index2);

        string seq1 = seqs[index1];
        string seq2 = seqs[index2];
        int start1 = (int)((seq1.length() - segmentLength)*drand48());
        int start2 = (int)((seq2.length() - segmentLength)*drand48());

        start_coords.push_back(start1);
        start_coords.push_back(start2);

        string s1 = seq1.substr(start1, segmentLength);
        s1 = add_errors(error_rate, s1);
        estList.add(id, make_id(id), s1);
    
        string s2 = seq2.substr(start2, segmentLength);
        s2 = add_errors(error_rate, s2);
        estList.add(id+1, make_id(id+1), s2);
    
        id += 2;
    }

    //*****************************************
    // Setup heursitics and distance functions

    //************
    // Set up d2
    char frame[5];
    sprintf(frame, "%d", window_length);

    char word[5];
    sprintf(word, "%d", word_length);

    //********************
    // Set up huristics
    HeuristicChain chain;
    chain.setESTList(&estList);
    chain.setupChain("tv");
    // Initialize the parameter manager    
    chain.getParamSetMgr()->setupParameters(tv_t, uv_u, uv_skip,
                                            tv_t, uv_u, uv_skip,
                                            tv_t, uv_u, uv_skip);
    // Initialize the chain and all underlying heuristics
    Heuristic* uv_tv = chain.getHeuristic("tv");
    assert( uv_tv != NULL );

    std::auto_ptr<ESTAnalyzer> d2(ESTAnalyzerFactory::create("d2"));    
    d2->setArgument("--frame", window_length);
    d2->setArgument("--word", word_length);
    d2->setArgument("--noNormalize", true);
    d2->setESTList(&estList);
    d2->setHeuristicChain(&chain);
    d2->initialize();

    std::auto_ptr<ESTAnalyzer> p2_d2(NULL);
    if (use_two_pass) {
        p2_d2.reset(ESTAnalyzerFactory::create("twoPassD2"));
        p2_d2->setArgument("--frame", window_length);
        p2_d2->setArgument("--word", word_length);
        p2_d2->setArgument("--noNormalize", true);
        p2_d2->setESTList(&estList);
        p2_d2->setHeuristicChain(&chain);
        p2_d2->initialize();
    }
  
    //*********************
    // Computations and ouput
    for (int i=0; i < id_overlap; i += 2) {
        d2->setReferenceEST(estList[i]);
        const float metric = d2->analyze(estList[i+1]);

        chain.setReferenceEST(estList.get(i));
        bool uv_tv_result = uv_tv->shouldAnalyze(estList.get(i+1));
        *out << start_coords[i] + segmentLength - start_coords[i+1]
             << "\t" << metric << "\t";

        if (use_two_pass) {
            p2_d2->setReferenceEST(estList[i]);
            const float metric2 = p2_d2->analyze(estList[i+1]);
            *out << metric2 << "\t";
	
            if (find_error && metric2 < metric) {
                cerr << estList[i]->getSequence()   << "\n";
                cerr << estList[i+1]->getSequence() << "\n";
                cerr << metric << "\t" << metric2   << "\n\n";
            }
        }

        *out << (uv_tv_result ? "TRUE" : "FALSE") << eol;
    }
  
    for (int i=id_overlap; i < id; i += 2) {
        d2->setReferenceEST(estList[i]);
        const float metric = d2->analyze(estList[i+1]);

        chain.setReferenceEST(estList[i]);
        bool uv_tv_result = uv_tv->shouldAnalyze(estList[i+1]);
        *out << 0 << "\t"  << metric << "\t";

        if (use_two_pass) {
            p2_d2->setReferenceEST(estList[i]);
            const float metric2 = p2_d2->analyze(estList[i+1]);
            *out << metric2 << "\t";
	
            if (find_error && metric2 < metric) {
                cerr << estList[i]->getSequence() << "\n";
                cerr << estList[i+1]->getSequence() << "\n";
                cerr << metric << "\t" << metric2 << "\n\n";
            }
        }

        *out << (uv_tv_result ? "TRUE" : "FALSE") << eol;
    }

    return 0;
}

  
void parse_args(int argc, char** argv) {
    int c;
    while ( (c = getopt(argc, argv, "U:S:T:r:o:f:s:w:W:x:t:m:n:e:N:2uhiHF")) != -1 ) {
        switch (c) {

            // Heuristic / distance function parameter selections
        case 'w' : window_length = atoi(optarg); use_two_pass = false; break;
        case 'W' : window_fraction = atof(optarg); use_two_pass = false; break;
        case 'x' : word_length = atoi(optarg); break;
        case 'U' : uv_u = atoi(optarg); break;
        case 'S' : uv_skip = atoi(optarg); break;
        case 'T' : tv_t = atoi(optarg); break;

            // Simulation paramteters
        case '2' : use_two_pass = false; break;
        case 's' : segmentLength = atoi(optarg); break;
        case 't' : num_trials = atoi(optarg); assert(num_trials >= 0); break;
        case 'm' : min_overlap = atoi(optarg); assert(min_overlap > 0); break;
        case 'n' : max_overlap = atoi(optarg); assert(max_overlap > 0); break;
        case 'e' : error_rate = (double)atof(optarg); assert(error_rate >= 0); break;
        case 'N' : N_rate = (double)atof(optarg); assert(N_rate >= 0); break;

            // Input parameters
        case 'f' : file = optarg; break;
        case 'r' : seed = atol(optarg); assert(seed >= 0); break;

            // Ouput parameters
        case 'i' : header = false; break;
        case 'u' : unixCRLF = true; break;
        case 'o' : output_file = optarg; break;
        case 'H' : header_only = true; break;

            // Helo
        case 'h' : help = true; break;
        case 'F' : find_error = true; break;
        case '?' : cout << argv[0] << ": Bad switch: -" << (char)c << endl; exit(1);
        }
    }
    if (window_fraction > 0) 
        window_length = (int) (segmentLength * window_fraction);

    if (max_overlap < 0) 
        max_overlap = segmentLength;
    assert(min_overlap < max_overlap);

    if (help) {
        cout << argv[0] << " parameters:\n\
\t Heuristic and distance metric parameters:\n\
\t\t-w: Set window length (default = 100); surpress two-pass\n\
\t\t-W: Set window length as a % of segment size (overrides -w); surpresses two-pass\n\
\t\t-x: Set word length (default = 6)\n\
\t\t-U: Set uv parameters u (default = 8)\n\
\t\t-S: Set uv parameters s (default = 8)\n\
\t\t-T: Set tv parameter t (default = 40)\n\
\tSimulation parameters:\n\
\t\t-2: Turn of two-pass calculations (on by default)\n\
\t\t-s: Set segment length (default = 500)\n\
\t\t-t: Set number of trials (defult = 1000)\n\
\t\t-m: Minimum overlap to be checked (default = 1)\n\
\t\t-n: Maximum overlap to to be checked (default = segmentLength)\n\
\t\t-e: Error rate (default = 0.01)\n\
\t\t-N: Probability of replacing a base with an N (default = 0)\
\t Input parameters:\n\
\t\t-f: File of genes (default = all_zf_cdnas.reduced.fa)\n\
\t\t-r: Set RNG seed (default = time + pid)\n\
\t Output parameters:\n\
\t\t-H: Print header only (default = false)\n\
\t\t-i: Supress header\n\
\t\t-u: End output with linux \\n instead of windows \\r\\n (default = false)\n\
\t\t-o: Output file (default = standadrd out)\n\
\t Other parameters:\n\
\t\t-h: Print help menu\n\
\t\t-F: Print out sequence pair (to stderr) if the X2pd2 < d2\n";
        exit(0);
    }
}
