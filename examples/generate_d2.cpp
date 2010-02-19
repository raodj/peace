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

#include <fstream>
#include <unistd.h>

#include "ClusterMakerFactory.h"
#include "ESTAnalyzerFactory.h"
#include "HeuristicFactory.h"
#include "HeuristicChain.h"

#include "ClusterMaker.h"
#include "ESTAnalyzer.h"
#include "Heuristic.h"
#include "EST.h"

#include "generate_d2.h"

#include <unistd.h>
#include <assert.h>

#include <time.h>

const int max_id_length = 1000;
const int max_string_length = 100000;

int segmentLength = 500;
int min_overlap = 1;
int max_overlap = 100;
int num_trials = 10;
int window_length = 100;
int word_length = 6;
string file = "all_zf_cdnas.reduced.fa";
string output_file = "";

double error_rate = 0.03;
bool unixCRLF = false;
long seed = -1;
bool help = false;
bool header = true;

string add_errors(double error_rate, string& s) {
  for (string::iterator i=s.begin(); i != s.end(); i++) {
    if (drand48() < error_rate) {
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
  cout << "#generate_d2";
  for (int i=0; i < argc; i++)
    cout << " " << argv[i] << " ";
  cout << endl;

  int c;
  while ( (c = getopt(argc, argv, "r:o:f:s:w:x:t:m:n:e:uhi")) != -1 ) {
      switch (c) {
      case 's' : segmentLength = atoi(optarg); break;
      case 'w' : window_length = atoi(optarg); break;
      case 'x' : word_length = atoi(optarg); break;
      case 't' : num_trials = atoi(optarg); assert(num_trials >= 0); break;
      case 'm' : min_overlap = atoi(optarg); assert(min_overlap > 0); break;
      case 'n' : max_overlap = atoi(optarg); assert(max_overlap > 0); break;
      case 'e' : error_rate = (double)atof(optarg); assert(error_rate >= 0); break;
      case 'i' : header = false; break;
      case 'u' : unixCRLF = true; break;
      case 'f' : file = optarg; break;
      case 'r' : seed = atol(optarg); assert(seed >= 0); break;
      case 'o' : output_file = optarg; break;
      case 'h' : help = true; break;
      case '?' : cout << argv[0] << ": Bad switch: -" << (char)c << endl; exit(1);
      }
  }
  assert(min_overlap < max_overlap);

  if (help) {
    cout << argv[0] << " parameters:\n\
\t-s: Set segment length (default = 500)\n\
\t-w: Set window length (default = 100)\n\
\t-x: Set word length (default = 6)\n\
\t-t: Set number of trials (defult = 1000)\n\
\t-m: Minimum overlap to be checked (default = 1)\n\
\t-n: Maximum overlap to to be checked (default = 100)\n\
\t-e: Error rate (default = 0.01)\n\
\t-i: Supress header\n\
\t-u: End output with linux \\n instead of windows \\r\\n (default = false)\n\
\t-f: File of genes (default = all_zf_cdnas.reduced.fa)\n\
\t-r: Set RNG seed (default = time + pid)\n\
\t-o: Output file (default = standadrd out)\n\
\t-h: Print help menu\n\
\n";
    exit(0);
  }


  ifstream fin(file.c_str());
  assert(fin);

  ostream* out = output_file=="" ? &cout : new ofstream(output_file.c_str());

  if (seed < 0) {
    seed = (long) time(NULL) + getpid();
  }
  srand48(seed);

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
      EST::create(id, make_id(id).c_str(), s1.c_str(), id);

      string s2 = seq.substr(start2, segmentLength);
      s2 = add_errors(error_rate, s2);
      EST::create(id+1, make_id(id+1).c_str(), s2.c_str(), id+1);

      id += 2;
    }
  }
  int id_overlap = id;


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
    EST::create(id, make_id(id).c_str(), s1.c_str(), id);
    
    string s2 = seq2.substr(start2, segmentLength);
    s2 = add_errors(error_rate, s2);
    EST::create(id+1, make_id(id+1).c_str(), s2.c_str(), id+1);
    
    id += 2;
  }

  char frame[5];
  sprintf(frame, "%d", window_length);

  char word[5];
  sprintf(word, "%d", word_length);


  std::auto_ptr<ESTAnalyzer> d2(ESTAnalyzerFactory::create("d2", 0, ""));
  char param0[]  = "generate_d2";   
  char param1[]  = "--frame";   
  char param2[]  = "--word";    
  char param3[]  = "--estFile", value3[] = "<none>";


  char* params[] = {param0, // First param (although not used) is needed as
                            // parser assumes it is exectuable name.
                    param1, frame, param2, word, // Set Window & Word size
                    param3, value3};  // Last parameter is mandatory

  //char* params[] = {param1, frame, param2, word, // Set Window & Word size
  //		    param3, value3};  // Last parameter is mandatory
  int paramCount = sizeof(params) / sizeof(char*);

  string eol = unixCRLF ? "\n" : "\r\n";

  HeuristicChain* chain = HeuristicChain::setupChain("uv-tv", 0, "");
  chain->initialize();
  chain->setReferenceEST(0); // Set reference EST using index value.
  Heuristic* uv = chain->getHeuristic("uv");
  Heuristic* tv = chain->getHeuristic("tv");

  d2->parseArguments(paramCount, params);
  d2->initialize();
  if (header) 
    *out << "overlap d2" << eol;

  for (int i=0; i < id_overlap; i += 2) {
      d2->setReferenceEST(i);
      const float metric = d2->analyze(i+1);

      chain->setReferenceEST(id);
      bool uv_result = uv->shouldAnalyze(id+1);
      bool tv_result = tv->shouldAnalyze(id+1);
      *out << start_coords[i] + segmentLength - start_coords[i+1] << " " << metric << " ";
      *out << (uv_result ? "TRUE" : "FALSE") << " ";
      *out << (tv_result ? "TRUE" : "FALSE") << eol;
  }
  
  for (int i=id_overlap; i < id; i += 2) {
      d2->setReferenceEST(i);
      const float metric = d2->analyze(i+1);

      chain->setReferenceEST(id);
      bool uv_result = uv->shouldAnalyze(id+1);
      bool tv_result = tv->shouldAnalyze(id+1);
      *out << 0 << " "  << metric << " ";
      *out << start_coords[i] + segmentLength - start_coords[i+1] << " " << metric << " ";
      cout << (uv_result ? "TRUE" : "FALSE") << " ";
      cout << (tv_result ? "TRUE" : "FALSE") << eol;
  }

  return 0;
}

  
