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

#include "ParameterSetManager.h"

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
double N_rate = 0.00;
bool unixCRLF = false;
long seed = -1;
bool help = false;
bool header = true;

string add_errors(double error_rate, string& s) {
  for (string::iterator i=s.begin(); i != s.end(); i++) {
    if (drand48() < N_rate) {
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
//   cout << "#generate_d2";
//   for (int i=0; i < argc; i++)
//     cout << " " << argv[i] << " ";
//   cout << endl;

  int c;
  while ( (c = getopt(argc, argv, "w:x:")) != -1 ) {
      switch (c) {
      case 'w' : window_length = atoi(optarg); break;
      case 'x' : word_length = atoi(optarg); break;
      }
  }

  if (help) {
    cout << argv[0] << " parameters:\n\
\t-w: Set window length (default = 100)\n\
\t-x: Set word length (default = 6)\n\
\t-h: Print help menu\n\
\n";
    exit(0);
  }

  string s1 = argv[optind];
  string s2 = argv[optind+1];

  EST::create(0, "string1", s1.c_str(), 0);
  EST::create(1, "string2", s2.c_str(), 1);

  char frame[5];
  sprintf(frame, "%d", window_length);

  char word[5];
  sprintf(word, "%d", word_length);


  std::auto_ptr<ESTAnalyzer> d2(ESTAnalyzerFactory::create("d2", 0, ""));
  std::auto_ptr<ESTAnalyzer> p2_d2(ESTAnalyzerFactory::create("twopassD2", 0, ""));
  ParameterSetManager::setupParameters();
  char param0[]  = "generate_d2";   
  char param1[]  = "--frame";   
  char param2[]  = "--word";    
  char param3[]  = "--estFile", value3[] = "<none>";
  char param4[]  = "--noNormalize";

  char* params[] = {param0, // First param (although not used) is needed as
                            // parser assumes it is exectuable name.
                    param1, frame, param2, word, // Set Window & Word size
                    param3, value3, param4};  // Last parameter is mandatory
  int paramCount = sizeof(params) / sizeof(char*);
  d2->parseArguments(paramCount, params);
  d2->initialize();
  p2_d2->parseArguments(paramCount, params);
  p2_d2->initialize();

  d2->setReferenceEST(0);
  const float metric = d2->analyze(1);

  p2_d2->setReferenceEST(0);
  const float metric2 = p2_d2->analyze(1);

  cout << metric << "\t" << metric2 << endl;

  return 0;
}

  
