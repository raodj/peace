#ifndef ALIGNMENT_ANALYZER_CPP
#define ALIGNMENT_ANALYZER_CPP

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
// Authors:   Dhananjai M. Rao          raodm@muohio.edu
//
//---------------------------------------------------------------------

#include <vector>
#include <algorithm>
#include "ArgParser.h"
#include "AlignmentAnalyzer.h"
#include "ESTList.h"

AlignmentAnalyzer::AlignmentAnalyzer()
    : ESTAnalyzer("align") {
}

AlignmentAnalyzer::~AlignmentAnalyzer() {
    // Nothing else to be done for now.
}

void
AlignmentAnalyzer::addCommandLineArguments(ArgParser& argParser) {
    // Let base class add common parameters.
    ESTAnalyzer::addCommandLineArguments(argParser);
    // Now setup arguments for this class
    const ArgParser::ArgRecord LocalArgsList[] = {
        {"--matchScore", "Score for matching bases in 2 reads",
         &matchScore, ArgParser::INTEGER},
        {"--mismatchScore", "Score for mismatching bases in 2 reads",
         &mismatchScore, ArgParser::INTEGER},
        {"--gapPenalty", "Penalty for a gap in one of the 2 reads",
         &gapPenalty, ArgParser::INTEGER},
        {"", "", NULL, ArgParser::INVALID}
    };
    // Add valid arguments to global arg parser
    argParser.addValidArguments(LocalArgsList);    
}

bool
AlignmentAnalyzer::initialize() {
    // Let the base class to its initializations
    ESTAnalyzer::initialize();
    // Setup our scoring matrix and look-up tables.
    if (baseCodes.empty()) {
        setScoringMatrix(matchScore, mismatchScore);
    }
    return true;
}

// Return the metric between the refEST and otherEST based on global
// alignment.
float
AlignmentAnalyzer::getMetric(const EST* otherEST) {
    const AlignResult result = getNWAlignment(refEST->getSequence(), otherEST->getSequence());
    float nwDistance = 1.0;
    if (result.score > 0) {
        // Normalize the nw-score into distance using read-length
        nwDistance = (result.str1.length() - result.score) /
            (result.str1.length() * 1.f);
    }
    /*
    std::cout << refEST->getInfo() << " vs " << otherEST->getInfo()
              << ", nw-distance = " << nwDistance << ", raw score = "
              << result.score << ", alignments:\n\t" << result.str1
              << "\n\t" << result.str2 << "\n\t" << refEST->getSequence()
              << "\n\t" << otherEST->getSequence() << std::endl;
    */
    return nwDistance;
}

void
AlignmentAnalyzer::setScoringMatrix(int match, int mismatch) {
    // The list of known bases and amino acids to look-up
    const std::string Codes = "ACDEQGHILKMFTNPRWYScVBX";
    // Initialize to all invalid values (i.e. -1)    
    baseCodes.resize(256, -1);
    // Set positional values in the base codes
    for (size_t i = 0; (i < Codes.size()); i++) {
        if (baseCodes[Codes[i]] != -1) {
            std::cerr << "Duplicate code found: " << Codes[i] << std::endl;
        }
        baseCodes[Codes[i]] = i;
    }

    // Now setup the scoring matrix with all values set to mismatch.
    const std::vector<int> Mismatches(Codes.size(), mismatch);
    scoreMatrix.resize(Codes.size(), Mismatches);
    // Now change the diagonals to matches
    for (size_t i = 0; (i < Codes.size()); i++) {
        scoreMatrix[i][i] = match;
    }
}

int
AlignmentAnalyzer::getNWScore(const std::string& s1,
                              const std::string& s2) const {
    int		r, c, rows, cols, tmp, ins, del, sub, score;
    rows = s1.length() + 1;
    cols = s2.length() + 1;
    int encodedBases1[rows];
    int encodedBases2[cols];

    if (rows < cols) {
        // goes columnwise
        int array[rows];

        // initiate first column
        array[0] = 0;
        for (r = 1; r < rows; r++) {
            array[r] = array[0] + gapPenalty * r;
            encodedBases1[r-1] = encodeBase(s1[r-1]);
        }
        for (r=1; r < cols; r++)
            encodedBases2[r-1] = encodeBase(s2[r-1]);

        // calculate the similarity matrix (keep current column only)
        for (c = 1; c < cols; c++) {
            // initiate first row (tmp hold values
            // that will be later moved to the array)
            tmp = array[0] + gapPenalty;
            int base2 = encodedBases2[c-1];
            for (r = 1; r < rows; r++) {
                int base1 = encodedBases1[r-1];
                ins = array[r] + gapPenalty;
                sub = array[r-1] + (this->scoreMatrix)[base1][base2];
                del = tmp + gapPenalty;
                
                // move the temp value to the array
                array[r-1] = tmp;
                
                // choose the greatest
                if (sub >= ins)
                    tmp = sub >= del ? sub : del;
                else
                    tmp = ins >= del ? ins : del;
            }
            
            // move the temp value to the array
            array[rows - 1] = tmp;
        }
        score = array[rows - 1];
    } else {
        // goes rowwise
        int array[cols];
        
        // initiate first row
        array[0] = 0;
        for (c = 1; c < cols; c++) {
            array[c] = array[0] + gapPenalty * c;
            encodedBases2[c-1] = encodeBase(s2[c-1]);
        }
        for (r=1; r < rows; r++)
            encodedBases1[r-1] = encodeBase(s1[r-1]);
        
        // calculate the similarity matrix (keep current row only)
        for (r = 1; r < rows; r++) {
            // initiate first column (tmp hold values
            // that will be later moved to the array)
            tmp = array[0] + gapPenalty;
            int base1 = encodedBases1[r-1];
            
            for (c = 1; c < cols; c++) {
                int base2 = encodedBases2[c-1];
                ins = tmp + gapPenalty;
                sub = array[c-1] + (this->scoreMatrix)[base1][base2];
                del = array[c] + gapPenalty;
                
                // move the temp value to the array
                array[c-1] = tmp;
                
                // choose the greatest
                if (sub >= ins)
                    tmp = sub >= del ? sub : del;
                else
                    tmp = ins >= del ? ins : del;
            }
            
            // move the temp value to the array
            array[cols - 1] = tmp;
        }
        score = array[cols - 1];
    }
    return score;
}

AlignmentAnalyzer::AlignResult
AlignmentAnalyzer::getNWAlignment(const std::string& s1,
                                  const std::string& s2) const {
    //calculate the alignment score and record the traceback path
    const int numOfRows = s1.length();
    const int numOfCols = s2.length();
    const std::vector<int> ZeroVec(numOfCols + 1);
    std::vector<std::vector<int>> trace(numOfRows + 1, ZeroVec), alignMatrix(trace);
    trace[0][0] = 0;
    
    //initialize the matrix
    for (int i = 1; i <= numOfCols; i++) {
        alignMatrix[0][i] = gapPenalty * i;
        trace[0][i] = 2; // west
    }
    for (int i = 1; i <= numOfRows; i++) {
        alignMatrix[i][0] = gapPenalty * i;
        trace[i][0] = 1; // north
    }
    
    std::vector<int> encodedBases1;
    std::vector<int> encodedBases2;
    for (int i = 0; i <= numOfRows - 1; i++) {
        encodedBases1.push_back(encodeBase(s1[i]));
    }
    for (int i = 0; i <= numOfCols - 1; i++) {
        encodedBases2.push_back(encodeBase(s2[i]));
    }

    // build the matrix row by row
    for (int i = 1; i <= numOfRows; i++) {
        for (int j = 1; j <= numOfCols; j++) {
            int flag = 1;
            // Initialize max to the first of the three terms (NORTH).
            int base1 = encodedBases1[i - 1];
            int base2 = encodedBases2[j - 1];
            int max = alignMatrix[i - 1][j] + gapPenalty;
            
            // See if the second term is larger (WEST).
            int west = alignMatrix[i][j - 1] + gapPenalty;
            if (max <= west) {
                max = west;
                flag = 2;
            }
            
            // See if the third term is the largest (NORTHWEST)
            int northwest = alignMatrix[i - 1][j - 1]
                + scoreMatrix[base1][base2];
            if (max <= northwest) {
                max = northwest;
                flag = 3;
            }
            
            alignMatrix[i][j] = max;
            trace[i][j] = flag;
        }
    }

    AlignResult result;    
    result.score = alignMatrix[numOfRows][numOfCols];

    //trace back and get the alignment strings
    std::string tStr1, tStr2;
    int row = numOfRows;
    int col = numOfCols;
    while ((row != 0) || (col != 0)) {
        int flag = trace[row][col];
        switch (flag) {
        case 1: //i-1, j, north
            tStr1.append(1, s1[row-1]);
            tStr2.append(1, '-');
            row = row - 1;
            break;
        case 2: //i, j-1, west
            tStr1.append(1, '-');
            tStr2.append(1, s2[col-1]);
            col = col - 1;
            break;
        case 3: //i-1, j-1, northwest
            tStr1.append(1, s1[row-1]);
            tStr2.append(1, s2[col-1]);
            row = row - 1;
            col = col - 1;
            break;
        }
    }
    
    //set str1 and str2 in result, they are reverse of tStr1 and tStr2
    std::reverse(tStr1.begin(), tStr1.end());
    std::reverse(tStr2.begin(), tStr2.end());
    result.str1 = tStr1;
    result.str2 = tStr2;

    return result;
}

AlignmentAnalyzer::AlignResult
AlignmentAnalyzer::getSWAlignment(const std::string& s1,
                                  const std::string& s2) const {
    AlignResult result;

    //calculate the alignment score and record the traceback path
    int numOfRows = s1.length();
    int numOfCols = s2.length();
    const std::vector<int> ZeroVec(numOfCols + 1);    
    std::vector<std::vector<int>> trace(numOfRows + 1, ZeroVec), alignMatrix(trace);
    trace[0][0] = 0;
    alignMatrix[0][0] = 0;
    
    //int* encodedBases2 = new int[numOfCols];
    std::vector<int> encodedBases2(numOfCols);
    for (int i = 1; i <= numOfCols; i++) {
        alignMatrix[0][i] = 0;
        trace[0][i] = 0; //start point
        encodedBases2[i-1] = encodeBase(s2[i-1]);
    }

    int maxScore = 0;
    int maxRow = 0;
    int maxCol = 0;

    //build the matrix row by row
    for (int i = 1; i <= numOfRows; i++) {
        int base1 = encodeBase(s1[i - 1]);
        alignMatrix[i][0] = 0;
        trace[i][0] = 0; //start point

        for (int j = 1; j <= numOfCols; j++) {
            int flag = 1;
            // Initialize max to the first of the three terms (NORTH).
            int base2 = encodedBases2[j - 1];
            int max = alignMatrix[i - 1][j] + this->gapPenalty;

            // See if the second term is larger (WEST).
            int west = alignMatrix[i][j - 1] + this->gapPenalty;
            if (max <= west) {
                max = west;
                flag = 2;
            }

            // See if the third term is the largest (NORTHWEST)
            int northwest = alignMatrix[i - 1][j - 1]
                + (this->scoreMatrix)[base1][base2];
            if (max <= northwest) {
                max = northwest;
                flag = 3;
            }

            if (max <= 0) {
                alignMatrix[i][j] = 0;
                trace[i][j] = 0; //start point
            } else {
                alignMatrix[i][j] = max;
                trace[i][j] = flag;
            }
            if (max > maxScore) {
                maxScore = max;
                maxRow = i;
                maxCol = j;
            }
        }
    }

    result.score = alignMatrix[maxRow][maxCol];

    //trace back and get the alignment strings
    std::string tStr1;
    std::string tStr2;
    int row = maxRow;
    int col = maxCol;
    int flag = trace[row][col];
    while (flag != 0) {
        switch (flag) {
        case 1: //i-1, j, north
            tStr1.append(1, s1[row-1]);
            tStr2.append(1, '-');
            row = row - 1;
            break;
        case 2: //i, j-1, west
            tStr1.append(1, '-');
            tStr2.append(1, s2[col-1]);
            col = col - 1;
            break;
        case 3: //i-1, j-1, northwest
            tStr1.append(1, s1[row-1]);
            tStr2.append(1, s2[col-1]);
            row = row - 1;
            col = col - 1;
            break;
        }
        flag = trace[row][col];
    }

    //set str1 and str2 in result, they are reverse of tStr1 and tStr2
    std::reverse(tStr1.begin(), tStr1.end());
    std::reverse(tStr2.begin(), tStr2.end());
    result.str1 = tStr1;
    result.str2 = tStr2;

    return result;
}

#endif
