#ifndef ALIGNMENT_ALGORITHM_H
#define ALIGNMENT_ALGORITHM_H

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
// Authors:   Jenna Zhang               zhangy9@muohio.edu
//            Dhananjai M. Rao          raodm@muohio.edu
//
//---------------------------------------------------------------------

#include <string>
#include <vector>

/**
   This class implements Global and local alignment algorithms. This
   class is meant to be used in the following manner:

   <ol>

   <li>First create an object with suitable scoring metrics specified
   in the constructor. The constructor sets up the internal scoring
   matrices used to seepdup  alignments.</li>

   <li>Then call the various alignment algorithms (with appropriate
   parameters) to perform different types of alignments supported by
   this class.</li>

   </ol>
*/
class AlignmentAlgorithm {
public:
    /** Instantiate alignment algorithm class with suitable scoring
        metrics.

        \param[in] match This value indicates the score to be used if
        a matching pair of bases (between the two different sequences
        being aligned) is found.

        \param[in] mismatch This value indicates the penalty score to
        be used if a mismatching pair of bases are being aligned.

        \param[in] gap The gap penalty value that indicates the
        penalty associated with opening up a gap during alignment.

		\param[in] aScoreDelPen The penalty value to be used for a
		nucleotide/base deletion occurrence when computing A-scores.

		\param[in] aScoreInsPen The penalty value to be used for a
		nucleotide/base insertion occurrence when computing A-scores.

		\param[in] aScoreSubPen The penalty value to be used for a
		nucleotide/base substitution occurrence when computing
		A-scores.
    */
    AlignmentAlgorithm(const int match = 1, const int mismatch = -1,
                       const int gap = -1,  const int aScoreDelPen = -15,
                       const int aScoreInsPen = -15,
                       const int aScoreSubPen = -5);
    
    /** The destructor.

        The destructor clean up any memory allocated by this class (if
        any).  Mostly the destructor is present to adhere to coding
        conventions.
    */
    ~AlignmentAlgorithm();

    /** Set scoring values to be used when computing the alignments.

        This method can be used to reset the scoring values to be used
        when computing alignments between two nucleotide sequences.
        The supplied values are used to recompute the internal scoring
        matrix that is used by this class to accelerate alignment
        algorithm(s).
        
        \param[in] match This value indicates the score to be used if
        a matching pair of bases (between the two different sequences
        being aligned) is found.
        
        \param[in] mismatch This value indicates the penalty score to
        be used if a mismatching pair of bases are being aligned.
        
        \param[in] gap The gap penalty value that indicates the
        penalty associated with opening up a gap during alignment.
    */
    void setAlignmentParams(const int match, const int mismatch, const int gap);

    /** Set/update the parameters used for computing a-scores.

		\param[in] aScoreDelPen The penalty value to be used for a
		nucleotide/base deletion occurrence when computing A-scores.

		\param[in] aScoreInsPen The penalty value to be used for a
		nucleotide/base insertion occurrence when computing A-scores.

		\param[in] aScoreSubPen The penalty value to be used for a
		nucleotide/base substitution occurrence when computing
		A-scores.

        \note These values will affect the next a-score computations
        performed by the pertinent method(s) in this class.
    */
    void setAScoreParams(const int aScoreDelPen,
                         const int aScoreInsPen,
                         const int aScoreSubPen);

    /** Obtain result from performing Needleman-Wunsch global
        alignment.

        This method can be used to perform global alignment on a given
        pair of nucleotide sequences and returns the results from
        performing the alignment as described below.  The global
        alignment is performed using the Needleman-Wunsch dynamic
        programming algorithm (see http://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm).
        
        \param[in] seq1 The first nucleotide sequence involved in the
        global alignment.  If applicable, this should preferably be
        the the reference gene/transcript.

        \param[in] seq2 The second nucleotide sequence involved in the
        global alignment.  If application, this should preferably be
        the assembler-generated contig.

        \param[in] bandWidth The number of nucleotides around the
        diagonal that need to be evaluated to determine alignment
        information.  Keeping this value small is critical for optimal
        performance and memory-footprint of this method.  However,
        this value cannot be too small to handle large gaps. Suggested
        value is about 20 for sequences that are expected to align
        well.  Pass-in -1 to force evaluation of full alignment
        matrix.
        
        \param[out] alignScore The final \i best alignment score
        generated by the algorithm as it tries to determine the best
        possible global alignment for the two given sequences. The
        scores ultimately depend on the scoring parameters set.

        \param[out] aScore The A-score values generated by the
        algorithm for the best possible global alignment for the two
        sequences.  Note that the first sequence is assumed to be the
        reference sequence for the computations.  The A-scores
        computing by this method are influenced by the scoring
        parameters (see setAScoreParams() method).  The A-score is
        computed using the formula:

        A-score = (2 * seq1_length) - (#aScoreDelPen * num_of_deletions) -
        (#aScoreInsPen * num_of_insertions) - (#aScoreSubPen * num_of_SNPs)
        
        \param[out] alignedSeq1 The final globally aligned result for
        the first sequence. This value is generated by this method.

        \param[out] alignedSeq1 The final globally aligned result for
        the second sequence. This value is generated by this method.
        
        \return This method returns the alignment score which is the
        same value that is stored in alignScore. The value is returned
        as a convenience to streamline coding (where this method is
        called/used).
    */
    int getNWAlignment(const std::string& seq1, const std::string& seq2,
                       int& alignScore, int& aScore, std::string& alignedSeq1,
                       std::string& alignedSeq2) const;

   int getNWAlignment(const std::string& seq1, const std::string& seq2,
					  const int bandwidth,
                       int& alignScore, int& aScore, std::string& alignedSeq1,
                       std::string& alignedSeq2) const;	

   /** Obtain result from performing Smith-Waterman local
          alignment.

          This method can be used to perform local alignment on a given
          pair of nucleotide sequences and returns the results from
          performing the alignment as described below.  The local
          alignment is performed using the Smith-Waterman dynamic
          programming algorithm (see http://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm).

          \param[in] seq1 The first nucleotide sequence involved in the
          global alignment.  If applicable, this should preferably be
          the the reference gene/transcript.

          \param[in] seq2 The second nucleotide sequence involved in the
          global alignment.  If application, this should preferably be
          the assembler-generated contig.

          \param[in] bandWidth The number of nucleotides around the
          diagonal that need to be evaluated to determine alignment
          information.  Keeping this value small is critical for optimal
          performance and memory-footprint of this method.  However,
          this value cannot be too small to handle large gaps. Suggested
          value is about 20 for sequences that are expected to align
          well.  Pass-in -1 to force evaluation of full alignment
          matrix.

          \param[out] alignScore The final \i best alignment score
          generated by the algorithm as it tries to determine the best
          possible global alignment for the two given sequences. The
          scores ultimately depend on the scoring parameters set.

          \param[out] aScore The A-score values generated by the
          algorithm for the best possible global alignment for the two
          sequences.  Note that the first sequence is assumed to be the
          reference sequence for the computations.  The A-scores
          computing by this method are influenced by the scoring
          parameters (see setAScoreParams() method).  The A-score is
          computed using the formula:

          A-score = (2 * seq1_length) - (#aScoreDelPen * num_of_deletions) -
          (#aScoreInsPen * num_of_insertions) - (#aScoreSubPen * num_of_SNPs)

          \param[out] alignedSeq1 The final globally aligned result for
          the first sequence. This value is generated by this method.

          \param[out] alignedSeq1 The final globally aligned result for
          the second sequence. This value is generated by this method.

          \return This method returns the alignment score which is the
          same value that is stored in alignScore. The value is returned
          as a convenience to streamline coding (where this method is
          called/used).
      */
   int getSWAlignment(const std::string& seq1, const std::string& seq2,
                          int& alignScore, int& aScore, std::string& alignedSeq1,
                          std::string& alignedSeq2) const;

protected:
    /** Get the encoding for a given nucleotide character.

        This is a convenience method that is used to obtain the
        encoding for a given nucleotide character.  The encoding is
        used to lookup scores in the scoreMatrix.

        \param[in] c The nucleotide character for which an encoding is
        to be returned.  The currently supported base-pair characters
        are: A, C, G, T, N, P, Y, R, W, S, K, M, D, V, H, B, X.  

        \note Passing-in unsupported characters will cause this method
        to crash and burn.
        
        \return The encoding for the specified base-pair character.
    */
    inline int encodeBase(const char c) const { return baseCode[(int) c]; }

    /** Helper method to compute values in the #scoreMatrix.

        This is an internal helper method that is called from the
        constructor or from the setAlignmentParams() method.  This
        method computes and stores appropriate values in the
        #scoreMatrix.  The #scoreMatrix is used to accelerate the
        various alignment algorithm routines in this class.

        \param[in] match This value indicates the score to be used if
        a matching pair of bases (between the two different sequences
        being aligned) is found.

        \param[in] mismatch This value indicates the penalty score to
        be used if a mismatching pair of bases are being aligned.
    */
    void computeScoreMatrix(const int match, const int mismatch);

	/** Convenience method to verify that the aScoreMatrix is
		symmetric.

		<p>There is a strong assumption that the aScoreMatrix (that is
		used to compute the a-scores) is symmetric to facilitate
		comparison. This method is invoked from the constructor to
		ensure this is indeed the case. This check is useful because
		the matrix is currently set by hand and this method ensures
		that there are no unforeseen typographic errors if-and-when
		the matrix is updated.</p>

		<p>This method checks to ensure that the matrix is symmetric
		and when asymmetry is detected it prints a suitable warning
		message to help troubleshoot the problem.</p>
		
		\note This method is invoked only when the \c
		DEVELOPER_ASSERTION compiler flag is enabled.

		\return This method always returns \c true.
	 */
	static bool verifyAScoreMatrixIsSymmetric();
	
private:
	/** The scoring matrix to make alignment go faster.

		This is a scoring matrix that provides rapid lookup of scores
		when comparing a pair of potentially matching nucleotide bases
		during alignment.  The rows and the columns in the matrix are
		setup in the order: A, C, G, T, N, P, Y, R, W, S, K, M, D, V,
		H, B, X.  This matrix is set via call to the
		computeAScoreMatrix() method.


		\note The order of values in this matrix must match the
		encoding associated with the base pairs as setup in the
		#baseCode array.
	*/
	int scoreMatrix[17][17];

	/** The scoring matrix to compute A-scores quickly

		This is a scoring matrix that provides rapid lookup of scores
		when comparing a pair of nucleotide bases during alignment.
		The rows and the columns in the matrix are setup in the order:
		A, C, G, T, N, P, Y, R, W, S, K, M, D, V, H, B, X.  This
		matrix is a constant matrix that is statically initialized.

		\note The order of values in this matrix must match the
		encoding associated with the base pairs as setup in the
		#baseCode array.
	*/
	static const int aScoreMatrix[17][17];
    
    /** Array to quickly lookup base-pair encodings for various base
        pairs.

        This array is setup in the constructor.  The array is used to
        quickly lookup the index values for various base pairs.  The
        index values are the values in to the two dimensions of the
        #scoreMatrix scoring matrix.  Consequently, the entries in
        this array are the organization of #scoreMatrix are tightly
        coupled (for improved performance).  For example, the \c
        baseCode['G'] would provide the row/column value (the default
        is \c 2) into the #scoreMatrix.
    */
    int baseCode[255];
    
	/** The gap penalty value.
        
		This instance variable contains the gap penalty value.  This
		value is setup in the constructor or via the
		setAlignmentParams() method.
	*/
	int gapPenalty;

    /** The penalty value to be used when computing a-scores when a
        deletion type difference is encountered.  The deletion type
        difference is encountered when aligning a reference
        gene/transcript with an assembler-generated contig.  This
        value is initialized in the constructor and can later on be
        changed via the setAScoreParams() method.
    */
    int aScoreDelPen;

    /** The penalty value to be used when computing a-scores when a
        substitution type difference is encountered.  Substitution
        type differences are detected when aligning a reference
        gene/transcript with an assembler-generated contig.  This
        value is initialized in the constructor and can later on be
        changed via the setAScoreParam() method.
    */
    int aScoreSubPen;

    /** The penalty value to be used when computing a-scores when a
        insertion type difference is encountered.  The insertion type
        difference is encountered when aligning a reference
        gene/transcript with an assembler-generated contig.  This
        value is initialized in the constructor and can later on be
        changed via the setAScoreParams() method.
    */
    int aScoreInsPen;
    
    /** A simple numeration for four directions.

        This enumeration was introduced to primarily make the
        AlignmentAlgorithm codes/methods more readable.
    */
    enum DIR {BAD_DIR, NORTH, WEST, NORTH_WEST};
};

#endif
