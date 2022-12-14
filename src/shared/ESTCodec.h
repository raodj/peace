#ifndef EST_CODEC_H
#define EST_CODEC_H

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

#include <functional>
#include "HashMap.h"

/** A helper class to serve as a EST enCOder-DECoder.

    <p>This class was introduced to try and centralize a bunch of EST
    encoding/decoding code that was spread out amonst multiple
    independent classes.  Most of the EST analysis algorithms try and
    encode the base pairs (\c A, \c T, \c C, and \c G) in ESTs into
    2-bits each (\c 00, \c 11, \c 10, and \c 01) to reduce memory
    footprint and to enable use of numerical operations.  The code to
    perform encoding and decoding was spread across multiple classes
    on a "as needed" basis.  However, more than three classes required
    pretty much the same code warranting the introduction of this
    class to minimize code-redundancy. <p>

    <p>There is one process-wide unique instance of this ESTCodec
    class.  A single instance is used to create commonly used tables
    to enable rapid CODEC operations.  These tables are created once,
    when the globally unique EST object is referenced.  The globally
    unique instance can be referenced via the getCodec() method. </p>
*/
class ESTCodec {
public:
    /** Obtain reference to process-wide unqiue instance of ESTCodec.

        This method must be used to obtain the process-wide unique
        instance of the ESTCodec object. The returned reference can be
        used to invoke other methods in this class.  Here is a typical
        usage:

        \code

        const ESTCodec& codec = ESTCodec::getCodec();
        int bitCodec = codec.encode('A');

        \endcode
    */
    inline static ESTCodec& getCodec() { return estCodec; }

    /** Obtain 2-bit code for a given base pair.

        This method can be used to obtain the 2-bit encoding for a
        given base pair (bp).  This method essentially translate base
        pairs (\c A, \c T, \c C, and \c G, both upper and lower case)
        into 2-bits codes (\c 00, \c 11, \c 10, and \c 01).

        \note In favor of speed, this method does not perform any
        special checks on the actual character in bp. It is the
        responsiblity of the caller to ensure that this method is
        invoked with appropriate parameter value.
        
        \param[in] bp The base pair character (both upper and lower
        cases are handled correctly) to be encoded.

        \return The 2-bit ecnoding for the given nucleotide.
    */
    inline static char encode(const char bp)
    { return charToInt[(int) bp]; }
    
    /** Obtain base pair character for a given 2-bit encoded base
        pair.

        This method can be used to obtain the actual base pair
        character for a given 2-bit encoding for a nucleotide obtained
        via a call to the ESTCodec::encode() method.  This method
        essentially translates 2-bits codes \c 00, \c 11, \c 10, and
        \c 01 into the base pair characters \c A, \c T, \c C, and \c G
        respectively.  For all other values this method returns '-'
        (hyphen character) as the encoding.

        \param[in] codedBP The 2-bit encoded nucleotide value to be
        converted to corresponding base pair character.

        \return This method returns the base character corresponding
        to the given 2-bit encoding.
    */
    inline static char decode(const unsigned int codedBP)
    { return (codedBP < 4 ? intToChar[codedBP] : '-'); }
    
    /** Obtain 2-bit \b complement code for a given base pair.

        This method can be used to obtain the \b complementary 2-bit
        encoding for a given base pair (bp).  This method essentially
        translate base pairs (\c A, \c T, \c C, and \c G, both upper
        and lower case) into 2-bits codes (\c 11, \c 00, \c 01, and \c
        10).

        \note In favor of speed, this method does not perform any
        special checks on the actual character in bp. It is the
        responsiblity of the caller to ensure that this method is
        invoked with appropriate parameter value.
        
        \param[in] bp The base pair character (both upper and lower
        cases are handled correctly) whose complementary encoding is
        required.

        \return The reverse-complementary ecnoding for the given
        nucleotide.
    */
    inline static char encode2rc(const char bp)
    { return charToIntComp[(int) bp]; }

    /** Obtain complementary nucleotide a given base pair.

        This method can be used to obtain the \b complementary
        nucleotide character for a given base pair (bp).  This method
        essentially translate bases \c A, \c T, \c C, and \c G, (both
        upper and lower case) into \c T, \c A, \c G, and \c C
        respectively.

        \note In favor of speed, this method does not perform any
        special checks on the actual character in bp. It is the
        responsiblity of the caller to ensure that this method is
        invoked with appropriate parameter value.
        
        \param[in] bp The base pair character (both upper and lower
        cases are handled correctly) whose complementary nucleotide is
        required.

        \return The complementary nucleotide corresponding to the
        given base pair.
    */
	inline static char getComp(const char bp) {
        return RCBases[(int) bp];
    }
    
    /** Obtain the reverse-complement for a given word.

        This method can be used to obtain the reverse-complement
        encoding for a given encoded word.  This method essentially
        translates a given encoded word to its reverse-complement
        representation.

        \note In favor of speed, this method does not perform any
        special checks on the word to be translated.  It is the
        responsiblity of the caller to ensure that this method is
        invoked with appropriate parameter value after the \c
        setRevCompTable method is invoked.
        
        \param[in] word The encoded word that must be translated to
        its corresponding reverse complement representation.

        \return The reverse-complement representation for a given
        word.
    */
    inline int encode2rc(const int word) const { return revCompTable[word]; }
    
    /** Set the reverse-complement translation table to be used whe
        next time encode2rc method is called.

        This method must be invoked to set the correct translation
        table to be used by the encode2rc(int) method.  If a
        translation table does not exist in the \c revCompTables, then
        a new reverse-complement table is created by the \c
        addRevCompTable method.
	
        \param[in] wordSize The number of base pairs in the word for
        which a reverse-complement translation table is to be created.
    */
    void setRevCompTable(const int wordSize);

    /** A functor to generate a encoded word (serves as a hash entry).
        
        This functor must be used to generate an encoded word from a
        "normal" (rather than reverse complement) fragment. This
        me thod handles 'n' entries in the EST in the following manner:

        <ul>

        <li>Whenever it encounters an 'n' base pair (that is created
        when ESTs are loaded in the class EST.cpp) it sets the
        ignoreMask to Mask (as Mask template parameter already tells
        us the number of bits we care about and that is used as an
        indicator of number hashes to ignore).</li>

        <li>The ignoreMask is shifted right dropping off the least
        significant bits each time this method is called. If all bits
        of ignoreMask are cleared then the ignore mask is zero and
        ignored.</li>

        <li>If the ignoreMask is non-zero, then the caller is expected
        not to use the hash returned from this method.  as this
        word/hash contains a 'n'. As bits get shifted the value for
        'n' drops off (when the ignoreMask is zero) and the hashes can
        actually be used by the caller.</li>

        </ul>

        \tparam Shift The number of bits by which the encoding for the
		given base pair must be shifted to the left. For example, when
		using a word of length 6 nt, this value would be 10.

        \tparam Mask The mask (with bits set to 1) that must be used
		to retain the signficiant values in the hash. For example,
		when using words of length 6, the Mask would be the binary
		<code>1111 1111 1111</code> or <code>0xfff</code>.

        \param[in] w The current hash value that has been computed
        thusfar.

        \param[in] bp The base pair ('A', 'T', 'C', 'G', or 'n') to be
        encoded by this method.

        \param[in,out] ignoreMask The ignore mask that is used to
        determine if the current hash/word has a 'n' entry in it and
        must be ignored.

        \return The hash for the current word being encoded.

        \note The caller must use the returned value for further
        operation only if the ignoreMask is zero.  Otherwise the
        encoder must be repeatedly called with subsequent bases until
        the ignoreMask is cleared by this method.
    */
    template <const int& Shift, const int& Mask>
    struct NormalEncoder : public std::binary_function<int, char, int> {
        /**
           The newer-version of this method that explicitly tracks and
           sets igoreMask to enable ignoring 'N' bases.
        */
        inline int operator()(const int w, const char bp, int& ignoreMask) const {
            // Setup the ignore mask if base pair is 'n' otherwise
            // drop off the lowest 2 bits.
            ignoreMask = (bp != 'N') ? (ignoreMask >> 2) : Mask;
            // Compute the hash for the word. Encoding for 'n' is 0.
            return ((w >> 2) | (ESTCodec::encode(bp) << Shift)) & Mask;
        }

        /**
           The earlier-version of this method does not track ignoring
           bases.  This method assumes that all 'N' characters have
           been suitably replaced with a random base letter.
        */
        inline int operator()(const int w, const char bp) const {
            // Compute the hash for the word.
            return ((w >> 2) | (ESTCodec::encode(bp) << Shift)) & Mask;
        }
    };

    /** A functor to generate a encoded word (serves as a hash entry).
        
        This functor must be used to generate an encoded word for
        a reverse complement (rather than normal) fragment. This
        method handles 'n' entries in the EST in the following
        manner:
        
        <ul>
        
        <li>Whenever it encounters an 'n' base pair (that is created
        when ESTs are loaded in the class EST.cpp) it sets the
        ignoreMask to Mask (as Mask template parameter already tells
        us the number of bits we care about and that is used as an
        indicator of number hashes to ignore).</li>
        
        <li>The ignoreMask is shifted right dropping off the least
        significant bits each time this method is called. If all bits
        of ignoreMask are cleared then the ignore mask is zero and
        ignored.</li>
      
        <li>If the ignoreMask is non-zero, then the caller is expected
        not to use the hash returned from this method.  as this
        word/hash contains a 'n'. As bits get shifted the value for
        'n' drops off (when the ignoreMask is zero) and the hashes can
        actually be used by the caller.</li>
        
        </ul>
        
        \tparam Shift The number of bits by which the encoding for the
        given base pair must be shifted. Currently this parameter is
        not used but it is present to provide a symmetric interface
        with the NormalEncoder.
            
        \tparam Mask The mask (with bits set to 1) that must be used
        to retain the signficiant values in the hash. For example,
        when using words of length 6, the Mask would be the binary
        <code>1111 1111 1111</code> or <code>0xfff</code>.
        
        \param[in] w The current hash value that has been computed
        thusfar.
        
        \param[in] bp The base pair ('A', 'T', 'C', 'G', or 'n') to be
        encoded by this method.
        
        \param[in,out] ignoreMask The ignore mask that is used to
        determine if the current hash/word has a 'n' entry in it and
        must be ignored.

        \note The caller must use the returned value for further
        operation only if the ignoreMask is zero.  Otherwise the
        encoder must be repeatedly called with subsequent bases until
        the ignoreMask is cleared by this method.
    */
    template <const int& Shift, const int& Mask>
    struct RevCompEncoder : public std::binary_function<int, char, int> {
        /**
           The newer-version of this method that explicitly tracks and
           sets igoreMask to enable ignoring 'N' bases.
        */        
        inline int operator()(const int w, const char bp, int& ignoreMask) const {
            // Setup the ignore mask if base pair is 'n' otherwise
            // drop off the lowest 2 bits.
            ignoreMask = (bp != 'N') ? (ignoreMask >> 2) : Mask;
            // Compute the hash for the word. 'n' becomes 0.
            return ((w << 2) | ESTCodec::encode2rc(bp)) & Mask;
        }
        /**
           The earlier-version of this method does not track ignoring
           bases.  This method assumes that all 'N' characters have
           been suitably replaced with a random base letter.
        */
        inline int operator()(const int w, const char bp) const {
            // Compute the hash for the word.
            return ((w << 2) | ESTCodec::encode2rc(bp)) & Mask;
        }
    };
    
    /** The destructor.
        
        The destructor frees up memory allocated to hold translation
        tables etc.  The destructor is called only once, when the
        process-wide unique instance is destroyed when program
        terminates.
    */
    ~ESTCodec();

protected:
    /** The constructor.

        The constructor is invoked only once when the process-wide
        unique static instance of the ESTCodec is created when the
        process starts.  The constructor initializes the CharToInt
        array that is used to translate base pairs (\c A, \c T, \c C,
        and \c G) into 2-bits codes (\c 00, \c 11, \c 10, and \c 01)
    */
    ESTCodec();

    /** Creates and adds a new reverse-complement translation table
        for the given word size.

        This method is a helper method that is invoked from the \c
        setRevCompTable whenever a new reverse-complement translation
        table is needed.  This method creates a reverse-complement
        table with 4<sup>wordSize</sup> entries.

        \param[in] wordSize The number of base pairs in the word for
        which a reverse-complement translation table is to be created.

        \return This method returns the newly created
        reverse-complement translation table.
    */
    int* addRevCompTable(const int wordSize);
    
private:
    /** A simple array to map characters \c A, \c T, \c C, and \c G to
        \c 0, \c 3, \c 2, and \c 1 respectively.
        
        This is a simple array of 255 entries that are used to convert
        the base pair encoding characters \c A, \c T, \c C, and \c G
        to \c 0, \c 3, \c 2, and \c 1 respectively.  This encoding is
        typically used to compute the hash as defined by various EST
        analysis algorithms.  This array is statically allocated. It
        is initialized in the constructor and is never changed during
        the life time of this class.

        \note This array is statically allocated to enable ready
        access from NormalEncoder and RevCompEncoder functors defined
        in this class.  Hopefully with static arrays, the compilers
        more readily optimize and inline the method calls to enocde
        and encode2rc.
    */
    static char charToInt[];

    /** A simple array to map numbers \c 0, \c 1, \c 2, and \c 3
        to characters \c A, \c C, \c G, and \c T respectively.
        
        This is a simple array of 4 entries that are used to convert
        the encoding for characters \c 0, \c 3, \c 2, and \c 1 to
        nucleotide characters \c A, \c T, \c C, and \c G respectively.
        This encoding is typically used to display characters and
        rebuild sequences.  This array is statically allocated. It is
        initialized in the constructor and is never changed during the
        life time of this class.
    */
    static char intToChar[];
    
    /** A simple array to map characters \c A, \c T, \c C, and \c G to
        complementary encodings \c 3, \c 0, \c 1, and \c 2
        respectively.
        
        This is a simple array of 255 entries that are used to convert
        the base pair encoding characters \c A, \c T, \c C, and \c G
        to \b complementary codes \c 3, \c 0, \c 1, and \c 2
        respectively.  This encoding is typically used to compute the
        hash as defined by various EST analysis algorithms.  This
        array is statically allocated. It is initialized in the
        constructor and is never changed during the life time of this
        class.

        \note This array is statically allocated to enable ready
        access from NormalEncoder and RevCompEncoder functors defined
        in this class.  Hopefully with static arrays, the compilers
        more readily optimize and inline the method calls to enocde
        and encode2rc.
    */
    static char charToIntComp[];

    /** A simple array to map characters \c A, \c T, \c C, and \c G to
        complementary encodings \c T, \c A, \c G, and \c C
        respectively.
        
        This is a simple array of 255 entries that are used to convert
        the base pair encoding characters \c A, \c T, \c C, and \c G
        to \b complementary codes \c T, \c A, \c G, and \c C
        respectively.  This array is statically allocated. It is
        initialized in the constructor and is never changed during the
        life time of this class.
    */
    static char RCBases[];
    
    /** A hash map that holds tables to aid in translating a given
        word to its reverse complement.

        <p>Converting a given encoded word (some fixed \em n number of
        base pairs, with each base pair encoded into 2-bits) to its
        reverse complement (that is, given the encoded sequence for \c
        attcggct it must be converted to the encoded sequence for \c
        agccgaat) needs to be computed as a part of EST analysis
        algorithms and heuristics. In order to enable rapid translation
        pre-populated tabes are used.</p>

        <p>However, the reverse-complement translation tables need to
        have entries corresponding to the size of words to be
        translated. Different algorithms use different word sizes
        (such as: 8 bps or 10 bps etc).  Accordingly, this hash_map is
        used to hold pre-computed reverse-complement translation
        tables. The key in the hash map is the word size. The
        translation tables contained in this hash map are used via the
        \c setRevCompTable method. If a reverse-complement entry does
        not exist, then a new entry is added by the addRevCompTable
        method. </p>
    */
    HashMap<int, int*> revCompTables;

    /** The reverse-complement translation table to be used by the
        encode2rc method.

        This array is set by the \c setRevCompTable method to refer to
        the reverse-complement translation table to translate words of
        given size to their corresponding reverse-complement
        encodings.
    */
    const int* revCompTable;
    
    /** The process-wide unique codec instance.
        
        This instance variable is a process-wide unique codec that is
        created when the process is started and is destroyed only when
        the process terminates.
    */
    static ESTCodec estCodec;
};

#endif
