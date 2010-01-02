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
    */
    inline static char encode(const char bp)
    { return charToInt[(int) bp]; }

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
    */
    inline static char encode2rc(const char bp)
    { return charToIntComp[(int) bp]; }


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
	
     */
    template <const int& Shift, const int& Mask>
    struct NormalEncoder : public std::binary_function<int, char, int> {
        inline int operator()(const int w, const char bp) const {
            return ((w >> 2) | (ESTCodec::encode(bp) << Shift)) & Mask;
        }
    };
    
    template <const int& Shift, const int& Mask>
    struct RevCompEncoder : public std::binary_function<int, char, int> {
        inline int operator()(const int w, const char bp) const {
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
