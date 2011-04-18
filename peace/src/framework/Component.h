#ifndef COMPONENT_H
#define COMPONENT_H

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

#include "ArgParser.h"

// Forward declarations to keep compiler fast and happy
class ArgParser;
class SubSystem;
class RuntimeContext;

/** A common base-class for a Component that constitutes a sub-system.

    This class serves as a common parent class for a component. In
    PEACE terminology, a <i>component</i> is a logical sub-unit of a
    sub-system.  In other words, a sub-system consists of a collection
    of components (and possibly sub-components).  This class provides
    some of the general and common methods that are typically
    supported by all components.  Furthermore, this class provides
    default implementations for the methods, thereby easing the
    development of derived components without deviating from the
    generic API.

    \note This class is not meant to be directly instantiated.
    Instead derived concrete component objects are created by various
    sub-systems once the initial set of arguments have been
    successfully parsed.
    
*/
class Component {
public:
    /** The destructor.

        The destructor is merely present to adhere to coding
        conventions.  Currently, this does not have any special tasks
        to perform.
    */
    virtual ~Component();

    /** Add the set of command line parameters for this component.
        
        This method is invoked by the SubSystem (that logically owns
        the component) when it is being requested for command line
        arguments.  This method is expected to add the various command
        line arguments that can be used to further customize the
        operations of this component.

        \note The default implementation does not add any arguments.
        
        \param[out] argParser The argument parser to which the command
        line arguments for this component are to be added.
    */
    virtual void addCommandLineArguments(ArgParser& argParser);

    /** Determine if initialization has been completed.

        This method must be used to detect if initialization of this
         object has been completed.  If initialization has been
         completed then further initialization is not necessary and
         can be avoided.  However, appropriate behavior and semantics
         of initialization is divested to the derived classes.
		 
		 \return This method returns \c true if initialization has been
		 completed.  Otherwise this method returns \c false.
    */
    inline bool isInitialized() const { return initializedFlag; }
	
    /** Initialize the operations of this component.

        This method is invoked when the SubSystem (that logically owns
        this component) is being initialized.  This method (along with
        helper methods) is expected to performs the necessary
        initialization operation required by this component.

        \note The default implementation in this class sets the
        Component::intializedFlag to \c true and returns \c true.
        
        \return This method returns \c true if initialization was
        successfully completed.  On errors it returns \c false.
    */
    virtual bool initialize();

    /** Invoked to permit the component to perform its core tasks.

        This method is invoked once all the sub-systems and components
        have been successfully initialized.  This method is expected
        to perform the core task for this component.

        \note The default implementation in the base class performs no
        operations and always returns zero.
        
        \return This method returns zero on success.  If errors occur
        during initialization then this method must return an non-zero
        error code.
    */
    virtual int run();
    
    /** Wind-up the operations of this component.

        This method is invoked when the SubSystem (that logically owns
        this component) is being finalized.  This method is expected
        to perform any cleanup operations (converse of operations
        performed in the initialize() method).
    */
    virtual void finalize() {}

    /** Get the name associated with this component.

        This method can be used to obtain the name associated with
        this component.  The name is the value that was set when this
        component was instantiated.

        \return The name associated with the string. There are no
        guarantees that this string is unique.
    */
    inline const std::string getName() const { return name; }

    /** Set the sub-system that logically contains this component.

        This method is invoked when this component is logically added
        to a sub-system.

        \param[in] subSystem The sub-system that logically owns this
        component.
    */
    void setSubSystem(SubSystem* subSystem);

    /** Obtain the runtime context associated with this component.

        This method provides a convenient mechanism to obtain the
        RuntimeContext object associated with this component.  The
        context is indirectly obtained from the SubSystem that
        logically contains this component.

        \return The RuntimeContext that contains references to the
        shared information in the current run of PEACE.
    */
	RuntimeContext* getContext() const;

	/** Set/change value for a string argument.

		This method can be used to set/change the value of a string
		argument (in the list of valid arguments) recognized by this
		argument parser.

		\param[in] arg The argument whose value is to be changed.  The
		argument string must exactly match one of the recognized
		values and the data type of the argument must be
		ArgParser::STRING.

		\param[in] value The new value to set for the specified
		argument (assuming it is valid).

		\return This method returns \c true if the specified argument
		(arg) was valid and its value was successfully updated.
		Otherwise this method returns \c false.
	*/
	bool setArgument(const std::string& arg, const std::string& value);

	/** Set/change value for a boolean argument.

		This method can be used to set/change the value of a boolean
		argument (in the list of valid arguments) recognized by this
		argument parser.

		\param[in] arg The argument whose value is to be changed.  The
		argument string must exactly match one of the recognized
		values and the data type of the argument must be a
		ArgParser::BOOLEAN.

		\param[in] value The new value to set for the specified
		argument (assuming it is valid).

		\return This method returns \c true if the specified argument
		(arg) was valid and its value was successfully updated.
		Otherwise this method returns \c false.
	*/
	bool setArgument(const std::string& arg, const bool value);

	/** Set/change value for an integer valued argument.

		This method can be used to set/change the value of an integer
		argument (in the list of valid arguments) recognized by this
		argument parser.

		\param[in] arg The argument whose value is to be changed.  The
		argument string must exactly match one of the recognized
		values and the data type of the argument must be a
		ArgParser::INTEGER.

		\param[in] value The new value to set for the specified
		argument (assuming it is valid).

		\return This method returns \c true if the specified argument
		(arg) was valid and its value was successfully updated.
		Otherwise this method returns \c false.
	*/
	bool setArgument(const std::string& arg, const int value);

	/** Set/change value for a floating valued argument.

		This method can be used to set/change the value of an integer
		argument (in the list of valid arguments) recognized by this
		argument parser.

		\param[in] arg The argument whose value is to be changed.  The
		argument string must exactly match one of the recognized
		values and the data type of the argument must be a
		ArgParser::FLOAT.

		\param[in] value The new value to set for the specified
		argument (assuming it is valid).

		\return This method returns \c true if the specified argument
		(arg) was valid and its value was successfully updated.
		Otherwise this method returns \c false.
	*/
	bool setArgument(const std::string& arg, const float value);

	/** Set/change value for a double valued argument.

		This method can be used to set/change the value of an double
		argument (in the list of valid arguments) recognized by this
		argument parser.

		\param[in] arg The argument whose value is to be changed.  The
		argument string must exactly match one of the recognized
		values and the data type of the argument must be a
		ArgParser::DOUBLE.

		\param[in] value The new value to set for the specified
		argument (assuming it is valid).

		\return This method returns \c true if the specified argument
		(arg) was valid and its value was successfully updated.
		Otherwise this method returns \c false.
	*/
	bool setArgument(const std::string& arg, const double value);

    /** Display the valid arguments along with brief description.

        This method can be used to display the set of valid arguments
        supported by this component (along with a brief description).

        \param[out] os The output stream to which the valid arguments
        are to be written.
    */
    void showArguments(std::ostream& os = std::cout);
    
protected: 
    /** The constructor.

        This is the only constructor for this class. The
        constructor(s) are not public to ensure that this class is not
        instantiated directly.  Instead one of the derived classes
        need to be instantiated.  Note that this class is a very
        generic interface-type class and does not provide any
        meanigful functionality.  Currently the constructor does not
        have any special tasks to perform and merely initializes
        instance variables to their default initial value.

		\param[in] name A generic name to be associated with this
		component.  This is just a debugging convenience than anything
		else.
    */
    Component(const std::string& name);

	/** Helper method to locate and change an argument defined by a
		derived class.

		This is a helper method that is used by the various forms of
		setArgument() method to update a argument.  This method first
		gathers the data in a temporary argument parser and then uses
		the corresponding method in the ArgParser to update the
		parameter values.

        \note This method is templatized so that it can handle the
		various value parameters elegantly.

        \param[in] arg The argument whose value is to be changed.  The
		argument string must exactly match one of the command-line
		arguments supported by the derived class.

        \param[in] argType The expected data type of the argument.

        \param[in] value The new value to be associated with the
        argument.

        \return This method returns \c true if the specified argument
		(arg) was valid and its value was successfully updated.
		Otherwise this method returns \c false.
	*/
	template <typename ValueType>
    bool setArgument(const std::string& arg, const ArgParser::ArgType argType,
                     const ValueType& value);

    /** Helper method to compute the start and ending indexes of the
        EST that this process owns.

        This method was introduced to keep the math and logic clutter
        involved in computing the list of owned ESTs out of the
        methods that use the information.  This method returns the
        range, such that: \c startIndex <= \em ownedESTidx < \c
        endIndex.
        
        \note This method must be invoked only after MPI::Intialize()
        has beeen called and the ESTs to be processed have be loaded
        to get a correct response.

        \param[in] estCount The number of cDNA fragments to be evenly
        divided amongst the various MPI parallel processes.
        
        \param[out] startIndex The starting (zero-based) index value
        of the contiguous range of ESTs that this process owns.

        \param[out] endIndex The ending (zero-based) index value of
        the contiguous range ESTs that this process owns.  The value
        returned in this parameter is \b not included in the range of
        values.
    */
    virtual void getOwnedESTidx(const int estCount, int& startIndex,
                                int& endIndex);
	
    /** Reference to the sub-system that owns this component.

        This instance variable holds a pointer to the sub-system that
        logically owns this component.  This value is set when the
        component is logically associated with a sub-system via the
        setSubSystem() API method.
    */
    SubSystem *subSystem;
    
private:
    /** A generic string to identify the component.

        This is a convenience string that is initialized to contain a
        string that can be used to identify this component. This is
        just a convenience.
    */
    const std::string name;

    /** Flag to detect if initialization has been completed.

        This flag is initialized to \c false in the constructor.  It
        is set to \c true in the initialize() method to indicate
        initialization has been completed.
    */
    bool initializedFlag;	
};

#endif
