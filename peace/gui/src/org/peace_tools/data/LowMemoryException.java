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

package org.peace_tools.data;


/**
 * A simple exception to report a low memory condition.
 *
 * This exception class provides a more pertinent exception to report
 * the condition where the JVM is low on memory and the user has chosen
 * not to proceed further with a file open operation.
 */
public class LowMemoryException extends Exception {
	/**
	 * The default constructor.
	 * 
	 * The default constructor to create an exception without any message.
	 * Constructs a new exception with null as its detail message. The cause
	 * is not initialized, and may subsequently be initialized by a call to
	 * {@link Throwable#initCause(Throwable)}. 
	 */
	public LowMemoryException() {
		// Nothing to be done here.
	}
	
	/**
	 * Constructor to create an exception with the given message.
	 *
	 * A convenience constructor to create an exception with a specific
	 * mesasge.
	 * 
	 * @message The message to be set for this exception.
	 */	
	public LowMemoryException(String message) {
		super(message);
	}
	
	/**
	 * A generated serialization UID to keep compiler happy.
	 */
	private static final long serialVersionUID = 5033361156599458567L;	
}
