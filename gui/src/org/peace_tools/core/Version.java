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

package org.peace_tools.core;

/**
 * This is a non-instantiable class that serves as a central point to
 * hold various version, copyright, and lisence information associated
 * with PEACE and its GUI.  
 */
public class Version {
	/** The current version information of the software.

	   This constant provides a convenient way to track and log the version
	   information of the software for tracking version changes.
	*/
	public static final String GUI_VERSION = "Version 0.95\n(Feb 19 2010)";
	
	/** A simple copyright message.

	   This constant provides a convenient way to display a consistent
	   copyright message in the software.
	*/
	public static final String COPYRIGHT  = 
		"Copyright (C) Miami University, 2009-";

	/** The standard disclaimer message.

	   This define aims to centralize the disclaimer for the software and
	   keeps the code clutter to a minimum in other source files.
	*/
	public static final String DISCLAIMER = 
	    "Miami University (MU), Oxford, OHIO makes no representations\n"    +
	    "or warranties about the suitability of the software, either \n"    +
	    "express or implied, including but not limited to the implied\n"    +
	    "warranties of merchantability, fitness for a particular purpose\n" + 
	    "or non-infringement.  MU shall not be liable for any damages\n"    +
	    "suffered by licensee as a result of using, result of using,\n"     +
	    "modifying or distributing this software or its derivatives.\n"     +
	    "\n\n"                                                              +
	    "By using or copying this Software, Licensee agrees to abide by\n"  +
	    "the intellectual property laws, and all other applicable laws\n"   +
	    "of the U.S., and the terms of this license.";

	/**
	 * Since this class is not meant to be instantiated, its constructor
	 * has been made private.
	 */
	private Version() {}
}
