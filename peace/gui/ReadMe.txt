      ___         ___           ___           ___           ___     
     /  /\       /  /\         /  /\         /  /\         /  /\    
    /  /::\     /  /:/_       /  /::\       /  /:/        /  /:/_   
   /  /:/\:\   /  /:/ /\     /  /:/\:\     /  /:/        /  /:/ /\  
  /  /:/~/:/  /  /:/ /:/_   /  /:/~/::\   /  /:/  ___   /  /:/ /:/_ 
 /__/:/ /:/  /__/:/ /:/ /\ /__/:/ /:/\:\ /__/:/  /  /\ /__/:/ /:/ /\
 \  \:\/:/   \  \:\/:/ /:/ \  \:\/:/__\/ \  \:\ /  /:/ \  \:\/:/ /:/
  \  \::/     \  \::/ /:/   \  \::/       \  \:\  /:/   \  \::/ /:/ 
   \  \:\      \  \:\/:/     \  \:\        \  \:\/:/     \  \:\/:/  
    \  \:\      \  \::/       \  \:\        \  \::/       \  \::/   
     \__\/       \__\/         \__\/         \__\/         \__\/    



          PEACE: Parallel EST Analyzer and Clustering Engine
             Copyright (c) Miami University,  Oxford, OHIO.
			 
              GPLed by Dhananjai M. Rao (raodm@muohio.edu)

------------------------------------------------------------------------

 This file is part of PEACE.
 
 PEACE is free software: you can redistribute it and/or modify it
 under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 PEACE is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with PEACE.  If not, see <http://www.gnu.org/licenses/>.
 
 Miami University makes no representations or warranties about the
 suitability of the software, either express or implied, including
 but not limited to the implied warranties of merchantability,
 fitness for a particular purpose, or non-infringement.  Miami
 University shall not be liable for any damages suffered by licensee
 as a result of using, result of using, modifying or distributing
 this software or its derivatives.

 By using or copying this Software, Licensee agrees to abide by the
 intellectual property laws, and all other applicable laws of the
 U.S., and the terms of GNU General Public License (version 3).

 Authors:   Dhananjai M. Rao              raodm@muohio.edu

------------------------------------------------------------------------

Package Organization:
---------------------

All the packages in the GUI subsystem for PEACE have been organized in
the following manner:

	- org.peace_tools: This is the top-level package.  This choice of
      package is meant to reflect the standard package naming
      convention recommended for Java packages. The choice mirrors the
      domain name for PEACE, namely: peace-tools.org.

	- org.peace_tools.workspace: This package contains the workspace
	  class hierarchy that maintains an in-memory representation of
	  the data loaded from a workspace XML file. In addition to
	  providing convenient access to the XML data, the files also help
	  in marshalling the in-memory data to XML format for persistence.

	- org.peace_tools.data: This package contains the in-memory
      classes that encapsulate the data associated with data files
      (EST files, MST files, Cluster files), jobs, and servers. This
      package essentially contains the "model" components -- that is,
      "model" as in the Model-View-Controller (MVC) pattern.

	- org.pace_tools.views: This package contains the GUI components
      that are specifically tailored to display the data in a suitable
      graphical format. These components constitute the "view" as in
      the Model-View-Controller (MVC) pattern. These views display the
      data stored in the "model" files in the org.peace_tools.data
      package.

	- org.peace_tools.core: This package contains the core files of
      the GUI. These files are essentially the controllers (as in the
      Model-View-Controller pattern) and perform other activities
      associated with managing the data and views in a workspace.
	
	- org.peace_tools.generic: This package contains generic/common
      files that are shared utility classes that are used by two or
      more classes from other packages. In addition, these files are
      "generic" in the sense that they are not really tied to PEACE
      GUI but can be readily used to build other GUIs as well.
