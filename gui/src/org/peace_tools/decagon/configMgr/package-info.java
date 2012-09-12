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
// Authors:   Dhananjai M. Rao              raodm@muohio.edu
//
//---------------------------------------------------------------------

/**
 * Provides all the classes and functionality associated with setting and
 * testing properties of DECAGON. 
 *   
 * Distributed Environment for Comparative Analysis of Genomic-Assemblers 
 * (DECAGON) is a framework designed to ease comparative analysis of 
 * assemblers. This package contains the set of classes that ease setup/
 * configuration of DECAGON. Specifically this package contains a 
 * wizard that guides the user through the following setups/configuration of 
 * DECAGON:
 *  
 * <ol>
 * 
 * <li>The path of MultiSim is set and validated.</li>
 * 
 * <li>Location of a shared space for sharing source transcripts, 
 * real and synthetic data sets, and results from various genomic 
 * assemblers.</li>
 * 
 * <li>Information for a DECAGON data base that is used to organize
 * and store results to ease further analysis and reporting</li>
 * 
 * </ol>
 * 
 * @since 0.96
 */
package org.peace_tools.decagon.configMgr;
