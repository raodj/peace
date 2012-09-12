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
 * Provides all the classes and functionality associated with integrating 
 * DECAGON into PEACE GUI.
 *   
 * Distributed Environment for Comparative Analysis of Genomic-Assemblers 
 * (DECAGON) is a framework designed to ease comparative analysis of 
 * assemblers. Specifically it provides the following main features:
 * 
 * <ul>
 * 
 * <li>Generate synthetic reads form a given set of genes using 
 * Metasim (http://ab.inf.uni-tuebingen.de/software/metasim/) a 
 * closed-source synthetic read generator. Metasim can generate
 * synthetic reads resembling various sequencing technologies 
 * such as: Sanger, Solid/454, Solexa/Illumina. Note that Metasim
 * is not distributed with PEACE because it is a closed source, non-GPL
 * software package that the user has to separately download and
 * install.</li>
 * 
 * <li>DECAGON can manage and operate off a shared storage space to
 * ease sharing of data sets for comparative analysis with various
 * team members.</li>
 * 
 * <li>DECAGON can store results from comparative analysis in a 
 * MySQL database. The database eases further analysis of the data
 * and eases publication of the data.</li>
 * 
 * </ul>
 * 
 * @since 0.95
 */
package org.peace_tools.decagon;
