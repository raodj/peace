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
 * Provides all the classes and functionality associated with the
 * Synthetic Dataset Generator (SDG).
 *   
 * Distributed Environment for Comparative Analysis of Genomic-Assemblers 
 * (DECAGON) is a framework designed to ease comparative analysis of 
 * assemblers. This package contains the sub-set of DECAGON classes that 
 * constitute the Synthetic DataSet Generator (SDG). The SDG has been
 * implemented as a wizard in PEACE-GUI. The wizard guides the user through 
 * the following setups in SDG to generate a synthetic data set:
 *  
 * <ol>
 * 
 * <li>The user selects the source data set that contains the genes from
 * where the synthetic data is to be created. This page also permits
 * the user to enter a description to be associated with the generated
 * data set</li>
 * 
 * <li>If the data set is to include Sanger sequences, then the user can
 * setup the parameters corresponding to this type of nucleotide sequences
 * via a suitable wizard page.</li>
 * 
 * <li>If the data set is to include 454 sequences, then the user can
 * setup the parameters corresponding to this type of nucleotide sequences
 * in a given wizard page.</li>
 * 
 * <li>If the data set is to include illumina sequences, then the user can
 * setup the parameters corresponding to this type of nucleotide sequences
 * in a given wizard page.</li>
 * 
 * </ol>
 * 
 * @since 0.97
 */
package org.peace_tools.decagon.sdg;
