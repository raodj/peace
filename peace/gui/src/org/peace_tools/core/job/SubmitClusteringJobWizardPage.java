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

package org.peace_tools.core.job;

import org.peace_tools.workspace.DataSet;
import org.peace_tools.workspace.FileEntry;
import org.peace_tools.workspace.GeneratedFileList;
import org.peace_tools.workspace.Job;

/**
 * <p>This class serves as the final page in a JobWizard. This page
 * provides the user with feedback as the job is submitted to be
 * run on the server. The page specifically deals with submitting
 * a clustering job. The base class performs all the core tasks.
 * This class essentially implements helper call back methods to
 * provide the necessary customized information to the base class.</p>
 * 
 * <p>Note that the process of submitting the job can take several
 * seconds (or even a minute if EST file to copied is large). 
 * Consequently the operations are run on a separate thread to ensure
 * that the GUI continues to provide feedback without appearing as if
 * it is hanging (it would appear to be hanging if the operations are
 * performed on the same thread).</p> 
 */
public class SubmitClusteringJobWizardPage extends GenericSubmitJobWizardPage { 
	/**
	 * The constructor. The constructor merely saves the reference
	 * to the wizard that owns this page. The core operations of
	 * setting up the GUI elements is delegated to the base class.
	 * 
	 * @param wizard The wizard that logically owns this page.
	 */
	public SubmitClusteringJobWizardPage(JobWizard wizard) {
		super(wizard, "Job Submission", 
			  "Submitting job for running on a server", true);
		this.wizard = wizard;
	}
	
	@Override
	protected Job getJob() {
		return wizard.createJobEntry();
	}

	@Override
	protected GeneratedFileList getGeneratedFiles() {
		return wizard.createGFL();
	}

	@Override
	protected FileEntry[] getFilesToCopy() {
		// Intentionally returning null as there are no
		// additional files to copy.
		return null;
	}

	@Override
	protected DataSet getDataSet() {
		return wizard.getDataSet();
	}

	@Override
	protected String getEASTWorkDir() {
		return "";
	}
	
	/**
	 * A reference to the wizard dialog that logically owns this
	 * page. This reference is used to enable and disable 
	 * buttons on this wizard appropriately.
	 */
	private final JobWizard wizard;
	
	/**
	 * A serialization UID to keep the compiler happy.
	 */
	private static final long serialVersionUID = 8643291107936550834L;
}
