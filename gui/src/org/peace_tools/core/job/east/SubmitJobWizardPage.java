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

package org.peace_tools.core.job.east;

import org.peace_tools.core.job.GenericSubmitJobWizardPage;
import org.peace_tools.workspace.DataSet;
import org.peace_tools.workspace.FileEntry;
import org.peace_tools.workspace.GeneratedFileList;
import org.peace_tools.workspace.Job;

/**
 * <p>This class serves as the final (one or two) page in a 
 * EASTJobWizard. The same class (with slightly different parameters
 * during instantiation) used to submit clustering as well as a
 * assembly job. This page provides the user with feedback as
 * a clustering or assembly job is submitted to be
 * run on the server. The base class performs all the core tasks.
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
public class SubmitJobWizardPage extends GenericSubmitJobWizardPage { 
	/**
	 * The constructor. The constructor merely saves the reference
	 * to the wizard that owns this page. The core operations of
	 * setting up the GUI elements is delegated to the base class.
	 * 
	 * @param wizard The wizard that logically owns this page.
	 * 
	 * @param clusteringJob Flag to indicate if this a clustering or
	 * assembly job submission page. If this flag is true then it is
	 * a clustering job submission page.
	 * 
	 * @param alsoClustering Flag to indicate if the wizard is also
	 * clustering. This flag is used to decide if the job must be
	 * started right away or if an assembly job must wait for a
	 * clustering job to finish first.
	 * 
	 * @param title The title to be set for this wizard page.
	 * 
	 * @param subTitle The subtitle to be set for this wizard page.

	 */
	public SubmitJobWizardPage(EASTJobWizard wizard, boolean clusteringJob, 
			boolean alsoClustering, String title, String subTitle) {
		super(wizard, title, subTitle, (clusteringJob || !alsoClustering)); 
		this.wizard        = wizard;
		this.clusteringJob = clusteringJob;
		this.alsoClustering= alsoClustering;
	}
	
	@Override
	protected Job getJob() {
		// If we are also clustering then to avoid user confusion in job IDs
		// and order in which entries are added to the workspace, we create
		// the clustering job entry first
		if (alsoClustering) {
			wizard.createJobEntry(false, true);
		}
		return wizard.createJobEntry(!clusteringJob, true);
	}

	@Override
	protected GeneratedFileList getGeneratedFiles() {
		return wizard.createGFL(!clusteringJob, true, getJob());
	}

	@Override
	protected FileEntry[] getFilesToCopy() {
		if (alsoClustering) {
			// Intentionally returning null as there are no
			// additional files to copy.
			return null;
		}
		// For direct assembly the MST file needs to be copied
		// to the server (the source EST file is handled separately
		// already).
		FileEntry[] filesToCopy = new FileEntry[1];
		filesToCopy[0] = wizard.getMSTFile();
		return filesToCopy;
	}

	@Override
	protected DataSet getDataSet() {
		return wizard.getDataSet();
	}

	@Override
	protected String getEASTWorkDir() {
		if (clusteringJob) {
			return wizard.getEASTWorkDir();
		}
		return "";
	}
	
	/**
	 * Call-back method to indicate that the thread for creating/submitting
	 * a job is done.
	 *
	 * This method overrides the default implementation in the base class 
	 * to change to the next wizard page, if this is another page after
	 * this one.
	 */
	@Override
	protected void done() {
		final int currPage = wizard.getCurrentPage();
		if (wizard.getPage(currPage + 1) != null) {
			// Force the wizard to skip to the next wizard page.
			wizard.changePage(currPage, currPage + 1);
		}
	}
	
	/**
	 * A reference to the wizard dialog that logically owns this
	 * page. This reference is used to enable and disable 
	 * buttons on this wizard appropriately.
	 */
	private final EASTJobWizard wizard;

	/**
	 * Flag to indicate if this page is meant to deal with a 
	 * clustering or an assembly job. If the flag is true then
	 * this instance is dealing with a clustering job.
	 */
	private final boolean clusteringJob;

	/**
	 * Flag to indicate if this page is working in a mode
	 * where the wizard has been launched for clustering+assembly.
	 * If this flag is false, then the wizard is meant to
	 * perform just assembly.
	 */
	private final boolean alsoClustering;

	/**
	 * A serialization UID to keep the compiler happy.
	 */
	private static final long serialVersionUID = 640853143689284383L;
}
