package org.peace_tools.core.job.east;

import java.awt.BorderLayout;

import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.border.EmptyBorder;

import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.Utilities;

/**
 * This class serves as the first interactive page in
 * a EAST Job wizard. This page is used for the following
 * purpose:
 * 
 * <ol>
 * 
 * <li>It is used by the EASTJobWizard to provide the user with
 * a brief description of the purpose of this wizard and a
 * high level summary of its operations.</li>
 * 
 * <li>It is used to obtain a brief description from the user
 * about this job.</li>
 * 
 * </ol>
 * 
 * The choice of operation in the above two options is decided by 
 * the flag value that is passed in via the constructor.
 */
public class OverviewPage extends GenericWizardPage { 	
	/**
	 * The constructor. The constructor sets up the various components
	 * on this wizard page. The components include: combo box to
	 * select data set, a text area for job description, and a check
	 * box for clustering. 
	 * 
	 * @param wizard The wizard that logically owns this page.
	 * 
	 * @param alsoCluster If this flag is true, then the wizard is expected to
	 * perform clustering and assembly. If the flag is false, then this wizard
	 * is launched just for an assembly job.
	 */
	public OverviewPage(EASTJobWizard wizard, boolean alsoCluster) {
		this.wizard = wizard;
		assert(this.wizard != null);
		// Setup the title(s) for this page and border
		setTitle("Overview", "Overview of tasks in this wizard.");
		
		JLabel message = new JLabel((alsoCluster ? CLUSTER_OVERVIEW_MSG : ASSEMBLY_OVERVIEW_MSG));
		add(message, BorderLayout.NORTH);
		setBorder(new EmptyBorder(20, 15, 10, 5));
		
		// Create panels with description of the job
		description = new JTextArea(3, 3);
		JScrollPane jsp = new JScrollPane(description);
		jsp.setMinimumSize(description.getPreferredSize());
		JComponent descBox = 
			Utilities.createLabeledComponents("Description for job(s):",
					"(This is for your reference & can be anything)", 0, 
					false, jsp);
		add(descBox, BorderLayout.CENTER);
	}

	/**
	 * Obtain the user-entered description for the job(s) being created.
	 * 
	 * @return The user-entered description for the job.
	 */
	protected String getDescription() {
		return description.getText();
	}
	
	/**
	 * Field to read and edit a brief description about the job.
	 * This information can be anything the user desires and is
	 * meaningful only to the user.
	 */
	private JTextArea description;
	
	/**
	 * A reference to the wizard dialog that logically owns this
	 * page. This reference is used to enable and disable 
	 * buttons on this wizard appropriately.
	 */
	private final EASTJobWizard wizard;
	
	/**
	 * A static overview message that is displayed in the first
	 * overview page displayed by this wizard to the user.
	 */
	private static final String ASSEMBLY_OVERVIEW_MSG = "<html>"+
		"This wizard guides you through the process of creating an EAST job for<br/>" +
		"assembling cDNA fragments using an existing Minimum Spanning Tree<br/>" +
		"(MST) and clustering solution.  The assembly job runs as a single<br/>" +
		"process on a given server on which EAST has been installed.  The time<br/>" +
		"taken by EAST to assemble <i>n</i> cDNA fragments is \u03F4" +
		"(<i>n</i>\u00B2).<br/><br/>" +
		"<b>Note</b>: If you would like to cluster (using default parameters)<br/>" +
		"and assemble please use the the other EAST assembly tool/wizard. For<br/>" +
		"clustering with custom parameters use the clustering tool directly." +
		"</html>";
	
	/**
	 * A static overview message that is displayed in the first
	 * overview page displayed by this wizard to the user. This message
	 * is display when the user is launching the wizard to perform
	 * clustering + assembly.
	 */
	private static final String CLUSTER_OVERVIEW_MSG = "<html>"+
		"This wizard guides you through the process of clustering (via PEACE)<br/>" +
		"and then assembly via EAST. Clustering can run in parallel but assembly<br/>" +
		"job runs as a single process on a given server.  The time taken by PEACE<br/>" +
		"cluster and EAST to assemble <i>n</i> cDNA fragments is \u03F4" +
		"(<i>n</i>\u00B2).<br/>" +
		"<b>This wizard does not permit you to customize</b><br/>" +
		"<b>clustering parameters.</b> Use the clustering wizard<br/>" +
		"to perform custom clustering.<br/><br/>" +
		"<b>Note</b>: You can perform assembly via EAST directly if an MST file has<br/>"+
		"already been created through a previous clustering job." +
		"</html>";

	
	/**
	 * A serialization UID to keep the compiler happy.
	 */
	private static final long serialVersionUID = -1041228009499130473L;	
}
