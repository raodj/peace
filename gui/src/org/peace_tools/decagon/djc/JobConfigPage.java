package org.peace_tools.decagon.djc;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;

import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.border.EmptyBorder;

import org.peace_tools.core.job.ServerPanel;
import org.peace_tools.decagon.jaxb.Job;
import org.peace_tools.decagon.jaxb.JobKind;
import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.Utilities;
import org.peace_tools.generic.WizardDialog;

/**
 * A wizard page that permits the user to select the server
 * on which a job is to be run. 
 * 
 * In addition to selecting the server, this page also permits
 * the user to setup additional properties for a job,
 * such as: number of compute nodes, CPUs/node, max memory,
 * and estimated runtime.</p>
 * 
 * <p>This wizard page is used by DECAGON's DecagonJobCreator
 * wizard to permit the user to configure properties for one 
 * or more jobs. In DECAGON, the user selects the server for
 * the first job and subsequent jobs run on the same server 
 * but can have different configurations for: number of 
 * compute nodes, CPUs/node, max memory, and estimated 
 * runtime.</p>
 */
public class JobConfigPage extends GenericWizardPage {
	/**
	 * Constructor to configure this page to select a server for 
	 * a given job. 
	 * 
	 * This constructor is primarily used by the DecagonJobCreator 
	 * to setup a wizard page for permitting the user to select a server
	 * for clustering and/or assembly. The constructor sets up the 
	 * various components on this wizard page. The components are 
	 * actually created via
	 * call to {@link ServerPanel#createServerPanel(boolean, boolean)}
	 * method.
	 * 
	 * @param wizard The wizard that logically owns this page. This
	 * parameter cannot be null.
	 * 
	 * @param job The DECAGON job object from where the information
	 * needed to appropriately initialize various Job parameters. This
	 * parameter cannot be null.
	 * 
	 * @param ssp The server selection page from where the currently
	 * selected server entry is to be obtained to display as the
	 * chosen server to be displayed. 
	 */
	public JobConfigPage(DecagonJobCreator wizard, Job job,
			ServerSelectionPage ssp) {
		this.wizard       = wizard;
		this.srvrSelPage  = ssp;
		this.job          = job;
		// The server panel that contains reused GUI widgets
		this.serverInfo   = new ServerPanel();
		
		// Setup the title(s) for this page and border
		setTitle(job.getName(), ((ssp == null) ?
				"Select server and configure job parameters" :
				"Configure job parameters (server is fixed)"));
		setBorder(new EmptyBorder(5, 5, 5, 5));
		
		// Setup a label to inform the user about the job for which
		// server information is being gathered.
		JPanel infoPanel = new JPanel(new BorderLayout(0, 5));
		infoPanel.add(new JLabel("Job details:", 
				Utilities.getIcon("images/16x16/Information.png"), JLabel.LEFT), BorderLayout.NORTH);
		Color bgColor = wizard.getMainFrame().getBackground();
		bgColor = new Color(bgColor.getRed(), bgColor.getGreen(), bgColor.getBlue() + 0xf);
		JComponent htmlInfo = wizard.getMainFrame().createHTMLComponent(job.getDescription(), bgColor);
		htmlInfo.setPreferredSize(new Dimension(200, 50));
		htmlInfo.setBorder(null);
		infoPanel.add(htmlInfo);
		
		// Create the generic server panel to select a common server.
		// If the job is the first job then choice of selecting a 
		// specific server entry is enabled/disabled.
		final boolean isParallelJob = JobKind.PARALLEL.equals(job.getType());
		JPanel srvrConfigPanel = serverInfo.createServerPanel(isParallelJob, false, null);
		// Layout the components in a suitable sub-panel.
		JPanel subPanel = new JPanel(new BorderLayout(0, 5));
		subPanel.add(infoPanel, BorderLayout.NORTH);
		subPanel.add(srvrConfigPanel, BorderLayout.CENTER);
		// Add the sub panel with all the GUI components to the wizard page.
		add(subPanel, BorderLayout.CENTER);
	}
	
	/**
	 * This method is called just before this page is to be displayed.
	 * 
	 * This page essentially updates the list of server entries
	 * displayed in the combo-box. If the {@link #srvrSelPage} is
	 * not null (which is typically the case) then this method 
	 * locks the server that is to be used for this job. 
	 * 
	 * @param dialog The wizard dialog that logically owns this page.
	 * Currently this parameter is not used.
	 * 
	 * @param currPage The zero-based logical sequence index of this
	 * wizard page in the overall wizard. Currently this parameter
	 * is not used.
	 *  
	 * @param prevPage The zero-based logical sequence index of the 
	 * previous page from where the user was directed to this page.
	 * Currently this parameter is not used by this method.
	 */
	@Override
	public void pageChanged(WizardDialog dialog, int currPage, int prevPage) {
		if (srvrSelPage != null) {
			// We need to narrow down to a single server option.
			serverInfo.lockSelectedServer(srvrSelPage.getSelectedServer());
		} else {
			// Let the server info update server list and associated information
			serverInfo.updateServerList(true, true, wizard.getDataSet());
		}
	}
	
	/**
	 * Obtain the panel that contains the GUI components and data 
	 * regarding the server and its configuration.
	 * 
	 * @return The server panel associated with this wizard page.
	 */
	protected ServerPanel getServerInfoPanel() {
		return this.serverInfo;
	}
	
	/**
	 * Obtain the DECAGON job associated with this configuration page.
	 * 
	 * This method essentially returns the DECAGON job for which this
	 * page is currently configured to obtain configuration 
	 * information. This method returns the job set when this
	 * object was instantiated.
	 * 
	 * @return The DECAGON job associated with this configuration
	 * page.
	 */
	protected Job getJob() {
		return job;
	}
	
	/**
	 * A reference to the wizard dialog that logically owns this
	 * page. This reference is used to enable and disable 
	 * buttons on this wizard appropriately.
	 */
	private final DecagonJobCreator wizard;
	
	/**
	 * A helper panel class that provides the necessary GUI features
	 * to interactively obtain server choice and job configuration
	 * information from the user.
	 */
	private final ServerPanel serverInfo;
	
	/**
	 * The DECAGON Job object for which server information is to be
	 * obtained. This information is used to appropriately setup 
	 * the GUI components on this wizard page and to update various
	 * values based on scaling factors supplied in the Job
	 * configuration. 
	 */
	private final Job job;
	
	/**
	 * This object is initialized in the constructor to point to the 
	 * page in the wizard in which the server to be used for running
	 * the pipeline was setup by the user. This object is never
	 * null.
	 */
	private final ServerSelectionPage srvrSelPage;

	/**
	 * A generated serialization GUID to keep the compiler happy. 
	 */
	private static final long serialVersionUID = -687750884685828099L;
}
