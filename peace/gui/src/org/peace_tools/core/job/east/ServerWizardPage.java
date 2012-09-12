package org.peace_tools.core.job.east;

import java.awt.BorderLayout;
import java.awt.Color;

import javax.swing.BorderFactory;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.border.EmptyBorder;

import org.peace_tools.core.job.ServerPanel;
import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.Utilities;
import org.peace_tools.generic.WizardDialog;

/**
 * <p>A wizard page that permits the user to select the server
 * on which a clustering or assembly job is to be run. In
 * addition to selecting the server, this page also permits
 * the user to setup additional properties for a job,
 * such as: number of compute nodes, CPUs/node, max memory,
 * and estimated runtime.</p>
 * 
 * <p>This wizard page may be used twice in the EAST job wizard. 
 * The first optional instance is used for specifying the
 * configuration for a clustering job. The second instance is
 * used for acquiring the assembly job configuration. Note that
 * if two instances are present, then this class provides
 * the API to tie the servers used by the two instances. That is,
 * when the user gets to select the server to be used for both
 * jobs on the first wizard page and the second wizard page
 * uses the same server (but permits configuration of other
 * job parameters).</p>
 */
public class ServerWizardPage extends GenericWizardPage {
	/**
	 * The constructor. The constructor sets up the various components
	 * on this wizard page. The components are actually created via
	 * call to {@link ServerPanel#createServerPanel(boolean, boolean)}
	 * method.
	 * 
	 * @param wizard The wizard that logically owns this page.
	 * 
	 * @param cluster If this flag is true that indicates that the
	 * server page must be configured and not for assembly. On the
	 * other hand, if this flag is false, then this page is configured
	 * for assembly.
	 * 
	 * @param mainSrvrPage An server wizard page that serves as the
	 * source of the server to be used for this page as well. If this
	 * parameter is not null, then the server option is locked based
	 * on the entry selected in this server page. If the parameter is
	 * null, then the user is permitted to select any  valid server
	 * from the set. 
	 */
	public ServerWizardPage(EASTJobWizard wizard, boolean clustering, 
			                ServerWizardPage mainSrvrPage) {
		this.wizard       = wizard;
		this.serverInfo   = new ServerPanel();
		this.mainSrvrPage = mainSrvrPage;
		
		// Setup the title(s) for this page and border
		if (clustering) {
			setTitle("Clustering Server", "Select server configuration for clustering");
		} else {
			setTitle("Assembly Server", "Select server configuration for assembly");
		}
		setBorder(new EmptyBorder(5, 5, 5, 5));
		
		// Setup a label to inform the user about the intent for this page
		JLabel infoLabel = new JLabel(clustering ? CLUSTERING_INFO_MSG : ASSEMBLY_INFO_MSG);
		infoLabel.setIcon(Utilities.getIcon("images/32x32/Note.png"));
		infoLabel.setBackground(new Color(0xe0, 0xe0, 0xff)); // pale blue
		infoLabel.setOpaque(true);
		infoLabel.setBorder(BorderFactory.createEtchedBorder());
		// Create the generic server panel to select a common server
		// and depending on the setup a configuration clustering or
		// just assembly.
		JPanel srvrConfigPanel = serverInfo.createServerPanel(clustering, false, null);
		// Layout the components in a suitable sub-panel.
		JPanel subPanel = new JPanel(new BorderLayout(0, 5));
		subPanel.add(infoLabel, BorderLayout.NORTH);
		subPanel.add(srvrConfigPanel, BorderLayout.CENTER);
		// Add the sub panel with all the GUI components to the wizard page.
		add(subPanel, BorderLayout.CENTER);
	}
	
	/**
	 * This method is called just before this page is to be displayed.
	 * This page essentially updates the list of server entries
	 * displayed in the combo-box. It also suggests a default MST file
	 * name to the user. However, the user can edit the default suggestion
	 * and set it to an appropriate value. In addition, it updates the memory
	 * requirement for running the job. 
	 */
	@Override
	public void pageChanged(WizardDialog dialog, int currPage, int prevPage) {
		if (mainSrvrPage != null) {
			// We need to narrow down to a single server option.
			serverInfo.lockSelectedServer(mainSrvrPage.getServerInfoPanel().getSelectedServer());
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
	 * A reference to the wizard dialog that logically owns this
	 * page. This reference is used to enable and disable 
	 * buttons on this wizard appropriately.
	 */
	private final EASTJobWizard wizard;
	
	/**
	 * A helper panel class that provides the necessary GUI features
	 * to interactively obtain server choice and job configuration
	 * information from the user.
	 */
	private final ServerPanel serverInfo;
	
	/**
	 * This object is initialized in the constructor (via a parameter
	 * from EASTJobWizard) to point to the first page in the wizard
	 * in which the common server configuration was setup by the
	 * user. If this value is not null, then this page sets its
	 * server options (presented by this object) down to the one 
	 * selected by the user in the {@link #mainSrvrPage}. 
	 */
	private final ServerWizardPage mainSrvrPage;
	
	/**
	 * A generic informational message that is displayed at the
	 * top of this wizard page to provide some additional information
	 * to the user when both clustering and assembly is requested.
	 */
	private static final String CLUSTERING_INFO_MSG = 
		"<html><b>The same server is used both for clustering &amp; assembly.</b><br/>" +
		"Configure parameters for parallel clustering (via PEACE) job.<br/>" +
		"Default clustering parameters will be used" +
		"</html>";

	/**
	 * A generated serialization GUID to keep the compiler happy. 
	 */
	private static final long serialVersionUID = -687750884685828099L;
	
	/**
	 * A generic informational message that is displayed at the
	 * top of this wizard page to provide some additional information
	 * to the user when only assembly is requested.
	 */
	private static final String ASSEMBLY_INFO_MSG = 
		"<html><b>The same server is used both for clustering &amp; assembly.</b><br/>" +
		"Configure parameters for serial assembly (via EAST) job.<br/>" +
		"Additional parameters for assembler are configured next." +
		"</html>";
}
