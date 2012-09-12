package org.peace_tools.decagon.djc;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.border.EmptyBorder;

import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.Utilities;
import org.peace_tools.generic.WizardDialog;
import org.peace_tools.workspace.Workspace;

/**
 * A wizard page to obtains various additional information from
 * the user.
 * 
 * This class serves as a wizard page in the DECAGON Job Creator.
 * This wizard page obtains some generic inputs regarding the
 * job from the user. The information gathered from the user
 * on this wizard-page includes:
 * 
 * <ul>
 * <li>The directory on the local machine where generated files
 * are to be stored.</li>
 * 
 * <li>A brief description about the job being run.</li>
 * 
 * </ul>
 * 
 */
public class AdditionalInfoPage extends GenericWizardPage implements ActionListener {
	/**
	 * A reference to the wizard dialog that logically owns this
	 * page. This reference is used to obtain information about
	 * the DADX file and the assembler pipeline.
	 */
	private final DecagonJobCreator wizard;
	
	/**
	 * The constructor.
	 * 
	 * The constructor initializes the wizard page and creates various
	 * input GUI components for 
	 * user are created later on (just before the page is actually
	 * displayed in the wizard) by the {@link #pageChanged(WizardDialog, int, int)}
	 * method.
	 * 
	 * @param wizard The DECAGON job creation wizard that logically owns
	 * this page.
	 */
	AdditionalInfoPage(DecagonJobCreator wizard) {
		// Save DADX helper object for later use
		this.wizard = wizard;
		// Set the title for this page.
		setTitle("Additional Inputs", "Provide additional information");
		setBorder(new EmptyBorder(5, 5, 5, 5));
		// Create the GUI input components
		storageDir   = new JTextField();
		JButton browse = Utilities.createButton("images/16x16/Browse.png",
				"Browse", "browse", this, 
				"Click to select storage directory interactively", true);
		// Place the storage text field and browse button horizontally
		// Create top-level container to hold the path and browse buttons.
		JPanel container = new JPanel();
		container.setLayout(new BoxLayout(container, BoxLayout.X_AXIS));
		container.add(storageDir); 
		container.add(Box.createHorizontalStrut(5));
		container.add(browse);
		
		JLabel genericMsg = new JLabel(DESCRIPTION_MSG, 
				Utilities.getIcon("images/32x32/Information.png"), JLabel.LEFT);
		JLabel nextStepsMsg = new JLabel(ADDITIONAL_INFORMATION_NEXT_STEPS_MSG, 
				Utilities.getIcon("images/32x32/Note.png"), JLabel.LEFT);
		// Create the text area for the user to enter information.
		description = new JTextArea(3, 40);
		JScrollPane jsp = new JScrollPane(description);
		// Pack all the components together.
		final JPanel subPanel = 
		Utilities.createLabeledComponents("Enter storage directory", 
				"Generated output files will be stored here",
				1, false, container, genericMsg, jsp, nextStepsMsg);
		add(subPanel);
	}

	/**
	 * Method to handle click on the browse button in this page.
	 * 
	 * This method is invoked as the action listener when the
	 * user clicks the "Browse.." button on this page. This method
	 * displays a file chooser dialog to permit the user to set
	 * the storage directory interactively. 
	 */
	@Override
	public void actionPerformed(ActionEvent e) {
		JFileChooser jfc = new JFileChooser();
		jfc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
		if (jfc.showOpenDialog(this) == JFileChooser.APPROVE_OPTION) {
			storageDir.setText(jfc.getSelectedFile().getAbsolutePath());
		}
	}

	/**
	 * This method is called just before this page is to be displayed
	 * in the  wizard.
	 * 
	 * This method populates the default location where the generated
	 * files for this pipeline are to be stored. The default location
	 * is suggested based on the data set selected for the job.
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
		long   jobCnt   = Workspace.get().getJobList().getSeqCounter();
		String subDir   = "Job" + jobCnt + "_" + (jobCnt + 
				wizard.getDADXHelper().getAssemblerDetails().getJob().size() - 1);
		String baseDir  = Workspace.get().getDirectory() + File.separator + subDir;
		// Ensure storage directory is new by adding a suffix if needed.
		int    suffix   = 0;
		String storeDir = baseDir + File.separator;
		do {
			File dir = new File(storeDir);
			if (dir.exists()) {
				suffix++;
				storeDir = baseDir + "_" + suffix + File.separator;
			} else {
				// The storage directory does not exists.
				break;
			}
		} while (suffix < 10);
		// Set up suggestion for storage directory
		storageDir.setText(storeDir);
	}

	/**
	 * Obtain the directory on local machine where generated files
	 * are to be copied.
	 * 
	 * This method must be used by various wizard pages to determine
	 * the directory where the various files generated by running the
	 * DADX pipeline are to be stored. This method returns correct information
	 * only after the corresponding wizard page has been displayed and
	 * the user has entered the description.
	 * 
	 * @return This directory on the local machine where the generated
	 * output files are to be copied. 
	 */
	public String getStorageDir() {
		return storageDir.getText();
	}
	
	/**
	 * Obtain the description entered by the user for this pipeline run.
	 * 
	 * This method can be used to obtain the description entered by the
	 * user for this pipeline run. This method returns correct information
	 * only after the corresponding wizard page has been displayed and
	 * the user has entered the description.
	 * 
	 * @return The description entered by the user.
	 */
	public String getDescription() {
		return description.getText();
	}
	
	/**
	 * A brief description for the set of jobs being created by 
	 * this wizard. This field is meant to be used by the user
	 * to provide brief description for the type of test being
	 * conducted. This description is set for all the jobs that
	 * are created by this wizard.
	 */
	private JTextArea description;

	/**
	 * A text field in which the path to store generated files
	 * is to be set by the user. All the output files generated
	 * by this assembler pipeline will be placed in this 
	 * directory.
	 */
	private JTextField storageDir;
	
	/**
	 * A static message that is displayed in the wizard page that
	 * is used to obtain a brief description about the job(s)
	 * and purpose of the test.
	 */
	private static final String DESCRIPTION_MSG = "<html>"+
		"Enter a brief description about the purpose of this test.<br/>" +
		"<font size=\"-2\">This information if for your future reference " +
		"and can be anything.</font>" +
		"</html>";
	
	/**
	 * A static message that is displayed in the additional 
	 * information wizard page that provides the user with details
	 * about forthcoming operations.
	 */
	private static final String ADDITIONAL_INFORMATION_NEXT_STEPS_MSG = "<html>"+
		"In the next page a summary of the jobs that will be created<br/>" +
		"to run this pipleline (on the select server) will be displayed.<br/> " +
		"Once you have verified the information subsquent steps will<br/>" +
		"create and submit the jobs and the first job will start running." +
		"</html>";
	
	/**
	 * A generated serialization UID to keep the compiler happy.
	 */
	private static final long serialVersionUID = 2983800963465222319L;
}