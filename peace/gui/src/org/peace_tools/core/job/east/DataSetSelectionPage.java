package org.peace_tools.core.job.east;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.Vector;

import javax.swing.Box;
import javax.swing.ImageIcon;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.ListCellRenderer;
import javax.swing.border.EmptyBorder;

import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.Utilities;
import org.peace_tools.generic.WizardDialog;
import org.peace_tools.workspace.DataSet;
import org.peace_tools.workspace.DataSet.DataFileType;
import org.peace_tools.workspace.FileEntry;
import org.peace_tools.workspace.FileEntry.FileEntryType;
import org.peace_tools.workspace.GeneratedFileList;
import org.peace_tools.workspace.Workspace;

/**
 * This class serves typically as the second interactive page in
 * a EAST Job wizard. This page is used in one of the following
 * two different ways:
 * 
 * <ol>
 * 
 * <li>It is used by the EASTJobWizard to list the various data sets
 * (that contain valid MST and clustering file) that can be used for
 * assembly by EAST.</li>
 * 
 * <li>It can be used by the PEACEEASTJobWizard class to list the
 * various valid data sets that can be used for clustering followed
 * by assembly.</li>
 * 
 * </ol>
 * 
 * The choice of operation in the above two options is decided by 
 * the flag value that is passed in via the constructor.
 */
public class DataSetSelectionPage extends GenericWizardPage 
implements ListCellRenderer, ActionListener {
	/**
	 * This is a simple internal class that is used to streamline the
	 * operations of this class. This class can contain either a DataSet
	 * object or a MSTData object as its sourceDataObject. It contains
	 * the HTML description that is built  
	 */
	 class Choice {
		/**
		 * The source workspace data object that is Choice object
		 * actually represents. This value is set when this object
		 * is created and is never changed. This object can either:
		 * 
		 * <ul>
		 * 
		 * 	<li>A DataSet object in cases where DataSetSectionPage is used
		 *      to select a data set for clustering and assembly.</li>
		 *      
		 *  <li>A MSTData object in cases where the DataSetSelectionPage is
		 *  used to select data for direct assembly.</li>
		 *  
		 * </ul>
		 */
		private final Object sourceDataObject;
		
		/**
		 * This string is created when the Choice object is constructed. This
		 * string contains an HTML fragment that used to display the 
		 * information associated with the choice in a pretty manner.
		 */
		private final String infoString;
				
		/**
		 * A boolean flag to indicate if this choice entry has a corresponding
		 * quality data file that can be used during assembly.
		 */
		private final boolean haveQualFile;
		
		/** 
		 * A simple integer value that can be used by a special renderer to
		 * select an appropriate icon to be associated with this entry in
		 * a list (or a combo box). Note that this is only a suggestion/hint
		 * based on the type of the {@link #sourceDataObject}. Valid values
		 * are in the range 0 (data set), 1 (MST data).
		 */
		private final int iconHint;
		
		/**
		 * Constructor to create a Choice entry with a DataSet object as the
		 * primary source object. This choice object logically represents the
		 * given data set and provides a convenient interface for remainder
		 * of the operations in this wizard page.
		 * 
		 * @param dataSet The data set to be used as the primary source of
		 * information in this choice object. This parameter cannot be null.
		 */
		public Choice(DataSet dataSet) {
			// Extract the information string to be set for this choice.
			File tmpData  = new File(dataSet.getPath());			
			String dsInfo = "<html><b>" + Utilities.trim(tmpData.getName(), 40) + "</b><br/>" +
				"<font size=\"-2\"><i>" + Utilities.trim(dataSet.getDescription(), 60) + "</i></font></html>";
			// Setup the various instance variables
			this.sourceDataObject = dataSet;
			this.infoString       = dsInfo;
			this.haveQualFile     = dataSet.getFileType().equals(DataFileType.SFF);
			this.iconHint         = 0;
		}
		
		/**
		 * Constructor to create a Choice entry with a MST file object as the
		 * primary source object. This choice object logically represents the
		 * given MST file and provides a convenient interface for remainder
		 * of the operations in this wizard page.
		 * 
		 * @param mst The MST data object to be used as the primary source of
		 * information in this choice object. This parameter cannot be null.
		 */
		public Choice(FileEntry mst) {
			// Extract the information string to be set for this choice.
			File mstFile  = new File(mst.getPath());
			File estFile  = new File(mst.getGFL().getDataSet().getPath());
			String dsInfo = "<html><b>" + Utilities.trim(mstFile.getName(), 40) + "</b><br/>" +
				"<font size=\"-2\"><b>EST File:</b>" + Utilities.trim(estFile.getName(), 40) + "<br/>" +
				"<font size=\"-2\"><i>" + Utilities.trim(mst.getDescription(), 60) + "</i></font></html>";
			// Setup the various instance variables
			this.sourceDataObject = mst;
			this.infoString       = dsInfo;
			this.haveQualFile     = mst.getGFL().getDataSet().getFileType().equals(DataFileType.SFF);
			this.iconHint         = 1;
		}

		/**
		 * Obtain the source data object represented by this choice object.
		 * 
		 * @return This method returns a DataSet object or a MSTData object. This
		 * object is the same object that was used to create this object.
		 */
		public Object getSourceDataObject() {
			return sourceDataObject;
		}

		/**
		 * Obtain the HTML formatted information string for this choice object.
		 * 
		 * @return A HTML formatted information string that provides information
		 * about this Choice.
		 */
		public String getInfoString() {
			return infoString;
		}

		/**
		 * Determine if the source object also has a quality file associated with it.
		 * 
		 * @return This method returns true if the source data object has a 
		 * quality file associated with it. Otherwise this method returns false.
		 */
		public boolean hasQualFile() {
			return haveQualFile;
		}

		/**
		 * Obtain a suggestion/hint for any icon to be associated with this 
		 * choice.
		 * 
		 * @return This method returns an integer that provides suggestions or 
		 * hints for the icon to be used (if any) when displaying this choice
		 * in the GUI. Valid return values are:
		 * 
		 *  <ol start="0">
		 *  
		 *   <li>Data Set</li>
		 *  
		 *   <li>MST data</li>
		 *  
		 *  </ol>
		 */
		public int getIconHint() {
			return iconHint;
		}
		
	}
	
	/**
	 * The constructor. The constructor sets up the various components
	 * on this wizard page. The components include: combo box to
	 * select data set, a text area for job description, and a check
	 * box for clustering. 
	 * 
	 * @param wizard The wizard that logically owns this page.
	 * 
	 * @param dataSetFlag If this flag is true, then this page presents
	 * the data sets in the workspace as primary inputs for assembly (
	 * for clustering+assembly). If this flag is false, then MST files
	 * are presented.
	 */
	public DataSetSelectionPage(EASTJobWizard wizard, boolean dataSetFlag) {
		this.wizard = wizard;
		assert(this.wizard != null);
		// Setup the title(s) for this page and border
		setTitle("DataSet to Assemble", 
				"Select the data set to be clustered & assembled");
		setBorder(new EmptyBorder(5, 5, 5, 5));
		// Create the combo-box with list of Choice objects
		Vector<Choice> comboBoxChoices = (dataSetFlag ? getDataSetsAsChoices() : 
			                                            getMSTsAsChoices());
		dataSetList = new JComboBox(comboBoxChoices);
		dataSetList.setRenderer(this);
		dataSetList.addActionListener(this);
		
		// Create panel with the combo box for the user to choose
		JComponent dataSetBox = 
			Utilities.createLabeledComponents("Select Data Set for Job:",
					"(The cDNA file in data set will be clustered & assembled)", 0, false, 
					dataSetList);
		
		// Create the quality file usage check box and associated information label.
		useQualData = new JCheckBox("Use quality data for assembly");
		Utilities.adjustFont(useQualData, 0, 8, 1);
		useQualDataInfo = new JLabel(QUAL_FILE_AVAILABLE, 
				Utilities.getIcon("images/32x32/Note.png"), 
				JLabel.LEFT);
		// Pack the check box and the associated informational label into a panel
		JPanel qualPanel = new JPanel(new BorderLayout(0, 2));
		qualPanel.add(useQualData, BorderLayout.NORTH);
		qualPanel.add(useQualDataInfo, BorderLayout.SOUTH);
		
		// Create the informational label.
		JLabel info = new JLabel(INFO_MSG, 
				Utilities.getIcon("images/32x32/Information.png"), 
				JLabel.LEFT);
		
		// Pack the input fields into a box
		JPanel subPanel = Utilities.createLabeledComponents(null, null, 0, true,
			info, Box.createVerticalStrut(10),
			dataSetBox, Box.createVerticalStrut(10),
			qualPanel);
		subPanel.setBorder(new EmptyBorder(5, 15, 10, 10));
		// Add the contents to this page
		add(subPanel, BorderLayout.CENTER);
		// Finally create the label to be used to display the options
		// and setup a custom renderer.
		cellRenderer = new JLabel();
	}

	/**
	 * Helper method to build a list of choices based on the data sets
	 * (with EST file and possibly a quality file) to be used in the
	 * combo box. This method is invoked only once from the constructor.
	 * It was primarily introduced to streamline the code in the 
	 * constructor.
	 * 
	 * @return A vector containing the various data sets wrapped in
	 * a Choice object.
	 */
	private Vector<Choice> getDataSetsAsChoices() {
		Vector<Choice> choiceList = new Vector<Choice>();
		for(DataSet ds: Workspace.get().getDataSets()) {
			Choice c = new Choice(ds);
			choiceList.add(c);
		}
		return choiceList;
	}

	/**
	 * Helper method to build a list of choices based on the MST files
	 * to be used in the combo box. This method is invoked only once 
	 * from the constructor.  It was primarily introduced to 
	 * streamline the code in the constructor.
	 * 
	 * @return A vector containing the various MST files in the data set
	 * wrapped in a Choice object.
	 */
	private Vector<Choice> getMSTsAsChoices() {
		Vector<Choice> choiceList = new Vector<Choice>();
		for(DataSet ds: Workspace.get().getDataSets()) {
			for(GeneratedFileList gfl: ds.getGflList()) {
				FileEntry mst = gfl.findEntry(FileEntryType.MST);
				if (mst != null) {
					Choice c = new Choice(mst);
					choiceList.add(c);
				}
			}
		}
		return choiceList;
	}
	
	/**
	 * This method is called by core Java GUI system whenever it
	 * needs to render a data set or MST entry for the user to
	 * choose from the {@link #dataSetList} combo box. This
	 * method appropriately updates the {@link #cellRenderer}
	 * label (created in the constructor) and return the
	 * label back. See JavaDoc on this method for additional
	 * details.
	 */
	@Override
	public Component getListCellRendererComponent(JList list, Object value,
			int index, boolean isSelected, boolean cellHasFocus) {
		// Set the default background color depending on selection
		cellRenderer.setOpaque(true);
        if (isSelected) {
            cellRenderer.setBackground(list.getSelectionBackground());
            cellRenderer.setForeground(list.getSelectionForeground());
        } else {
        	cellRenderer.setBackground(Color.white);
        	cellRenderer.setForeground(list.getForeground());
        }

        // Set the text to be rendered based on the current choice.
        Choice currChoice = (Choice) value;
        if (currChoice != null) {
        	cellRenderer.setText(currChoice.getInfoString());
        	//Set the icon to make things look pretty
        	final String IconNames[] = {"NewCluster", "NewMST"};
        	ImageIcon icon = Utilities.getIcon("images/24x24/" + 
        			IconNames[currChoice.getIconHint()] + ".png");
        	cellRenderer.setIcon(icon);
        }
        return cellRenderer;
	}
	
	/**
	 * Obtain the DataSet of MSTData object that has been selected 
	 * by the user.
	 * 
	 * This method is used by the EASTJobWizard to obtain the
	 * source entry selected by the user. 
	 * 
	 * @return The DataSet or MSTData object corresponding to the 
	 * selection made by the user. If the user does not have a 
	 * valid selection to make, then this method returns null.
	 */
	protected Object getSelection() {
		int selIdx = dataSetList.getSelectedIndex();
		if (selIdx == -1) {
			return null;
		}
		Choice currChoice = (Choice) dataSetList.getSelectedItem();
		return currChoice.getSourceDataObject();
	}
	
	/**
	 * This method is called whenever the user selects a 
	 * data set or MST from the {@link #dataSetList} combo box.
	 * This method appropriately enables/disables the
	 * {@link #useQualData} check box and updates the
	 * informational label {@link #useQualDataInfo} (that is
	 * displayed below the check box).
	 * 
	 * @param arg0 The action event associated with a call.
	 * Currently this event is ignored (and can be null).
	 */
	@Override
	public void actionPerformed(ActionEvent arg0) {
		Choice currDS = (Choice) dataSetList.getSelectedItem();
		boolean haveQualFile = currDS.hasQualFile();
		// Enable/disable check box depending on whether we have
		// quality file or not
		useQualData.setEnabled(haveQualFile);
		useQualData.setSelected(haveQualFile);
		// Update the info-label to provide the user with additional
		// information
		useQualDataInfo.setText(haveQualFile ? QUAL_FILE_AVAILABLE : 
			QUAL_FILE_NOT_AVAILABLE);
		// Update tool tip based on selection?
	}
	
	/**
	 * Method to setup this page just before it is displayed.
	 * 
	 * This method is invoked by the core WizardDialog panel just
	 * before this page is displayed. This method updates the
	 * check box status based on the current selection in the
	 * combo box {@link #dataSetList}.
	 */
	@Override
	public void pageChanged(WizardDialog dialog, int currPage, int prevPage) {
		int selIdx = Math.max(dataSetList.getSelectedIndex(), 0);
		if (dataSetList.getItemCount() > selIdx) {
			dataSetList.setSelectedIndex(selIdx);
		}
	}
	
	/**
	 * Method to setup this page just before it is displayed.
	 * 
	 * This method is invoked by the core WizardDialog panel just
	 * before this page is displayed. This method updates the
	 * check box status based on the current selection in the
	 * combo box {@link #dataSetList}.
	 */
	@Override
	public boolean pageChanging(WizardDialog dialog, int currPage, int prevPage) {
		if (currPage < prevPage) {
			return true;
		}
		// Before letting the user move on, ensure that the files
		// we are looking for are valid.
		Object selection = this.getSelection();
		if (selection == null) {
			// No valid selection? This is not right.
			return false;
		}
		// Extract the path for the files we care about.
		String dataSetPath = null;
		String mstFilePath = null;
		if (selection instanceof DataSet) {
			dataSetPath = ((DataSet) selection).getPath();
		} else {
			FileEntry mst = ((FileEntry) selection);
			mstFilePath = mst.getPath();
			dataSetPath = mst.getGFL().getDataSet().getPath();
		}
		// In all cases data set file must exist.
		File cDNAFile = new File(dataSetPath);
		if (!cDNAFile.exists() || !cDNAFile.canRead()) {
			return false;
		}
		// Check MST file if set
		if (mstFilePath != null) {
			File mstFile = new File(mstFilePath);
			if (!mstFile.exists() || !mstFile.canRead()) {
				
			}
		}
		return false;
	}

	/**
	 * A reference to the wizard dialog that logically owns this
	 * page. This reference is used to enable and disable 
	 * buttons on this wizard appropriately.
	 */
	private final EASTJobWizard wizard;
	
	/**
	 * A combo-box to select the data set to be used for this job.
	 * This field is created and populated in the constructor.
	 */
	private JComboBox dataSetList;
	
	/**
	 * This label is used by this class to render the Choice objects
	 * presented by this class to the user. This cell renderer is created
	 * once (in the constructor) and reused each time the 
	 * {@link #getListCellRendererComponent(JList, Object, int, boolean, boolean)}
	 * method is called (by JComboBox)
	 */
	private JLabel cellRenderer;
	
	/**
	 * A check box that provides a Boolean value to indicate if the quality
	 * file associated with a given data set should be used for assembly.
	 * This check box is enabled and disabled (depending on whether the
	 * data set has a quality file associated with it) so that the user
	 * can select or deselect) the option.
	 */
	private JCheckBox useQualData;
	
	/**
	 * A simple HTML formatted string that provides some basic information
	 * to help the user about as to why quality data option usage is
	 * currently unavailable.
	 */
	private static final String QUAL_FILE_NOT_AVAILABLE =
		"<html><font size=\"-2\">" +
		"This option is currently disabled because a suitable quality<br/>" +
		"file is not available for this data set. You need to add a<br/>" +
		"suitable quality file to this data set to enable this option." +
		"</font></html>";

	/**
	 * A simple HTML formatted string that provides some basic information
	 * to inform the user  as to how the quality data is going to be 
	 * used purely for assembly purposes.
	 */
	private static final String QUAL_FILE_AVAILABLE =
		"<html><font size=\"-2\">" +
		"If the above box is checked, then the quality file associated<br/>" +
		"with this data set will be used during assembly to improve the<br/>" +
		"overall quality the contigs generated after assembly." +
		"</font></html>";

	/**
	 * A label that is used to display information to the user as to why 
	 * the {@link #useQualData} check box is enabled or disabled. The
	 * text in this label is changed between {@link #QUAL_FILE_AVAILABLE}
	 * and {@link #QUAL_FILE_NOT_AVAILABLE} depending on the currently
	 * selected Choice in the {@link #dataSetList}
	 */
	private JLabel useQualDataInfo;
	
	/**
	 * A generic informational message that is displayed at the
	 * top of this wizard page to provide some additional information
	 * to the user.
	 */
	private static final String INFO_MSG = 
		"<html>Select the data set that contains the cDNA file to be<br>" +
		"clustered & assembled. Subsequent wizard pages<br>" +
		"will permit setting assembly parameters and server <br>"+
		"configuration for clustering.</html>";
	
	/**
	 * A serialization UID to keep the compiler happy.
	 */
	private static final long serialVersionUID = -1437405859971218347L;
}
