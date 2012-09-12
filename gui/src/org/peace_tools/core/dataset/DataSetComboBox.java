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

package org.peace_tools.core.dataset;

import java.awt.Color;
import java.awt.Component;
import java.io.File;
import java.util.Vector;

import javax.swing.ImageIcon;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.ListCellRenderer;

import org.peace_tools.generic.Utilities;
import org.peace_tools.workspace.DataSet;
import org.peace_tools.workspace.DataSet.DataFileType;
import org.peace_tools.workspace.FileEntry;
import org.peace_tools.workspace.FileEntry.FileEntryType;
import org.peace_tools.workspace.GeneratedFileList;
import org.peace_tools.workspace.Workspace;

/**
 * A GUI sub-component used permit user to select a given type of
 * data set from the workspace.
 * 
 * This component was primarily introduced to permit various wizards
 * to share common functionality of permitting a user to select
 * a data set entry in the workspace. This GUI component provides a 
 * customized JComboBox that can is populated with data set entries. 
 * In addition, this class also provides a custom cell renderer that
 * provides a more detailed information about the data set entries.
 * 
 */
public class DataSetComboBox extends JComboBox implements ListCellRenderer {
	/**
	 * This is a simple internal class that is used to streamline the
	 * operations of this class. This class can contain either a DataSet
	 * object or a MSTData object as its sourceDataObject. It contains
	 * the HTML description that is built  
	 */
	 public class Choice {
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
			String dsInfo = "<html><b>" + Utilities.trim(tmpData.getName(), 30) + "</b><br/>" +
				"<font size=\"-2\">" + Utilities.trim(tmpData.getPath(), 30) + "</font><br/>" +
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
			String dsInfo = "<html><b>" + Utilities.trim(mstFile.getName(), 30) + "</b> " +
				"[<font size=\"-2\">" + Utilities.trim(mstFile.getPath(), 30) + "</font>]<br/>" +
				"<font size=\"-2\"><b>EST File:</b>" + Utilities.trim(estFile.getName(), 30) + 
				" [<font size=\"-2\">" + Utilities.trim(estFile.getPath(), 30) + "</font>]<br/>" +
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
	 * The constructor. The constructor sets up the entries in the
	 * combo-box from the workspace.
	 * 
	 * @param wizard The wizard that logically owns this page.
	 * 
	 * @param dataSetFlag If this flag is true, then this page presents
	 * the data sets in the workspace as primary inputs for assembly (
	 * for clustering+assembly). If this flag is false, then MST files
	 * are presented.
	 * 
	 * @param fastaOnly If this flag is true then only data sets that
	 * are associated with FASTA files are populated in the combo-box.
	 */
	public DataSetComboBox(final boolean dataSetFlag, 
			final boolean fastaOnly) {
		// Create the list of Choice objects
		Vector<Choice> comboBoxChoices = (dataSetFlag ? getDataSetsAsChoices(fastaOnly) : 
			                                            getMSTsAsChoices());
		// Populate the choices in our combo box.
		for(Choice entry: comboBoxChoices) {
			addItem(entry);
		}
		// Finally create the label to be used to display the options
		// and setup a custom renderer.
		cellRenderer = new JLabel();
		// Setup this class as the cell renderer
		setRenderer(this);
	}

	/**
	 * Helper method to build a list of choices based on the data sets
	 * (with EST file and possibly a quality file) to be used in the
	 * combo-box. This method is invoked only once from the constructor.
	 * It was primarily introduced to streamline the code in the 
	 * constructor.
	 * 
	 * @param fastaOnly If this flag is true then only data sets that
	 * are associated with FASTA files are returned.
	 * 
	 * @return A vector containing the various data sets wrapped in
	 * a Choice object.
	 */
	private Vector<Choice> getDataSetsAsChoices(final boolean fastaOnly) {
		Vector<Choice> choiceList = new Vector<Choice>();
		for(DataSet ds: Workspace.get().getDataSets()) {
			if (!ds.isGood()) {
				// Not a useful/valid choice.
				continue;
			}
			if ((fastaOnly) && (!ds.getFileType().equals(DataFileType.FASTA))) {
				// We care only for data sets associated with FASTA files.
				continue;
			}
			// Found valid choice to be added to the resulting list.
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
				if ((mst != null) && (mst.isGood())) {
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
	 * Obtain the DataSet or MSTData object that has been selected 
	 * by the user.
	 * 
	 * This method must be used to obtain the source entry selected
	 * by the user in this combo-box. 
	 * 
	 * @return The DataSet or MSTData object corresponding to the 
	 * selection made by the user. If the user does not have a 
	 * valid selection to make, then this method returns null.
	 */
	public Object getSelection() {
		int selIdx = getSelectedIndex();
		if (selIdx == -1) {
			return null;
		}
		Choice currChoice = (Choice) getSelectedItem();
		return currChoice.getSourceDataObject();
	}
	
	/**
	 * This label is used by this class to render the Choice objects
	 * presented by this class to the user. This cell renderer is created
	 * once (in the constructor) and reused each time the 
	 * {@link #getListCellRendererComponent(JList, Object, int, boolean, boolean)}
	 * method is called (by JComboBox)
	 */
	private final JLabel cellRenderer;
	
	/**
	 * A serialization UID to keep the compiler happy.
	 */
	private static final long serialVersionUID = -3482454971391881407L;
}
