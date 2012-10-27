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

package org.peace_tools.core.job.baton;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JTextArea;
import javax.swing.SpinnerNumberModel;
import javax.swing.border.EmptyBorder;

import org.peace_tools.generic.GenericWizardPage;
import org.peace_tools.generic.Utilities;
import org.peace_tools.generic.WizardDialog;
import org.peace_tools.workspace.DataSet;
import org.peace_tools.workspace.Server;
import org.peace_tools.workspace.Workspace;

class BatonHeadWizardPage extends GenericWizardPage 
implements ActionListener{
	
		public BatonHeadWizardPage(BatonJobWizard job){
			this.wizard = job;
			assert(this.wizard != null);
			// Setup the title(s) for this page and border
			setTitle("K-mer selection", 
					"Select number of nucleotides comprising the BATON head.");
			setBorder(new EmptyBorder(5, 5, 5, 5));
			// Create the combo-box with list of data sets.
			kmers = new JSpinner(new SpinnerNumberModel(2,2,3,1));
			Utilities.adjustDimension(kmers, 50, 4);
		
			JLabel title = new JLabel("Select k-mers: ");
			kmers.setBackground(Color.white);
					
			// Create panel with the combo box for the user to choose
			
					
			// Create the informational label.
			JLabel info = new JLabel(INFO_MSG, 
					Utilities.getIcon("images/32x32/Information.png"), 
					JLabel.LEFT);
			
			// Pack the input fields into a box
			JPanel subPanel = Utilities.createLabeledComponents("Select k-mer value for this job", ("(Number of nucleotides comprising the BATON head)"), 0, true,
				kmers);
			subPanel.setBorder(new EmptyBorder(5, 15, 10, 10));
			// Add the contents to this page
			add(subPanel, BorderLayout.CENTER);
			
	
		}
		
		/**
		 * This method returns the title set for this wizard page. 
		 * 
		 * This method is invoked by the core wizard dialog whenever it
		 * needs to display the title for this page.
		 * 
		 * @return The title set for this page.
		 */
		@Override
		public String getTitle() {
			return title;
		}
		
		/**
		 * Method to return the component for this wizard page.
		 * 
		 * @return This method simply returns this to display this 
		 * component as the overview page.
		 */
		@Override
		public Component getPage() {
			return this;
		}
		
		/**
		 * This method returns the sub-title set for this web page. 
		 * 
		 * This method is invoked by the core wizard dialog whenever it
		 * needs to display the sub-title for this page.
		 * 
		 * @return The sub-title set for this page.
		 */
		@Override
		public String getSubTitle() {
			return subTitle;
		}
		
		@Override
		public void actionPerformed(ActionEvent arg0) {
			// TODO Auto-generated method stub
			
		}
		
		/**
		 * Obtain information about the k-mer amount to be displayed in 
		 * the summary page. 
		 */
		protected String getSummary(String indent) {
			
			return indent + "K-mers: " + kmers.getValue().toString() +  "\n";
		}
		
		protected int getKmers(){
			return Integer.parseInt(kmers.getValue().toString());
		}
		
		/**
		 * A reference to the wizard dialog that logically owns this
		 * page. This reference is used to enable and disable 
		 * buttons on this wizard appropriately.
		 */
		private final WizardDialog wizard;
		
		/**
		 * A combo-box to select the data set to be used for this job.
		 * This field is created and populated in the constructor.
		 */
		private JSpinner kmers;

		/**
		 * Field to read and edit a brief description about the job.
		 * This information can be anything the user desires and is
		 * meaningful only to the user.
		 */
		private JTextArea description;
	
		
		/**
		 * A generic informational message that is displayed at the
		 * top of this wizard page to provide some additional information
		 * to the user.
		 */
		private static final String INFO_MSG = 
			"<html>Select the data set that contains the EST file to be<br>" +
			"processed by this job. Subsequent wizard pages will<br>" +
			"permit setting up additional information for the job.";
		
		/**
		/**
		 * The title for this page 
		 */
		private final String title = "K-mer selection"; 
		/**
		 * The sub-title for this page			
		 */
		private final String subTitle = "Select number of nucleotides comprising the BATON head.";
		/**
		 * A serialization UID to keep the compiler happy.
		 */
		private static final long serialVersionUID = -1193172876744546759L;
}
	
