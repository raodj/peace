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

package org.peace_tools.generic;

import java.awt.BorderLayout;
import java.awt.Component;

import javax.swing.JPanel;


/**
 * This class serves as a generic wizard page.
 * 
 *  This page provides a generic mechanism to display a component
 *  in a wizard. Typically, this page is used for displaying some
 *  information such as introductory message to provide a 
 *  overview of the activities involved with a given wizard. 
 *  This page is rather simple in the sense that it simply displays 
 *  a component. It does not perform any other operation and readily
 *  lets the user switch the wizard page to the next page.
 *  
 *  <p><b>Note:</b>This page uses a BorderLayout by default. You can add any 
 *  component to this page or even change the default layout.</p> 
 */
public class GenericWizardPage extends JPanel implements WizardPage {
	/**
	 * The default constructor.
	 * 
	 * The default constructor merely sets up the this panel with a
	 * border layout. The message to be displayed is set later on
	 * through a call to the setMessage() method.
	 */
	public GenericWizardPage() {
		super(new BorderLayout(0, 0));
	}
	
	/**
	 * This method must be used to set the title and subtitle for
	 * this panel.
	 * 
	 * @param title The title to be set for this wizard page.
	 * 
	 * @param subTitle The subtitle to be set for this wizard page.
	 */
	public void setTitle(String title, String subTitle) {
		this.title    = title;
		this.subTitle = subTitle;
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

	/**
	 * This method returns the title set for this web page. 
	 * 
	 * This method is invoked by the core wizard dialog whenever it
	 * needs to display the title for this page.
	 * 
	 * @return The title set for this page.
	 */
	@Override
	public String getTitle() {
		// TODO Auto-generated method stub
		return title;
	}

	/**
	 * Method to setup this page just before it is displayed.
	 * 
	 * This method is invoked by the core WizardDialog panel just
	 * before this page is displayed. This method does not have
	 * any special setup to do and is empty.
	 */
	@Override
	public void pageChanged(WizardDialog dialog, int currPage, int prevPage) {
	}

	/**
	 * Method to veto a page change (if needed).
	 * 
	 * This method is invoked by the core WizardDialog panel just
	 * before the user changes from this page. This method does not
	 * veto such a change and always returns true.
	 * 
	 * @return This method always returns true to indicate that the
	 * user can navigate off from this page.
	 */
	@Override
	public boolean pageChanging(WizardDialog dialog, int currPage, int nextPage) {
		return true;
	}
	
	/**
	 * The title to be set for this wizard page. This value is set
	 * via the setTitle method and is returned via the getTitle
	 * method.
	 */
	private String title;
	
	/**
	 * The subtitle to be set for this wizard page. This value is set
	 * via the setTitle method and is returned via the getSubTitle
	 * method.
	 */
	private String subTitle;

	/**
	 * The generated serial UID for serializing this class.
	 */
	private static final long serialVersionUID = -9111663142182703629L;
}
