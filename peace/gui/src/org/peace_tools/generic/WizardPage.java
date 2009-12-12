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
// Authors:   Dhananjai M. Rao          raodm@muohio.edu
//
//---------------------------------------------------------------------

package org.peace_tools.generic;

import java.awt.Component;

/**
 * A simple interface that is implemented by all wizard pages so  
 * that they can be displayed in a wizard dialog. This interface
 * provides the hand shake between a wizard dialog and a wizard
 * page.
 */
public interface WizardPage {
	/**
	 * Obtain the actual component associated with this wizard page.
	 * 
	 * This method must be implemented to return the actual graphics
	 * component associated with this page. This is typically a panel
	 * that contains all the dialog elements associated with this page.
	 * 
	 * @return The actual graphics component to be displayed for this
	 * page.
	 */
	Component getPage();
	
	/**
	 * This method must be implemented to return a title for this
	 * wizard page. The title is displayed in the wizard dialog's
	 * title pane. In addition the title is used to display the
	 * list of pages in the overview sequence.
	 * 
	 * @return A non-null title for this page. The title typically is
	 * 2-3 short words that describe overall functionality of the
	 * wizard page.
	 */
	public String getTitle();
	
	/**
	 * This method must be implemented to return a sub-title for this
	 * wizard page. The sub-title is displayed in the wizard dialog's
	 * title pane just below the title. 
	 * 
	 * @return An optional sub-title for the wizard page. If the wizard
	 * page does not have a sub-title then this method can return null.
	 */
	public String getSubTitle();

	/**
	 * This method is invoked to indicate indicate to the currently
	 * visible page that the user is moving away from the current
	 * page by clicking either on the "Previous" or the "Next" button.
	 * 
	 * @param dialog The wizard dialog that is invoking this method.
	 * 
	 * @param currPage The logical (zero-based) index of the current
	 * page in the sequence of pages in the wizard.
	 * 
	 * @param nextPage The logical (zero-based) index of the next
	 * page that is going to be displayed.
	 * 
	 * @return If this method returns true, then the switching is
	 * completed. Otherwise, the switching is canceled and the
	 * current page continues to be displayed.
	 */
	public boolean pageChanging(WizardDialog dialog, int currPage, 
			int nextPage);
	
	/**
	 * This method is invoked on a wizard page just before it is
	 * going to be displayed. This method can be used by the
	 * wizard page to update the information displayed by it
	 * or enable/disable some of the buttons in the parent
	 * 
	 * @note This method must call wizardDialog.setButtonStatus()
	 * with suitable parameters.
	 * 
	 * @param dialog The wizard dialog that is invoking this method.
	 * 
	 * @param currPage The logical (zero-based) index of the current
	 * page in the sequence of pages in the wizard that is going to
	 * be displayed.
	 * 
	 * @param prevPage The logical (zero-based) index of the previous
	 * page that was displayed. For the first page this value is
	 * -1.
	 */
	public void pageChanged(WizardDialog dialog, int currPage, int prevPage);
}
