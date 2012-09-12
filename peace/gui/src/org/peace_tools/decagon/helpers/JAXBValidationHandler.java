//--------------------------------------------------------------------
//
// This DECAGON file is part of PEACE.
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

package org.peace_tools.decagon.helpers;

import javax.xml.bind.ValidationEvent;
import javax.xml.bind.ValidationEventHandler;
import javax.xml.bind.ValidationEventLocator;

import org.peace_tools.generic.UserLog;
import org.peace_tools.generic.Log.LogLevel;

/**
 * A simple class to gather and report validation errors when reading parameters
 * from an XML file.
 *
 * This class is used by the {@link ParameterSetGUIHelper#unmarshal(String)}
 * method to intercept XML-Schema validation events when unmarshalling (reading)
 * data from an parameters (XML) description file. This handler gathers data
 * about zero or more validation errors so that it can be suitably reported
 * back to the user.
 */
public class JAXBValidationHandler implements ValidationEventHandler {
	/**
	 * An internal string builder that is used to track validation
	 * error messages reported to this listener. Errors are added to
	 * this builder in the {@link #handleEvent(ValidationEvent)} method.
	 * If no errors occur then this buffer will be empty.
	 */
	private final StringBuilder errorLog = new StringBuilder();
	
	/**
	 * A format string that is appropriately filled-in to generate the
	 * necessary exception message.
	 */
	private static final String ERROR_FORMAT = "Event [Severity: %d]: %s.\n" +
		"Linked Exception: %s\n." + 
		"Encountered at Line %d, Column %d, Offset %d\n" +
		"in URL: %s (Object: %s, Node: %s)\n";
	
	/**
	 * Method to receive validation errors from JAXB schema validator.
	 * 
	 * This method is called by the internal JAXB schema validator whenever
	 * a Schema violation is discovered. This method adds the necessary
	 * information about the error to the {@link #errorLog}.
	 * 
	 * @param event The validation event that contains information about
	 * the error that was detected.
	 * 
	 * @return This method always returns true to indicate that the
	 * validator can proceed further to report any other issues.
	 */
	@Override
	public boolean handleEvent(ValidationEvent event) {
		// Extract any information about linked exception.
		final String lnex = (event.getLinkedException() != null ? 
				event.getLinkedException().getMessage() : "n/a");
		// Extract information about location to streamline code here.		
		final ValidationEventLocator location = event.getLocator();
		final String url  = (location.getURL()    == null ? "n/a" : location.getURL().toString());
		final String obj  = (location.getObject() == null ? "n/a" : location.getObject().toString());
		final String node = (location.getNode()   == null ? "n/a" : location.getNode().toString());
		// Generate error message.
		final String errorMessage =
			String.format(ERROR_FORMAT, event.getSeverity(), 
					event.getMessage(), lnex,
					location.getLineNumber(), location.getColumnNumber(),
					location.getOffset(), url, obj, node); 
		errorLog.append(errorMessage);
		UserLog.log(LogLevel.ERROR, "DECAGON", errorMessage);
		return true;
	}
	
	/**
	 * Obtain the errors (if any) tracked by this validation handler.
	 * 
	 * @return The errors tracked by this validation handler. If no errors
	 * were reported then this method returns an empty string.
	 * @return
	 */
	public String getErrorLog() {
		return errorLog.toString();
	}
}
