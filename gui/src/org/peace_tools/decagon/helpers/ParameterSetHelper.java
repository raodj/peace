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

package org.peace_tools.decagon.helpers;

import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

import javax.xml.XMLConstants;
import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Unmarshaller;
import javax.xml.transform.stream.StreamSource;
import javax.xml.validation.Schema;
import javax.xml.validation.SchemaFactory;

import org.peace_tools.decagon.DecagonVariables;
import org.peace_tools.decagon.jaxb.CmdLineKind;
import org.peace_tools.decagon.jaxb.ComparatorKind;
import org.peace_tools.decagon.jaxb.Condition;
import org.peace_tools.decagon.jaxb.ConnectorKind;
import org.peace_tools.decagon.jaxb.Parameter;
import org.peace_tools.decagon.jaxb.ParameterKind;
import org.peace_tools.decagon.jaxb.ParameterSet;
import org.peace_tools.generic.ProgrammerLog;
import org.peace_tools.generic.Utilities;

/**
 * Convenience class to help read parameters from an XML file and process
 * conditions.
 * 
 * This is a top-level class that can be used as the primary interface 
 * into this package. This class aids the remainder of DECAGON components
 * with the following operations:
 * 
 * <ul>
 * 
 * <li>Unmarshal (load/read) a parameter set from a given XML description
 * file. The XML file must be compatible with the ParameterSet.xsd file.</li>
 * 
 * <li>Manage DECAGON variables to be used for processing conditions and
 * parameterized values and command-line arguments.</li>  
 * 
 * <li>Process conditions associated with a parameter and generate overall 
 * boolean result of processing conditions.</li>
 * 
 * <li>Substitution of parameterized values and generation of command-line
 * arguments.</li>
 * 
 * </ul>
 */
public class ParameterSetHelper {
	/**
	 * The set of parameters that are currently being handled by this
	 * parameter set.
	 */
	protected ParameterSet paramSet;
	
	/**
	 * The set of DECAGON variables that are used to process parameter
	 * sets and conditions.
	 */
	protected DecagonVariables decVar;
	
	/**
	 * The default constructor for this class.
	 * 
	 * The constructor is relatively straightforward and merely initializes
	 * the various instance variables to their default initial value. The
	 * constructor instantiates a standard set of DECAGON variables and
	 * holds a reference to it in the {@link #decVar} instance variable.
	 * The helper is instantiated without any parameters to process. The
	 * parameters can be loaded from a suitable XML file via call to
	 * {@link #unmarshal(String)} method.
	 */
	public ParameterSetHelper() {
		decVar = new DecagonVariables();
	}
	
	/**
	 * The default and only constructor for this class.
	 * 
	 * The constructor is relatively straightforward and merely initializes
	 * the various instance variables to their default initial value.
	 * The constructor instantiates a standard set of DECAGON variables and
	 * holds a reference to it in the {@link #decVar} instance variable.
	 * 
	 * @param paramSet The set of parameters that are to be managed
	 * by this object. This parameter cannot be null.
	 */
	public ParameterSetHelper(ParameterSet paramSet) {
		this.decVar   = new DecagonVariables();
		this.paramSet = paramSet;
	}
	
	/**
	 * Load parameters from an XML configuration file.
	 * 
	 * This method provides a convenient mechanism to load a set of
	 * parameters from a XML file. The XML file must be compatible
	 * with the ParameterSet.xsd schema file. Existing parameters
	 * (if any) are replaced with the newly loaded parameters.
	 * 
	 * @param paramsFileName The name of the parameter file from where
	 * the XML description for the parameters is to be loaded. The
	 * file is opened via call to {@link Utilities#getStream(String)}
	 * method.
	 * 
	 * @throws Exception This method throws different types of exceptions
	 * depending of the type of error that occurs when loading the
	 * parameters from the XML description file.
	 * 
	 * @see Utilities#getStream(String)
	 */
	public void unmarshal(final String paramsFileName) throws Exception {
		// Load the schema for parameters for schema validation.
		SchemaFactory sf = SchemaFactory.newInstance(XMLConstants.W3C_XML_SCHEMA_NS_URI);
		StreamSource schemaFile = new StreamSource(Utilities.getStream("decagon/ParameterSet.xsd"));
		Schema paramSetSchema = sf.newSchema(schemaFile);

		// create JAXB context and instantiate unmarshaller
		JAXBContext context = JAXBContext.newInstance(ParameterSet.class);
		Unmarshaller um = context.createUnmarshaller();
		um.setSchema(paramSetSchema);
		// Setup handler to intercept validation errors.
		JAXBValidationHandler pvh = new JAXBValidationHandler();
		um.setEventHandler(pvh);
		// Get the unmarshaller to load the data for us.
		final InputStream xml = Utilities.getStream(paramsFileName);
		paramSet = (ParameterSet) um.unmarshal(xml);
		// Check and throw exception if any validation errors occur
		final String validationErrors = pvh.getErrorLog();
		if (validationErrors.length() != 0) {
			// Validation errors occurred
			throw new JAXBException(validationErrors);
		}
	}
	
	/**
	 * Obtain the list of parameters being used by this class.
	 * 
	 * @return The full list of parameters used by this class.
	 */
	public List<Parameter> getParameters() {
		return this.paramSet.getParameter();
	}
	
	/**
	 * Obtain the set of variables currently being used for processing
	 * the parameters associated with this class.
	 * 
	 * @return The set of variables being used with the set of parameters
	 * for this class.
	 */
	public DecagonVariables getVariables() {
		return this.decVar;
	}
	
	/**
	 * Find the parameter object with the given name.
	 * 
	 * This method iterates of the set of parameters (that this object
	 * is managing) to find the parameter with the given name.
	 * 
	 * @param paramName The name of the parameter to search for. The
	 * name is specified as an attribute of a parameter. Example:
	 * 
	 * <code>
	 *     <parameter name="SangerFixed">
	 *         <!-- Elements here are not shown for brevity -->
	 *     </parameter>
	 * </code>
	 * 
	 * @return The parameter object with the given name. If the parameter
	 * object was not found, then this method returns null.
	 */
	public Parameter getParameter(final String name) {
		for(Parameter param: paramSet.getParameter()) {
			if (name.equals(param.getName())) {
				return param;
			}
		}
		return null;
	}

	/**
	 * Method to translate variables to their values.
	 * 
	 * This method must be used to translate variables (if any) specified
	 * in {@link Parameter#getCmdLine()} and {@link Parameter#getValue()}
	 * elements to corresponding values. For example, assume variables
	 * <code>VAR1</code> and <code>VAR2</code> have values <code>Test</code>
	 * and <code>0.5</code> as values. Furthermore assumed variable
	 * <code>VAR3</code> is undefined. Then given the string
	 * <code>--dir %VAR1% --thres%% %VAR2% --code %VAR3%</code> 
	 * this method returns the string 
	 * <code>--dir Test --thres% 0.5 --code</code>. The current
	 * variables defined in {@link #decVar} are used.
	 * 
	 * @param strWithVars The string that may contain variable definitions
	 * to be processed by this method. This parameter cannot be null 
	 * but can be an empty string.
	 * 
	 * @return A string with values of variables suitably replaced.
	 */
	public String processVariables(final String strWithVars) {
		StringBuilder result = new StringBuilder(strWithVars.length());
		int currPos          = 0;
		int perStartPos      = -1;
		while ((perStartPos = strWithVars.indexOf('%', currPos)) != -1) {
			// Found a '%' sign. Everything before it is not a variable.
			result.append(strWithVars.substring(currPos, perStartPos));
			// Find the ending '%' sign next.
			final int perEndPos = strWithVars.indexOf('%', perStartPos + 1);
			if (perEndPos > -1) {
				// Found a pair of '%' signs. Everything between them is
				// assumed to be a variable. To be substituted.
				final String var = strWithVars.substring(perStartPos + 1, perEndPos);
				if (var.length() == 0) {
					// This is a '%%' use which translates to a '%'
					result.append('%');
				} else if (decVar.isDefined(var)) {
					// Append the variable's actual value.
					result.append(decVar.getStrValue(var));
				}
				// Update the current position after the second '%' sign.
				currPos = perEndPos + 1;
			} else {
				// Did not find trailing '%' sign This is a bad case.
				// Just break out of the loop to avoid infinite loops
				break;
			}
		}
		// All the trailing characters after the last variable need
		// to be appended to the final result.
		if (currPos < strWithVars.length()) {
			result.append(strWithVars.substring(currPos));
		}
		// Return the string with variables replaced with their values.
		return result.toString();
	}
	
	/**
	 * Obtain value for a specific boolean parameter given parameter name.
	 * 
	 * This method provides a convenient mechanism to obtain the 
	 * value for a {@link ParameterKind#BOOLEAN} parameter as a
	 * Boolean. The value by default is stored as a string. This method
	 * converts the string to a Boolean and returns the Boolean value.
	 * Variables referenced in the value are suitably replaced
	 * via call to {@link #processVariables(String)} method.
	 * 
	 * @param name The name of the boolean parameter whose value is
	 * to be returned.
	 * 
	 * @return The boolean value set for the parameter (assuming the
	 * parameter was found and its kind was {@link ParameterKind#BOOLEAN}).
	 * If the parameter was not found or if its kind was not 
	 * {@link ParameterKind#BOOLEAN} this method returns null.
	 */
	public Boolean getParameterValueAsBoolean(final String name) {
		final Parameter param = getParameter(name);
		return getParameterValueAsBoolean(param);
	}

	/**
	 * Obtain value for a specific boolean parameter given parameter name.
	 * 
	 * This method provides a convenient mechanism to obtain the 
	 * value for a {@link ParameterKind#BOOLEAN} parameter as a
	 * Boolean. The value by default is stored as a string. This method
	 * converts the string to a Boolean and returns the Boolean value.
	 * Variables referenced in the value are suitably replaced
	 * via call to {@link #processVariables(String)} method.
	 * 
	 * @param param The boolean parameter whose value is to be converted and
	 * returned.
	 * 
	 * @return The boolean value set for the parameter (assuming the
	 * parameter was found and its kind was {@link ParameterKind#BOOLEAN}).
	 * If the parameter was not of {@link ParameterKind#BOOLEAN} kind then
	 * this method returns null.
	 */
	public Boolean getParameterValueAsBoolean(final Parameter param) {
		if ((param != null) && (ParameterKind.BOOLEAN.equals(param.getKind()))) {
			final String actualValue = processVariables(param.getValue());
			return Boolean.valueOf(actualValue);
		}
		return null;
	}

	/**
	 * Obtain value for a specific integer or double parameter 
	 * given parameter name.
	 * 
	 * This method provides a convenient mechanism to obtain the 
	 * value for a {@link ParameterKind#INTEGER} or {@link ParameterKind#DOUBLE}
	 * parameter as a Number. The value by default is stored as a string. 
	 * This method converts the string to an Integer or Double (depending
	 * on the parameter's kind) and returns the Number value.
	 * Variables referenced in the value are suitably replaced
	 * via call to {@link #processVariables(String)} method.
	 * 
	 * @param name The name of the integer or double parameter whose value is
	 * to be returned.
	 * 
	 * @return A number that represents the integer or double value for the
	 * parameter, if the parameter was found and its kind was either
	 * {@link ParameterKind#INTEGER} or {@link ParameterKind#DOUBLE}. Otherwise
	 * this method returns null.
	 */
	public Number getParameterValueAsNumber(final String name) {
		final Parameter param = getParameter(name);
		return getParameterValueAsNumber(param);
	}

	/**
	 * Obtain value for a specific integer or double parameter. 
	 * 
	 * This method provides a convenient mechanism to obtain the 
	 * value for a {@link ParameterKind#INTEGER} or {@link ParameterKind#DOUBLE}
	 * parameter as a Number. The value by default is stored as a string. 
	 * This method converts the string to an Integer or Double (depending
	 * on the parameter's kind) and returns the Number value.
	 * Variables referenced in the value are suitably replaced
	 * via call to {@link #processVariables(String)} method.
	 * 
	 * @param param The integer or double parameter whose value is
	 * to be returned.
	 * 
	 * @return A number that represents the integer or double value for the
	 * parameter, if the parameter was either {@link ParameterKind#INTEGER} 
	 * or {@link ParameterKind#DOUBLE}. Otherwise this method returns null.
	 */
	public Number getParameterValueAsNumber(final Parameter param) {
		if ((param != null) && (ParameterKind.INTEGER.equals(param.getKind()))) {
			final String actualValue = processVariables(param.getValue());
			return Integer.valueOf(actualValue);
		} else if ((param != null) && (ParameterKind.DOUBLE.equals(param.getKind()))) {
			final String actualValue = processVariables(param.getValue());
			return Double.valueOf(actualValue);
		}
		return null;
	}

	/**
	 * Convenience method to generate a command-line string from parameters.
	 * 
	 * This is a convenience method that can be used to generate a command-line
	 * string from the given set of parameters. This method
	 * calls the {@link #getCmdLineArguments()} method to obtain
	 * the list of command-line arguments and then converts it to
	 * a single string using {@link Utilities#toString(List, String)}
	 * method a single space as the delimiter between consecutive
	 * arguments.
	 * 
	 * @return The parameters and their values suitably converted to
	 * an equivalent command-line argument list.
	 */
	public String getCmdLine() {
		return Utilities.toString(getCmdLineArguments(true), " ");
	}
	
	/**
	 * Method to obtain the command-line arguments from parameters.
	 * 
	 * This is a convenience method that can be used to generate a 
	 * list of command-line arguments from the given set of parameters.
	 * Conditional parameters (parameters that have conditions
	 * associated with them) are included only if the conditions
	 * validate to true. Otherwise they are not included in the 
	 * resulting list of command-line arguments. Parameters of 
	 * kind {@value CmdLineKind#IGNORE} are not included 
	 * in the command-line. Variables in the command-line and/or
	 * value are translated to their actual values via call 
	 * to {@link #processVariables(String)} method.
	 * 
	 * @param addVars If this flag is true then variables that are
	 * added to the command-line are also added to {@link #decVar} 
	 * so that their values can be used for processing subsequent
	 * conditions and variables. Moreover, variables that are not
	 * added (because condition evaluates to false) those variables
	 * are moved from {@link #decVar} impacting subsequent 
	 * condition processing.
	 * 
	 * @return The parameters and their values suitably converted to
	 * an equivalent command-line argument list.
	 */
	public List<String> getCmdLineArguments(final boolean addVars) {
		List<String> argList  = new ArrayList<String>(20);
		for(Parameter param: paramSet.getParameter()) {
			// Get helper method to add parameter to the argList after
			// performing necessary checks.
			final String actualValue = addToCmdLine(param, argList);
			// Appropriately add/remove entries if the variable is to be used.
			if (actualValue == null) {
				decVar.removeVariable(param.getName());
			} else {
				decVar.addVariable(param.getName(), actualValue);
			}
		}
		return argList;
	}

	/**
	 * Convenience method to process a parameter and suitably add an entry
	 * to be used as command-line argument.
	 * 
	 * This method provides a convenience interface to convert a parameter
	 * entry to a suitable command-line argument. Conditional parameters 
	 * (that is, parameters that have conditions associated with them) are 
	 * included only if the conditions validate to true. Otherwise they are
	 * not included in the resulting list of command-line arguments. 
	 * Parameters of kind {@value CmdLineKind#IGNORE} are not included 
	 * in the command-line. Variables in the command-line and/or
	 * value are translated to their actual values via call 
	 * to {@link #processVariables(String)} method.
	 * 
	 * @param param The parameter entry to be processed. This value 
	 * cannot be null.
	 * 
	 * @param argList The list to which suitable entries are to be added.
	 * This method may add zero, one, or two entries to this list depending
	 * on the type of the parameter. This parameter cannot be null.
	 *  
	 * @return This method returns the actual value added. The actual
	 * value added is null if no entries were added (because of conditions
	 * or because the parameter is tagged to be ignored).
	 */
	public String addToCmdLine(final Parameter param, List<String> argList) {
		// Check and process condition on this parameter (if any)
		if (!checkCondition(param.getCondition())) {
			// Condition check resulted in false. This parameter is ignored.
			return null;
		}
		// Process the parameter and obtain it's actual value for
		// further use after the switch-statement.
		String actualValue = (param.getValue() != null ? 
				processVariables(param.getValue()) : null);
		switch (param.getCmdLineKind()) {
		case IMPLICIT_VALUE:
			if (ParameterKind.BOOLEAN.equals(param.getKind())) {
				if (Boolean.valueOf(actualValue)) {
					argList.add(processVariables(param.getCmdLine()));
				}
			}
			break;
		case EXPLICIT_VALUE:
			argList.add(processVariables(param.getCmdLine()));
			argList.add(actualValue);
			break;
		case VALUE_ONLY:
			argList.add(actualValue);
			break;
		case IGNORE:
			// This parameter is to be ignored.
			actualValue = null;
			break;
		default:
			ProgrammerLog.log("Unhandled command-line kind" + param.getCmdLineKind());
		}
		return actualValue;
	}
	
	/**
	 * Method to process a given condition in a depth-first manner.
	 * 
	 * This method is a convenience method that can be used to process
	 * a single condition in a depth-first manner. The depth-first
	 * operation occurs because the sub-conditions are processed first
	 * prior to processing the given condition. This  method 
	 * operates in the following manner:
	 * 
	 * <ol>
	 * 	<li>The sub-conditions are processed first via a call to
	 *  {@link #checkCondition(List)} helper method.</li>
	 *  <li>The given condition is processed based on the
	 *  comparator specified in {@link Condition#comp} and
	 *  depending on whether numerical or string values have
	 *  been specified in the condition.</li>
	 *  <li>The resulting condition and result from sub-conditions
	 *  (if any) are combined using the connector for the condition
	 *  to generate the final boolean result.</li> 
	 * </ol>
	 * 
	 * @param cond The condition to be evaluated by this method.
	 * This value must not be null. The variables in the condition
	 * are evaluated using the current values in {@link #decVar}.
	 * 
	 * @return This method returns true if the evaluating the
	 * condition results in true. Otherwise it returns false.
	 */
	public boolean checkCondition(final Condition cond) {
		// Check and process any sub-conditions specified. Sub-conditions
		// are flagged as true if the list if empty.
		final boolean subCondRes = checkCondition(cond.getCondition());
		boolean condRes          = false; // Will be updated below based on comparator
		// Check this current condition.
		final boolean isDef      = decVar.isDefined(cond.getVar());
		// Handle various comparators for the given condition.
		if (ComparatorKind.DEF.equals(cond.getComp())) {
			condRes = isDef;
		} else if (ComparatorKind.NDEF.equals(cond.getComp())) {
			condRes = !isDef;
		} else if (isDef) { // Have variable for further checks.
			// For rest of the comparators we need to use either the numeric
			// or string values.
			final String strValue = decVar.getStrValue(cond.getVar());
			final Double numValue = decVar.getNumValue(cond.getVar());
			// Check if numeric or string comparison is to be used and compute
			// comparison for final checks below.
			final boolean useNum  = (cond.getNumValue() != null) && (numValue != null);
			final int     numComp = (useNum ? numValue.compareTo(cond.getNumValue()) : -1);
			final boolean useStr  = (cond.getStrValue() != null) && (strValue != null);
			final int     strComp = (useStr ? strValue.compareTo(cond.getStrValue()) : -1);
			// If both values are not defined then this condition degrades to false
			if (useNum || useStr) {
				// Compute result for this condition based on comparator
				switch (cond.getComp()) {
				case EQ : condRes = (useNum ? (numComp == 0) : (strComp == 0));
				break;
				case NEQ: condRes = (useNum ? (numComp != 0) : (strComp != 0));
				break;
				case GT : condRes = (useNum ? (numComp > 0)  : (strComp > 0));
				break;
				case GTE: condRes = (useNum ? (numComp >= 0) : (strComp >= 0));
				break;
				case LT : condRes = (useNum ? (numComp < 0)  : (strComp <= 0));
				break;
				case LTE: condRes = (useNum ? (numComp <= 0) : (strComp <= 0));
				}
			}
		}
		// Now combine the results for sub-conditions with this condition
		// if we have sub-conditions to deal with.
		boolean result = condRes; // The final result to be returned.
		if (!cond.getCondition().isEmpty()) {
			switch (cond.getConnector()) {
			case AND: result = condRes && subCondRes;
			break;
			case OR : result = condRes || subCondRes;
			}
		}
		return result;
	}
	
	/**
	 * Method to process conditions in a given list.
	 * 
	 * This method can be used to process conditions associated with a
	 * given parameter. The list of conditions associated with the 
	 * parameter must be passed-in as the parameter. The current
	 * set of variables in the {@link #decVar} is used for validating
	 * the conditions. This method does not add any entries to
	 * {@link #decVar}. This method operates in the following manner:
	 * 
	 * <ol>
	 * <li>If the condition list passed-in is empty then this method
	 * exits immediately returning true.</li>
	 * <li>Each condition in the list is processed in a depth-first
	 * manner via call to {@link #checkCondition(Condition)} method.</li>
	 * <li>The resulting boolean value is combined based on the 
	 * {@link Condition#connector} specified for each condition to
	 * determine the overall boolean result for the list.</li>
	 * </ol>
	 * 
	 * @param condList The list of conditions to be processed by
	 * this method.
	 * 
	 * @return The overall boolean result (true or false) that 
	 * indicates if the condition results in a true or false value.
	 */
	public boolean checkCondition(final List<Condition> condList) {
		if (condList.isEmpty()) {
			// No conditions to be processed. The default return is true.
			return true;
		}
		// Process each condition in the list and accumulate the result
		boolean result = checkCondition(condList.get(0));
		// Process remaining conditions and process results based on connectors.
		ConnectorKind connector = condList.get(0).getConnector();
		for(int i = 1; (i < condList.size()); i++) {
			// Compute boolean result for the condition.
			final boolean condResult = checkCondition(condList.get(i));
			// Compute overall result thus far.
			switch (connector) {
			case AND:
				result = (result && condResult);
				break;
			case OR:
				result = (result || condResult);
			}
			// Update the connector in preparation for next condition (if any)
			connector = condList.get(i).getConnector();
		}
		// Return overall result back to the caller.
		return result;
	}
}
