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

import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.FontMetrics;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JFormattedTextField;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JTextField;
import javax.swing.SpinnerNumberModel;
import javax.swing.UIManager;
import javax.swing.border.Border;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.text.BadLocationException;
import javax.swing.text.NumberFormatter;

import org.peace_tools.core.SummaryWriter;
import org.peace_tools.decagon.jaxb.Parameter;
import org.peace_tools.decagon.jaxb.ParameterKind;
import org.peace_tools.generic.BorderWithComponent;
import org.peace_tools.generic.Utilities;

/**
 * Convenience class to help read parameters from an XML file and generate
 * GUI components for them.
 * 
 * This is a top-level class that forms the primary interface into this
 * package. This class aids the remainder of DECAGON GUI components
 * with the following operations:
 * 
 * <ul>
 * 
 * <li>Unmarshal (load/read) a parameter set from a given XML description
 * file. The XML file must be compatible with the ParameterSet.xsd file.</li>
 * 
 * <li>Create GUI components for the various parameters specified in
 * the parameter set.</li>
 * 
 * </ul>
 */
public class ParameterSetGUIHelper extends ParameterSetHelper {
	/**
	 * A helper class to copy text from a text field to a parameter.
	 * 
	 * This class is an internal helper class that is added to the
	 * JTextField components created by this class. This class
	 * intercepts notifications about changes to the text in the 
	 * GUI-element and updates the value of the corresponding
	 * parameter.
	 *  
	 */
	private class TextParameterUpdater implements DocumentListener {
		/**
		 * The parameter whose value is to be updated whenever the
		 * text associated with this parameter is changed by the user.
		 */
		private final Parameter param;
		
		/**
		 * The only constructor for this class.
		 * 
		 * @param param The parameter whose value should be updated 
		 * whenever a notification about text changes is received by
		 * this class. This parameter cannot be null.
		 */
		TextParameterUpdater(final Parameter param) {
			this.param = param;
		}
		
		/**
		 * Helper method to centralize update of parameter.
		 * 
		 * This is a convenience method that is called from various methods
		 * in this class to appropriately set the parameter value.
		 * 
		 * @param de The document event from where the text is to be 
		 * obtained and parameter is to be updated.
		 */
		private void updateParam(DocumentEvent de) {
			try {
				param.setValue(de.getDocument().getText(0, de.getDocument().getLength()));
			} catch (BadLocationException e) {
				e.printStackTrace();
			}
		}
		
		@Override
		public void insertUpdate(DocumentEvent de) {
			updateParam(de);
		}

		@Override
		public void removeUpdate(DocumentEvent de) {
			updateParam(de);
		}

		@Override
		public void changedUpdate(DocumentEvent de) {
			updateParam(de);
		}
	}

	/**
	 * Convenience method to generate a panel with all the parameters
	 * in it.
	 * 
	 * <p>This is a convenience method that can be used to obtain a JPanel
	 * with all the parameters organized in it. The JPanel returned
	 * by this will use the TwoColumnLayout layout manager to organize
	 * the parameters in a two column format. The returned panel may be 
	 * wrapped in a JScrollPane to handle parameter sets that may 
	 * have a large number of parameters.</p>
	 * 
	 * <p>Note that any conditions associated with parameters are 
	 * evaluated based on the current set of variables defined. As
	 * variables are processed, their values are updated or variables
	 * are undefined based on results from conditions.</p>
	 * 
	 * @param panel An optional panel to which the GUI elements are to
	 * be added. The panel must use a TwoColumnLayout as its layout manager.
	 * If this parameter is null then a JPanel is created by this method
	 * and returned. 
	 * 
	 * @return A JPanel containing the GUI elements that will permit
	 * the user to appropriately update the various parameters.
	 */
	public JPanel createGUI(JPanel panel) {
		if (panel == null) {
			panel = new JPanel(new TwoColumnLayout(7, 10));
		}
		for(Parameter param : paramSet.getParameter()) {
			if (checkCondition(param.getCondition()) && !param.isHidden()) {
				createGUI(param, panel);
			}
		}
		return panel;
	}
	
	/**
	 * A convenience method to create the GUI layout for a given parameter.
	 * 
	 * This is a convenience method that can be used to create the
	 * GUI layout for a given parameter. This method adds the
	 * components corresponding to the given parameter to the 
	 * given panel. This method assumes that the panel uses a 
	 * TwoColumnLayout as its layout manager.
	 * 
	 * @param param The parameter for which the GUI components are to be
	 * created and added to the panel.
	 * 
	 * @param panel The panel to which the GUI components for the given
	 * parameter are to be added. This method assumes that the panel uses a 
	 * TwoColumnLayout as its layout manager.
	 */
	public void createGUI(final Parameter param, final JPanel panel) {
		switch(param.getKind()) {
		case BOOLEAN:
			createBooleanInput(param, panel);
			break;
		case INTEGER:
			createNumericSpinnerInput(param, panel);
			break;
		case DOUBLE:
			createDoubleInput(param, panel);
			break;
		case STRING:
			createStringInput(param, panel);
			break;
		case FILE:
		case DIRECTORY:
			createFileInput(param, panel);
			break;
		case CHOICE:
			createChoiceInput(param, panel);
			break;
		case DUMMY:
			this.addToContainer(panel, TwoColumnLayout.FIRST_COLUMN, param, null);
			break;
		default:
			System.err.println("Unhandled parameter type: " + param.getKind());
		}
	}
	
	/**
	 * Helper method to create GUI components for boolean parameters.
	 * 
	 * This method is invoked from the {@link #createGUI(Parameter, JPanel)}
	 * method to create the GUI components for boolean parameters. This
	 * method performs the following task:
	 * 
	 * <ol>
	 * 
	 * <li>If the parameter is flagged as a controller then it creates a
	 * check-box and adds it as the border component for the given panel.
	 * Checking/un-checking the check-box will cause the components to be
	 * enabled/disabled.</li>
	 * 
	 * <li>If the parameter is not a controller this method creates a 
	 * combo box with "yes" and "no" as options. It then adds the description
	 * and the combo-box to the first and second column of the given panel
	 * via call to {@link #addToContainer(JPanel, String, Parameter, JComponent)}
	 * method.
	 * 
	 * </ol>
	 * 
	 * @param param The boolean parameter for which the GUI components
	 * are to be created by this method. This parameter cannot be null
	 * and must be of kind {@link ParameterKind#BOOLEAN}.
	 * 
	 * @param panel The panel to which the GUI components for the given
	 * parameter are to be added. This method assumes that the panel uses a 
	 * TwoColumnLayout as its layout manager.
	 */
	private void createBooleanInput(final Parameter param, final JPanel panel) {
		if (param.isController()) {
			// Setup a check-box as the border for the panel.
			final JCheckBox controller = new JCheckBox(param.getDescription(), true);
			controller.addActionListener(new ActionListener() {
				@Override
				public void actionPerformed(ActionEvent ae) {
					Utilities.setEnabled(panel, controller.isSelected(), controller);
					param.setValue(Boolean.toString(controller.isSelected()));
				}				
			});
			final Border inputBorder = new BorderWithComponent(controller, panel, 
					BorderFactory.createEtchedBorder());
			if (panel.getBorder() == null) {
				panel.setBorder(BorderFactory.createCompoundBorder(inputBorder,
						BorderFactory.createEmptyBorder(5, 10, 5, 10)));
			} else {
				panel.setBorder(BorderFactory.createCompoundBorder(inputBorder,
						panel.getBorder()));
			}
			return;
		}
		
		
		// Setup action listeners to be used to with radio buttons
		ActionListener al = new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent ae) {
				param.setValue(ae.getActionCommand());
			}
		};
		// Extract the default to be set to streamline code below
		final Boolean boolVal  = getParameterValueAsBoolean(param);
		final boolean defValue = (boolVal != null ? boolVal : false);
		// In other cases we create radio buttons with "Yes" and "No" as
		// the two options.
		JRadioButton yes = new JRadioButton("Yes", defValue ? true : false);
		yes.setActionCommand("true");
		yes.addActionListener(al);
		JRadioButton no  = new JRadioButton("No", !defValue ? true : false);
		yes.setActionCommand("false");
		yes.addActionListener(al);
		// Create button group to ensure only Yes or No is selected
		ButtonGroup  group = new ButtonGroup();
		group.add(yes);
		group.add(no);
		// Layout the radio buttons.
		JPanel yesNoPanel = new JPanel(new FlowLayout(FlowLayout.LEFT, 5, 0));
		yesNoPanel.add(yes);
		yesNoPanel.add(no);
		// Add input to top-level container and return back to caller
		addToContainer(panel, TwoColumnLayout.SECOND_COLUMN, param, yesNoPanel);
	}
	
	/**
	 * Helper method to create GUI components for integer or double parameters.
	 * 
	 * This method is invoked from the {@link #createGUI(Parameter, JPanel)}
	 * method to create the GUI components for integer  
	 * ({@link ParameterKind#INTEGER}) parameters. This method creates
	 * a JSpinner to obtain integer inputs. The size of the JSpinner is
	 * set to display at least 10 digits. The spinner is initialized to
	 * the default value (if any) specified for the parameter. 
	 * 
	 * It then adds the description and the spinner to the first 
	 * and second column of the given panel via call to
	 *  {@link #addToContainer(JPanel, String, Parameter, JComponent)}
	 * method.
	 * 
	 * <p><b>NOTE:</b>This method is invoked from {@link #createDoubleInput(Parameter, JPanel)}
	 * if the double value has a minimum or maximum value set.</p>
	 * 
	 * @param param The integer parameter for which the GUI components
	 * are to be created by this method. This parameter cannot be null
	 * and must be of kind {@link ParameterKind#INTEGER}.
	 * 
	 * @param panel The panel to which the GUI components for the given
	 * parameter are to be added. This method assumes that the panel uses a 
	 * TwoColumnLayout as its layout manager.
	 */
	private void createNumericSpinnerInput(final Parameter param, final JPanel panel) {
		SpinnerNumberModel model = new SpinnerNumberModel();
		// Setup min and max values if specified.
		if (param.getMax() != null) {
			model.setMaximum(param.getMax());
		}
		if (param.getMin() != null) {
			model.setMinimum(param.getMin());
		}
		// Setup default value if specified.
		final Number numVal = getParameterValueAsNumber(param);
		if (numVal != null) {
			model.setValue(numVal);
			if (ParameterKind.DOUBLE.equals(param.getKind())) {
				model.setStepSize(new Double(0.1));
			}
		}
		// Create the spinner
		final JSpinner spinner = new JSpinner(model);
		spinner.addChangeListener(new ChangeListener() {
			@Override
			public void stateChanged(ChangeEvent ce) {
				param.setValue(spinner.getValue().toString());
			}
		});
		spinner.setName(param.getCmdLine());
		// Setup preferred size of the spinner to display at least 10 digits
		final FontMetrics fm    = spinner.getFontMetrics(spinner.getFont());
		final Dimension prefDim = spinner.getPreferredSize(); // preferred dimension
		prefDim.width           = fm.stringWidth("000000000000"); // preferred width
		spinner.setPreferredSize(prefDim);
		spinner.setMaximumSize(prefDim);
		// Add the description and input component to the panel.
		addToContainer(panel, TwoColumnLayout.SECOND_COLUMN, param, spinner);
	}
	
	/**
	 * Helper method to create GUI components for double parameters.
	 * 
	 * This method is invoked from the {@link #createGUI(Parameter, JPanel)}
	 * method to create the GUI components for double  
	 * ({@link ParameterKind#DOUBLE}) parameters. This method creates
	 * a JFormattedTextField to obtain double inputs. The size of the 
	 * input field is set to display at least 10 digits. The text field is
	 * initialized to the default value (if any) specified for the parameter. 
	 * 
	 * It then adds the description and the text field to the first 
	 * and second column of the given panel via call to
	 *  {@link #addToContainer(JPanel, String, Parameter, JComponent)}
	 * method.
	 * 
	 * @param param The double parameter for which the GUI components
	 * are to be created by this method. This parameter cannot be null
	 * and must be of kind {@link ParameterKind#DOUBLE}.
	 * 
	 * @param panel The panel to which the GUI components for the given
	 * parameter are to be added. This method assumes that the panel uses a 
	 * TwoColumnLayout as its layout manager.
	 */
	private void createDoubleInput(final Parameter param, final JPanel panel) {
		// If we have a minimum or maximum value then prefer a spinner
		if ((param.getMax() != null) || (param.getMin() != null)) {
			createNumericSpinnerInput(param, panel);
			return;
		}
		// Create a formatted text input to enter only numbers.
		NumberFormatter formatter = new NumberFormatter();
		JFormattedTextField input = new JFormattedTextField(formatter);
		input.setName(param.getCmdLine());
		input.setColumns(10);
		// Add two pixels to height to make the field look nicer and
		// lock its maximum size to make layouts look nicer.
		Dimension prefDim = input.getPreferredSize();
		prefDim.height   += 2;
		input.setMaximumSize(prefDim);
		// Setup default value (if specified).
		final Number numVal = getParameterValueAsNumber(param);
		if (numVal != null) {
			input.setValue(numVal);
		}
		input.addPropertyChangeListener("value", new PropertyChangeListener() {
			@Override
			public void propertyChange(PropertyChangeEvent pce) {
				if (pce.getNewValue() != null) {
					param.setValue(pce.getNewValue().toString());
				}
			}			
		});
		// Add input in top-level container and return back to caller
		addToContainer(panel, TwoColumnLayout.SECOND_COLUMN, param, input);
	}

	/**
	 * Helper method to create GUI components for string parameters.
	 * 
	 * This method is invoked from the {@link #createGUI(Parameter, JPanel)}
	 * method to create the GUI components for string  
	 * ({@link ParameterKind#STRING}) parameters. This method creates
	 * a JTextField to obtain double inputs. The size of the 
	 * input field is set to display at least 10 characters. The 
	 * text in the text field is initialized to the default value
	 * (if any) specified for the parameter. 
	 * 
	 * It then adds the description and the text field to two
	 * consecutive rows such that they occupy the whole row. 
	 * 
	 * @param param The string parameter for which the GUI components
	 * are to be created by this method. This parameter cannot be null
	 * and must be of kind {@link ParameterKind#STRING}.
	 * 
	 * @param panel The panel to which the GUI components for the given
	 * parameter are to be added. This method assumes that the panel uses a 
	 * TwoColumnLayout as its layout manager.
	 */
	private void createStringInput(final Parameter param, final JPanel panel) {
		final JTextField input = new JTextField(10);
		if (param.getValue() != null) {
			input.setText(processVariables(param.getValue()));	
		}
		// Fix maximum height of the text field to handle resizing nicely
		Dimension prefDim = input.getPreferredSize();
		prefDim.height   += 2; // Add couple of pixels to make input look nice
		Dimension maxDim  = input.getMaximumSize();
		maxDim.height     = prefDim.height; // Fix maximum height.
		input.setMaximumSize(maxDim); // Set the dimensions
		input.setPreferredSize(prefDim);
		input.getDocument().addDocumentListener(new TextParameterUpdater(param));
		// Add input in top-level container and return back to caller
		addToContainer(panel, TwoColumnLayout.FIRST_COLUMN, param, input);
	}
	
	/**
	 * Helper method to create GUI components for file or directory
	 * parameters.
	 * 
	 * This method is invoked from the {@link #createGUI(Parameter, JPanel)}
	 * method to create the GUI components for file 
	 * ({@link ParameterKind#FILE}) or directory 
	 * ({@link ParameterKind#DIRECTORY}) parameters. This method creates
	 * a JTextField and a "Browse" button that the user can use to 
	 * provide the necessary information. The 
	 * text in the text field is initialized to the default file name
	 * or directory (if any) specified for the parameter. 
	 * 
	 * It then adds the description and the text field to two
	 * consecutive rows such that they occupy the whole row. This
	 * is of course accomplished via call to 
	 * {@link #addToContainer(JPanel, String, Parameter, JComponent)} method. 
	 * 
	 * @param param The string parameter for which the GUI components
	 * are to be created by this method. This parameter cannot be null
	 * and must be of kind {@link ParameterKind#STRING}.
	 * 
	 * @param panel The panel to which the GUI components for the given
	 * parameter are to be added. This method assumes that the panel uses a 
	 * TwoColumnLayout as its layout manager.
	 */
	private void createFileInput(final Parameter param, final JPanel panel) {
		// Create the text box into which the user can enter a file
		// or directory path.
		final JTextField path = new JTextField(10);
		if (param.getValue() != null) {
			path.setText(processVariables(param.getValue()));	
		}
		// Create the browse button to go with this input.
		final JButton browse = new JButton("Browse...");
		browse.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent ae) {
				final boolean directory  = ParameterKind.BOOLEAN.equals(param.getKind());
				JFileChooser jfc = new JFileChooser();
				jfc.setFileSelectionMode((directory ? JFileChooser.DIRECTORIES_ONLY : JFileChooser.FILES_ONLY));
				if (jfc.showOpenDialog(path) == JFileChooser.APPROVE_OPTION) {
					path.setText(jfc.getSelectedFile().getAbsolutePath());
				}
			}			
		});
		// Create top-level container to hold the path and browse buttons.
		JPanel container = new JPanel();
		container.setLayout(new BoxLayout(container, BoxLayout.X_AXIS));
		container.add(path);
		container.add(Box.createHorizontalStrut(5));
		container.add(browse);
		// Add compound container in top-level panel
		addToContainer(panel, TwoColumnLayout.FIRST_COLUMN, param, container);
	}

	/**
	 * Helper method to create GUI components for choice parameters.
	 * 
	 * This method is invoked from the {@link #createGUI(Parameter, JPanel)}
	 * method to create the GUI components for integer  
	 * ({@link ParameterKind#CHOICE}) parameters. This method creates
	 * a JComboBox to permit user to select appropriate choice.
	 * 
	 * It then adds the description and the combo-box to the first 
	 * and second column of the given panel via call to
	 *  {@link #addToContainer(JPanel, String, Parameter, JComponent)}
	 * method.
	 * 
	 * @param param The choice parameter for which the GUI components
	 * are to be created by this method. This parameter cannot be null
	 * and must be of kind {@link ParameterKind#CHOICE}.
	 * 
	 * @param panel The panel to which the GUI components for the given
	 * parameter are to be added. This method assumes that the panel uses a 
	 * TwoColumnLayout as its layout manager.
	 */
	private void createChoiceInput(final Parameter param, final JPanel panel) {
		// Extract the comma-separated list of choices.
		final String actualVal = processVariables(param.getValue());
		String[] nameValueList = actualVal.split("(:|,)");
		// Create separate lists for name and value to ease combo-box operations
		final String choiceList[] = new String[nameValueList.length / 2];
		final String valueList[]  = new String[nameValueList.length / 2];
		for(int srcIdx = 0, destIdx = 0; (destIdx < choiceList.length); srcIdx += 2, destIdx++) {
			choiceList[destIdx] = nameValueList[srcIdx];
			valueList [destIdx] = nameValueList[srcIdx + 1];
		}
		// Create the combo-box
		final JComboBox input = new JComboBox(choiceList);
		input.setMaximumSize(input.getPreferredSize());
		// Setup initial value.
		param.setValue(valueList[0]);
		// Setup action listener to handle user's choices.
		input.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent ae) {
				String value = valueList[input.getSelectedIndex()];
				param.setValue(value);
			}
		});
		// Add the description and input component to the panel.
		addToContainer(panel, TwoColumnLayout.SECOND_COLUMN, param, input);
	}
	
	/**
	 * Helper method to add description and input component to a given panel.
	 * 
	 * This helper method is invoked from various methods in this class
	 * to add the description (for a given parameter) and a given GUI
	 * input component to a given panel. This method is used to:
	 * 
	 * <ol>
	 * 
	 * <li>Add description to first column and input component to second
	 * column by specifying the inputColumn parameter to 
	 * {@link TwoColumnLayout#SECOND_COLUMN}.</li>
	 * 
	 * <li>Add description to occupy a full row and input component to
	 * occupy a full second row by specifying the inputColumn parameter to
	 * {@link TwoColumnLayout#FIRST_COLUMN}.
	 * 
	 * </ol>
	 * 
	 * @param container The panel to which the description (in from of a
	 * JLabel) and the input are to be added. This method assumes that
	 * the JPanel uses a TwoColumnLayout as its layout manager.
	 * 
	 * @param inputColumn The column to which the input component is to be
	 * added. The description is always added to 
	 * {@link TwoColumnLayout#FIRST_COLUMN}.
	 * 
	 * @param param The parameter whose description is to be added to the
	 * given panel.
	 * 
	 * @param input The GUI component to be added to the panel. This value
	 * can be null.
	 */
	private void addToContainer(final JPanel container, final String inputColumn, 
			final Parameter param, final JComponent input) {
		// Create a JLabel that contains the description for the parameter
		final String desc  = "<html>" + param.getDescription().trim() + "</html>";
		JLabel description = new JLabel(desc);
		description.setToolTipText(param.getCmdLine());
		// Add to the top-level container to layout description and component.
		container.add(description, TwoColumnLayout.FIRST_COLUMN);
		if (input != null) {
			input.setToolTipText(desc);
			container.add(input, inputColumn);
		}
	}
	
	/**
	 * Convenience method to summary information about parameters.
	 * 
	 * This method can be used to summarize information about the 
	 * various parameters currently defined in this object.
	 * Currently, this method simply lists all parameters and their
	 * values in the same order in which they are listed in the 
	 * DADX file.
	 * 
	 * @param sw The summary writer to be used to generate the summary
	 * of the variables.
	 * 
	 * @param subSummary If this flag is true then the summary information
	 * is generated using {@link SummaryWriter#addSubSummary(String, String, String)}
	 * method. If this flag is false, then the summary information is
	 * generated using {@link SummaryWriter#addSummary(String, String, String)}
	 * method.
	 */
	public void summarize(SummaryWriter sw, boolean subSummary) {
		for(Parameter param: getParameters()) {
			if (param.isHidden() || param.isController()) {
				// Skip this parameter.
				continue;
			}
			// Get the description to be displayed while giving 
			// preference to summary description (if specified).
			String description = param.getSummary();
			if ((description == null) || (description.isEmpty())) {
				description = param.getDescription();
			}
			if (subSummary) {
				sw.addSubSummary(description, param.getValue(), "");
			} else {
				sw.addSummary(description, param.getValue(), "");
			}
		}
	}
	
	
	public static void main(String args[]) throws Exception {
		/*
		ParameterSet ps = new ParameterSet();
		ps.setDescription("Test description");
		Parameter param = new Parameter();
		ps.getParameter().add(param);
		
		JAXBContext context = JAXBContext.newInstance(ParameterSet.class);
		Marshaller m = context.createMarshaller();
		m.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, Boolean.TRUE);
		
		ObjectFactory of = new ObjectFactory();

		Writer w = new FileWriter("SimpleTest.xml");
		m.marshal(of.createParameterSet(ps), w);
		*/
		// Turn off metal's use of bold fonts
		UIManager.put("swing.boldMetal", Boolean.FALSE);
		ParameterSetGUIHelper psgh = new ParameterSetGUIHelper();
		psgh.unmarshal("/decagon/MetaSimIllumina.xml");
		System.out.println("Command-line = " + psgh.getCmdLine());
		JFrame frame = new JFrame("Test");
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		JPanel inputPanel = psgh.createGUI(null);
		inputPanel.invalidate();
		JScrollPane jsp = new JScrollPane(inputPanel);
		jsp.setBorder(BorderFactory.createEmptyBorder(5, 10, 5, 10));
		frame.add(jsp);
		frame.pack();
		//Utilities.setEnabled(inputPanel, false);
		frame.setVisible(true);		
	}
}
