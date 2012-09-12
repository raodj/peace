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
import java.awt.FlowLayout;
import java.awt.GridLayout;
import java.awt.Rectangle;
import java.awt.Window;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.List;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JScrollPane;
import javax.swing.SwingUtilities;
import javax.swing.SwingWorker;
import javax.swing.border.EmptyBorder;

/**
 * A convenience class for running background tasks with updates in a 
 * a custom foreground dialog box.
 * 
 * This is a convenience class that has been designed to run lengthy
 * GUI-interacting tasks in a dedicated background thread. The background
 * operations of this class are achieved via the SwingWorker class. 
 */
public class BackgroundTask extends SwingWorker<Throwable, String> 
implements PropertyChangeListener, ActionListener {

	/**
	 * The interface to be implemented by the user to perform the actual
	 * task in the background. 
	 */
	public interface UserTask {
		/**
		 * This method is invoked from a background swing thread to
		 * permit a long running task to operate.
		 *  
		 * @param bTask The helper class that displays the dialog and
		 * updates various GUI components based on the progress reported
		 * by this method.
		 * 
		 * @exception Exception This method can throw an exception on errors.
		 */
		public void run(BackgroundTask bTask) throws Exception;

		/**
		 * This method is invoked from event dispatch thread to
		 * indicate background task is done.
		 *  
		 * @param bTask The helper class that displays the dialog and
		 * updates various GUI components based on the progress reported
		 * by this method.
		 */
		public void done(BackgroundTask bTask);
		
		/**
		 * This method is invoked after all operations have completed and
		 * the dialog (displayed by this class) has been closed.
		 * 
		 * This is the last method that is invoked on the user task once
		 * the dialog (displayed by this class) has been closed. This method
		 * can be used by derived classes to perform any finalization 
		 * operations.
		 * 
		 * @param bTask The background task that displays the dialog (and
		 * updates various GUI) that is invoking this method. 
		 */
		public void dialogClosed(BackgroundTask bTask);
		
	}

	/**
	 * The constructor.
	 * 
	 * The constructor initializes the various GUI components to their
	 * default initial value and lay's out the various GUI components
	 * in a dialog to be displayed by this method.
	 * 
	 * @param task The user task to be run in the background. This parameter
	 * cannot be null.
	 *  
	 * @param parent The parent GUI component which is creating this
	 * background task.
	 *  
	 * @param dialogTitle The title for the dialog that displays the 
	 * operations being performed by the background thread.
	 * 
	 * @param infoMessage The informational message providing summary 
	 * information about the operations being performed. This can be 
	 * simple JLabel or a more complex JPanel containing various GUI
	 * components.
	 * 
	 * @param mainProgBar If this flag is true, then a main progress bar
	 * is displayed. The number of steps in the progress bar matches
	 * the number of sub-steps specified.
	 * 
	 * @param subSteps An array indicating the number of sub-steps to be
	 * performed by this class. This array cannot be null but can be an
	 * empty array in which case no sub-steps are displayed and the 
	 * main-progress bar (if present) remains in indeterminate mode.
	 * 
	 * @param showStepsInProgBar If this flag is set to true then the
	 * current step being processed is shown in the the main progress
	 * bar (in addition to being highlighted). 
	 * 
	 * @param scrollSteps If this flag is true and simpeStepStyle is
	 * false, then the sub-steps are placed within a scroll pane. This
	 * is handy when several (more than 6) of sub-steps are to be displayed.
	 * 
	 * @param simpleStepStyle If this flag is true then a simple style 
	 * for showing one step at a time is used. The current step being
	 * processed is shown under the main message. If this flag is 
	 * false then all steps are shown in a vertical list and as 
	 * progress is made each step get's a check mark before it indicating
	 * its successful completion.
	 * 
	 * @param numCols If the simpleStepStyle flag is false then 
	 * individual steps are shown in a grid. This value determines
	 * the number of columns in the grid. This value must be at least 1.
	 *  
	 * @param showLogs If this flag is true then a JTextArea with 
	 * log messages is included in the dialog created by this class.
	 * 
	 * @param showButtons If this flag is true then "OK" and optionally a
	 * "Cancel" button are displayed in the dialog box shown by this
	 * class. If this flag is false then the buttons are not displayed at all.
	 * 
	 * @param enableCancelButton If this flag is true then a "Cancel" button is
	 * enabled. This is meaningful only if the showButtons options is true.
	 * 
	 * @param parameter A convenience place holder that can be used to 
	 * pass parameters to the background process.
	 */
	public BackgroundTask(UserTask task, JComponent infoMessage, boolean mainProgBar, 
			String[] subSteps, boolean showStepsInProgBar, boolean scrollSteps,
			boolean simpleStepStyle, int numCols, boolean showLogs, 
			boolean showButtons, boolean enableCancelButton, Object parameter) {
		// Save the task.
		this.task      = task;
		this.parameter = parameter;
		this.result    = null; // initial no exceptions
		// Create the various components in a panel to be used in the dialog.
		topProgBar = new JProgressBar();
		topProgBar.setIndeterminate(true);
		topProgBar.setString("Starting...");
		topProgBar.setMaximum(subSteps.length);
		if (showStepsInProgBar) {
			topProgBar.setStringPainted(true);
		}
		// Create labels for the sub-steps
		final JPanel stepPanel = createSubSetps(subSteps, simpleStepStyle, numCols);
		// Create scrolling log area displayed in this panel.
		logs = new ProcessOutputDisplay();
		JScrollPane jsp = new JScrollPane(logs.getTextPane());
		// Create the "OK" and "Cancel" buttons
		okButton = Utilities.createButton("images/16x16/Ok.png", "Close", 
				"ok", this, "Close the dialog box", false);
		cancelButton = Utilities.createButton("images/16x16/Cancel.png", "Cancel", 
				"cancel", this, "Cancel background operations if possible", enableCancelButton);
		// Create the panel with "OK" and "Cancel" buttons at the bottom.
		JPanel buttonPanel = new JPanel(new FlowLayout(FlowLayout.CENTER, 10, 3));
		if (showButtons) {
			// Add buttons to display them appropriately
			buttonPanel.add(okButton);
			buttonPanel.add(cancelButton);
		}

		// Now create the top-level panel to hold the various components
		// we have created.
		topPanel = new JPanel();
		topPanel.setLayout(new BoxLayout(topPanel, BoxLayout.Y_AXIS));
		//topPanel.add(Box.createVerticalStrut(5));
		topPanel.add(infoMessage); 
		// For simple step-style the sub-steps go here.
		if (simpleStepStyle && (subSteps.length > 0)) {
			topPanel.add(Box.createVerticalStrut(5));
			topPanel.add(stepPanel);
		}
		if (mainProgBar) {
			topPanel.add(Box.createVerticalStrut(5));
			topPanel.add(topProgBar);
		}
		if (!simpleStepStyle && (subSteps.length > 0)) {
			topPanel.add(Box.createVerticalStrut(5));
			if (scrollSteps) {
				// Create a scroll pane to scroll the sub-steps
				subStepScrollPane = new JScrollPane(stepPanel);
				subStepScrollPane.setBorder(null);
				topPanel.add(subStepScrollPane);
			} else {
				// No scroll pane to be used in this situation
				subStepScrollPane = null;
				topPanel.add(stepPanel);
			}
		} else {
			// No scroll pane to be used in this situation
			subStepScrollPane = null;
		}
		if (showLogs) {
			topPanel.add(Box.createVerticalStrut(8));
			topPanel.add(new JLabel("<html><b>Logs:</b></html>", 
					Utilities.getIcon("images/16x16/UserLog.png"), JLabel.LEFT));
			topPanel.add(Box.createVerticalStrut(2));
			topPanel.add(jsp);
		}
		if (showButtons) {
			topPanel.add(Box.createVerticalStrut(8));
			topPanel.add(buttonPanel);
		}
		topPanel.setBorder(new EmptyBorder(10, 5, 5, 10));
		// To make layout look good, make all components have left alignment.
		for(Component child: topPanel.getComponents()) {
			if (child instanceof JComponent) {
				JComponent jc = (JComponent) child;
				jc.setAlignmentX(JComponent.LEFT_ALIGNMENT);
			}
		}
		
		// Setup property listener on parent class
		addPropertyChangeListener(this);
	}

	/**
	 * Helper method to create a list of sub-steps and place them in a top-level
	 * panel.
	 * 
	 * This is a helper method that was introduced to streamline the code in
	 * the constructor. This method creates a label for each sub-step 
	 * associated with the background task and adds them to a top-level panel.
	 * It returns the top-level panel from this method.
	 * 
	 * @param subSteps The array of sub-steps passed-in by the user for which
	 * corresponding labels are to be created.
	 * 
	 * @param simpleStyle Flag to indicate if individual steps are to be
	 * shown or just the current step is to be shown. If this parameter is
	 * true then only the current step is shown (it is called simple style).
	 * 
	 * @param numCols If simpleStyle is false, then individual steps are
	 * shown in a grid layout. This parameter determines the number of
	 * columns in the grid layout.
	 * 
	 * @return A JPanel containing the labels (if any) that were created
	 * by this method. Note that this method always returns a panel.
	 */
	private JPanel createSubSetps(final String[] subSteps, 
			final boolean simpleStyle, int numCols) {
		subStepLabels    = new JLabel[subSteps.length];
		JPanel stepPanel = new JPanel();
		
		if ((subSteps.length > 0) && (!simpleStyle)) {
			final int numRows = (subSteps.length / numCols) + 
				((subSteps.length % numCols > 0) ? 1 : 0);
			stepPanel.setLayout(new GridLayout(numRows, numCols, 5, 5));
			stepPanel.setBorder(new EmptyBorder(1, 10, 1, 10));
		} else {
			// Simple style.
			simpleStepLabel = new JLabel(" ");
			stepPanel.setLayout(new BorderLayout(0, 0));
			stepPanel.add(simpleStepLabel, BorderLayout.CENTER);
		}
		
		for(int step = 0; (step < subSteps.length); step++) {
			subStepLabels[step] = new JLabel(subSteps[step],
				Utilities.getIcon("images/16x16/Box.png"), JLabel.LEFT);
			subStepLabels[step].setVerticalAlignment(JLabel.TOP);
			subStepLabels[step].setOpaque(true);
			if (!simpleStyle) {
				stepPanel.add(subStepLabels[step]);
			}
		}
		return stepPanel;
	}

	/**
	 * Obtain the parameter value set for this class.
	 * 
	 * This method returns a convenience place-holder object that can
	 * be used to pass-in a value to the background thread.
	 * 
	 * @return The parameter (if any) set when this class was created.
	 * This method can return null if the parameter passed-in was null.
	 */
	public Object getParameter() {
		return parameter;
	}

	/**
	 * Cut an entry in the log displayed by this class.
	 * 
	 * This method can be used to cut a log in the log panel displayed
	 * in this class.
	 * 
	 * @param logEntry The log entry to be appended to the logs. This
	 * parameter cannot be null.
	 */
	public void log(final String logEntry) {
		publish(logEntry);
	}
	
	/**
	 * Cut an entry in the log displayed by this class with a specific
	 * style.
	 * 
	 * This method can be used to cut a log in the log panel displayed
	 * in this class.
	 * 
	 * @param logEntry The log entry to be appended to the logs. This
	 * parameter cannot be null.
	 * 
	 * @param style The style to be used for this log entry.
	 */
	public void log(final String logEntry, ProcessOutputDisplay.StyleKind style) {
		logs.addLog(logEntry, style);
	}

	/**
	 * Obtain the helper component used to display logs.
	 * 
	 * @return The helper class that is used to display the
	 * logs.
	 */
	public ProcessOutputDisplay getLog() {
		return logs;
	}
	
	/**
	 * Method to display dialog and start background process.
	 * 
	 * This method displays the background dialog box and starts the 
	 * background operations.
	 * 
	 * @param showDialog If this flag is true then this method creates
	 * a modal dialog to display the {@link #topPanel}. If this parameter
	 * is false, then a modal dialog is not created and it is the caller's
	 * responsibility to ensure that the {@link #topPanel} has been
	 * suitably displayed in a different dialog prior to calling this
	 * method.
	 * 
	 * @param parent The parent component to be used for creating
	 * the modal dialog box if showDialog parameter is true.
	 * 
	 * @param title The title for the dialog.
	 * 
	 * @param icon The icon to be used for the dialog title bar.
	 */
	public void start(final boolean showDialog, final Component parent, 
			final String title, final ImageIcon icon) {
		if (showDialog) {
			// Finally, create our top-level dialog and set various properties.
			Window owner = (parent instanceof Window) ? (Window) parent : SwingUtilities.getWindowAncestor(parent);
			dialog = new JDialog(owner, title);
			dialog.setIconImage(icon.getImage());
			dialog.setContentPane(topPanel);
			dialog.setDefaultCloseOperation(JDialog.DO_NOTHING_ON_CLOSE);
			dialog.setResizable(false);
			dialog.setModal(true);
			dialog.setAlwaysOnTop(true);
			dialog.pack();
			dialog.setLocationRelativeTo(dialog.getParent());
		}
		// Start background job.
		execute();
		// The modal dialog should be made visible last!
		if (showDialog) {
			dialog.setVisible(true);
		}
	}

	/**
	 * Indicate a step of progress has been done.
	 * 
	 * This method is a convenience method that can be used to indicate
	 * that a step of the operation has successfully completed.
	 */
	public void updateProgress() {
		setProgress(getProgress() + 1);
	}

	/**
	 * The resulting exception (if any) at the end of running the 
	 * background task. This value is meaningful only after
	 * the background task is completed. The same value can 
	 * also be obtained via call to {@link #get()} method.
	 * 
	 * @return The value returned my this method is meaningful 
	 * only after the background task has completed. This method
	 * is convenience method to be used in {@link UserTask#done(BackgroundTask)}
	 * method only.
	 */
	public Throwable getResult() {
		return result;
	}

	/**
	 * Convenience method to display a message to the user.
	 * 
	 * This is a convenience method that can be used to display either a
	 * success message or a failure message. The choice is made 
	 * based on the value of {@link #result}. If {@link #result} is null,
	 * then successMsg is displayed. Otherwise the error message is
	 * displayed along with the exception details.
	 * 
	 * <p>NOTE: The use of this method is meaningful only after the
	 * background task has completed its operations.</p>
	 * 
	 * @param successMsg The success message to be displayed to the user
	 * if {@link #result} is null.
	 * 
	 * @param errorMsg The error message to be displayed to the user
	 * if {@link #result} is not null.
	 */
	public void showMessage(final String successMsg, final String errorMsg) {
		// Display the exception if one occurred.
		JComponent msgPanel = null;
		if (result == null) {
			msgPanel = new JLabel(successMsg);
		} else {
			final String htmlErrMsg = "<html>" + errorMsg + "<br>"
				+ "<dl><dt><b>Error:</b></dt><dd>" + Utilities.wrapStringToHTML(result.getMessage(), 65)
				+ "</dd></dl></html>";
			msgPanel = Utilities.collapsedMessage(htmlErrMsg, Utilities.toString(result));
		}
		// Display the message
		JOptionPane.showMessageDialog(dialog, msgPanel, "Result",
				(result == null) ? JOptionPane.INFORMATION_MESSAGE : JOptionPane.ERROR_MESSAGE);
	}
	
	@Override
	protected Throwable doInBackground() throws Exception {
		try {
			// Let derived class do its actual job.
			task.run(this);
		} catch (Throwable err) {
			result = err;
			// Log the error for convenience.
			logs.addLog("\n", ProcessOutputDisplay.StyleKind.STDOUT);
			logs.addLog(Utilities.toString(err), ProcessOutputDisplay.StyleKind.STDERR);
		}
		return result;
	}

	@Override
	protected void process(List<String> logEntries) {
		for(String entry: logEntries) {
			logs.addLog(entry, ProcessOutputDisplay.StyleKind.STDOUT);
		}
	}

	@Override
	public void propertyChange(PropertyChangeEvent pce) {
		final String prop = pce.getPropertyName();
		int stepNum       = -1;

		if ("state".equals(prop)) {
			if (StateValue.STARTED.equals(pce.getNewValue())) {
				if (subStepLabels.length > 0) {
					stepNum = 0;
				}
			} else if (StateValue.DONE.equals(pce.getNewValue())) {
				okButton.setEnabled(true);
				task.done(this); // Let derived class do more operations
			}
		} else if ("progress".equals(prop)) {
			stepNum = (Integer) pce.getNewValue();
		}
		// Update progress bar settings.
		if (stepNum != -1) {
			topProgBar.setIndeterminate(false);
			topProgBar.setValue(stepNum);
			if (stepNum < subStepLabels.length) {
				topProgBar.setString(subStepLabels[stepNum].getText());
				subStepLabels[stepNum].setIcon(Utilities.getIcon("images/16x16/RArrow.png"));
				if (subStepScrollPane != null) {
					Rectangle rect = subStepLabels[stepNum].getBounds();
					subStepLabels[stepNum].scrollRectToVisible(rect);
				}
			} else {
				topProgBar.setString("Done!");
			}
		}
		// Update status icons associated with the previous steps (if any)
		if ((subStepLabels != null) && (subStepLabels.length > 0)) {
			for(int step = 0; (step < stepNum); step++) {
				subStepLabels[step].setIcon(Utilities.getIcon("images/16x16/CheckedBox.png"));
			}
		}
		// Update simple label with current step in progress
		if (simpleStepLabel != null) {
			if ((stepNum >= 0) && (stepNum < subStepLabels.length)) { 
				simpleStepLabel.setText(subStepLabels[stepNum].getText());
			}
		}
	}

	/**
	 * Helper method to dispose the dialog.
	 * 
	 * <p>This is a convenience method that can be used to dispose the dialog
	 * displayed by this class. This method is typically used when this
	 * dialog is created without any buttons in it. This method is 
	 * typically used from the {@link UserTask#done(BackgroundTask)} method.</p>
	 * 
	 * <b>NOTE</b>: This method calls {@link UserTask#dialogClosed(BackgroundTask)}
	 * method.
	 */
	public void disposeDialog() {
		dialog.setVisible(false);
		dialog.dispose();
		// Let the caller know this class is all done.
		task.dialogClosed(this);
	}
	
	@Override
	public void actionPerformed(ActionEvent ae) {
		final String cmd = ae.getActionCommand();
		if ("ok".equals(cmd)) {
			disposeDialog();
		} else if ("cancel".equals(cmd)) {
			// Suggest to base class to cancel operation if possible.
			cancel(false);
		}
	}

	/**
	 * Obtain the dialog being used by this component.
	 * 
	 * @return The dialog being used by this component to display
	 * progress information about the task running in the background.
	 * If a dialog has not been created then this method returns
	 * null.
	 */
	public JDialog getDialog() {
		return dialog;
	}
	
	/**
	 * Obtain the top-level with various GUI components.
	 * 
	 * @return The top-panel with various GUI components.
	 */
	public JPanel getTopPanel() {
		return topPanel;
	}
	
	/**
	 * Get the scroll pane that is used to display sub-steps.
	 * 
	 * This method can be used to obtain the scroll pane that is being used
	 * to display a long list of sub-steps. A scroll pane is used only
	 * if this object has been so configured when the object was
	 * instantiated. 
	 * 
	 * @return The scroll pane that is used to display the list of sub-steps
	 * in {@link #subStepLabels}. This method can return null if a scroll pane
	 * is not being used. 
	 */
	public JScrollPane getSubStepScrollPane() {
		return subStepScrollPane;
	}
	
	/**
	 * The top-level dialog box displayed by this component. The dialog
	 * box is created in the {@link #start()} method if the user
	 * requests it. Otherwise this value is null (in which case it is
	 * the caller's responsibility to have the topPanel appropriately
	 * displayed in another modal dialog).
	 */
	private JDialog dialog;

	/**
	 * The top-level panel that is used to hold the various components
	 * created by this class. The {@link #dialog} has this panel
	 * as its content pane.
	 */
	private final JPanel topPanel;
	
	/**
	 * The top-level progress bar that is displayed just below the top-level
	 * label in this class.
	 */
	private final JProgressBar topProgBar;

	/**
	 * The labels that indicate the sub-steps in the dialog box and
	 * the current sub-step being performed. This list is created
	 * immaterial of whether the  user configures this class to display
	 * each step separately in the constructor. However, it is not
	 * always displayed (in cases where the {@link #simpleStepLabel} is used to
	 * display one step at a time.
	 */
	JLabel[] subStepLabels;

	/**
	 * The logs that are displayed in the the dialog box created 
	 * by this class. The logs are based on a styled document that
	 * permits different types of logs (that are highlighted in 
	 * slightly different ways) to be displayed.
	 */
	final ProcessOutputDisplay logs;

	/**
	 * The OK button that is displayed at the bottom of this dialog box.
	 */
	final JButton okButton;

	/**
	 * The Cancel button that is displayed at the bottom of this dialog box.
	 */
	final JButton cancelButton;

	/**
	 * A label that is used when multiple steps are shown one at a time.
	 * This label is used only if the user configures to display one
	 * step at a time. Otherwise this object is null and the {@link #subStepLabels}
	 * are used to display status of all steps.
	 */
	 JLabel simpleStepLabel;
	
	/**
	 * A convenience place holder that can be used to pass in any
	 * parameters to the background thread.
	 */
	private final Object parameter;

	/**
	 * The resulting exception (if any) at the end of running the 
	 * background task. This value is meaningful only after
	 * the background task is completed. The same value can 
	 * also be obtained via call to {@link #get()} method. This value
	 * is set in the {@link #doInBackground()} method.
	 */
	private Throwable result;
	
	/**
	 * The actual task to be performed by this method. The task is 
	 * essentially a Runnable object whose {@link Runnable#run()} method
	 * is invoked by this class in a background thread.
	 */
	private final UserTask task;
	
	
	/**
	 * A scroll pane that can contains sub-steps if the user configured
	 * to use the a scroll pane for labels. This scroll pane is
	 * created in the constructor.
	 */
	private final JScrollPane subStepScrollPane;	
}
