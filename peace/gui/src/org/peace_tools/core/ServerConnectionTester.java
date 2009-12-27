package org.peace_tools.core;

import java.awt.BorderLayout;
import java.awt.Cursor;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;
import javax.swing.border.EmptyBorder;

import org.peace_tools.core.session.ServerSession;
import org.peace_tools.core.session.SessionFactory;
import org.peace_tools.generic.Utilities;
import org.peace_tools.workspace.Server;

public class ServerConnectionTester implements ActionListener, Runnable {

	/**
	 * The constructor.
	 * 
	 * The constructor coordinates the task of creating the
	 * informative dialog. The actual task of creating the
	 * tester thread is performed by the start() method. 
	 * 
	 * @param mainFrame The main frame that logically owns this
	 * class.
	 * 
	 * @param server The server with which a connection is to be 
	 * established.
	 */
	public ServerConnectionTester(MainFrame mainFrame, Server server) {
		this.server = server;
		dialog      = createDialog(mainFrame);
		dialog.pack();
		dialog.setLocationRelativeTo(mainFrame);
	}
	
	/**
	 * Starts the background thread to test connection.
	 * 
	 * This method must be used to start up the background
	 * connection thread that attempts to establish a connection
	 * with the remote server.
	 */
	public void start() {
		// Create the background connect thread.
		connThread  = new Thread(this);
		connThread.setDaemon(true);
		connThread.start();
	}
	
	/**
	 * The background thread method that performs the connection test.
	 * 
	 * This method is invoked from a separate background thread. This
	 * method is responsible for performing the core operations
	 * pertaining to testing the connection to the server.
	 */
	@Override
	public void run() {
		// First display the dialog.
		SwingUtilities.invokeLater(new Runnable() {
			@Override
			public void run() {
				progressBar.setIndeterminate(true);
				dialog.setVisible(true);
			}
		});
		
		try {
			// Create a suitable server session first.
			ServerSession session = SessionFactory.createSession(dialog, server);
			// Setup the purpose for this session.
			session.setPurpose(String.format(PURPOSE_MSG, server.getName()));
			// Establish the connection.
			session.connect();
			// OK! connection succeeded. Run a simple command
			String outputs[] = {"", ""};
			int exitCode = session.exec("set", outputs);
			if (exitCode != 0) {
				throw new Exception("Connection was established to " +
						server.getName() + "\nbut could not even run a simple " +
						"set command to determine enviornment.\nThis is a problem " +
						"with the server!");
			}
			// So far so good. Wind down the connection
			session.disconnect();
			// Update the status of the server
			server.setStatus(Server.ServerStatusType.GOOD);
			// Let the user know about the success along with environment
			String info = "Successfully connected to " + server.getName();
			JPanel fullInfo = Utilities.collapsedMessage(info, outputs[0]);
			JOptionPane.showMessageDialog(dialog, fullInfo, 
					"Connection successful", JOptionPane.INFORMATION_MESSAGE);
		} catch (Exception exp) {
			progressBar.setIndeterminate(false);
			String errMsg  = String.format(ERROR_MSG, server.getName());
			String expInfo = Utilities.toString(exp);
			JPanel info    = Utilities.collapsedMessage(errMsg, expInfo);
			JOptionPane.showMessageDialog(dialog, info, 
					"Can't connect to server", JOptionPane.ERROR_MESSAGE);
			// Update the status of the server to be invalid.
			server.setStatus(Server.ServerStatusType.CONNECT_FAILED);
		} finally {
			// Hide the dialog.
			dialog.setVisible(false);
			dialog.dispose();
		}
	}
	
	/**
	 * Helper method to create the dialog to be displayed.
	 * 
	 * This method is a helper method that was introduced to
	 * streamline the constructor. It creates a JDialog containing
	 * a message, a progress bar, and an "Interrupt" button. This
	 * dialog is displayed just before the connection with the
	 * server is attempted. 
	 * 
	 * @param parent The frame that logically owns this dialog.
	 * 
	 * @return The dialog to be displayed to the user.
	 */
	private JDialog createDialog(JFrame parent) {
		// Create the message and progress bar to be displayed.
		String info = String.format(INFO_MSG, server.getName());
		JLabel msg = new JLabel(info, UIManager.getIcon("OptionPane.informationIcon"), JLabel.LEFT);
		progressBar = new JProgressBar(JProgressBar.HORIZONTAL);
		JPanel infoPanel = new JPanel(new BorderLayout(5, 10));
		infoPanel.add(msg, BorderLayout.NORTH);
		infoPanel.add(progressBar, BorderLayout.SOUTH);
		
		// Create the interrupt button and its panel.
		JButton interrupt = Utilities.createButton(null, "Interrupt", 
				"interrupt", this, 
				"Interrupt the process of trying to connect to the server", true);
		// Wrap button in a panel to ensure it does not become too long.
		JPanel buttonPanel = new JPanel(new FlowLayout(FlowLayout.CENTER));
		buttonPanel.add(interrupt);
		
		// Add message, progress bar, and button to a content pane
		JPanel contentPane = new JPanel(new BorderLayout(0, 0));
		contentPane.setBorder(new EmptyBorder(20, 20, 20, 20));
		contentPane.add(infoPanel, BorderLayout.NORTH);
		contentPane.add(buttonPanel, BorderLayout.SOUTH);
		
		// Finally Create the connection tester dialog.
		JDialog dialog = new JDialog(parent, "Server Connection Test", true);
		dialog.setIconImage(Utilities.getIcon("images/16x16/ServerConnect.png").getImage());
		dialog.setContentPane(contentPane);
		dialog.setDefaultCloseOperation(JDialog.DO_NOTHING_ON_CLOSE);
		dialog.setResizable(false);
		// Setup the cursor to be the wait cursor.
		dialog.setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
		// Return dialog for further use.
		return dialog;
	}
	
	/**
	 * Handle user clicking the "Interrupt" button.
	 * 
	 * This action listener handles the case when the user clicks
	 * on the interrupt button. This method essentially interrupts
	 * the background thread that is attempting to connect to the
	 * remote server.
	 * 
	 * @param event The action event associated with the button click.
	 * This event is not really used.
	 */
	@Override
	public void actionPerformed(ActionEvent event) {
		if (connThread != null) {
			connThread.interrupt();
		}
	}
	
	/**
	 * The server to which this class is attempting to establish is
	 * a connection to verify that the GUI can actually connect to
	 * the server.
	 */
	private final Server server;
	
	/**
	 * The roving (indeterminate mode) progress bar to illustrate 
	 * that some work is being done (and to ensure that the GUI
	 * is not hanging).
	 */
	private JProgressBar progressBar;
	
	/**
	 * The background thread that is attempting to connect to the
	 * remote server. This thread is created in the constructor.
	 */
	private Thread connThread;
	
	/**
	 * The dialog that is being used to display some message and
	 * progress information to the user. This dialog is created
	 * in the constructor but is displayed and disposed in the
	 * background thread.
	 */
	private JDialog dialog;
	
	/**
	 * A simple message that is formatted (to fill in server name)
	 * and displayed to the user informing him of the progress
	 * being made.
	 */
	private static final String INFO_MSG = "<html>" +
			"Attempting to establish an independent connection<br>" +
			"to server: %s<br><br>" +
			"<i>Please wait...</i></html>";
	
	/**
	 * An error message that is formatted (to fill in server name)
	 * and displayed to the user informing him about the problem
	 * that occurred when attempting to connect to the server.
	 */
	private static final String ERROR_MSG = "<html>" +
			"Unable to establish an independent connection<br>" +
			"to server: %s<br>" +
			"because an error Occured. Refer to the detailed information<br>" +
			"below for additional information about the error that occured." +
			"</html>";
	
	/**
	 * A simple message that is formatted (to fill in server name)
	 * and set as the purpose for the session/connection.
	 */
	private static final String PURPOSE_MSG = "<html>" +
			"Attempting to connect to the<br>" +
			"server: %s<br>" +
			"to verify connectability status." +
			"</html>";
}
