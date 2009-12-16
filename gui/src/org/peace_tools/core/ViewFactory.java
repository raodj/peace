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

package org.peace_tools.core;

import java.awt.Color;
import java.awt.Component;
import java.awt.Container;
import java.awt.Cursor;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.HashMap;

import javax.swing.Icon;
import javax.swing.JComponent;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.ProgressMonitorInputStream;
import javax.swing.SwingUtilities;

import org.peace_tools.data.ClusterFile;
import org.peace_tools.data.ClusterTreeTableModel;
import org.peace_tools.data.DataStore;
import org.peace_tools.data.ESTList;
import org.peace_tools.data.MST;
import org.peace_tools.data.MSTTreeModel;
import org.peace_tools.generic.ProgrammerLog;
import org.peace_tools.generic.Utilities;
import org.peace_tools.generic.dndTabs.DnDTabListener;
import org.peace_tools.generic.dndTabs.DnDTabPos;
import org.peace_tools.generic.dndTabs.DnDTabbedPane;
import org.peace_tools.views.ClusterSummaryView;
import org.peace_tools.views.ClusterTreeTableView;
import org.peace_tools.views.DataSetFileListView;
import org.peace_tools.views.DataSetTreeView;
import org.peace_tools.views.JobListView;
import org.peace_tools.views.MSTFileView;
import org.peace_tools.views.ProgrammerLogPane;
import org.peace_tools.views.ServerListView;
import org.peace_tools.views.UserLogPane;
import org.peace_tools.workspace.DataSet;
import org.peace_tools.workspace.MSTClusterData;
import org.peace_tools.workspace.MSTData;

/**
 * A factory to help with creation of views.
 * 
 * This class provides a convenience mechanism to create and add a 
 * view to a given main frame. This class is a base class that is
 * not meant to be instantiated. Instead ask MainFrame for a valid
 * ViewFactory by calling mainFrame.getViewFactory() method.
 */
public abstract class ViewFactory implements DnDTabListener {
	/**
	 * Set of enumeration values to identify the different types of
	 * views that this factory can create. 
	 */
	public enum ViewType { EST_FILE, MST_FILE, CLUSTER_FILE, 
		CLUSTER_SUMMARY, DATASET_TREE, DATASET_FILE, JOB_LIST, SERVER_LIST, 
		USER_LOGS, PROGRAMMER_LOGS, TEXT_VEIW, DEFAULT_VIEW}

	/**
	 * The constructor. 
	 * 
	 * The constructor merely saves the reference to the owning 
	 * main frame for future reference.
	 * 
	 * @param mf The main frame that logically owns this view factory.
	 */
	protected ViewFactory(MainFrame mf) {
		this.mainFrame = mf;
		assert(mainFrame != null);
		// Register this view factory as a listener to receive 
		// notification on tab deletions.
		DnDTabbedPane.addListener(this);
	}

	/**
	 * Method to create the standard set of static views.
	 * 
	 * This is an interface method that is used to create the set
	 * of static views associated with this view factory. This
	 * method is typically called only once from the main method.
	 * However, calling this method multiple times on the same view
	 * factory has no major side effects.
	 */
	public synchronized void createStaticViews() {
		if (staticViews == null) {
			staticViews   = new HashMap<ViewType, JComponent>();
			staticViewPos = new HashMap<ViewType, DnDTabPos>();
			createViews();
			createLogViews();
		}
	}

	/**
	 * Method to display a static view that is currently not visible.
	 * 
	 * This method can be used to re-display the view that is currently
	 * not visible. If a view is already visible then this method
	 * does not have any side effect. 
	 * 
	 * @param view The type of view to be displayed by this method. 
	 */
	public void showStaticView(ViewType view) {
		JComponent staticView = staticViews.get(view);
		if (staticView == null) {
			// The view is not available or is already visible.
			return;
		}
		if (staticView.isVisible()) {
			// The view is already visible. Simply emphasize it.
			this.emphasize(staticView);
		}
		// Make the view visible and restore its position.
		DnDTabPos tabPos = staticViewPos.get(view);
		if (tabPos != null) {
			// Found a valid component and position. Restore the view.
			staticView.setVisible(true);
			tabPos.restorePosition(mainFrame.getDesktop(), staticView);
		}
	}
	
	/**
	 * This is a helper method that invoked from the constructor to
	 * create the programmer and user log panes. This method was
	 * introduced to keep the code clutter down in the constructor.
	 */
	private void createLogViews() {
		DnDTabbedPane centerPane = mainFrame.getCenterPane();
		// Create and setup the user log and programmer log.
		UserLogPane ulp = new UserLogPane();
		ProgrammerLogPane plp = new ProgrammerLogPane();
		// Add the programmer and user log in drag panes.
		DnDTabbedPane logArea = 
			centerPane.createSplitPane("User Logs", 
					Utilities.getIcon("images/16x16/UserLog.png"), ulp, 
					DnDTabbedPane.Location.BOTTOM, 0.8, 575);
		// Add the programmer log behind the user pane.
		logArea.createSplitPane("Programmer Logs",
				Utilities.getIcon("images/16x16/ProgLog.png"), plp, 
				DnDTabbedPane.Location.CENTER);
		// Save all the views we created for future references
		staticViews.put(ViewType.USER_LOGS, ulp);
		staticViews.put(ViewType.PROGRAMMER_LOGS, plp);
	}

	/**
	 * This is a helper method that invoked from the constructor to
	 * create the standard views (tabs that display various information
	 * from the work space) . This method was mainly introduced to keep 
	 * the code clutter down in the constructor.
	 */
	private void createViews() {
		DnDTabbedPane centerPane = mainFrame.getCenterPane();
		// Create the hierarchical data set view and add it to the left
		// of the center pane.
		DataSetTreeView dstv = new DataSetTreeView(mainFrame);
		// Wrap the tree view in a scroll pane to let it scroll.
		JScrollPane jsp = new JScrollPane(dstv);
		jsp.setBorder(null);
		DnDTabbedPane leftPane = 
			centerPane.createSplitPane("Hierarchy", 
					Utilities.getIcon("images/16x16/DataSet.png"), jsp,
					DnDTabbedPane.Location.LEFT, 0.2, 200);
		// Add the FileListView
		DataSetFileListView dsflv = new DataSetFileListView(mainFrame);
		leftPane.createSplitPane("File List", 
				Utilities.getIcon("images/16x16/FileView.png"), dsflv, 
				DnDTabbedPane.Location.CENTER);
		// Create and add job list to the bottom of the hierarchy
		JobListView jlv = new JobListView();
		DnDTabbedPane leftBotPane = 
			leftPane.createSplitPane("Jobs", 
					Utilities.getIcon("images/16x16/Job.png"), 
					jlv, DnDTabbedPane.Location.BOTTOM, 0.5, 350);
		// Create server list and add that below the list of jobs
		ServerListView slv = new ServerListView();
		leftBotPane.createSplitPane("Servers",
				Utilities.getIcon("images/16x16/Server.png"),
				slv, DnDTabbedPane.Location.BOTTOM, 0.5, 200);
		// Save all the views we created for future references
		staticViews.put(ViewType.DATASET_TREE, dstv);
		staticViews.put(ViewType.DATASET_FILE, dsflv);
		staticViews.put(ViewType.JOB_LIST, jlv);
		staticViews.put(ViewType.SERVER_LIST, slv);
	}

	/**
	 * Helper method to create a default view for a given data type.
	 * 
	 * @param estFileName The EST file name from where the FASTA data is
	 * to be loaded.
	 * 
	 * @param dataFileName The name of the file where the additional data for
	 * this view is stored.
	 * 
	 * @param viewType The type of view to be displayed.
	 * @return A valid view if successfully created.
	 * @throws Exception This method throws various exceptions when errors
	 * occur during view creation.
	 */
	private JComponent createDefaultView(String estFileName, 
			String dataFileName, ViewType viewType) throws Exception {
		// Check and load the FASTA file
		JComponent view = null;
		ESTList ests = null;
		if (estFileName != null) {
			ests = DataStore.get().getFASTA(estFileName, mainFrame);
		}
		// Check and load the data file depending on the view type
		switch (viewType.ordinal()) {
		case 1:
			MST mst = DataStore.get().getMSTData(dataFileName, mainFrame);
			MSTTreeModel mstModel = new MSTTreeModel(mst, ests);
			view = new MSTFileView(mainFrame, mstModel);
			break;

		case 2:
			ClusterFile ct = DataStore.get().getClusterData(dataFileName, mainFrame);
			ClusterTreeTableModel model = new ClusterTreeTableModel(ct, ests); 
			view = new ClusterTreeTableView(model, mainFrame);
			break;

		case 0:
		default:
		}
		return view;
	}

	/**
	 * Helper method to create a textual view for a given file.
	 * 
	 * @param fileName The text file whose contents is to be loaded
	 * and displayed
	 * 
	 * @return A valid view if successfully created.
	 * @throws Exception This method throws various exceptions when errors
	 * occur during view creation.
	 */
	private JComponent createTextView(String fileName) throws Exception {
		// Check and load the text file
		InputStream fis = new FileInputStream(fileName);
		// Wrap the input stream in progress monitor
		fis = new ProgressMonitorInputStream(mainFrame, 
				"Please wait while the file is loaded...\n" +
				"(file: " + fileName + ")", fis);
		// Load the data from file into a big string
		String data = Utilities.readFullStream(fis);
		// Put the information in a jtext area
		JTextArea display = new JTextArea();
		display.setText(data);
		display.setEditable(false);
		display.setCaretPosition(0);
		// Put text area in a scroll pane and return it
		return new JScrollPane(display);
	}

	/**
	 * Method to create a graph type of view for a given entry.
	 * 
	 * This method performs the actual loading process via  a
	 * background thread so that the GUI does not become unresponsive
	 * as the data is loaded.
	 * 
	 * @param wsEntry The entry object that provides the necessary
	 * view 
	 */
	public void createSummaryView(final MSTClusterData cluster) {
		final String dataFileName = cluster.getPath();
		ViewType viewType   = ViewType.CLUSTER_SUMMARY;
		// Check and handle duplicate file name & view type.
		String viewSignature = dataFileName + "_" + viewType;
		if (views.get(viewSignature) != null) {
			// The view already exists. Nothing further to be done.
			emphasize(views.get(viewSignature).get(0));
			return;
		}
		Thread creator = new Thread(new Runnable() {
			@Override
			public void run() {
				try {
					// Set cursor to indicate longer operations
					mainFrame.getContentPane().setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
					ClusterFile ct = DataStore.get().getClusterData(dataFileName, mainFrame);
					// Load the ESt file as well if needed.
					String estFileName = cluster.getDataSet().getPath();
					ESTList ests = DataStore.get().getFASTA(estFileName, mainFrame);
					// Create the view.
					final ClusterSummaryView csv = new ClusterSummaryView(mainFrame, ct, ests);
					// Add the view to the center panel in the main frame.
					SwingUtilities.invokeAndWait(new Runnable() {
						@Override
						public void run() {
							addView(csv, dataFileName, ViewType.CLUSTER_SUMMARY);
						}
					});
				} catch (Exception exp) {
					// Reset cursor type
					mainFrame.getContentPane().setCursor(Cursor.getDefaultCursor());
					JPanel msg = Utilities.collapsedMessage(CANT_MAKE_VIEW, 
							Utilities.toString(exp));
					JOptionPane.showMessageDialog(mainFrame, msg, 
							"Error creating view", 
							JOptionPane.ERROR_MESSAGE);
				} finally {
					// Reset cursor type
					mainFrame.getContentPane().setCursor(Cursor.getDefaultCursor());
				}

			}			
		});
		// Start the view creation in background.
		creator.start();
	}

	/**
	 * Utility method to add a newly created view to display and internal tables.
	 * 
	 * This is a refactored utility method that is used to add a view to the
	 * internal tables and to the main frame.
	 * 
	 * @param view The newly created view to be added.
	 * @param dataFileName The data file from where the view was loaded.
	 * @param viewType The type of view that has been created.
	 */
	private void addView(JComponent view, String dataFileName, ViewType viewType) {
		String viewSignature = dataFileName + "_" + viewType;
		// Set name for view so that it is easy to locate the view
		// again at a later date.
		view.setName(viewSignature);
		// Add the view to the center panel in the main frame.
		File file = new File(dataFileName);
		mainFrame.getCenterPane().createSplitPane(file.getName(), 
				TabIcons[viewType.ordinal()], view, 
				DnDTabbedPane.Location.CENTER);
		mainFrame.getCenterPane().setSelectedComponent(view);
		// Add the view to the list in a cache.
		ArrayList<JComponent> list = views.get(viewSignature);
		if (list == null) {
			list = new ArrayList<JComponent>();
			views.put(viewSignature, list);
		}
		// Add view to the list.
		list.add(view);		
	}

	/**
	 * Method to ease creation of a specific view. 
	 * 
	 * This method is the preferred approach to creating and adding
	 * a view to the MainFrame. A view is a specific display of a 
	 * given data item. This method creates and adds a view as
	 * needed.
	 * 
	 * @param wsEntry The work space entry for which a view must be
	 * created (if one does not exist).
	 * @param duplicate If this flag is true, then a duplicate view is created
	 * even if a view already exists.
	 * @param text If this flag is true, then a textual view of the data is
	 * created.
	 */
	public void createView(Object wsEntry, boolean duplicate, final boolean textView) {
		String estFileName   = null;
		String dataFileName  = null;
		ViewType viewType    = ViewType.EST_FILE;
		// Setup the above variables based on entry type
		if (wsEntry instanceof DataSet) {
			DataSet ds   = (DataSet) wsEntry;
			dataFileName = ds.getPath();
			estFileName  = dataFileName;
		} else if (wsEntry instanceof MSTData) {
			MSTData mst  = (MSTData) wsEntry;
			estFileName  = mst.getDataSet().getPath();
			dataFileName = mst.getPath();
			viewType     = ViewType.MST_FILE;
		} else if (wsEntry instanceof MSTClusterData) {
			MSTClusterData cluster = (MSTClusterData) wsEntry;
			estFileName            = cluster.getDataSet().getPath();
			dataFileName           = cluster.getPath();
			viewType               = ViewType.CLUSTER_FILE;
		}
		// Check and handle duplicate file name & view type.
		String viewSignature = dataFileName + "_" + (textView ? ViewType.TEXT_VEIW : viewType);
		if ((!duplicate) &&  (views.get(viewSignature) != null)) {
			// The view already exists. Nothing further to be done.
			emphasize(views.get(viewSignature).get(0));
			return;
		}
		// Create the views via  a background thread to ensure that the
		// gui remains interactive. First create a bunch of final variables
		// so that the values are accessible inside the thread below.
		final String finalDataFileName = dataFileName;
		final String finalEstFileName  = estFileName;
		final ViewType finalViewType   = viewType;
		Thread viewCreator = new Thread(new Runnable() {
			@Override
			public void run() {
				try {
					// Set cursor to indicate longer operations
					mainFrame.getContentPane().setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
					// Check and create suitable view.
					JComponent view = null;
					if (textView) {
						view = createTextView(finalDataFileName);
					} else {
						view = createDefaultView(finalEstFileName, finalDataFileName, finalViewType);
					}
					// Now we should have a valid view.
					if (view == null) {
						throw new IllegalArgumentException("View was null!");
					}
					// Add the newly created view from the swing thread.
					final JComponent finalView = view;
					SwingUtilities.invokeAndWait(new Runnable() {
						@Override
						public void run() {
							addView(finalView, finalDataFileName, finalViewType);	
						}
					});
				} catch (Exception exp) {
					// Reset cursor type
					mainFrame.getContentPane().setCursor(Cursor.getDefaultCursor());
					final JPanel msg = Utilities.collapsedMessage(CANT_MAKE_VIEW, 
							Utilities.toString(exp));
					JOptionPane.showMessageDialog(mainFrame, msg, 
							"Error creating view", 
							JOptionPane.ERROR_MESSAGE);
				} finally {
					// Reset cursor type
					mainFrame.getContentPane().setCursor(Cursor.getDefaultCursor());
				}
			}
		});
		// Start the view creation thread.
		viewCreator.start();
	}

	/**
	 * A helper method to create thread to make a tab blink to emphasize
	 * its presence. This method is called when a view needs to be 
	 * emphasized to draw the attention of the user.  This method 
	 * creates a new daemon thread that changes the background of the
	 * component frequently drawing the attention of the user.
	 * 
	 * @param component The component to which the user's attention is
	 * to be drawn.
	 */
	public void emphasize(final JComponent component) {
		// Check to see if this component is inside a tab pane.
		Container p = component.getParent();
		while ((p != null) && !(p instanceof DnDTabbedPane)) {
			p = p.getParent();
		}
		// Save the tab information for further use below.
		final DnDTabbedPane tabParent = (DnDTabbedPane) p;
		final int index = (tabParent == null) ? - 1 :
			tabParent.indexOfComponent(component);
		// First ensure the emphasized page is visible.
		if ((tabParent != null) && (index != -1)) {
			tabParent.setSelectedIndex(index);
		}
		// Create the thread to change colors rapidly
		Thread blinker = new Thread(new Runnable() {
			@Override
			public void run() {
				Color normal  = component.getBackground();
				Color lighter = Color.white;
				Color darker  = component.getBackground().darker();
				try {
					for(int i = 0; (i < 5); i++) {
						Color bgColor = (i % 2 == 0) ? lighter : darker;
						if ((tabParent != null) && (index != -1)) {
							tabParent.setEnabledAt(index, (i % 2) == 0);
							tabParent.repaint();
						} else {
							component.setBackground(bgColor);
							component.repaint();
						}
						Thread.sleep(150);
					}
				} catch (InterruptedException ie) {
					ProgrammerLog.log(ie);
				} finally {
					if ((tabParent != null) && (index != -1)) {
						tabParent.setEnabledAt(index, true);
					} else {
						component.setBackground(normal);
					}
				}
			}
		});
		blinker.setDaemon(true);
		blinker.start();
	}

	/**
	 * Method to intercept events generated by DnDTabbedPane to 
	 * report removal of views.
	 *
	 * This method removes the view from the list of views in the
	 * views hash map. This method is called from the DnDTabbedPane
	 * 
	 * 
	 * @param ce The component event to be processed.
	 */
	@Override
	public void tabDeleted(Component component) {
		// Locate its entry in the hash map.
		ArrayList<JComponent> viewList = views.get(component.getName());
		if (viewList != null) {
			viewList.remove(component);
			if (viewList.size() == 0) {
				views.remove(component.getName());
			}
		} else {
			// If this is one of the static views that don't get really
			// removed but just hidden, then save its position.
			for(ViewType view: ViewType.values()) { 
				if (staticViews.get(view) == component) {
					// Found a valid component. Save its position.
					DnDTabPos tabPos = new DnDTabPos();
					tabPos.savePosition(mainFrame.getDesktop(), component);
					// Save it in reference map for future reference
					staticViewPos.put(view, tabPos);
					// Hide the component
					component.setVisible(false);
				}
			}
		}
	}

	/**
	 * A convenient reference to the main frame that logically owns
	 * this view factory. This reference is used to add views to
	 * appropriate locations in the main frame.
	 */
	private final MainFrame mainFrame;

	/**
	 * A static message that is to be displayed if a view could
	 * not be opened.
	 */
	private static final String CANT_MAKE_VIEW = 
		"<html>Unable to open a view for the specified item.<br>" +
		"Refer to the details below for additional information.</html>";

	/**
	 * A hash map that maintains the list of views that are currently
	 * visible in this frame. Entries are added and removed by the 
	 * ViewFactory class.
	 */
	private HashMap<String, ArrayList<JComponent> > views = 
		new HashMap<String, ArrayList<JComponent> >();

	/**
	 * The list of static views associated with any main frame.
	 * The static views are created only once for a given main frame.
	 * Each view is represented by a unique enumeration value.
	 */
	private HashMap<ViewType, JComponent > staticViews = null;

	/**
	 * The last saved positions for the various static views. This
	 * hash map saves the last saved positions for the various static
	 * views. This enables restoring static views (assuming they are
	 * visible) at appropriate positions.
	 */
	private HashMap<ViewType, DnDTabPos> staticViewPos = null;
	
	/**
	 * The array of icons that are displayed in a the tabs along
	 * with titles for each view. The icons are arranged in the
	 * same order in which the views are enumerated so that the
	 * ordinal values of enumerations can be used to obtain
	 * icons for tab titles.
	 */
	private static final Icon TabIcons[] = { 
		Utilities.getIcon("images/16x16/EST.png"),
		Utilities.getIcon("images/16x16/MST.png"),
		Utilities.getIcon("images/16x16/Cluster.png"),
		Utilities.getIcon("images/16x16/ClusterSummary.png")
	};
}
