package org.peace_tools.views;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.ArrayList;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.Icon;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JToggleButton;
import javax.swing.JToolBar;
import javax.swing.JTree;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeCellRenderer;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.MutableTreeNode;

import org.peace_tools.data.ClusterFile;
import org.peace_tools.data.DataStore;
import org.peace_tools.data.ESTList;
import org.peace_tools.generic.CustomBorder;
import org.peace_tools.generic.IconTreeNode;
import org.peace_tools.generic.Utilities;
import org.peace_tools.workspace.DataSet;
import org.peace_tools.workspace.FWAnalyzer;
import org.peace_tools.workspace.Filter;
import org.peace_tools.workspace.Heuristic;
import org.peace_tools.workspace.Job;
import org.peace_tools.workspace.JobSummary;
import org.peace_tools.workspace.MSTClusterData;
import org.peace_tools.workspace.MSTData;
import org.peace_tools.workspace.Param;
import org.peace_tools.workspace.Server;
import org.peace_tools.workspace.Workspace;

/**
 * A helper (utility) class to assemble various properties about 
 * Jobs, MST file, and Cluster file in the form of a JTree.
 * 
 * <p>This class is a helper class that is used by different views in
 * this package to present the user with a consistent set of information
 * about the file (may it be a MST file or Cluster file) being 
 * viewed.</p>
 * 
 * <p><b>Note</b>: All the API methods exposed by this class are static
 * methods and they can be directly invoked without instantiating an
 * object.</p>
 */
public class PropertiesTreeMaker {
	/**
	 * Build a property tree with MST file as the main entry.
	 * 
	 * This method can be used to build a properties tree with the MST
	 * entry as the main entry. This method calls the various internal
	 * helper methods in this class to build the complete three.
	 * 
	 * @param mst The MST workspace entry for which a detailed properties tree
	 * is to be built by this method.
	 * 
	 * @param parent The parent component to be used to display dialog boxes
	 * on errors.
	 * 
	 * @return The detailed properties tree for the given MST entry
	 */
	public static JTree makeProperties(MSTData mst, Component parent) {
		// First create the root element for the tree.
		DefaultMutableTreeNode root = makeMSTProperties(mst, true, true, parent);
		// Add job information if available.
		root.add(makeJobProperties(mst.getJobSummary()));
		// Create the tree model to be populated with necessary information.
		DefaultTreeModel model = new DefaultTreeModel(root);
		JTree propertiesTree = new JTree(model);
		// Setup properties for the JTree
		propertiesTree.setDragEnabled(false);
		propertiesTree.setEditable(false);
		propertiesTree.setRootVisible(true);
		// Set renderer to display icons correctly
		addPropertiesRenderer(propertiesTree);
		// Return the Jtree to the user.
		return propertiesTree;
	}
	
	/**
	 * Method to compute and setup summary information.
	 * 
	 * This method is invoked from the constructor to create and populate the 
	 * summary information tab. The summary information is computed once when
	 * the view is created. After that the summary information is built using
	 * the PropertiesTreeMaker class.
	 * 
	 * @return The scroll pane that contains the summary information.
	 */
	public static JSplitPane createPropertiesLayout(String title, JTree summaryInfo, 
			Component mainView, JToolBar toolbar, int toolPosition) {
		// Wrap the properties tree into a scroll pane to enable scrolling.
		final JScrollPane summaryPane = new JScrollPane(summaryInfo);
		summaryPane.setPreferredSize(new Dimension(150, 200));
		summaryPane.setMinimumSize(new Dimension(50, 50));
		// Set a header with information to make it look nice.
		JLabel heading = new JLabel("<html><b>" + title + "</b></html>",
				Utilities.getIcon("images/16x16/Information.png"), JLabel.LEFT);
		heading.setBorder(BorderFactory.createEmptyBorder(2, 5, 2, 5));
		summaryPane.setColumnHeaderView(heading);
		summaryPane.setViewportBorder(new CustomBorder("dsss"));
		// Now wrap the main view and the summaryPane into a suitable split pane
		final JSplitPane contentPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, summaryPane, mainView);
		contentPane.setDividerSize(3);
		contentPane.setOneTouchExpandable(false);
		contentPane.setBorder(null);
		contentPane.setDividerLocation(150);
		// Now add a suitable button to the tool bar, if tool bar is valid.
		if (toolbar != null) {
			// Add button to show or hide summary information
			toolbar.add(Box.createHorizontalStrut(10), toolPosition);
			JToggleButton summaryButton = 
				new JToggleButton("Summary ...", 
						Utilities.getIcon("images/16x16/Information.png"), true);
			summaryButton.setToolTipText("Show/Hide summary information about this cluster file");
			toolbar.add(summaryButton, toolPosition + 1);
			// Finally add a action listener to handle clicking this button.
			summaryButton.addActionListener(new ActionListener() {
				@Override
				public void actionPerformed(ActionEvent event) {
					// Show or hide the summary pane depending on current settings.
					// split/un-split the content pane depending on visibility
					if (!summaryPane.isVisible()) {
						summaryPane.setVisible(true);
						contentPane.setLeftComponent(summaryPane);
						int dividerPos = summaryPane.getPreferredSize().width;
						if (dividerPos >= contentPane.getWidth() - 50) {
							dividerPos = contentPane.getWidth() - 50;
						}
						contentPane.setDividerLocation(dividerPos);
					} else {
						// Save the current size as preferred size for future 
						summaryPane.setPreferredSize(summaryPane.getSize());
						contentPane.setLeftComponent(null);
						summaryPane.setVisible(false);
					}
					contentPane.revalidate();
					contentPane.repaint();
				}
			});
			// Now finally return the content pane for the user to work with
			return contentPane;
		}
		// Return the final content pane for the user to work with.
		return contentPane;
	}
	
	/**
	 * Build a property tree with cluster file as the main entry.
	 * 
	 * This method can be used to build a properties tree with the cluster
	 * entry as the main entry. This method calls the various internal
	 * helper methods in this class to build the complete three.
	 * 
	 * @param cluster The cluster workspace entry for which a detailed properties tree
	 * is to be built by this method.
	 * 
	 * @param parent The final parent component that owns the JTree being created.
	 * This parent is used to report error messages that may occur when building
	 * properties.
	 * 
	 * @return The detailed properties tree for the given cluster entry.
	 */
	public static JTree makeProperties(MSTClusterData cluster, Component parent) {
		// First create the root element for the tree.
		DefaultMutableTreeNode root = new IconTreeNode(cluster, CustomIcons[CLUSTER_ICON]);
		// Add information about the clustering file.
		root.add(new DefaultMutableTreeNode("InternalID: " + cluster.getID()));
		root.add(new DefaultMutableTreeNode("Description: " + cluster.getDescription()));
		// Add information about the physical cluster file on disk on local machine
		root.add(makeFileProperties(cluster.getPath(), "Cluster File Properties"));
		// Add clustering summary information.
		root.add(makeSummaryProperties(cluster, parent));
		// Add information about the source FASTA file.
		root.add(makeDataSetProperties(cluster.getDataSet(), parent));
		// Add MST information if available.
		MSTData mst = Workspace.get().getMSTData(cluster.getJobSummary().getJobID());
		if (mst != null) {
			root.add(makeMSTProperties(mst, false, true, parent));
		} else {
			// MST information could not be located. Add dummy node.
			DefaultMutableTreeNode mstNode = new IconTreeNode("MST (data unavailable", 
															  CustomIcons[ERROR_ICON]);
			root.add(new DefaultMutableTreeNode("Internal ID: " + cluster.getMSTID()));
			root.add(new DefaultMutableTreeNode("Note: MST entry may have been deleted by user"));
			root.add(mstNode);
		}
		// Add job information if available.
		root.add(makeJobProperties(cluster.getJobSummary()));
		// Create the tree model to be populated with necessary information.
		DefaultTreeModel model = new DefaultTreeModel(root);
		JTree propertiesTree = new JTree(model);
		// Setup properties for the JTree
		propertiesTree.setDragEnabled(false);
		propertiesTree.setEditable(false);
		propertiesTree.setRootVisible(true);
		// Set renderer to display icons correctly
		addPropertiesRenderer(propertiesTree);
		// Return the Jtree to the user.
		return propertiesTree;
	}

	/**
	 * Helper method to create a summary node for a given clustering.
	 * 
	 * This is an internal helper method that is used to build the properties
	 * node for a given cluster entry/
	 * 
	 * @param cluster The cluster workspace entry for which the properties sub-tree 
	 * is to be built
	 * 
	 * @return A sub-tree with necessary statistical summary information to be included 
	 * in the main properties tree at an appropriate location.
	 */
	private static MutableTreeNode makeSummaryProperties(MSTClusterData cluster, Component parent) {
		// If the clustering file exists then let's also make summary information.
		DefaultMutableTreeNode summary = null;
		try {
			ClusterFile clsFile = DataStore.get().getClusterData(cluster.getPath(), parent);
			if (clsFile != null) {
				double[] clsStats = clsFile.getRoot().computeStatistics();
				summary = new IconTreeNode("Cluster Statistics", CustomIcons[STATS_ICON]);
				summary.add(new DefaultMutableTreeNode("Cluster count: " + clsStats[0]));
				summary.add(new DefaultMutableTreeNode("Smallest cluster size: " + clsStats[1]));
				summary.add(new DefaultMutableTreeNode("Largest cluster size: " + clsStats[2]));
				summary.add(new DefaultMutableTreeNode("Average cluster size: " + clsStats[3]));
				summary.add(new DefaultMutableTreeNode("Cluster size SD: " + clsStats[4]));
			}
		} catch (Exception e) {
			JPanel msg = Utilities.collapsedMessage("Unable to compute clustering statistics", 
					Utilities.toString(e));
			JOptionPane.showMessageDialog(parent, msg, "Error summarizing clustering file", 
					JOptionPane.ERROR_MESSAGE);
		}
		// Create dummy summary node if needed.
		if (summary == null) {
			summary = new DefaultMutableTreeNode("Summary Statistics not available");
		}
		return summary;
	}

	/**
	 * Helper method to create a properties node for a given job summary node.
	 * 
	 * This is an internal helper method that is used to build the properties
	 * node for a given job summary node. This method prefers to use the actual
	 * Job information if available. If the complete Job information could not
	 * be found, then the summary information in the job summary is used instead.
	 * 
	 * @param summary The Job summary entry for which the properties sub-tree 
	 * is to be built
	 * 
	 * @return A sub-tree with necessary Job properties information to be included 
	 * in the main properties tree at an appropriate location.
	 */
	private static DefaultMutableTreeNode makeJobProperties(JobSummary summary) {
		// Check if a fully job entry is available.
		final String jobID = summary.getJobID();
		final Job job = Workspace.get().getJobList().getjob(jobID);
		if (job != null) {
			// Found full job information. Use it
			return makeJobProperties(job);
		}
		// Could not find full job information. Populate some information
		// using the summary data available to us.
		// Create the top-level sub-tree node to hold file information
		DefaultMutableTreeNode jobNode = new IconTreeNode("Job (only summary available)", 
														  CustomIcons[JOB_ICON]);
		jobNode.add(new DefaultMutableTreeNode("Note: Full job information is not available."));
		jobNode.add(new DefaultMutableTreeNode("Internal ID: " + summary.getJobID()));
		// Add server information (if available) through helper method.
		jobNode.add(makeServerProperties(summary.getServerID()));
		// Add other information we have
		jobNode.add(new DefaultMutableTreeNode("Total CPUs used: " + summary.getCPUs()));
		jobNode.add(new DefaultMutableTreeNode("Heuristics (summary): " + summary.getHeuristicsSummary()));
		jobNode.add(new DefaultMutableTreeNode("Filters (summary): " + summary.getFiltersSummary()));
		// Return job node for further use
		return jobNode;
	}
	
	/**
	 * Helper method to create a properties node for a given MST entry.
	 * 
	 * This is an internal helper method that is used to build the properties
	 * node for a given MST entry.
	 * 
	 * @param mst The MST workspace entry for which the properties sub-tree 
	 * is to be built
	 * 
	 * @param needFileInfo If this flag is true then this method adds a sub-tree
	 * with properties about the physical MST file.
	 *  
	 * @param parent The parent component to be used to display dialog boxes
	 * on errors.
	 * 
	 * @return A sub-tree with necessary MST properties information to be included 
	 * in the main properties tree at an appropriate location.
	 */
	private static DefaultMutableTreeNode makeMSTProperties(MSTData mst,
			boolean needDataSetInfo, boolean needFileInfo, Component parent) {
		// Create the top-level sub-tree node to hold file information
		DefaultMutableTreeNode mstNode = new IconTreeNode(mst, CustomIcons[MST_ICON]);
		// Add data set properties if needed.
		if (needDataSetInfo) {
			// Get the helper method to populate information about data set
			mstNode.add(makeDataSetProperties(mst.getDataSet(), parent));
		}
		// Add ID and user supplied description next.
		mstNode.add(new DefaultMutableTreeNode("Internal ID: " + mst.getID()));
		mstNode.add(new DefaultMutableTreeNode("Description: " + mst.getDescription()));
		// Set the type of analyzer that we are working with
		mstNode.add(new DefaultMutableTreeNode("Type: " + mst.getType()));
		// Add sub-tree about the type of analyzer being using.
		mstNode.add(makeAnalyzerProperties(mst.getFWAnalyzer()));
		// Add file information if requested
		if (needFileInfo) {
			// Get the helper method to populate physical information about
			// the MST file under the root.
			mstNode.add(makeFileProperties(mst.getPath(), "MST File Properties"));
		}
		// Return the MST information back to the caller
		return mstNode;
	}

	/**
	 * Helper method to create a properties node for a given MST entry.
	 * 
	 * This is an internal helper method that is used to build the properties
	 * node for a given MST entry.
	 * 
	 * @param mst The MST workspace entry for which the properties sub-tree 
	 * is to be built
	 * 
	 * @param needFileInfo If this flag is true then this method adds a sub-tree
	 * with properties about the physical MST file.
	 *  
	 * @return A sub-tree with necessary MST properties information to be included 
	 * in the main properties tree at an appropriate location.
	 */
	private static DefaultMutableTreeNode makeJobProperties(Job job) { 
		// Create the top-level sub-tree node to hold file information
		DefaultMutableTreeNode jobNode = new IconTreeNode(job, CustomIcons[JOB_ICON]);
		// Add the core information about the job first.
		jobNode.add(new DefaultMutableTreeNode("Internal ID: " + job.getJobID()));
		jobNode.add(new DefaultMutableTreeNode("Description: " + job.getDescription()));
		// Add server information next.
		jobNode.add(makeServerProperties(job.getServerID()));
		// Add job configuration parameters. 
		jobNode.add(new DefaultMutableTreeNode("Path (on server): " + job.getPath()));
		jobNode.add(new DefaultMutableTreeNode("#Nodes: " + job.getNodes()));
		jobNode.add(new DefaultMutableTreeNode("CPUs per Node: " + job.getCPUsPerNode()));
		jobNode.add(new DefaultMutableTreeNode("Max requested memory (MB): " + job.getMemory()));
		jobNode.add(new DefaultMutableTreeNode("Max requested runtime (hours): " + job.getMaxRunTime()));
		// Now add information about heuristics based on the analyzer being used.
		MSTData mstEntry = Workspace.get().getMSTData(job.getJobID());
		boolean autoHeuristics = ((mstEntry != null) && 
				(mstEntry.getFWAnalyzer().getType().equals(FWAnalyzer.FWAnalyzerType.TWOPASSD2)));
		// Create node for heuristics.
		DefaultMutableTreeNode heuristicsNode = new IconTreeNode("Heuristics", CustomIcons[HEURISTICS_ICON]);
		// Add heuristics only if autoHeuristics is not enabled.
		if (!autoHeuristics) {
			for(Heuristic heuristic: job.getHeuristicList()) {
				heuristicsNode.add(makeHeuristicProperties(heuristic));
			}
		} else {
			heuristicsNode.add(new DefaultMutableTreeNode("Automatically configured u/v and t/v heuristics"));
		}
		jobNode.add(heuristicsNode);
		// Add information about filters
		DefaultMutableTreeNode filtersNode = new IconTreeNode("Filters", CustomIcons[FILTERS_ICON]);
		for(Filter filter: job.getFilterList()) {
			filtersNode.add(makeFilterProperties(filter));
		}
		jobNode.add(filtersNode);
		
		// Return the MST information back to the caller
		return jobNode;
	}

	/**
	 * Helper method to create a properties node for a given server entry.
	 * 
	 * This is an internal helper method that is used to build the properties
	 * node for a given server entry.
	 * 
	 * @param serverID The internal string identifier that uniquely identifies a
	 * server entry in the workspace.
	 * 
	 * @return A sub-tree with necessary server properties information to be included 
	 * in the main properties tree at an appropriate location.
	 */
	private static MutableTreeNode makeServerProperties(String serverID) {
		Server server = Workspace.get().getServerList().getServer(serverID);
		DefaultMutableTreeNode srvrNode = null;
		if (server != null) {
			// Found valid server node. Populate necessary server information.
			srvrNode = new IconTreeNode("Server: " + server.getName(), CustomIcons[SERVER_ICON]);
			srvrNode.add(new DefaultMutableTreeNode("Internal ID: " + server.getID()));
			srvrNode.add(new DefaultMutableTreeNode("Description: " + server.getDescription()));
			srvrNode.add(new DefaultMutableTreeNode("Install Path:" + server.getInstallPath()));
			srvrNode.add(new DefaultMutableTreeNode("Login ID:" + server.getUserID()));
		} else {
			// Did not find server node. Create dummy entry.
			srvrNode = new IconTreeNode("Server entry not found", CustomIcons[ERROR_ICON]);
			srvrNode.add(new DefaultMutableTreeNode("Server ID: " + serverID));
		}
		// Let the caller have the sub-tree node.
		return srvrNode;
	}
	
	/**
	 * Helper method to create a properties node for a given data set entry.
	 * 
	 * This is an internal helper method that is used to build the properties
	 * node for a given data set entry.
	 * 
	 * @param dataSet The data set workspace entry for which the properties sub-tree 
	 * is to be built
	 * 
	 * @return A sub-tree with necessary data set properties information to be included 
	 * in the main properties tree at an appropriate location.
	 */
	private static MutableTreeNode makeHeuristicProperties(Heuristic heuristic) {
		// Create the top-level sub-tree node to hold file information
		DefaultMutableTreeNode hurNode = new DefaultMutableTreeNode("Heuristic: " + heuristic.getName());
		// Add parameter information to the heuristic as well.
		final ArrayList<Param> params = heuristic.getParameters();
		for(Param parameter: params) {
			hurNode.add(new DefaultMutableTreeNode(parameter.getName() + ": " + parameter.getValue()));
		}
		// Return the heuristic information back to the caller
		return hurNode;
	}

	/**
	 * Helper method to create a properties node for a given data set entry.
	 * 
	 * This is an internal helper method that is used to build the properties
	 * node for a given data set entry.
	 * 
	 * @param filter The filter workspace entry for which the properties sub-tree 
	 * is to be built
	 * 
	 * @return A sub-tree with necessary data set properties information to be included 
	 * in the main properties tree at an appropriate location.
	 */
	private static MutableTreeNode makeFilterProperties(Filter filter) {
		// Create the top-level sub-tree node to hold file information
		DefaultMutableTreeNode hurNode = new DefaultMutableTreeNode("Filter: " + filter.getName());
		// Add parameter information to the filter as well.
		final ArrayList<Param> params = filter.getParameters();
		for(Param parameter: params) {
			hurNode.add(new DefaultMutableTreeNode(parameter.getName() + ": " + parameter.getValue()));
		}
		// Return the heuristic information back to the caller
		return hurNode;
	}

	/**
	 * Helper method to create a properties node for a given data set entry.
	 * 
	 * This is an internal helper method that is used to build the properties
	 * node for a given data set entry.
	 * 
	 * @param dataSet The data set workspace entry for which the properties sub-tree 
	 * is to be built
	 * 
	 * @return A sub-tree with necessary data set properties information to be included 
	 * in the main properties tree at an appropriate location.
	 */
	private static MutableTreeNode makeDataSetProperties(DataSet dataSet, Component parent) {
		// Create the top-level sub-tree node to hold file information
		DefaultMutableTreeNode dsetNode = new IconTreeNode(dataSet, CustomIcons[DATASET_ICON]);
		// Add ID and user supplied description next.
		dsetNode.add(new DefaultMutableTreeNode("Description: " + dataSet.getDescription()));
		// Get the helper method to populate physical information about
		// the FASTA file for this data set.
		dsetNode.add(makeFileProperties(dataSet.getPath(), "FASTA File Properties"));
		// If the FASTA file exists then let's also make summary information.
		DefaultMutableTreeNode summary = null;
		try {
			ESTList estList = DataStore.get().getFASTA(dataSet.getPath(), parent);
			if (estList != null) {
				double[] fastaStats = estList.computeStatistics();
				summary = new DefaultMutableTreeNode("Summary Statistics");
				summary.add(new DefaultMutableTreeNode("Fragment count: " + fastaStats[0]));
				summary.add(new DefaultMutableTreeNode("Shortest fragment: " + fastaStats[1]));
				summary.add(new DefaultMutableTreeNode("Longest fragment: " + fastaStats[2]));
				summary.add(new DefaultMutableTreeNode("Average fragment length: " + fastaStats[3]));
				summary.add(new DefaultMutableTreeNode("SD in length: " + fastaStats[4]));
			}
		} catch (Exception e) {
			JPanel msg = Utilities.collapsedMessage("Unable to compute summary statistics for FASTA file.", 
					Utilities.toString(e));
			JOptionPane.showMessageDialog(parent, msg, "Error summarizing FASTA file", 
					JOptionPane.ERROR_MESSAGE);
		}
		// Create dummy summary node if needed.
		if (summary == null) {
			summary = new DefaultMutableTreeNode("Summary Statistics not available");
		}
		dsetNode.add(summary);
		// Return the data set information back to the caller
		return dsetNode;
	}

	/**
	 * Helper method to create a properties node for a given analyzer entry.
	 * 
	 * This is an internal helper method that is used to build the properties
	 * node for a given MST entry.
	 * 
	 * @param analyzer The analyzer workspace entry for which the properties sub-tree 
	 * is to be built
	 * 
	 * @return A sub-tree with necessary analyzer properties information to be included 
	 * in the main properties tree at an appropriate location.
	 */
	private static MutableTreeNode makeAnalyzerProperties(FWAnalyzer analyzer) {
		// Create the top-level sub-tree node to hold file information
		DefaultMutableTreeNode fwNode = new DefaultMutableTreeNode(analyzer);
		// Add ID and user supplied description next.
		fwNode.add(new DefaultMutableTreeNode("Type: " + analyzer.getType()));
		fwNode.add(new DefaultMutableTreeNode("Metric: " + (analyzer.isDistance() ? "Distance" : "Similarity")));
		// Add word and frame/window size appropriately
		if (analyzer.getType().equals(FWAnalyzer.FWAnalyzerType.TWOPASSD2)) {
			fwNode.add(new DefaultMutableTreeNode("Frame (aka Window) size: Auto/Adaptive"));
			fwNode.add(new DefaultMutableTreeNode("Word size: Auto/Adaptive"));
		} else {
			fwNode.add(new DefaultMutableTreeNode("Frame (aka Window) size: " + analyzer.getFrameSize()));
			fwNode.add(new DefaultMutableTreeNode("Word size: " + analyzer.getWordSize()));
		}
		// Add cache size and other cache information.
		fwNode.add(new DefaultMutableTreeNode("Cache type: " + analyzer.getCacheType()));
		fwNode.add(new DefaultMutableTreeNode("Cache size: " + analyzer.getCacheSize()));
		// Return the MST information back to the caller
		return fwNode;
	}
	
	/**
	 * Helper method to create a properties node for a given file.
	 * 
	 * This is an internal helper method that is used to build the properties
	 * node for a given file.
	 * 
	 * @param fileName The file name for which the properties node is to be built
	 * 
	 * @param subTreeTitle The title string to be set for the sub-tree being built
	 * by this method. This string is usually in the form <code>"MST File Properties"</code>.
	 *  
	 * @return A sub-tree to be included in the main properties tree at an
	 * appropriate location.
	 */
	private static MutableTreeNode makeFileProperties(String fileName, String subTreeTitle) {
		//  Create a file object to populate file information.
		final File fileInfo = new File(fileName);
		
		// Create the top-level sub-tree node to hold file information
		DefaultMutableTreeNode fileNode = new IconTreeNode("File properties", CustomIcons[FILE_ICON]);
		// Next add full path information to the fileNode.
		fileNode.add(new DefaultMutableTreeNode("Name: " + fileInfo.getName()));
		fileNode.add(new DefaultMutableTreeNode("Path: " + fileInfo.getAbsolutePath()));
		// Add additional information depending on file status
		if (fileInfo.exists()) {
			// Add additional information to the path if file exists.
			fileNode.add(new DefaultMutableTreeNode("Size: " + fileInfo.length() + " bytes"));
			// Indicate if file is read only or not.
			fileNode.add(new DefaultMutableTreeNode("Read only: " + !fileInfo.canWrite()));
		} else {
			fileNode.add(new IconTreeNode("File not found", CustomIcons[ERROR_ICON]));
		}
		// Return the sub-tree with file information back to the caller
		return fileNode;
	}
	
	/**
	 * The default and only constructor.
	 * 
	 * The constructor is private because the methods in this class
	 * are static and they are meant to be used directly.
	 */
	private PropertiesTreeMaker() {
		// Nothing to be done here.
	}
	
	/**
	 * Method to set a custom cell renderer to handle icons better.
	 * 
	 * This method sets a custom tree cell renderer that essentially 
	 * provides a better representative set of Icons to make the 
	 * overall display look a bit prettier.
	 * 
	 * @param tree The tree to which a custom properties renderer is
	 * to be added.
	 */
	private static void addPropertiesRenderer(JTree tree) {
		tree.setCellRenderer(new DefaultTreeCellRenderer() {
			private static final long serialVersionUID = 2130899045423674944L;
			@Override
	        public Component getTreeCellRendererComponent(JTree tree,
	        		Object value, boolean sel, boolean expanded,
	                boolean leaf, int row, boolean hasFocus) {
	        	// Let the base class do the standard processing.
	            super.getTreeCellRendererComponent(tree, value, sel,
	                            expanded, leaf, row, hasFocus);
	            if (value instanceof DefaultMutableTreeNode) {
	            	// No default icons.
	            	setIcon(null);
	            }
	            if (value instanceof IconTreeNode) {
	            	// A custom icon is to be used in this case.
	            	IconTreeNode itn = (IconTreeNode) value;
	            	setIcon(itn.getIcon());
	            	// Make text bold
	            	setText("<html><b>" + getText() + "</b></html>");
	            }
	            return this;
	        }			
		});
	}

	/**
	 * The static set of icons that are repeatedly used by
	 * this class as the user-defined data for selected nodes in the
	 * properties tree created by this class.
	 */
    private static final Icon CustomIcons[] = {
    	Utilities.getIcon("images/16x16/MST.png"),
    	Utilities.getIcon("images/16x16/Cluster.png"),
    	Utilities.getIcon("images/16x16/DataSet.png"),
    	Utilities.getIcon("images/16x16/Job.png"),
    	Utilities.getIcon("images/16x16/File.png"),
    	Utilities.getIcon("images/16x16/Server.png"),
    	Utilities.getIcon("images/16x16/Heuristic.png"),
    	Utilities.getIcon("images/16x16/Filter.png"),
    	Utilities.getIcon("images/16x16/Error.png"),
    	Utilities.getIcon("images/16x16/ClusterSummary.png")
    };
    
    /**
     * A static constant to refer to the MST icon position in the
     * {@link #CustomIcons} array. The constant is used to make
     * the code more readable.
     */
    private static final int MST_ICON = 0;

    /**
     * A static constant to refer to the Cluster icon position in the
     * {@link #CustomIcons} array. The constant is used to make
     * the code more readable.
     */
    private static final int CLUSTER_ICON = 1;

    /**
     * A static constant to refer to the data set icon position in the
     * {@link #CustomIcons} array. The constant is used to make
     * the code more readable.
     */
    private static final int DATASET_ICON = 2;

    /**
     * A static constant to refer to the job icon position in the
     * {@link #CustomIcons} array. The constant is used to make
     * the code more readable.
     */
    private static final int JOB_ICON = 3;

    /**
     * A static constant to refer to the file icon position in the
     * {@link #CustomIcons} array. The constant is used to make
     * the code more readable.
     */
    private static final int FILE_ICON = 4;

    /**
     * A static constant to refer to the server icon position in the
     * {@link #CustomIcons} array. The constant is used to make
     * the code more readable.
     */
    private static final int SERVER_ICON = 5;

    /**
     * A static constant to refer to the heuristics icon position in the
     * {@link #CustomIcons} array. The constant is used to make
     * the code more readable.
     */
    private static final int HEURISTICS_ICON = 6;

    /**
     * A static constant to refer to the filters icon position in the
     * {@link #CustomIcons} array. The constant is used to make
     * the code more readable.
     */
    private static final int FILTERS_ICON = 7;

    /**
     * A static constant to refer to the error icon position in the
     * {@link #CustomIcons} array. The constant is used to make
     * the code more readable.
     */
    private static final int ERROR_ICON = 8;
    
    /**
     * A static constant to refer to the statistics icon position in the
     * {@link #CustomIcons} array. The constant is used to make
     * the code more readable.
     */
    private static final int STATS_ICON = 9;
}
