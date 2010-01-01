/**
 * Contains the classes that provide data for the various views
 * in the GUI.
 * 
 * <p>
 * These classes constitute the "data" objects as per the terminology
 * using the Model-View-Controller (MVC) design pattern. The classes
 * in the {@link org.peace_tools.views} package use these classes to 
 * provide the data for the views. The classes in the {@link org.peace_tools.core}
 * package (particularly {@link org.peace_tools.core.ViewFactory}) 
 * serve as the "controller" to establish the link between the "data"
 * and "view" objects. The classes in this package are 
 * <b>expected to be dependent only on the classes in the workspace 
 * package</b>.</p>
 * 
 * <p>
 * The classes in this package serve as a bridge between the in-memory
 * representation of persistent data in a Workspace (represented by the
 * set of classes in the {@link org.peace_tools.workspace} package) and 
 * other data (such as EST or clustering information). These class enable
 * reusing the data set hierarchy maintained by the 
 * {@link org.peace_tools.workspace.Workspace} object and to display
 * it in various views such as JTable, JTree, etc. In addition, the 
 * classes also acts to monitor work space entries and update data 
 * set views.
 * </p>
 * 
 * @since 0.9
 */
package org.peace_tools.data;
