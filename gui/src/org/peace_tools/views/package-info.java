/**
 * Contains the classes that provide views for the various data
 * objects in the GUI.
 * 
 * <p>
 * These classes constitute the "view" objects as per the terminology
 * using the Model-View-Controller (MVC) design pattern. The classes
 * in the {@link org.peace_tools.data} package use these classes to 
 * provide the data for the views. The classes in the 
 * {@link org.peace_tools.core}
 * package (particularly {@link org.peace_tools.core.ViewFactory}) 
 * serve as the "controller" to establish the link between the "data"
 * and "view" objects. The classes in this package are 
 * can be dependent on the classes in other packages. There is some
 * mutual dependence between this package and the 
 * {@link org.peace_tools.core} package.</p>
 * 
 * <p>
 * The classes in this package provide the actual graphical representation
 * of persistent data in a Workspace (represented by the set of classes 
 * in the {@link org.peace_tools.workspace} package) and 
 * other data (such as EST or clustering information). These class utilize
 * the data set hierarchy maintained by the {@link org.peace_tools.data} package
 * classes to display information it in various views such as JTable, JTree, etc. 
 * </p>
 *
 * <p><b>Note:</b> Typically a single data object can be associated with
 * multiple views, possibly in multiple frames.</p>
 * 
 * @since 0.9
 */
package org.peace_tools.views;
