/**
 * Contains the classes that provide a custom overlap view of
 * clusters.
 * 
 * <p>
 * These classes constitute the "view" objects as per the terminology
 * using the Model-View-Controller (MVC) design pattern. The classes
 * in the {@link org.peace_tools.data} package are used to 
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
 * This package contains a overlap view that provides a visual 
 * representation of the various fragments and clusters that resulted
 * from a clustering operation. The view provides a convenient way to 
 * identify super clusters and possibly other interesting artifacts visually.
 * </p>
 *
 * @since 0.95
 */
package org.peace_tools.views.overlap;
