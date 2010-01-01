/**
 * Contains the classes associated with the server wizard.
 * 
 * <p>
 * This package contains the classes associated with the server
 * wizard. The server wizard guides the user through the process of
 * adding a new server entry to the current work space.
 * </p> 
 *
 * <p>Note that this wizard does not perform the actual installation
 * of the PEACE C++ clustering engine on a server. That task is 
 * delegated to the {@link org.peace_tools.core.PEACEInstaller} class.
 * This wizard actually creates an instance of the installer.</p>
 * 
 * <p>This wizard has been designed to be consistent with the core
 * infrastructure provided by 
 * {@link org.peace_tools.generic.WizardDialog}.
 * </p>
 *
 * @see org.peace_tools.core.PEACEInstaller
 * @see org.peace_tools.generic.WizardDialog
 * @see org.peace_tools.generic.WizardPage
 * 
 * @since 0.9
 */
package org.peace_tools.core.server;
