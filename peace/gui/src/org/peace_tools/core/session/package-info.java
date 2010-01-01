/**
 * Contains the classes associated with providing a session on a local
 * or remote server.
 * 
 * <p>
 * This package contains the classes associated with providing a 
 * session on a local or remote server. The sessions are packaged
 * separately to provide a more rigorous API by ensuring that
 * sessions are created only via the  
 * {@link org.peace_tools.core.session.SessionFactory}.
 * </p> 
 *
 * <p>Note that the remote server sessions are established via 
 * SSH (Secure Shell Protocol). The core SSH protocol is implemented
 * by using a slightly modified (to fix a couple of bugs) the Ganymede
 * SSH library (see <a href="http://www.ganymed.ethz.ch/ssh2/"> 
 * http://www.ganymed.ethz.ch/ssh2</a>). </p>
 * 
 * <p>This wizard has been designed to be consistent with the core
 * infrastructure provided by 
 * {@link org.peace_tools.generic.WizardDialog}.
 * </p>
 *
 * @since 0.9
 */
package org.peace_tools.core.session;
