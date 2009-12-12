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
// Authors:   Dhananjai M. Rao          raodm@muohio.edu
//
//---------------------------------------------------------------------

package org.peace_tools.generic;

import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.LayoutManager;
import javax.swing.ImageIcon;
import javax.swing.JPanel;

/**
 * A Custom panel that can have an Image for its Background
 * 
 * <p>
 * This class provides a CustomPanel for fancier display. The CustomPanel
 * extends the Java Swing JPanel class and customizes its functionality by
 * provding additional option to use a Image (PNG/JPEG/Gif) as the background
 * for the JPanel. The graphic formats (such as GIF/JPEG etc.) supported will
 * depend on the Java version being used.This image provides a fancy display on
 * top of which the other GUI objects class can be conveniently displayed.
 * 
 * <p>
 * This class provides all the functionality of the JPanel class. Refer to the
 * API documentation on the JPanel class for additional information on the
 * methods that can be invoked on this class.
 * </p>
 * 
 * <p>
 * The image to be used for the background must be set through a suitable call
 * to the setImage() method. By default this class has no image.
 * </p>
 * 
 * @see javax.swing.JPanel
 * 
 * @see javax.swing.ImageIcon
 */
public class CustomPanel extends JPanel {
	/**
	 * Creates a default CustomPanel Object.
	 * 
	 * <p>
	 * This constructor can be used to create a default CustomPanel object. Once
	 * the CustomPanel object has been created, its properties can be further
	 * modified through suitable methods.
	 * </p>
	 * 
	 * <p>
	 * By default the CustomPanel has the same settings as the default JPanel
	 * (super) object. In addition, the default background image is set to null,
	 * indicating that there is no background image set. In this case, the
	 * background color set for the panel is automatically used by the super
	 * class.
	 * </p>
	 */
	public CustomPanel() {
		backgroundImage = null;
	}

	/**
	 * Creates a default CustomPanel Object.
	 * 
	 * <p>
	 * This constructor can be used to create a default CustomPanel object. Once
	 * the CustomPanel object has been created, its properties can be further
	 * modified through suitable methods.
	 * </p>
	 * 
	 * <p>
	 * By default the CustomPanel has the same settings as the default JPanel
	 * (super) object. In addtion, the default background image is set to null,
	 * indicating that there is no background image set. In this case, the
	 * background color set for the panel is automatically used by the super
	 * class.
	 * </p>
	 * 
	 * @param layout The layout manager to be set for this panel.
	 */
	public CustomPanel(LayoutManager layout) {
		super(layout);
		backgroundImage = null;
		useImageSize = true;
	}
	
	/**
	 * Set the image to be used for the background
	 * 
	 * This method can be used to set the image that must be used as the
	 * background by the CustomPanel object. The file name (with necessary path)
	 * must be passed in as the parameter. This method automatically loads the
	 * image from the file and utilizes it for painting the background. The
	 * image can be a GIF or JPEG. Other formats may be supported but
	 * compatibility with a given format depends on the Java version being used.
	 * 
	 * <p>If the image is loaded successfully from the file, this method returns
	 * true. In addtion, it updates the size of the panel to reflect the
	 * dimensions of the image so that the entire image can be displayed. If an
	 * error occurs this method returns false and the CustomContainer reverts
	 * back to its earlier image (if any).</p>
	 * 
	 * @param fileName The name of the image file to be used
	 * 
	 * @return This method returns true on success and false on failure.
	 */
	public boolean setImage(String fileName) {
		// Call to getIcon blocks until the whole image is loaded.
		backgroundImage = Utilities.getIcon(fileName);
		if (backgroundImage != null) {
			setSize(backgroundImage.getIconWidth(), backgroundImage.getIconHeight());
			return true;
		}
		// When control drops here that means the file name was bad.
		return false;
	}
	
	/**
	 * Method to set if the panel should use the background image size
	 * as the preferred size for the panel.
	 * 
	 * @param useImageSize If the parameter is true, then the panel reports
	 * the background image size as the preferred image size.
	 */
	public void setUseImageSize(boolean useImageSize) {
		this.useImageSize = useImageSize;
	}
	
	/**
	 * Obtain the image associated with this panel (if any).
	 * 
	 * @return The image associated with this panel. This method returns
	 * null if an image has not been set.
	 */
	ImageIcon getImage() { return backgroundImage; }

	/**
	 * Paint the background image as and when necessary
	 * 
	 * This method is automatically invoked by the core Java AWT runtime engine
	 * whenever the background for the CustomPanel needs to be repainted. The
	 * Graphics object that must be used to paint the panel is passed in as the
	 * parameter by the core AWT engine.
	 * 
	 * <p>
	 * This method paints the background image, if a valid background image is
	 * available. Otherwise it simply drops control to the super JPanel class.
	 * Refer to the documentation on the paint() method on the JPanel class for
	 * information on the working of the JPanel.</p>
	 * 
	 * @param g The graphics object to be used for painting.
	 * 
	 * @see javax.swing.JPanel
	 * 
	 * @see java.awt.Graphics
	 */
	public void paintComponent(Graphics g) {
		// Let super do its thing first.
		super.paintComponent(g);
		// Draw background image on top if we have one.
		if (backgroundImage != null) {
			// We have a valid background image to work with...
			Dimension mySize = getSize();
			// Figure out where on the panel we start drawing
			int panelX = Math.min(0, mySize.width - backgroundImage.getIconWidth());
			// Figure out the starting x point on the image.
			int imgX = Math.min(0, backgroundImage.getIconWidth() - mySize.width);
			// Draw image right aligned
			g.drawImage(backgroundImage.getImage(), panelX, 0, 
					mySize.width, backgroundImage.getIconHeight(),
					imgX, 0, backgroundImage.getIconWidth(), 
					backgroundImage.getIconHeight(),
					getBackground(), backgroundImage.getImageObserver());
		}
	}

	/**
	 * Provide feedback on the preferred size of this panel
	 * 
	 * This method overrides the default implementation of the
	 * getPreferredSize() method in the super class. This method basically
	 * returns the current size of the CustomPanel. If an image has been
	 * specified for this panel, this method basically returns the size of the
	 * image (unless the size has been modified after the image was set).
	 * 
	 * @return The preferred size for this CustomPanel
	 * 
	 * @see javax.swing.JPanel
	 */
	public Dimension getPreferredSize() {
		if ((backgroundImage != null) && (useImageSize)) {
			return getSize();
		}
		return super.getPreferredSize();
	}
	
	/**
	 * Flag to indicate if the size of the image should be used
	 * as the preferred size for this panel.
	 */
	private boolean useImageSize;
	
	/**
	 * The image used for the background
	 * 
	 * This member object is used to maintain the image that is used as the
	 * background for this panel.
	 */
	private ImageIcon backgroundImage;
	
	/**
	 * The serial version ID associated with this class for serialization.
	 */
	private static final long serialVersionUID = 1L;
}
