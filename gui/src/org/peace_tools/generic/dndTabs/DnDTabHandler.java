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

package org.peace_tools.generic.dndTabs;

import java.awt.Color;
import java.awt.Component;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.Polygon;
import java.awt.Toolkit;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.Transferable;
import java.awt.datatransfer.UnsupportedFlavorException;
import java.awt.dnd.DnDConstants;
import java.awt.dnd.DragGestureEvent;
import java.awt.dnd.DragGestureListener;
import java.awt.dnd.DragSource;
import java.awt.dnd.DragSourceContext;
import java.awt.dnd.DragSourceDragEvent;
import java.awt.dnd.DragSourceDropEvent;
import java.awt.dnd.DragSourceEvent;
import java.awt.dnd.DragSourceListener;
import java.awt.dnd.DragSourceMotionListener;
import java.awt.dnd.DropTarget;
import java.awt.dnd.DropTargetDragEvent;
import java.awt.dnd.DropTargetDropEvent;
import java.awt.dnd.DropTargetEvent;
import java.awt.dnd.DropTargetListener;
import java.awt.image.BufferedImage;

import javax.swing.Icon;

/**
 * This class provides the necessary infrastructure for handling Drag-and-Drop
 * features of a Tab from one TabbedPane to another. In other words, this class
 * enables a user to move tabs from one tab to another within the same
 * application.
 * 
 * This implementation has been adapted from the code examples provided in the
 * Java Swing text book.
 */
public class DnDTabHandler implements DragSourceListener, DragGestureListener,
		DropTargetListener, DragSourceMotionListener {
	/**
	 * The only DataFlavor that is currently supported by all of the DnD related
	 * classes and methods in this class. The flavor class essentially defines
	 * the type of data being dragged-n-dropped. In our case we are just
	 * dragging a generic AWT component (which implies it could be any thing
	 * from a simple panel all the way to a complex display).
	 */
	public static final DataFlavor TAB_FLAVOR = new DataFlavor(Component.class,
			"javax.awt.Component");

	/**
	 * This member object hold a reference to the JTabbedPane with which this
	 * DnD handler is associated. This value is set in the constructor and is
	 * never changed during the life time of this class.
	 */
	private final DnDTabbedPane tabPane;

	/**
	 * The drag source object associated with this DnD handler. The drag source
	 * is needed to install a gesture handler and to initiate a drag operation.
	 */
	private DragSource dragSource;

	/**
	 * The tab (JComponent) to be moved via this DnD handler. This value has
	 * been made static so that multiple instances of the DnD handler can access
	 * this value.
	 */
	private transient static Component componentMoved;

	/**
	 * The name of the tab currently being moved. This value is set when the
	 * drag operation commences.
	 */
	private transient static String tabName;

	/**
	 * The icon for the tab currently being moved. This value is set when the
	 * drag operation commences.
	 */
	private transient static Icon tabIcon;

	/**
	 * This instance variable is filled in by the source pane when a drag
	 * operation is started and reset to null once the drag operation is
	 * completed. The context is used to change the cursor depending on its
	 * location in the target drop panel.
	 */
	private transient static DragSourceContext dragContext;

	/**
	 * The 4 block arrow cursors that are internally used by this DnD handler.
	 * These cursors are generated by a static helper method when the first
	 * DnDHandler class is created.
	 */
	private transient static Cursor ArrowCursors[] = null;

	/**
	 * The constructor. The DnDTabbedPane whith which this DnD Handler is
	 * associated must be passed in as the parameter. The constructor creates a
	 * new DragSource and DropTarget for this tree to handle dnd operations.
	 * 
	 * @param tab
	 *            The DnDTabbedPane this DnD handler is associated.
	 */
	public DnDTabHandler(DnDTabbedPane tab) {
		tabPane = tab;
		dragSource = new DragSource();
		// Setup gesture recognizer to permit the owning
		// DndTab to be moved/dragged. When the user tries to
		// drag the tab, methods in this class will be
		// invoked to handle various cases
		dragSource.createDefaultDragGestureRecognizer(tabPane,
				DnDConstants.ACTION_MOVE, this);
		// Add a drop target with the DndTabbedPane to
		// indicate we will accept tabs dropped on us from a
		// DnD operation.
		tabPane.setDropTarget(new DropTarget(tabPane, DnDConstants.ACTION_MOVE,
				this));
		// Add the motion listener to provide visual cues as
		// to where the dragged tab is going to dock.
		dragSource.addDragSourceMotionListener(this);
		dragContext = null;
		// Create static, shared cursors if needed.
		createCursors(tabPane);
	}

	/**
	 * This method is a part of the DragSourceListener interface. This method is
	 * called whenever a drag operation is terminated immaterial of wether it
	 * was successful or not. If the drag operation was indeed successful, then
	 * this method removes the dragged panel from the tabPane (and the panel is
	 * added to another tab by the drop listener).
	 * 
	 * @param dsde
	 *            The source drop event to be processed.
	 */
	@Override
	public void dragDropEnd(DragSourceDropEvent dsde) {
		dragContext = null;
		if (!dsde.getDropSuccess() || (componentMoved == null)) {
			// The dnd was not successful. Nothing further to do.
			return;
		}

		// Remove the tab from ourselves.
		tabPane.removeTab(componentMoved);
	}

	/**
	 * This method is part of the DragSourceMotionListener method. This method
	 * is invoked whenever the mouse is moved during a drag operation. This
	 * method currently just saves the context information for further use.
	 * 
	 */
	@Override
	public void dragMouseMoved(DragSourceDragEvent dsde) {
		dragContext = dsde.getDragSourceContext();
	}

	/**
	 * This method is a part of the DragSourceListener interface. This method is
	 * invoked whenever the user drags a tab onto another location that will
	 * accept the object. In this case, this method changes the cursor to let
	 * the user know that they are a valid drop point.
	 * 
	 * @param dsde
	 *            The source drag event to be processed.
	 */
	@Override
	public void dragEnter(DragSourceDragEvent dsde) {
		// System.out.println("Target action = " + dsde.getTargetActions());
		// dsde.getDragSourceContext().setCursor(DragSource.DefaultMoveDrop);
	}

	/**
	 * This method is a part of the DragSourceListener interface. This method is
	 * invoked whenever the user drags a tab associated with this DnD handler
	 * onto another location that will accept the object. In this case, this
	 * method simply changes the cursor to let the user know that they are at a
	 * valid drop point.
	 * 
	 * @param dsde
	 *            The source drag event to be processed.
	 */
	@Override
	public void dragOver(DragSourceDragEvent dsde) {
		// dsde.getDragSourceContext().setCursor(DragSource.DefaultMoveDrop);
	}

	/**
	 * This method is a part of the DragSourceListener interface. Currently this
	 * method is not used.
	 * 
	 * @param dsde
	 *            This parameter is not used.
	 */
	@Override
	public void dropActionChanged(DragSourceDragEvent dsde) {
		dragEnter(dsde);
	}

	/**
	 * This method is a part of the DragSourceListener interface. This method is
	 * invoked whenever the user drags a tab associated with this DnD handler
	 * onto a location that will <b>not</b> accept the component. In this case,
	 * this method simply changes the cursor to let the user know that they are
	 * at an invalid drop point.
	 * 
	 * @param dse
	 *            The source drag event to be processed.
	 */
	@Override
	public void dragExit(DragSourceEvent dse) {
		dse.getDragSourceContext().setCursor(DragSource.DefaultMoveNoDrop);
	}

	/**
	 * This method is a part of the DragGestureListener interface. This method
	 * is invoked whenever the user starts to drag an object in the tab
	 * associated with this dnd Handler. This method creates a suitable
	 * transferable object and initiates the DnD operation via the drag source.
	 * 
	 * @param dge
	 *            The drag gesture event to be processed.
	 */
	@Override
	public void dragGestureRecognized(DragGestureEvent dge) {
		if (dge.getDragAction() != DnDConstants.ACTION_MOVE) {
			// Sorry, only move is supported.
			return;
		}
		if ((tabPane.getComponentCount() < 1)
				|| (tabPane.getSelectedComponent() == null)) {
			// No component is selected or highlighted. Nothing
			// further to do here.
			return;
		}
		componentMoved = tabPane.getSelectedComponent();
		tabName = tabPane.getTitleAt(tabPane.getSelectedIndex());
		tabIcon = tabPane.getIconAt(tabPane.getSelectedIndex());

		TransferableTab transferable = new TransferableTab();
		dragSource.startDrag(dge, DragSource.DefaultMoveNoDrop, transferable,
				this);
	}

	/**
	 * This method is part of the DropTargetListener. This method is invoked
	 * whenever the user drags a tab onto the DnDTabbedPane object associated
	 * with this class. This method checks to see the location of the cursor
	 * with respect to the bounds of the tabbed node and decides on the accept
	 * action to be performed. If it cannot accept the tab then it indicates it
	 * will <b>not</b> accept it.
	 * 
	 * @param dtde
	 *            The drop target event to be processed.
	 */
	@Override
	public void dragEnter(DropTargetDragEvent dtde) {
		// Check and change the cursor suitably.
		if (dtde.getDropAction() == DnDConstants.ACTION_MOVE) {
			Point location = dtde.getLocation();
			DnDTabbedPane.Location position = getLocation(location.x,
					location.y);
			if (dragContext != null) {
				dragContext.setCursor(ArrowCursors[position.ordinal()]);
				if (!position.equals(tabPane.getDockCue())) {
					tabPane.setDockCue(position);
					tabPane.repaint();
				}
			}
		} else {
			// We don't accept anything other than moves.
			dtde.rejectDrag();
		}
	}

	/**
	 * Helper method used to determine where a given x and y coordinates
	 * (typically the mouse pointer positions) is with respect to the
	 * DnDTabbedPane that owns this handler.
	 * 
	 * @note The return value from this method is meaningful only when the mouse
	 *       pointer is within the bounds of the owning tabbed pane.
	 * 
	 * @param x
	 *            The x-coordinate to test.
	 * @param y
	 *            The y-coordinate to test
	 * @return This method returns one of the five predetermined location
	 *         constants (left, right, top, bottom, or center) depending on the
	 *         position of the mouse.
	 * 
	 */
	protected DnDTabbedPane.Location getLocation(int x, int y) {
		// The number of pixels at the edge of the window that
		// defines the gutter which will cause the window to be
		// placed adjacent rather than within the the dnd tabbed
		// pane.
		final int GutterSize = 25;
		// Check if the
		if (x < GutterSize) {
			return DnDTabbedPane.Location.LEFT;
		} else if (x > (tabPane.getWidth() - GutterSize)) {
			return DnDTabbedPane.Location.RIGHT;
		} else if (y < 20) {
			return DnDTabbedPane.Location.TOP;
		} else if (y > tabPane.getHeight() - GutterSize) {
			return DnDTabbedPane.Location.BOTTOM;
		}
		return DnDTabbedPane.Location.CENTER;
	}

	/**
	 * This method is part of the DropTargetListener. It simply calls the
	 * dragEnter() method to perform the same action that dragEnter() method
	 * does.
	 */
	public void dragOver(DropTargetDragEvent dtde) {
		dragEnter(dtde);
	}

	/**
	 * This method is part of the DropTargetListener. It simply calls the
	 * dargEnter() method to perfrom the same action as it does.
	 * 
	 * @param dtde
	 *            This parameter is not used.
	 */
	public void dropActionChanged(DropTargetDragEvent dtde) {
		dragEnter(dtde);
	}

	/**
	 * This method is part of the DropTargetListener. It currently resets the
	 * drag cue position in the tab so that docking cue box is no longer drawn
	 * by the tab.
	 * 
	 * @param dte
	 *            This parameter is not used.
	 */
	public void dragExit(DropTargetEvent dte) {
		tabPane.setDockCue(null);
		tabPane.repaint();
	}

	/**
	 * This method is part of the DropTargetListener interface. It is invoked
	 * whenever the user drops a Component onto the DnDTabbePane associated with
	 * this handler. This method adds the Component transferred to the
	 * DnDTabbedPane and highlights it.
	 * 
	 * @param dtde
	 *            The drop event to be processed.
	 */
	public void drop(DropTargetDropEvent dtde) {
		// First reset any dock cues on our tab
		tabPane.setDockCue(null);
		if ((componentMoved == null)
				|| (dtde.getDropAction() != DnDConstants.ACTION_MOVE)) {
			// This drop is not acceptable. sorry.
			dtde.rejectDrop();
			return;
		}
		// Now, depending on the location of the mouse, we need
		// add to the tab or create new split panes.
		Point mousePos = dtde.getLocation();
		DnDTabbedPane.Location location = getLocation(mousePos.x, mousePos.y);
		if ((location.equals(DnDTabbedPane.Location.CENTER))
				&& (tabPane == componentMoved.getParent())) {
			// The component was dropped on ourselves. Guess the user
			// want's to shuffle some tabs around.
			tabPane.remove(componentMoved);
			tabPane.addTab(tabName, tabIcon, componentMoved);
			final int index = tabPane.indexOfComponent(componentMoved);
			tabPane.setTabComponentAt(index, new DnDTabButton(tabPane, index));
			tabPane.setSelectedComponent(componentMoved);
			// Report that a move was not really done because the tab did
			// not really move. This will ensure that the tab is not
			// removed from the source by the dragDropEnd() method.
			dtde.dropComplete(false);
		} else {
			// The tab is to be added to the edges or to the center
			// Do that via the DnDTabbedPane method.
			tabPane.createSplitPane(tabName, tabIcon, componentMoved, location);
			// The drop was completed.
			dtde.dropComplete(true);
		}
	}

	/**
	 * This inner class provides a wrapper that impelements Transferable
	 * interface to simply provide a place holder for DnD operations. The actual
	 * component to be moved is stored as a static member in the TabDnDHandler
	 * class. Refer to the documentation on the Transferable interface for
	 * details on the API requirements for each one of the methods in this
	 * class.
	 * 
	 */
	public class TransferableTab implements Transferable {
		/**
		 * An array of flavours that is supported by this Transferable node.
		 * Currently this node supports only one data flavor for DnD operations,
		 * which is just a Component object that indicates the component to be
		 * moved.
		 */
		private transient DataFlavor flavors[] = { TAB_FLAVOR };

		/**
		 * The only constructor for creating a transferable tab. Currently this
		 * class does not require any other information as this is just a empty
		 * place holder.
		 */
		public TransferableTab() {
		}

		/**
		 * This method is a part of the Transferable interface. As per API
		 * reqirements this method returns the set of data flavours supported by
		 * this Transferable object. Currently this object supports only 1 data
		 * flavor.
		 * 
		 * @return The array of data flavors supported by this object.
		 */
		public synchronized DataFlavor[] getTransferDataFlavors() {
			return flavors;
		}

		/**
		 * This method is a part of the Transferable interface. As per API
		 * reqirements this method checks to see if the specified DataFlavor has
		 * the same class as a Component. If so, this method return's true.
		 * Otherwise it return's false.
		 * 
		 * @param flavor
		 *            The DataFlavor to be checked for compatibility.
		 * 
		 * @return This method returns true if flavor.getRepresentationClass()
		 *         is the same as Component.class.
		 */
		public boolean isDataFlavorSupported(DataFlavor flavor) {
			return (flavor.getRepresentationClass() == Component.class);
		}

		/**
		 * This method is a part of the Transferable interface. As per API
		 * reqirements this method first checks to ensure that the DataFlavor is
		 * supported (via isDataFlavorSupported() method) and if so, it returns
		 * the static componentMoved member of the TabDnDHandler. If the flavor
		 * is not supported, it throws an UnsupportedFlavorException.
		 * 
		 * @param flavor
		 *            The DataFlavor corresponding to which an object must be
		 *            obtained from this Transferable class and returned to the
		 *            caller. Currently, only TAB_FLAVOR is supported.
		 * 
		 * @return The Component object contained in the TabDnDHandler class.
		 * 
		 * @exception UnsupportedFlavorException
		 *                This method throws an exception if the DataFlavor is
		 *                not compatible (that flavor.getRepresentationClass()
		 *                is not the same as TreePath.class).
		 */
		public synchronized Object getTransferData(DataFlavor flavor)
				throws UnsupportedFlavorException {
			if (isDataFlavorSupported(flavor)) {
				return DnDTabHandler.componentMoved;
			} else {
				throw new UnsupportedFlavorException(flavor);
			}
		}
	}

	/**
	 * This is a helper method that is used to create some of the custom cursors
	 * used by this class.
	 * 
	 * @param tabPane
	 *            The tabbed pane to operate on
	 */
	private static synchronized void createCursors(DnDTabbedPane tabPane) {
		if (ArrowCursors != null) {
			// Arrow cursors have already been created.
			return;
		}
		Toolkit toolKit = tabPane.getToolkit();
		Dimension cSize = toolKit.getBestCursorSize(16, 16);
		// Create points for a polygon representing a right block arrow.
		int halfW = cSize.width / 2;
		int halfH = cSize.height / 2;
		int quatH = (cSize.height / 2) - (cSize.height / 4);
		int quat3H = (cSize.height / 2) + (cSize.height / 4);

		int[] pointsX = { 0, halfW, halfW, cSize.width, halfW, halfW, 0 };
		int[] pointsY = { quatH, quatH, 0, halfH, cSize.height, quat3H, quat3H };
		int[] hotX = { cSize.width - 1, halfW, 0, halfW };
		int[] hotY = { halfH, cSize.height - 1, halfH, 0 };

		Polygon blockArrow = new Polygon(pointsX, pointsY, pointsX.length);

		// Create the cursors by rotating each one by 90 degrees.
		Cursor tempCursors[] = new Cursor[4];
		for (int curID = 0; (curID < 4); curID++) {
			BufferedImage cursorImage = new BufferedImage(cSize.width,
					cSize.height, BufferedImage.TYPE_4BYTE_ABGR_PRE);
			Graphics2D g = (Graphics2D) cursorImage.getGraphics();
			g.rotate(Math.toRadians(curID * 90), cSize.width / 2,
					cSize.height / 2);
			g.setColor(Color.BLACK);
			g.fillPolygon(blockArrow);
			tempCursors[curID] = toolKit.createCustomCursor(cursorImage,
					new Point(hotX[curID], hotY[curID]), "Cursor"
							+ (curID * 90));
		}
		// Copy data from temp cursors to ArrowCursors in the order
		// of DndTabbedPane.Location enumeration which makes it easy
		// to use them.
		ArrowCursors = new Cursor[5];
		ArrowCursors[DnDTabbedPane.Location.TOP.ordinal()] = tempCursors[3];
		ArrowCursors[DnDTabbedPane.Location.BOTTOM.ordinal()] = tempCursors[1];
		ArrowCursors[DnDTabbedPane.Location.LEFT.ordinal()] = tempCursors[2];
		ArrowCursors[DnDTabbedPane.Location.RIGHT.ordinal()] = tempCursors[0];
		ArrowCursors[DnDTabbedPane.Location.CENTER.ordinal()] = DragSource.DefaultMoveDrop;
	}
}
