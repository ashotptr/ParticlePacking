/*
 * Open Source Physics software is free software as described near the bottom of this code file.
 *
 * For additional information and documentation on Open Source Physics please see:
 * <http://www.opensourcephysics.org/>
 */

package org.opensourcephysics.tools;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.Toolkit;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.StringSelection;
import java.awt.datatransfer.Transferable;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusAdapter;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.font.FontRenderContext;
import java.awt.geom.Rectangle2D;
import java.io.File;
import java.net.URL;
import java.util.ArrayList;
import java.util.HashMap;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.Icon;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.JTextField;
import javax.swing.JTextPane;
import javax.swing.JToolBar;
import javax.swing.JTree;
import javax.swing.SwingConstants;
import javax.swing.ToolTipManager;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.event.HyperlinkEvent;
import javax.swing.event.HyperlinkListener;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;
import javax.swing.filechooser.FileFilter;
import javax.swing.text.Document;
import javax.swing.tree.DefaultTreeCellRenderer;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.TreePath;
import javax.swing.tree.TreeSelectionModel;

import org.opensourcephysics.controls.XML;
import org.opensourcephysics.controls.XMLControl;
import org.opensourcephysics.controls.XMLControlElement;
import org.opensourcephysics.display.OSPRuntime;

/**
 * This is a JPanel that displays an OSP LibraryCollection in a tree.
 *
 * @author Douglas Brown
 */
public class LibraryTreePanel extends JPanel {
	
  protected static Color lightRed = new Color(255, 180, 200);
  protected static Color darkRed = new Color(220, 0, 0);
  protected static Color lightGreen = new Color(100, 200, 100);
  protected static Color defaultForeground;
  protected static Icon openFileIcon;
  protected static HyperlinkListener hyperlinkListener;
  protected static JFileChooser chooser;
  protected static FileFilter htmlFilter, folderFilter;
  protected static HashMap<URL, HTMLPane> htmlPanesByURL
  		= new HashMap<URL, HTMLPane>();
  protected static HashMap<LibraryTreeNode, HTMLPane> htmlPanesByNode
			= new HashMap<LibraryTreeNode, HTMLPane>();
  
	static {
    String imageFile = "/org/opensourcephysics/resources/tools/images/open.gif";        //$NON-NLS-1$
    openFileIcon = new ImageIcon(LibraryTreePanel.class.getResource(imageFile));
    hyperlinkListener = new HyperlinkListener() {
      public void hyperlinkUpdate(HyperlinkEvent e) {
        if(e.getEventType()==HyperlinkEvent.EventType.ACTIVATED) {
          try {
            if(!org.opensourcephysics.desktop.OSPDesktop.browse(e.getURL().toURI())) {
              // try the old way
              org.opensourcephysics.desktop.ostermiller.Browser.init();
              org.opensourcephysics.desktop.ostermiller.Browser.displayURL(e.getURL().toString());
            }
          } catch(Exception ex) {}
        }
      }
    };
	}
	
  // instance fields
	protected LibraryBrowser browser;
	protected LibraryCollection collection; // the collection displayed by this panel
	protected String pathToCollection;
  protected LibraryTreeNode rootNode;
  protected DefaultTreeModel treeModel;
  protected JTree tree;
  protected JScrollPane treeScroller = new JScrollPane();
  protected JScrollPane htmlScroller = new JScrollPane();
  protected JToolBar editorbar;
  protected JButton cutButton, copyButton, pasteButton;
  protected JButton addCollectionButton, addResourceButton;
  protected JButton moveUpButton, moveDownButton;
  protected Box editorPanel, fileBox;
  protected JPanel displayPanel;
  protected HTMLPane emptyHTMLPane;
  protected JSplitPane splitPane;
  protected JTextField nameField, htmlField, basePathField, targetField;
  protected JLabel nameLabel, htmlLabel, basePathLabel, targetLabel;
  protected JLabel typeLabel, typeField;
  protected JButton openHTMLButton, openBasePathButton, openFileButton;
  protected ArrayList<JLabel> labels = new ArrayList<JLabel>();
  protected JPopupMenu popup;
  protected MouseAdapter treeMouseListener, convertPathMouseListener;
  protected TreeSelectionListener treeSelectionListener;
  protected Action cutAction, copyAction, pasteAction;
  protected Action addCollectionAction, addResourceAction;
  protected Action moveUpAction, moveDownAction;
  protected XMLControl pasteControl;
  protected boolean editing, isCollectionChanged;
  protected XMLControl revertControl;
  protected int typeFieldWidth;
  protected String command;

  /**
   * Constructs an empty LibraryTreePanel.
   * 
   * @param browser the LibraryBrowser that will display this panel
   */
  public LibraryTreePanel(LibraryBrowser browser) {
    super(new BorderLayout());
    this.browser = browser;
    createGUI();
  }

  /**
   * Sets the collection displayed in the tree.
   * 
   * @param collection the collection
   * @param path the collection xml file path
   * @param editable true if the collection is user-editable
   */
  public void setCollection(LibraryCollection collection, String path, boolean editable) {
    this.collection = collection;
    pathToCollection = path;
    // clean up existing tree, if any
    if (tree!=null) {
      tree.removeTreeSelectionListener(treeSelectionListener);
      tree.removeMouseListener(treeMouseListener);    	
    }
    setEditing(false);
    if (collection!=null) {
	    // create new tree
		  rootNode = new LibraryTreeNode(collection, this);
		  rootNode.setEditable(editable);
	    createTree(rootNode);
	    tree.setSelectionRow(0);
	    splitPane.setDividerLocation(treeScroller.getPreferredSize().width);
	    isCollectionChanged = false;
	  }
  }

  /**
   * Gets the collection displayed in the tree.
   * 
   * @return the collection
   */
  public LibraryCollection getCollection() {
    return collection;
  }

  /**
   * Gets the selected node.
   * 
   * @return the selected node, or null if none
   */
  protected LibraryTreeNode getSelectedNode() {
  	return (LibraryTreeNode) tree.getLastSelectedPathComponent();
  }
  
  /**
   * Sets the selected node.
   *
   * @param node the node to select
   */
  protected void setSelectedNode(LibraryTreeNode node) {
    if(node!=null) {
      tree.setSelectionPath(new TreePath(node.getPath()));
    }
    else
    	tree.setSelectionPath(null);
  }

  /**
   * Sets the editing state.
   *
   * @param edit true to start editing, false to stop
   */
  protected void setEditing(boolean edit) {
  	editing = edit;
  	if (edit) {      		
  		displayPanel.add(editorPanel, BorderLayout.NORTH);
      add(editorbar, BorderLayout.NORTH);
      showInfo(getSelectedNode());
    }
  	else {
  		displayPanel.remove(editorPanel);
      remove(editorbar);
  	}  	
  	validate();
  	if (editing) {
  		revertControl = new XMLControlElement(collection);
  	}
  }

  /**
   * Returns true if the collection is editable.
   *
   * @return true if editable
   */
  protected boolean isEditable() {
  	boolean editable = rootNode!=null && rootNode.isEditable();
  	if (editable && !pathToCollection.startsWith("http:")) { //$NON-NLS-1$
  		File file = new File(pathToCollection);
  		editable = !file.exists() || file.canWrite();
  	}
  	return editable;
  }
  
  /**
   * Gets the editing state.
   *
   * @return true if editing
   */
  protected boolean isEditing() {
  	return editing;
  }
  
  /**
   * Displays the resource data for the specified node.
   *
   * @param node the LibraryTreeNode
   */
  protected void showInfo(LibraryTreeNode node) {
  	if (node==null) {
  		htmlScroller.setViewportView(emptyHTMLPane);
    	nameField.setText(null);
    	typeField.setText(" "); //$NON-NLS-1$
    	basePathField.setText(null);
    	htmlField.setText(null);
    	targetField.setText(null);
    	nameField.setBackground(Color.white);
    	basePathField.setBackground(Color.white);
    	htmlField.setBackground(Color.white);
    	targetField.setBackground(Color.white);
    	nameField.setEnabled(false);
    	basePathField.setEnabled(false);
    	htmlField.setEnabled(false);
    	targetField.setEnabled(false);
    	typeField.setEnabled(false);
    	nameLabel.setEnabled(false);
    	htmlLabel.setEnabled(false);
    	basePathLabel.setEnabled(false);
    	targetLabel.setEnabled(false);
    	typeLabel.setEnabled(false);
    	openHTMLButton.setEnabled(false);
    	openBasePathButton.setEnabled(false);
    	openFileButton.setEnabled(false);
    	return;
  	}
  	// show node data
  	HTMLPane htmlPane = getHTMLPane(node);
		htmlPane.setCaretPosition(0);
		try {
			htmlScroller.setViewportView(htmlPane);
		} catch (Exception e) {}
		
  	String path = node.isRoot()? pathToCollection: node.getAbsoluteTarget();
  	if (!browser.commandField.getText().equals(path)) {
  		browser.commandField.setText(path);
  		
      // check to see if resource is available
  		if (path!=null && path.startsWith("http://www.compadre.org/osp/")) //$NON-NLS-1$
      	path = "http://www.compadre.org/osp/"; //$NON-NLS-1$
  		Resource res = ResourceLoader.getResourceZipURLsOK(path);
      if (path!=null && res==null) {
      	browser.commandField.setForeground(darkRed);
//      	browser.commandField.setBackground(lightRed);
      }
      else {
	  		browser.commandField.setBackground(Color.white);
	    	browser.commandField.setForeground(defaultForeground);
      }
  		browser.commandField.setCaretPosition(0);
  		if (node.isRoot()) browser.commandButton.setEnabled(false);
  	}
  	// show editor data if editing
  	if (isEditing()) {
	  	if (!nameField.getText().equals(node.getName())) {
	  		nameField.setText(node.getName());
	  		nameField.setCaretPosition(0);
	  	}
	  	String base = basePathField.hasFocus()? node.record.getBasePath(): node.getBasePath();
	  	if (!basePathField.getText().equals(base)) {
	  		basePathField.setText(base);
	  		basePathField.setCaretPosition(0);
	  	}
	  	if (!htmlField.getText().equals(node.record.getHTMLPath())) {
	  		htmlField.setText(node.record.getHTMLPath());
	  		htmlField.setCaretPosition(0);
	  	}
	  	boolean isValidHTML = true;
	  	if (!"".equals(node.record.getHTMLPath())) { //$NON-NLS-1$
	  		isValidHTML = node.getHTMLURL()!=null;
	  	}
			htmlField.setForeground(isValidHTML? defaultForeground: darkRed);
			htmlField.setBackground(Color.white);
//			htmlField.setBackground(isValidHTML? Color.white: lightRed);
	  	if (!targetField.getText().equals(node.getTarget())) {
	  		targetField.setText(node.getTarget());
	  		targetField.setCaretPosition(0);
	  	}
	  	boolean isValidTarget = true;
	  	if (node.getTarget()!=null) {
	  		isValidTarget = node.getTargetURL()!=null;
	  	}
	  	targetField.setForeground(isValidTarget? defaultForeground: darkRed);
	  	targetField.setBackground(Color.white);
//	  	targetField.setBackground(isValidTarget? Color.white: lightRed);
	  	String type = node.record.getType();
	  	type = ToolsRes.getString("LibraryResource.Type."+type); //$NON-NLS-1$
	  	typeField.setText(type);
	  	boolean hasBasePath = !"".equals(node.record.getBasePath()); //$NON-NLS-1$
	  	boolean isCollection = node.record instanceof LibraryCollection;
	  	nameField.setEnabled(true);
	  	basePathField.setEnabled(true);
	  	htmlField.setEnabled(true);
	  	typeField.setEnabled(true);
	  	targetField.setEnabled(!isCollection);
	  	nameLabel.setEnabled(true);
	  	htmlLabel.setEnabled(true);
	  	basePathLabel.setEnabled(true);
	  	targetLabel.setEnabled(!isCollection);
	  	typeLabel.setEnabled(true);
	  	openHTMLButton.setEnabled(true);
	  	openBasePathButton.setEnabled(true);
	  	openFileButton.setEnabled(!isCollection);
	  	basePathField.setForeground(hasBasePath || basePathField.hasFocus()? 
	  			defaultForeground: lightGreen);
	    nameField.setBackground(Color.white);
	    basePathField.setBackground(Color.white);
	    
	    // set tooltips
			path = htmlField.getText();
			if (!path.equals(XML.getPathRelativeTo(path, base))) {
				htmlField.setToolTipText(ToolsRes.getString("LibraryTreePanel.Tooltip.Relative")); //$NON-NLS-1$
			}
			else if (!path.equals(XML.getResolvedPath(path, base))) {
				htmlField.setToolTipText(ToolsRes.getString("LibraryTreePanel.Tooltip.Absolute")); //$NON-NLS-1$
			}
			else htmlField.setToolTipText(null);
			path = targetField.getText();
			if (!path.equals(XML.getPathRelativeTo(path, base))) {
				targetField.setToolTipText(ToolsRes.getString("LibraryTreePanel.Tooltip.Relative")); //$NON-NLS-1$
			}
			else if (!path.equals(XML.getResolvedPath(path, base))) {
				targetField.setToolTipText(ToolsRes.getString("LibraryTreePanel.Tooltip.Absolute")); //$NON-NLS-1$
			}
			else targetField.setToolTipText(null);
  	}
  	repaint();
  }
  
  /**
   * Gets the HTMLPane that describes a given tree node.
   * 
   * @param node the node
   * @return the HTMLPane
   */
  protected HTMLPane getHTMLPane(final LibraryTreeNode node) {
  	HTMLPane htmlPane = null;
  	final URL url = node.getHTMLURL();
  	if (url!=null) {
	  	htmlPane = htmlPanesByURL.get(url);
	  	if (htmlPane==null) {
	  		htmlPane = new HTMLPane();
	  		htmlPane.setText("<h2>"+node+"</h2>"); //$NON-NLS-1$ //$NON-NLS-2$
	  		htmlPanesByURL.put(url, htmlPane);
	  		Runnable runner = new Runnable() {
	  			public void run() {
	  				HTMLPane htmlPane = htmlPanesByURL.get(url);
	      		try {
	      			htmlPane.setPage(url);
	    			} catch (Exception ex) {}
	  			}
	  		};
	  		new Thread(runner).start();
	  	}
	  	else {
	  		if (!url.equals(htmlPane.getPage())) {
		  		Runnable runner = new Runnable() {
		  			public void run() {
		  				HTMLPane htmlPane = htmlPanesByURL.get(url);
			    		try {
				  			htmlPane.getDocument().putProperty(Document.StreamDescriptionProperty, null);
				  			htmlPane.setPage(url);
							} catch (Exception ex) {}	  			
		  			}
		  		};
		  		new Thread(runner).start();
	  		}
	  	}
	  	return htmlPane;
  	}
  	htmlPane = htmlPanesByNode.get(node);
  	if (htmlPane==null) {
  		htmlPane = new HTMLPane();
  		htmlPanesByNode.put(node, htmlPane);
  		String htmlStr = node.getHTMLString();  		
  		htmlPane.setText(htmlStr);
  	}
  	return htmlPane;
  }
  
  /**
   * Creates the GUI and listeners.
   */
  protected void createGUI() {
    // create popup menu
    popup = new JPopupMenu();
    // create actions
		addCollectionAction = new AbstractAction() {
	    public void actionPerformed(ActionEvent e) {
	    	LibraryTreeNode node = getSelectedNode();
    		htmlPanesByNode.remove(node);
	    	LibraryCollection collection = (LibraryCollection)node.record;
      	LibraryCollection newCollection = new LibraryCollection(null);
      	collection.addResource(newCollection);
      	LibraryTreeNode newNode = new LibraryTreeNode(newCollection, LibraryTreePanel.this);
      	insertChildAt(newNode, node, node.getChildCount());
      	collectionChanged();
	    }
		};
		addResourceAction = new AbstractAction() {
	    public void actionPerformed(ActionEvent e) {
	    	LibraryTreeNode node = getSelectedNode();
    		htmlPanesByNode.remove(node);
	    	LibraryCollection collection = (LibraryCollection)node.record;
	    	LibraryResource record = new LibraryResource(null);
	    	collection.addResource(record);
	    	LibraryTreeNode newNode = new LibraryTreeNode(record, LibraryTreePanel.this);
	    	insertChildAt(newNode, node, node.getChildCount());
      	collectionChanged();
	    }
		};
		copyAction = new AbstractAction() {
		  public void actionPerformed(ActionEvent e) {
      	LibraryTreeNode node = getSelectedNode();
      	if (node!=null) {
	        XMLControl control = new XMLControlElement(node.record);
	        String target = XML.forwardSlash(node.getTarget());
	        // set base path if target not absolute and treePanel not editing
	        if(!isEditing() && !target.startsWith("/") && target.indexOf(":")==-1) { //$NON-NLS-1$ //$NON-NLS-2$
		        control.setValue("base_path", node.getBasePath()); //$NON-NLS-1$
	        }
	        StringSelection data = new StringSelection(control.toXML());
	        Clipboard clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
	        clipboard.setContents(data, data);
	        enableButtons();
		    }
		  }		
		};
	  cutAction = new AbstractAction() {
		  public void actionPerformed(ActionEvent e) {
	    	LibraryTreeNode node = getSelectedNode();
	    	if (node!=null) {
	    		copyAction.actionPerformed(null);
	    		removeNode(node);
	        enableButtons();
	      	collectionChanged();
	    	}
		  }
		};
		pasteAction = new AbstractAction() {
		  public void actionPerformed(ActionEvent e) {
      	LibraryTreeNode parent = getSelectedNode();
      	if (parent==null || !(parent.record instanceof LibraryCollection))
      		return;      	
		    try {
		      Clipboard clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
		      Transferable data = clipboard.getContents(null);
		      String dataString = (String) data.getTransferData(DataFlavor.stringFlavor);
		      if(dataString!=null) {
		        XMLControlElement control = new XMLControlElement();
		        control.readXML(dataString);
		        if (LibraryResource.class.isAssignableFrom(control.getObjectClass())) {
		      		htmlPanesByNode.remove(parent);
		        	LibraryResource record =(LibraryResource) control.loadObject(null);
		      		LibraryCollection collection = (LibraryCollection)parent.record;
	          	collection.addResource(record);
	          	LibraryTreeNode newNode = new LibraryTreeNode(record, LibraryTreePanel.this);
	          	insertChildAt(newNode, parent, parent.getChildCount());
	          	collectionChanged();
		        }
		      }
		    } catch(Exception ex) {}
		    
        enableButtons();
		  }		
		};
		moveUpAction = new AbstractAction() {
		  public void actionPerformed(ActionEvent e) {
      	LibraryTreeNode node = getSelectedNode();
      	if (node!=null) {
      		LibraryTreeNode parent = (LibraryTreeNode) node.getParent();
	        if(parent!=null) {
		        int i = parent.getIndex(node);
		        if(i>0) {
		      		htmlPanesByNode.remove(parent);
		        	treeModel.removeNodeFromParent(node);
		        	treeModel.insertNodeInto(node, parent, i-1);
		        	LibraryCollection collection = (LibraryCollection)parent.record;
		        	collection.removeResource(node.record);
		        	collection.insertResource(node.record, i-1);
		        	setSelectedNode(node);
		        	enableButtons();
		        	collectionChanged();
		        }
	        }
		    }
		  }		
		};
		moveDownAction = new AbstractAction() {
		  public void actionPerformed(ActionEvent e) {
      	LibraryTreeNode node = getSelectedNode();
      	if (node!=null) {
      		LibraryTreeNode parent = (LibraryTreeNode) node.getParent();
	        if(parent!=null) {
	          int i = parent.getIndex(node);
	          int end = parent.getChildCount();
		        if(i<end-1) {
		      		htmlPanesByNode.remove(parent);
		        	treeModel.removeNodeFromParent(node);
		        	treeModel.insertNodeInto(node, parent, i+1);
		        	LibraryCollection collection = (LibraryCollection)parent.record;
		        	collection.removeResource(node.record);
		        	collection.insertResource(node.record, i+1);
		        	setSelectedNode(node);
		        	enableButtons();
		        	collectionChanged();
		        }
	        }
		    }
		  }		
		};

		convertPathMouseListener = new MouseAdapter() {
      public void mouseClicked(MouseEvent e) {
      	EntryField field = (EntryField)e.getSource();
      	String path = field.getText();
      	if ("".equals(path)) return; //$NON-NLS-1$
      	final boolean isTarget = field==targetField;
      	
        final LibraryTreeNode node = getSelectedNode();
        if (node !=null && OSPRuntime.isPopupTrigger(e)) {        	
          String base = node.getBasePath();
          if ("".equals(base)) return; //$NON-NLS-1$
        	JPopupMenu popup = new JPopupMenu();
      		final String relPath = XML.getPathRelativeTo(path, base);
      		final String absPath = XML.getResolvedPath(path, base);
      		if (!path.equals(relPath)) {
            JMenuItem item = new JMenuItem(ToolsRes.getString("LibraryTreePanel.MenuItem.SetToRelative")); //$NON-NLS-1$
            popup.add(item);
            item.addActionListener(new ActionListener() {
              public void actionPerformed(ActionEvent e) {
              	if (isTarget)
              		node.setTarget(relPath);
              	else node.setHTMLPath(relPath);
              }
            });
      		}
      		else if (!path.equals(absPath)) {
            JMenuItem item = new JMenuItem(ToolsRes.getString("LibraryTreePanel.MenuItem.SetToAbsolute")); //$NON-NLS-1$
            popup.add(item);
            item.addActionListener(new ActionListener() {
              public void actionPerformed(ActionEvent e) {
              	if (isTarget)
              		node.setTarget(absPath);
              	else node.setHTMLPath(absPath);             	
              }
            });
      		}
        	if (popup.getComponentCount()>0)
        		popup.show(field, e.getX()+2, e.getY()+2);
        }
      }
    };

    // create tree listeners
    treeSelectionListener = new TreeSelectionListener() {
      public void valueChanged(TreeSelectionEvent e) {
	    	LibraryTreeNode node = getSelectedNode();
      	showInfo(node);
        enableButtons();
      }
    };
    treeMouseListener = new MouseAdapter() {
      public void mouseClicked(MouseEvent e) {
        // select node and show popup menu
        TreePath path = tree.getPathForLocation(e.getX(), e.getY());
        if(path==null) {
          return;
        }
        tree.setSelectionPath(path);
        LibraryTreeNode node = (LibraryTreeNode) tree.getLastSelectedPathComponent();
        if (OSPRuntime.isPopupTrigger(e)) {
          getPopup(node).show(tree, e.getX(), e.getY()+8);
        }
        else if (e.getClickCount()==2 && node.getTarget()!=null) {
        	firePropertyChange("target", null, node); //$NON-NLS-1$
        }
      }
    };
    // create toolbar and buttons
    addCollectionButton = new JButton();
    addCollectionButton.setOpaque(false);
    addCollectionButton.setBorder(LibraryBrowser.buttonBorder);
    addCollectionButton.addActionListener(addCollectionAction);
    addResourceButton = new JButton(addResourceAction);
    addResourceButton.setOpaque(false);
    addResourceButton.setBorder(LibraryBrowser.buttonBorder);
    copyButton = new JButton(copyAction);
    copyButton.setOpaque(false);
    copyButton.setBorder(LibraryBrowser.buttonBorder);
    cutButton = new JButton(cutAction);
    cutButton.setOpaque(false);
    cutButton.setBorder(LibraryBrowser.buttonBorder);
    pasteButton = new JButton(pasteAction);
    pasteButton.setOpaque(false);
    pasteButton.setBorder(LibraryBrowser.buttonBorder);
    moveUpButton = new JButton(moveUpAction);
    moveUpButton.setOpaque(false);
    moveUpButton.setBorder(LibraryBrowser.buttonBorder);
    moveDownButton = new JButton(moveDownAction);
    moveDownButton.setOpaque(false);
    moveDownButton.setBorder(LibraryBrowser.buttonBorder);
    
    editorbar = new JToolBar();
    editorbar.setFloatable(false);
    editorbar.setBorder(BorderFactory.createEmptyBorder(1, 0, 1, 0));
    editorbar.add(addResourceButton);
    editorbar.add(addCollectionButton);
    editorbar.addSeparator();
    editorbar.add(copyButton);
    editorbar.add(cutButton);
    editorbar.add(pasteButton);
    editorbar.addSeparator();
    editorbar.add(moveUpButton);
    editorbar.add(moveDownButton);
    
    // create default html pane
    emptyHTMLPane = new HTMLPane();    
    // create display panel for right side of split pane
    displayPanel = new JPanel(new BorderLayout());
    displayPanel.add(htmlScroller, BorderLayout.CENTER);
    
    // create split pane
    splitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, treeScroller, displayPanel);
    add(splitPane, BorderLayout.CENTER); 
    treeScroller.setPreferredSize(new Dimension(240, 400));
    
    // create editorPanel and components
    editorPanel = Box.createVerticalBox();

    nameField = new EntryField();
    nameField.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
      	LibraryTreeNode node = getSelectedNode();
      	if (node!=null) {
      		htmlPanesByNode.remove(node);
      		LibraryTreeNode parent = (LibraryTreeNode)node.getParent();
      		if (parent!=null) htmlPanesByNode.remove(parent);
      		node.setName(nameField.getText());
      		if (node.isRoot()) {
      			browser.refreshTabTitle(pathToCollection, collection);
      		}
      	}
      }
    });
    
    // create typeField
    typeField = new JLabel(" ") {  //$NON-NLS-1$
    	public Dimension getPreferredSize() {
    		Dimension dim = nameField.getPreferredSize();
    		dim.width = typeFieldWidth;
    		return dim;
    	}
    };
    typeField.setBorder(nameField.getBorder());
    typeField.setBackground(nameField.getBackground());
    typeField.setFont(nameField.getFont());
    typeField.setOpaque(true);
    typeField.addMouseListener(new MouseAdapter() {
    	public void mousePressed(MouseEvent e) {
      	final LibraryTreeNode node = getSelectedNode();
      	if (node!=null && !(node.record instanceof LibraryCollection)) {
        	JPopupMenu popup = new JPopupMenu();
        	ActionListener typeListener = new ActionListener() {
      		  public void actionPerformed(ActionEvent e) {
      		  	String type = e.getActionCommand();
      		  	if (!type.equals(node.record.getType())) {
            		htmlPanesByNode.remove(node);
            		LibraryTreeNode parent = (LibraryTreeNode)node.getParent();
            		if (parent!=null) htmlPanesByNode.remove(parent);
      		  		node.setType(type);
      		  		type = ToolsRes.getString("LibraryResource.Type."+node.record.getType()); //$NON-NLS-1$
      		  		typeField.setText(type);
    	        	collectionChanged();
    	        	showInfo(node);
      		  	}
            }
          };
        	for (String next: LibraryResource.resourceTypes) {
            JMenuItem item = new JMenuItem(ToolsRes.getString("LibraryResource.Type."+next)) { //$NON-NLS-1$
            	public Dimension getPreferredSize() {
            		Dimension dim = typeField.getPreferredSize();
            		dim.width -= 2;
            		return dim;
            	}
            };
            popup.add(item);
            item.addActionListener(typeListener);
            item.setActionCommand(next);
        	}
          popup.show(typeField, 0, typeField.getHeight());      		
      	}
    	}
    });

    htmlField = new EntryField();
    htmlField.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
      	LibraryTreeNode node = getSelectedNode();
      	if (node!=null) { 
      		node.setHTMLPath(htmlField.getText());
      	}
      }
    });
    htmlField.addMouseListener(convertPathMouseListener);
    
    openHTMLButton = new JButton(openFileIcon);
    openHTMLButton.setBorder(BorderFactory.createEmptyBorder(1, 1, 1, 2));
    openHTMLButton.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
      	LibraryTreeNode node = getSelectedNode();
      	if (node!=null) { 
          int result = JFileChooser.CANCEL_OPTION;
          JFileChooser chooser = getFileChooser();
          chooser.setDialogTitle(null);
          chooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
          chooser.setAcceptAllFileFilterUsed(true);
          chooser.addChoosableFileFilter(htmlFilter);
      	  result = chooser.showOpenDialog(LibraryTreePanel.this);
          File file = chooser.getSelectedFile();
      	  chooser.removeChoosableFileFilter(htmlFilter); 
          chooser.setSelectedFile(new File(""));  //$NON-NLS-1$
          if(result==JFileChooser.APPROVE_OPTION) {
            browser.library.chooserDir = chooser.getCurrentDirectory().toString();
        		if (file!=null) {
        			String path = XML.forwardSlash(file.getAbsolutePath());
              String base = node.getBasePath();
              if (!"".equals(base)) { //$NON-NLS-1$
              	path = XML.getPathRelativeTo(path, base);
              }
        			node.setHTMLPath(path);
        		}
          }
      	}
      }
    });
    
    basePathField = new EntryField();
    basePathField.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
      	// special handling to prevent setting base path if same as inherited
      	LibraryTreeNode node = getSelectedNode();
      	if (basePathField.getBackground()!=Color.yellow) {
      		return;
      	}
      	if (node!=null) {
      		if (!basePathField.getText().equals(node.record.getBasePath())) {
        		htmlPanesByNode.remove(node);
        		LibraryTreeNode parent = (LibraryTreeNode)node.getParent();
        		if (parent!=null) htmlPanesByNode.remove(parent);
        		node.setBasePath(basePathField.getText());
          	collectionChanged();
        	}
      	}
      }
    });
    basePathField.addFocusListener(new FocusAdapter() {
      public void focusGained(FocusEvent e) {
      	LibraryTreeNode node = getSelectedNode();
      	if ("".equals(node.record.getBasePath())) { //$NON-NLS-1$
      		basePathField.setText(null);
        	basePathField.setForeground(htmlField.getForeground());
        	basePathField.setBackground(Color.white);
      	}
      }
    });
    openBasePathButton = new JButton(openFileIcon);
    openBasePathButton.setBorder(BorderFactory.createEmptyBorder(1, 1, 1, 2));
    openBasePathButton.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
      	LibraryTreeNode node = getSelectedNode();
      	if (node!=null) { 
          int result = JFileChooser.CANCEL_OPTION;
          JFileChooser chooser = getFileChooser();
          chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
          chooser.setAcceptAllFileFilterUsed(false);
          chooser.addChoosableFileFilter(folderFilter);
          chooser.setDialogTitle(ToolsRes.getString("LibraryTreePanel.FileChooser.Title.Base")); //$NON-NLS-1$
      	  result = chooser.showDialog(LibraryTreePanel.this, 
      	  		ToolsRes.getString("LibraryTreePanel.FileChooser.Button.Select")); //$NON-NLS-1$
          File file = chooser.getSelectedFile();
      	  chooser.removeChoosableFileFilter(folderFilter); 
          chooser.setSelectedFile(new File(""));  //$NON-NLS-1$
          if(result==JFileChooser.APPROVE_OPTION) {
            browser.library.chooserDir = chooser.getCurrentDirectory().toString();
        		if (file!=null) {
          		htmlPanesByNode.remove(node);
          		LibraryTreeNode parent = (LibraryTreeNode)node.getParent();
          		if (parent!=null) htmlPanesByNode.remove(parent);
        			node.setBasePath(XML.forwardSlash(file.getAbsolutePath()));
        			collectionChanged();
        			showInfo(node);
        		}
          }
      	}
      }
    });

    targetField = new EntryField();
    targetField.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
      	LibraryTreeNode node = getSelectedNode();
      	if (node!=null && !(node.record instanceof LibraryCollection)) {
      		node.setTarget(targetField.getText());
      	}
      }
    });
    targetField.addMouseListener(convertPathMouseListener);
    
    openFileButton = new JButton(openFileIcon);
    openFileButton.setBorder(BorderFactory.createEmptyBorder(1, 1, 1, 2));
    openFileButton.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
      	LibraryTreeNode node = getSelectedNode();
      	if (node!=null) {
          int result = JFileChooser.CANCEL_OPTION;
          JFileChooser chooser = getFileChooser();
          chooser.setDialogTitle(null);
          chooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
          chooser.setAcceptAllFileFilterUsed(true);
      	  result = chooser.showOpenDialog(LibraryTreePanel.this);
          File file = chooser.getSelectedFile();
          chooser.setSelectedFile(new File(""));  //$NON-NLS-1$
          if(result==JFileChooser.APPROVE_OPTION) {
            browser.library.chooserDir = chooser.getCurrentDirectory().toString();
        		if (file!=null) {
        			String path = XML.forwardSlash(file.getAbsolutePath());
              String base = node.getBasePath();
              if (!"".equals(base)) { //$NON-NLS-1$
              	path = XML.getPathRelativeTo(path, base);
              }
          		node.setTarget(path);
        		}
          }
      	}
      }
    });

    nameLabel = new JLabel();
    nameLabel.setBorder(BorderFactory.createEmptyBorder(2, 2, 2, 2));
    nameLabel.setHorizontalAlignment(SwingConstants.TRAILING);
    typeLabel = new JLabel();
    typeLabel.setBorder(BorderFactory.createEmptyBorder(2, 6, 2, 2));
    typeLabel.setHorizontalAlignment(SwingConstants.TRAILING);
    htmlLabel = new JLabel();
    htmlLabel.setBorder(BorderFactory.createEmptyBorder(2, 2, 2, 2));
    htmlLabel.setHorizontalAlignment(SwingConstants.TRAILING);
    basePathLabel = new JLabel();
    basePathLabel.setBorder(BorderFactory.createEmptyBorder(2, 2, 2, 2));
    basePathLabel.setHorizontalAlignment(SwingConstants.TRAILING);
    targetLabel = new JLabel();
    targetLabel.setBorder(BorderFactory.createEmptyBorder(2, 2, 2, 2));
    targetLabel.setHorizontalAlignment(SwingConstants.TRAILING);
    labels.add(nameLabel);
    labels.add(htmlLabel);
    labels.add(basePathLabel);
    labels.add(targetLabel);

    Box box = Box.createHorizontalBox();
    box.setBorder(BorderFactory.createEmptyBorder(2, 0, 0, 2));
    box.add(nameLabel);
    box.add(nameField);    
    box.add(typeLabel);    
    box.add(typeField);    
    editorPanel.add(box);
    
    box = Box.createHorizontalBox();
    box.setBorder(BorderFactory.createEmptyBorder(2, 0, 0, 2));
    box.add(htmlLabel);
    box.add(htmlField);    
    box.add(openHTMLButton);    
    editorPanel.add(box);

    box = Box.createHorizontalBox();
    box.setBorder(BorderFactory.createEmptyBorder(2, 0, 0, 2));
    box.add(basePathLabel);
    box.add(basePathField);    
    box.add(openBasePathButton);    
    editorPanel.add(box);

    fileBox = Box.createHorizontalBox();
    fileBox.setBorder(BorderFactory.createEmptyBorder(2, 0, 2, 1));
    fileBox.add(targetLabel);
    fileBox.add(targetField);    
    fileBox.add(openFileButton);    
    editorPanel.add(fileBox);
  }
  
  /**
   * Refreshes the GUI including locale-dependent resource strings.
   */
  protected void refreshGUI() {
  	// set button and label text
  	addCollectionButton.setText(ToolsRes.getString("LibraryTreePanel.Button.AddCollection")); //$NON-NLS-1$
  	addResourceButton.setText(ToolsRes.getString("LibraryTreePanel.Button.AddResource")); //$NON-NLS-1$
  	copyButton.setText(ToolsRes.getString("LibraryTreePanel.Button.Copy")); //$NON-NLS-1$
  	cutButton.setText(ToolsRes.getString("LibraryTreePanel.Button.Cut")); //$NON-NLS-1$
  	pasteButton.setText(ToolsRes.getString("LibraryTreePanel.Button.Paste")); //$NON-NLS-1$
  	moveUpButton.setText(ToolsRes.getString("LibraryTreePanel.Button.Up")); //$NON-NLS-1$
  	moveDownButton.setText(ToolsRes.getString("LibraryTreePanel.Button.Down")); //$NON-NLS-1$
  	nameLabel.setText(ToolsRes.getString("LibraryTreePanel.Label.Name")); //$NON-NLS-1$
  	typeLabel.setText(ToolsRes.getString("LibraryTreePanel.Label.Type")); //$NON-NLS-1$
  	htmlLabel.setText(ToolsRes.getString("LibraryTreePanel.Label.HTML")); //$NON-NLS-1$
  	basePathLabel.setText(ToolsRes.getString("LibraryTreePanel.Label.BasePath")); //$NON-NLS-1$
  	targetLabel.setText(ToolsRes.getString("LibraryTreePanel.Label.TargetFile")); //$NON-NLS-1$

		// adjust size of labels so they right-align
    int w = 0;
    Font font = nameLabel.getFont();
    FontRenderContext frc = new FontRenderContext(null, false, false); 
    for(JLabel next: labels) {
      Rectangle2D rect = font.getStringBounds(next.getText()+" ", frc); //$NON-NLS-1$
      w = Math.max(w, (int) rect.getWidth()+4);
    }
    Dimension labelSize = new Dimension(w, 20);
    for(JLabel next: labels) {
      next.setPreferredSize(labelSize);
      next.setMinimumSize(labelSize);
    }
    // determine required size of type label
    typeFieldWidth = 0;
    for (String next: LibraryResource.resourceTypes) {
    	next = ToolsRes.getString("LibraryResource.Type."+next); //$NON-NLS-1$
      Rectangle2D rect = font.getStringBounds(next+" ", frc); //$NON-NLS-1$
      typeFieldWidth = Math.max(typeFieldWidth, (int) rect.getWidth()+24);    	
    }
  }

  /**
   * Enables/disables buttons based on selected node and clipboard state.
   */
  protected void enableButtons() {
  	LibraryTreeNode node = getSelectedNode();
  	boolean nodeIsCollection = node!=null && node.record instanceof LibraryCollection;
  	addCollectionButton.setEnabled(nodeIsCollection);
  	addResourceButton.setEnabled(nodeIsCollection);
  	copyButton.setEnabled(node!=null);
  	cutButton.setEnabled(node!=null);
  	pasteButton.setEnabled(nodeIsCollection && isClipboardPastable());
		boolean canMoveUp = false, canMoveDown = false;
		if (node!=null && node.getParent()!=null) {
			LibraryTreeNode parent = (LibraryTreeNode) node.getParent();
      int i = parent.getIndex(node);
      canMoveUp = i>0;
      canMoveDown = i<parent.getChildCount()-1;
    }
  	moveUpButton.setEnabled(canMoveUp);
  	moveDownButton.setEnabled(canMoveDown);
  }
  
  /**
   * Discards collection edits and reverts to the previous state.
   */
  protected void revert() {
  	if (revertControl!=null) {
  		// copy revertControl to ensure new library resources
  		revertControl = new XMLControlElement(revertControl); 
  		LibraryCollection collection = (LibraryCollection)revertControl.loadObject(null);
  		isCollectionChanged = false;
  		setCollection(collection, pathToCollection, rootNode.isEditable());
  	}
  }

  /**
   * Creates the tree.
   * 
   * @param root the root node
   */
  protected void createTree(LibraryTreeNode root) {
    treeModel = new DefaultTreeModel(root);
    tree = new JTree(treeModel);
    root.createChildNodes();
    tree.setCellRenderer(new LibraryTreeNodeRenderer());
    tree.getSelectionModel().setSelectionMode(TreeSelectionModel.SINGLE_TREE_SELECTION);
    ToolTipManager.sharedInstance().registerComponent(tree);
    // listen for tree selections and display the contents
    tree.addTreeSelectionListener(treeSelectionListener);
    // listen for mouse events to display node info and inform propertyChangeListeners
    tree.addMouseListener(treeMouseListener);
    // put tree in scroller
    treeScroller.setViewportView(tree);
  }
  
  /**
   * Determines if the clipboard can be pasted.
   * 
   * @return true if the clipboard contains a LibraryTreeNode XMLControl string
   */
  protected boolean isClipboardPastable() {
  	pasteControl = null;
    Clipboard clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
    Transferable data = clipboard.getContents(null);
    String dataString = null;
		try {
			dataString = (String)data.getTransferData(DataFlavor.stringFlavor);
		} catch (Exception e) {} 
		if(dataString!=null) {
      XMLControlElement control = new XMLControlElement();
      control.readXML(dataString);
      if (LibraryResource.class.isAssignableFrom(control.getObjectClass())) {
      	pasteControl = control;
      	return true;
      }
    }
    return false;
  }
  
  /**
   * Returns a popup menu with items appropriate for a given tree node.
   * 
   * @param node the node
   * @return the popup menu
   */
  protected JPopupMenu getPopup(final LibraryTreeNode node) {
  	popup.removeAll(); 
  	if (!isEditing()) {
   		JMenuItem item = new JMenuItem(ToolsRes.getString("LibraryTreePanel.Button.Copy")); //$NON-NLS-1$
      popup.add(item);
      item.addActionListener(copyAction);
      return popup;
  	}
  	boolean isCollection = node.record instanceof LibraryCollection;
		boolean canMoveUp = false, canMoveDown = false;
		if (node.getParent()!=null) {
			LibraryTreeNode parent = (LibraryTreeNode) node.getParent();
      int i = parent.getIndex(node);
      canMoveUp = i>0;
      canMoveDown = i<parent.getChildCount()-1;
    }
		if (isCollection) {
  		// add resource
      JMenuItem item = new JMenuItem(ToolsRes.getString("LibraryTreePanel.Button.AddResource")); //$NON-NLS-1$
      popup.add(item);
      item.addActionListener(addResourceAction);
  		// add collection
      item = new JMenuItem(ToolsRes.getString("LibraryTreePanel.Button.AddCollection")); //$NON-NLS-1$
      popup.add(item);
      item.addActionListener(addCollectionAction);
      popup.addSeparator();			
		}
 		JMenuItem item = new JMenuItem(ToolsRes.getString("LibraryTreePanel.Button.Copy")); //$NON-NLS-1$
    popup.add(item);
    item.addActionListener(copyAction);
 		item = new JMenuItem(ToolsRes.getString("LibraryTreePanel.Button.Cut")); //$NON-NLS-1$
    popup.add(item);
    item.addActionListener(cutAction);
		if (isCollection && isClipboardPastable()) {
   		item = new JMenuItem(ToolsRes.getString("LibraryTreePanel.Button.Paste")); //$NON-NLS-1$
      popup.add(item);
      item.addActionListener(pasteAction);  		
		}
    if (canMoveUp || canMoveDown) {
      popup.addSeparator();
      if (canMoveUp) {
	   		item = new JMenuItem(ToolsRes.getString("LibraryTreePanel.Button.Up")); //$NON-NLS-1$
	      popup.add(item);
	      item.addActionListener(moveUpAction);
      }
      if (canMoveDown) {
	   		item = new JMenuItem(ToolsRes.getString("LibraryTreePanel.Button.Down")); //$NON-NLS-1$
	      popup.add(item);
	      item.addActionListener(moveDownAction);
      }
    }
    return popup;
  }

  
  /**
   * Inserts a child into a parent node at a specified index.
   *
   * @param child the child node
   * @param parent the parent node
   * @param index the index
   * @return true if added
   */
  protected boolean insertChildAt(LibraryTreeNode child, LibraryTreeNode parent, int index) {
  	if (tree==null || parent.getChildCount()<index) return false;
  	DefaultTreeModel model = (DefaultTreeModel)tree.getModel();
  	model.insertNodeInto(child, parent, index);
  	TreePath path = new TreePath(child.getPath());
  	tree.scrollPathToVisible(path);
    tree.setSelectionPath(path);
    return true;
  }

  /**
   * Removes a given tree node.
   *
   * @param node the node
   */
  protected void removeNode(LibraryTreeNode node) {
  	if (rootNode==null || node==rootNode) return;
  	LibraryTreeNode parent = (LibraryTreeNode)node.getParent();
		htmlPanesByNode.remove(parent);
		htmlPanesByNode.remove(node);
  	LibraryCollection collection = (LibraryCollection)parent.record;
  	collection.removeResource(node.record);
  	DefaultTreeModel model = (DefaultTreeModel)tree.getModel();
		model.removeNodeFromParent(node);
  	TreePath path = new TreePath(parent.getPath());
  	tree.scrollPathToVisible(path);
    tree.setSelectionPath(path);
  }

  /**
   * Called whenever the collection changes due to a user edit.
   */
  protected void collectionChanged() {
  	isCollectionChanged = true;
		firePropertyChange("collection_edit", null, null); //$NON-NLS-1$
  }
  
  /**
   * Determines if the collection has been changed since the last save.
   * @return true if changed
   */
  protected boolean isCollectionChanged() {
  	return isEditable() && isCollectionChanged;
  }
  
  /**
   * Saves the current collection.
   */
  protected void save() {
  	if (isEditable()) {
			XMLControl control = new XMLControlElement(collection);
			control.write(pathToCollection);
			isCollectionChanged = false;
		}
  }

  /**
   * Gives the user an opportunity to save changes.
   * @param name the name of the collection
   * @return <code>false</code> if the user cancels, otherwise <code>true</code>
   */
  protected boolean saveChanges(String name) {
    if (!isCollectionChanged()) return true;
    if (org.opensourcephysics.display.OSPRuntime.applet!=null) return true;
    int i = JOptionPane.showConfirmDialog(this,
         ToolsRes.getString("LibraryBrowser.Dialog.SaveChanges.Message") + " \"" + name + "\"?", //$NON-NLS-1$ //$NON-NLS-2$ //$NON-NLS-3$
         ToolsRes.getString("LibraryBrowser.Dialog.SaveChanges.Title"), //$NON-NLS-1$
         JOptionPane.YES_NO_CANCEL_OPTION,
         JOptionPane.QUESTION_MESSAGE);
    if (i==JOptionPane.CLOSED_OPTION || i==JOptionPane.CANCEL_OPTION) {
      return false;
    }
    if (i == JOptionPane.YES_OPTION) {
      save();
    }
    else revert();
    return true;
  }

  /**
   * A JTextPane that displays html pages for LibraryTreeNodes.
   */
  protected static class HTMLPane extends JTextPane {
  	
  	public HTMLPane() {
      setEditable(false);
      setFocusable(false);
      setContentType("text/html"); //$NON-NLS-1$
      addHyperlinkListener(hyperlinkListener);  
  	}
  	
    public void paintComponent(Graphics g) {
      if(OSPRuntime.antiAliasText) {
        Graphics2D g2 = (Graphics2D) g;
        RenderingHints rh = g2.getRenderingHints();
        rh.put(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
        rh.put(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
      }
      super.paintComponent(g);
    }

  }

  
  /**
   * A JTextField for editing LibraryTreeNode data.
   */
  protected static class EntryField extends JTextField {
  	
  	EntryField() {
  		getDocument().putProperty("parent", this); //$NON-NLS-1$
      addFocusListener(focusListener);
      addActionListener(actionListener);
      getDocument().addDocumentListener(documentListener);
  	}
  	
  	static DocumentListener documentListener = new DocumentListener() {   
      public void insertUpdate(DocumentEvent e) {
      	EntryField field = (EntryField)e.getDocument().getProperty("parent"); //$NON-NLS-1$
      	field.setBackground(Color.yellow);
      }
      public void removeUpdate(DocumentEvent e) {
      	EntryField field = (EntryField)e.getDocument().getProperty("parent"); //$NON-NLS-1$
      	field.setBackground(Color.yellow);
      }
			public void changedUpdate(DocumentEvent e) {}
  	};
  	
    static FocusListener focusListener = new FocusAdapter() {
      public void focusGained(FocusEvent e) {
      	EntryField field = (EntryField)e.getSource();
      	field.selectAll();
	      field.setBackground(Color.white);
      }
      public void focusLost(FocusEvent e) {
      	EntryField field = (EntryField)e.getSource();
      	ActionListener action = field.getActionListeners()[0];
      	action.actionPerformed(null);
	      field.setBackground(Color.white);
      }
    };

  	static ActionListener actionListener = new ActionListener() {
      public void actionPerformed(ActionEvent e) {
      	if (e==null) return;
      	EntryField field = (EntryField)e.getSource();
        field.setBackground(Color.white);
      }
    };
    
  }

  /**
   * A tree node renderer to render LibraryTreeNodes.
   */
  private class LibraryTreeNodeRenderer extends DefaultTreeCellRenderer {
    public Component getTreeCellRendererComponent(JTree tree, Object value, boolean sel, boolean expanded, boolean leaf, int row, boolean hasFocus) {
      super.getTreeCellRendererComponent(tree, value, sel, expanded, leaf, row, hasFocus);
      LibraryTreeNode node = (LibraryTreeNode) value;
      Icon icon = node.record.getIcon();
      Color c = getForeground();
      if (node.record instanceof LibraryCollection) {
      	setToolTipText(node.isRoot()? pathToCollection: null);
      	icon = expanded? getOpenIcon(): getClosedIcon();
      	if (node.getTarget()!=null)
      		c = Color.red;
      }
      else {
      	setToolTipText(node.getAbsoluteTarget());
      }
      setIcon(icon!=null? icon: LibraryResource.unknownIcon);
      setForeground(c);
      return this;
    }
  }
  
//______________________________  static methods ___________________________
  
  /**
   * Gets a shared file chooser.
   * 
   * @return the file chooser
   */
  protected static JFileChooser getFileChooser() {
  	if (chooser==null) {
  		String chooserDir = LibraryBrowser.getBrowser().library.chooserDir;
      chooser = (chooserDir==null)? 
      		new JFileChooser() : new JFileChooser(new File(chooserDir));
      htmlFilter = new FileFilter() {
        // accept directories and html files
        public boolean accept(File f) {
          if(f==null) return false;
          if(f.isDirectory()) return true;
          String ext = XML.getExtension(f.getName());
          String[] accept = new String[] {"html", "htm"}; //$NON-NLS-1$ //$NON-NLS-2$
          for (String next: accept) {
            if (next.equals(ext)) return true;
          }
          return false;
        }
        public String getDescription() {
          return ToolsRes.getString("LibraryTreePanel.HTMLFileFilter.Description"); //$NON-NLS-1$
        } 	     	
      };
      folderFilter = new FileFilter() {
        // accept directories only
        public boolean accept(File f) {
          if(f!=null && f.isDirectory()) return true;
          return false;
        }
        public String getDescription() {
          return ToolsRes.getString("LibraryTreePanel.FolderFileFilter.Description"); //$NON-NLS-1$
        } 	     	
      };
  	}
  	return chooser;
  }

}

/*
 * Open Source Physics software is free software; you can redistribute
 * it and/or modify it under the terms of the GNU General Public License (GPL) as
 * published by the Free Software Foundation; either version 2 of the License,
 * or(at your option) any later version.

 * Code that uses any portion of the code in the org.opensourcephysics package
 * or any subpackage (subdirectory) of this package must must also be be released
 * under the GNU GPL license.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston MA 02111-1307 USA
 * or view the license online at http://www.gnu.org/copyleft/gpl.html
 *
 * Copyright (c) 2007  The Open Source Physics project
 *                     http://www.opensourcephysics.org
 */
