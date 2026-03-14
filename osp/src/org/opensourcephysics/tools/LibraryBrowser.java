/*
 * Open Source Physics software is free software as described near the bottom of this code file.
 *
 * For additional information and documentation on Open Source Physics please see:
 * <http://www.opensourcephysics.org/>
 */

package org.opensourcephysics.tools;

import java.awt.*;
import java.awt.event.*;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.net.URL;

import javax.swing.*;
import javax.swing.border.Border;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;

import org.opensourcephysics.controls.OSPLog;
import org.opensourcephysics.controls.XML;
import org.opensourcephysics.controls.XMLControl;
import org.opensourcephysics.controls.XMLControlElement;
import org.opensourcephysics.desktop.OSPDesktop;
import org.opensourcephysics.display.DisplayRes;
import org.opensourcephysics.display.GUIUtils;
import org.opensourcephysics.display.OSPRuntime;
import org.opensourcephysics.display.TextFrame;
import org.opensourcephysics.tools.LibraryCollection;
import org.opensourcephysics.tools.LibraryTreePanel;

/**
 * A GUI for browsing OSP digital library collections.
 *
 * @author Douglas Brown
 */
public class LibraryBrowser extends JPanel {
	
	// static fields
	private static LibraryBrowser browser;
	protected static Border buttonBorder;
	protected static String libraryFileName = "my_library.xml"; //$NON-NLS-1$
	protected static String collectionFileName = "my_collection.xml"; //$NON-NLS-1$
  protected static String helpName = "library_browser_help.html";                                 //$NON-NLS-1$
//  protected static String helpBase = "http://www.opensourcephysics.org/online_help/tools/"; //$NON-NLS-1$
  protected static String helpBase = "http://www.cabrillo.edu/~dbrown/OSP/html/"; //$NON-NLS-1$
  protected static String htmlAboutPath = "org/opensourcephysics/resources/tools/html/library_browser.html"; //$NON-NLS-1$
  protected static boolean webConnected;
  protected static JFrame frame;
  protected static JDialog externalDialog;
  protected static JMenuBar menubar;
  protected static Icon expandIcon, contractIcon, heavyExpandIcon, heavyContractIcon;
  
	static {
    buttonBorder = BorderFactory.createEtchedBorder();
    Border space = BorderFactory.createEmptyBorder(1,2,2,2);
    buttonBorder = BorderFactory.createCompoundBorder(buttonBorder, space);
    space = BorderFactory.createEmptyBorder(0,1,0,1);
    buttonBorder = BorderFactory.createCompoundBorder(space, buttonBorder);
    menubar = new JMenuBar();
    String imageFile = "/org/opensourcephysics/resources/tools/images/expand.png";        //$NON-NLS-1$
    expandIcon = new ImageIcon(LibraryTreePanel.class.getResource(imageFile));
    imageFile = "/org/opensourcephysics/resources/tools/images/contract.png";        //$NON-NLS-1$
    contractIcon = new ImageIcon(LibraryTreePanel.class.getResource(imageFile));
    imageFile = "/org/opensourcephysics/resources/tools/images/expand_bold.png";        //$NON-NLS-1$
    heavyExpandIcon = new ImageIcon(LibraryTreePanel.class.getResource(imageFile));
    imageFile = "/org/opensourcephysics/resources/tools/images/contract_bold.png";        //$NON-NLS-1$
    heavyContractIcon = new ImageIcon(LibraryTreePanel.class.getResource(imageFile));
	}
	
	// instance fields
  protected Library library = new Library();
  protected String libraryPath;
  protected JToolBar toolbar;
  protected Action commandAction;
  protected JLabel commandLabel;
  protected JTextField commandField;
  protected JMenu fileMenu, collectionsMenu, helpMenu;
  protected JMenuItem newItem, openItem, saveItem, saveAsItem, exportItem,
  		exitItem, deleteItem, libraryManagerItem, aboutItem, logItem, helpItem;
  protected JButton commandButton, editButton;
  protected ActionListener loadCollectionAction;
  protected boolean exitOnClose;
  protected JTabbedPane tabbedPane;
  protected JScrollPane htmlScroller;
  protected PropertyChangeListener treePanelListener;
  protected boolean keyPressed, textChanged;
  protected TextFrame helpFrame;
  protected JTextPane htmlAboutPane;

	/**
	 * Gets the shared singleton browser.
	 * 
	 * @return the shared LibraryBrowser
	 */
  public static LibraryBrowser getBrowser() {
  	if (browser==null) {
  		browser = getBrowser(null);
  	}
  	return browser;
  }

  
  /**
	 * Gets the shared singleton browser in a JDialog or, if none, in a shared JFrame.
	 * 
	 * @param dialog a JDialog (if null, browser is returned in a JFrame)
	 * @return the shared LibraryBrowser
	 */
  public static LibraryBrowser getBrowser(JDialog dialog) {
  	boolean newFrame = false;
  	if (frame==null && dialog==null) {
  		newFrame = true;
  		frame = new JFrame();
  	}
  	LibraryBrowser.externalDialog = dialog;
 	
  	if (browser==null) {
  		String userHome = System.getProperty("user.home").replace('\\', '/'); //$NON-NLS-1$
  		String ospFolder = OSPRuntime.isWindows()? "/My Documents/OSP/": "/Documents/OSP/"; //$NON-NLS-1$ //$NON-NLS-2$
    	String libraryPath = userHome+ospFolder+libraryFileName;
      File libraryFile = new File(libraryPath);
    	// create new library if none exists
      boolean libraryExists = libraryFile.exists();
      if (!libraryExists) {
      	String collectionPath = userHome+ospFolder+collectionFileName;      	
  			File collectionFile = new File(collectionPath);
      	// create new collection if none exists
        if (!collectionFile.exists()) {
          String name = ToolsRes.getString("LibraryCollection.Name.Local"); //$NON-NLS-1$
    			LibraryCollection collection = new LibraryCollection(name);
    			String base = XML.getDirectoryPath(collectionPath);
    			collection.setBasePath(XML.forwardSlash(base));
    			// save new collection
    			XMLControl control = new XMLControlElement(collection);
    			control.write(collectionPath);
        }
        Library library = new Library();
        String name = ToolsRes.getString("LibraryCollection.Name.Local"); //$NON-NLS-1$
        library.addCollection(collectionPath, name);
        library.save(libraryPath);
      }
  		browser = new LibraryBrowser(libraryPath);

      LibraryTreePanel treePanel = browser.getSelectedTreePanel();
  		if (treePanel!=null) {
  			treePanel.setSelectedNode(treePanel.rootNode);
  			treePanel.showInfo(treePanel.rootNode);
  		}
  		OSPLog.getOSPLog(); // instantiate log in case of exceptions, etc 
  	}
  	
  	
  	browser.setTitle(ToolsRes.getString("LibraryBrowser.Title")); //$NON-NLS-1$
    if (externalDialog!=null) {
    	externalDialog.setContentPane(browser);
    	externalDialog.setJMenuBar(menubar);
    	externalDialog.addWindowListener(new WindowAdapter() {
	      public void windowClosing(WindowEvent e) {
	        browser.exit();
	      }
	    });
    	externalDialog.pack();
    }
    else {
	  	frame.setContentPane(browser);
    	frame.setJMenuBar(menubar);
	    frame.setDefaultCloseOperation(WindowConstants.DO_NOTHING_ON_CLOSE);
	    // add window listener to exit
	    frame.addWindowListener(new WindowAdapter() {
	      public void windowClosing(WindowEvent e) {
	        browser.exit();
	      }
	    });
	    try {
	      java.net.URL url = LibraryBrowser.class.getResource(OSPRuntime.OSP_ICON_FILE);
	      ImageIcon icon = new ImageIcon(url);
	      frame.setIconImage(icon.getImage());
	    } catch(Exception ex) {} 
    	frame.pack();
    	if (newFrame) {
        // center on screen
        Dimension dim = Toolkit.getDefaultToolkit().getScreenSize();
        int x = (dim.width - frame.getBounds().width) / 2;
        int y = (dim.height - frame.getBounds().height) / 2;
        frame.setLocation(x, y);    		
    	}
    }
  	
  	return browser;
  }
  
  /**
   * Returns true if connected to the web.
   * 
   * @return true if web connected
   */
	public static boolean isWebConnected() {
		return webConnected;
	}
	
  /**
   * Imports a library with a specified path.
   * 
   * @param path the path to the Library xml file
   */
	public void importLibrary(final String path) {		
		Runnable runner = new Runnable() {
			public void run() {
		  	library.importLibrary(path);		
        rebuildCollectionsMenu();
			}
		};
		new Thread(runner).start();
	}
	
  /**
   * Adds an OSP-sponsored library with a specified path.
   * 
   * @param path the path to the Library xml file
   */
	public void addOSPLibrary(final String path) {		
		Runnable runner = new Runnable() {
			public void run() {
		  	library.addOSPLibrary(path);		
        rebuildCollectionsMenu();
			}
		};
		new Thread(runner).start();
	}
	
  /**
   * Adds a ComPADRE collection with a specified path.
   * 
   * @param path the ComPADRE query
   */
	public void addComPADRECollection(String path) {
		library.addComPADRECollection(path, LibraryComPADRE.getCollectionName(path));
	}
	
  /**
   * Rebuilds the collection menu.
   */
	public void rebuildCollectionsMenu() {
		library.rebuildCollectionsMenu(this);
	}
	
//____________________ private and protected methods ____________________________

  /**
   * Private constructor to prevent instantiation except for singleton.
   * 
   * @param libraryPath the path to a Library xml file
   */
  private LibraryBrowser(String libraryPath) {
  	super(new BorderLayout());
  	this.libraryPath = libraryPath;
  	library.browser = this;
    createGUI();
    refreshGUI();
		rebuildCollectionsMenu();
    editButton.requestFocusInWindow();
    ToolsRes.addPropertyChangeListener("locale", new PropertyChangeListener() { //$NON-NLS-1$
      public void propertyChange(PropertyChangeEvent e) {
        refreshGUI();
        rebuildCollectionsMenu();
        if (library.manager!=null)
        	library.manager.refreshGUI();
        LibraryTreePanel.htmlPanesByNode.clear();
        LibraryTreePanel treePanel = getSelectedTreePanel();
        if (treePanel!=null)
        	treePanel.showInfo(treePanel.getSelectedNode());
      }
    });
  }
  
  /**
   * Gets the selected LibraryTreePanel, if any.
   * 
   * @return the selected treePanel, or null if none
   */
  protected LibraryTreePanel getSelectedTreePanel() {
  	return (LibraryTreePanel)tabbedPane.getSelectedComponent();
  }
  
  /**
   * Gets the LibraryTreePanel at a specified tab index.
   * 
   * @param index the tab index
   * @return the treePanel
   */
  protected LibraryTreePanel getTreePanel(int index) {  	
  	return (LibraryTreePanel)tabbedPane.getComponentAt(index);
  }

  /**
   * Gets the title of the tab associated with a given path.
   * 
   * @param path the collection path
   * @return the tab title
   */
  protected String getTabTitle(String path) {
  	int i = getTabIndexFromPath(path);
  	return i>-1? getTabTitle(i): null;
  }

  /**
   * Gets the title of the tab at a given index.
   * 
   * @param index the tab index
   * @return the tab title
   */
  protected String getTabTitle(int index) { 
  	String title = tabbedPane.getTitleAt(index);
		if (title.endsWith("*")) //$NON-NLS-1$
			title = title.substring(0, title.length()-1);
		return title;
  }

  /**
   * Gets the index of the tab associated with a given path.
   * 
   * @param path the collection path
   * @return the tab index
   */
  protected int getTabIndexFromPath(String path) { 
  	for (int i=0; i<tabbedPane.getTabCount(); i++) {
  		LibraryTreePanel next = getTreePanel(i);
  		if (next.pathToCollection.equals(path)) 
  			return i;
  	}
  	return -1;
  }

  /**
   * Gets the index of the tab associated with a given path.
   * 
   * @param path the collection path
   * @return the tab index
   */
  protected int getTabIndexFromTitle(String title) { 
  	for (int i=0; i<tabbedPane.getTabCount(); i++) {
  		String next = tabbedPane.getTitleAt(i);
  		if (next.equals(title)) 
  			return i;
  	}
  	return -1;
  }

  /**
   * Loads the collection with a given path.
   * 
   * @param path the collection path
   */
  protected void loadCollection(String path) {
    // select tab if collection is already loaded
    int i = getTabIndexFromPath(path);
    if (i>-1) {
    	tabbedPane.setSelectedIndex(i);
    	return;
    }
    // otherwise add new tab
  	int n = addTab(path, false);
  	if (n>-1)
  		tabbedPane.setSelectedIndex(n);
  }
  
  /**
   * Loads the collection with a given path.
   * 
   * @param path the collection path
   * @return the loaded collection
   */
  protected LibraryCollection loadXML(String path) {
    LibraryCollection collection = null;  	
    // look first for comPADRE collections
    if (LibraryComPADRE.isComPADREPath(path)) {
//    	boolean primarySubjectOnly = LibraryComPADRE.isPrimarySubjectOnly(path);
//    	path = LibraryComPADRE.getCollectionPath(path, primarySubjectOnly);
    	collection = LibraryComPADRE.getCollection(path);
    }
    else {
    	XMLControlElement control = new XMLControlElement();
    	control.read(path);
    	if (!control.failedToRead() && LibraryCollection.class==control.getObjectClass()) {
      	collection = (LibraryCollection)control.loadObject(null);
    	}    	
    }  	
  	return collection;
  }
  
  /**
   * Adds a tab displaying a collection with a given path.
   * 
   * @param path the collection path
   * @param addToMyLibrary true to automatically add the collection to my library
   * @return the tab index added
   */
  protected int addTab(String path, boolean addToMyLibrary) {
  	if (path==null) return -1;
  	if (!webConnected && path.startsWith("http:")) { //$NON-NLS-1$
  		JOptionPane.showMessageDialog(this, 
  				ToolsRes.getString("LibraryBrowser.Dialog.ServerUnavailable.Message"), //$NON-NLS-1$
  				ToolsRes.getString("LibraryBrowser.Dialog.ServerUnavailable.Title"), //$NON-NLS-1$
  				JOptionPane.WARNING_MESSAGE);  		
  		return -1;
  	}
  	
  	setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
  	LibraryCollection collection = loadXML(path);  	
  	
  	if (collection!=null) {
  		LibraryTreePanel treePanel = createLibraryTreePanel();
  		boolean isLocal = !path.startsWith("http:"); //$NON-NLS-1$
  		treePanel.setCollection(collection, path, isLocal);
  		tabbedPane.addTab("", treePanel); //$NON-NLS-1$
  		refreshTabTitle(path, collection);
  		int i = tabbedPane.getTabCount()-1;
  		tabbedPane.setToolTipTextAt(i, path);
  		if (addToMyLibrary) {
  			addToCollections(treePanel.pathToCollection);
  		}
  		refreshGUI();
    	setCursor(Cursor.getDefaultCursor());
  		return i;
    }
  	setCursor(Cursor.getDefaultCursor());
  	String s = ToolsRes.getString("LibraryBrowser.Dialog.CollectionNotFound.Message"); //$NON-NLS-1$
  	JOptionPane.showMessageDialog(this, 
  			s+":\n"+path, //$NON-NLS-1$
				ToolsRes.getString("LibraryBrowser.Dialog.CollectionNotFound.Title"), //$NON-NLS-1$
				JOptionPane.WARNING_MESSAGE);  		
  	return -1;
  }
  
  /**
   * Refreshes the title of a tab based on the properties of a LibraryCollection
   * and the path associated with that collection.
   * 
   * @param path the collection path
   * @param collection the LibraryCollection itself
   */
  protected void refreshTabTitle(String path, LibraryCollection collection) {
  	int n = getTabIndexFromPath(path);
  	if (n==-1) return;
  	// title is name of collection or, if unnamed, file or server name
  	String title = collection.getName();
    if (title.equals("")) { //$NON-NLS-1$
    	String basePath = XML.getDirectoryPath(path);
    	if (basePath.startsWith("http:")) { //$NON-NLS-1$
    		basePath = basePath.substring(5);
    		while (basePath.startsWith("/")) { //$NON-NLS-1$
    			basePath = basePath.substring(1);
    		}

    		int i = basePath.indexOf("/"); //$NON-NLS-1$
    		if (i>-1)
    			basePath = basePath.substring(0, i);
    		title = basePath;
    	}
    	else { 
    		title = XML.getName(path);    	
    	}
    }
    Component c = tabbedPane.getTabComponentAt(n);
    if (path.contains(LibraryComPADRE.TRACKER_SERVER_TREE) && c==null) {
    	final TabTitle tabTitle = new TabTitle(contractIcon, heavyContractIcon);
	  	tabTitle.iconLabel.setToolTipText(ToolsRes.getString("LibraryBrowser.Tooltip.Contract")); //$NON-NLS-1$
    	Action action = new AbstractAction() {
  		  public void actionPerformed(ActionEvent e) {
  		  	boolean primaryOnly = tabTitle.normalIcon==contractIcon;
  		  	int index = getTabIndexFromTitle(tabTitle.titleLabel.getText());
  		  	if (index>-1) {
	  		  	tabTitle.setIcons(primaryOnly? expandIcon: contractIcon, primaryOnly? heavyExpandIcon: heavyContractIcon);
	  		  	tabTitle.iconLabel.setToolTipText(primaryOnly? ToolsRes.getString("LibraryBrowser.Tooltip.Expand"): //$NON-NLS-1$
	  		  		ToolsRes.getString("LibraryBrowser.Tooltip.Contract")); //$NON-NLS-1$
	  	  		LibraryTreePanel treePanel = getTreePanel(index);
	  		  	String query = LibraryComPADRE.getCollectionPath(treePanel.pathToCollection, primaryOnly);
	  		  	LibraryCollection collection = LibraryComPADRE.getCollection(query);
	  		  	treePanel.setCollection(collection, query, false);
  		  	}
  		  }
    	};
    	tabTitle.setAction(action);
    	tabbedPane.setTabComponentAt(n, tabTitle);
    }
    tabbedPane.setTitleAt(n, title);
		collection.displayName = title;
  }

  /**
   * Creates the visible components of this panel.
   */
  protected void createGUI() {
    setPreferredSize(new Dimension(800, 450));

    loadCollectionAction = new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        loadCollection(e.getActionCommand());
      }
    };
    
    // create command action, label, field and button
		commandAction = new AbstractAction() {
		  public void actionPerformed(ActionEvent e) {
		  	if (e==null) return;
      	commandField.setBackground(Color.white);
      	commandField.setForeground(LibraryTreePanel.defaultForeground);
		  	if (!commandButton.isEnabled()) return;
        String path = commandField.getText().trim();
        if (path.equals("")) return; //$NON-NLS-1$
				path = XML.forwardSlash(path);
        path = ResourceLoader.getNonURIPath(path);
        Resource res = null;
    		String xmlPath = path;
        
        // if path has no extension, look for xml file with same name
        if (!path.startsWith("http://www.compadre.org/OSP/") //$NON-NLS-1$
        		&& XML.getExtension(path)==null) {
      		while (xmlPath.endsWith("/")) //$NON-NLS-1$
      			xmlPath = xmlPath.substring(0, xmlPath.length()-1);
      		if (!xmlPath.equals("")) { //$NON-NLS-1$
      			String name = XML.getName(xmlPath);
      			xmlPath += "/"+name+".xml"; //$NON-NLS-1$ //$NON-NLS-2$
            res = ResourceLoader.getResource(xmlPath);
      		}
        }
        
        if (res!=null)
        	path = xmlPath;
        else 
        	res = ResourceLoader.getResourceZipURLsOK(path);
        
        if (res==null) {
        	commandField.setForeground(LibraryTreePanel.darkRed);
//        	commandField.setBackground(LibraryTreePanel.lightRed);
        	return;
        }
        
  			XMLControl control = new XMLControlElement(path);
  			if (!control.failedToRead()) { 
          // see if path points to a collection
  				if (control.getObjectClass()==LibraryCollection.class) {
	      		loadCollection(path);
	      		refreshGUI();
	      		LibraryTreePanel treePanel = getSelectedTreePanel();
	      		if (treePanel!=null && treePanel.pathToCollection.equals(path)) {
		      		treePanel.setSelectedNode(treePanel.rootNode);
		      		commandField.setBackground(Color.white);
		      		commandField.repaint();
	      		}
	      		return;
	  			}
//    			// see if path points to a library
//  				if (Library.class==control.getObjectClass()) {
//  					Library newLibrary = (Library)control.loadObject(null);
//  					if (library.importLibrary(path, newLibrary)) {
//  						library.rebuildCollectionsMenu(browser);
//  		        JOptionPane.showMessageDialog(LibraryBrowser.this, 
//  		        		ToolsRes.getString("LibraryBrowser.Dialog.LibraryAdded.Message") //$NON-NLS-1$
//  		        		+":\n\""+newLibrary.getName()+"\""); //$NON-NLS-1$ //$NON-NLS-2$
//  					}
//  				}
   			}
       
  			// send command
		  	LibraryResource record = null;
        LibraryTreePanel treePanel = getSelectedTreePanel();
        if (treePanel!=null && treePanel.getSelectedNode()!=null) {
        	record = treePanel.getSelectedNode().record.getClone();
        	record.setBasePath(treePanel.getSelectedNode().getBasePath());
        }
	    	else {
	    		record = new LibraryResource(""); //$NON-NLS-1$
	    		record.setTarget(path);
	    	}
	    	LibraryBrowser.this.firePropertyChange("target", null, record); //$NON-NLS-1$
		  }		
		};
    commandLabel = new JLabel();
    commandLabel.setAlignmentX(CENTER_ALIGNMENT);
    commandLabel.setBorder(BorderFactory.createEmptyBorder(0, 3, 0, 2));
    commandField = new JTextField();
    LibraryTreePanel.defaultForeground = commandField.getForeground();
    commandField.addActionListener(commandAction);
    commandField.getDocument().addDocumentListener(new DocumentListener() {   
      public void insertUpdate(DocumentEvent e) {
      	String text = commandField.getText();
        commandButton.setEnabled(!"".equals(text)); //$NON-NLS-1$
        textChanged = keyPressed;
        LibraryTreePanel treePanel = getSelectedTreePanel();
        if (treePanel!=null) {
        	treePanel.command = text;
        	LibraryTreeNode node = treePanel.getSelectedNode();
        	if (node!=null && node.isRoot() && treePanel.pathToCollection.equals(text))
        		commandButton.setEnabled(false);
        }
        else {
        	commandField.setBackground(Color.yellow);
        	commandField.setForeground(LibraryTreePanel.defaultForeground);
        }
      }
      public void removeUpdate(DocumentEvent e) {
        commandButton.setEnabled(!"".equals(commandField.getText())); //$NON-NLS-1$
        textChanged = keyPressed;
        LibraryTreePanel treePanel = getSelectedTreePanel();
        if (treePanel!=null) {
        	treePanel.command = commandField.getText();
        }
        else {
        	commandField.setBackground(Color.yellow);
        	commandField.setForeground(LibraryTreePanel.defaultForeground);
        }
      }
			public void changedUpdate(DocumentEvent e) {}
  	});
    commandField.addKeyListener(new KeyAdapter() {
      public void keyPressed(KeyEvent e) {
        keyPressed = true;
      }
      public void keyReleased(KeyEvent e) {
        LibraryTreePanel treePanel = getSelectedTreePanel();
        if (treePanel!=null && textChanged && e.getKeyCode()!=KeyEvent.VK_ENTER) {
        	commandField.setBackground(Color.yellow);
        	commandField.setForeground(LibraryTreePanel.defaultForeground);
          treePanel.setSelectedNode(null);
        }
        textChanged = keyPressed = false;
      }
    });
    commandField.addFocusListener(new FocusAdapter() {
      public void focusGained(FocusEvent e) {
      	commandField.selectAll();
      }
    });

    commandButton = new JButton(commandAction);
    commandButton.setOpaque(false);
    commandButton.setBorder(buttonBorder);
    

    
    tabbedPane = new JTabbedPane(SwingConstants.TOP) {
    	@Override
    	public void setTitleAt(int i, String title) {
  			super.setTitleAt(i, title);
    		Component c = tabbedPane.getTabComponentAt(i);
    		if (c!=null) {
    			TabTitle tabTitle = (TabTitle)c;
    			tabTitle.setTitle(title);
    		}
    	}
    };
    tabbedPane.addChangeListener(new ChangeListener() {
      public void stateChanged(ChangeEvent e) {
      	refreshGUI();
      	LibraryTreePanel treePanel = getSelectedTreePanel();
      	if (treePanel!=null) {
	        LibraryTreeNode node = treePanel.getSelectedNode();
	        if (node!=null) {
	        	String path = node.isRoot()? treePanel.pathToCollection: node.getAbsoluteTarget();
	        	commandField.setText(path);
//		        treePanel.showInfo(node);
	        }
	        else {
    	  		commandField.setText(treePanel.command);
    	  		commandField.setCaretPosition(0);
	        }
      	}
        commandField.setBackground(Color.white);
      	commandField.setForeground(LibraryTreePanel.defaultForeground);
      	if (library.manager!=null && library.manager.isVisible())
      		library.manager.refreshGUI();
      }
    });
    tabbedPane.addMouseListener(new MouseAdapter() {
      public void mousePressed(MouseEvent e) {
        if(OSPRuntime.isPopupTrigger(e)) {
          // make popup and add items
          JPopupMenu popup = new JPopupMenu();
          // close this tab
          JMenuItem item = new JMenuItem(ToolsRes.getString("MenuItem.Close")); //$NON-NLS-1$
          item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
              int i = tabbedPane.getSelectedIndex();
              LibraryTreePanel treePanel = getTreePanel(i);
              if (!treePanel.saveChanges(getTabTitle(i))) return;
              tabbedPane.removeTabAt(i);
              refreshGUI();
            }
          });
          popup.add(item);
          // add tab to Collections menu
          if (!library.containsPath(getSelectedTreePanel().pathToCollection, false)) {
	          item = new JMenuItem(ToolsRes.getString("LibraryBrowser.MenuItem.AddToLibrary")); //$NON-NLS-1$
	          item.addActionListener(new ActionListener() {
	            public void actionPerformed(ActionEvent e) {
	              LibraryTreePanel treePanel = getSelectedTreePanel();
	              addToCollections(treePanel.pathToCollection);
	           }
	          });
	          popup.addSeparator();
	          popup.add(item);
          }
          popup.show(tabbedPane, e.getX(), e.getY()+8);
        }
      }
    });
    
    // create property change listener for treePanels
    treePanelListener = new PropertyChangeListener() {
  		public void propertyChange(PropertyChangeEvent e) {
  			String propertyName = e.getPropertyName();
  			if (propertyName.equals("collection_edit")) { //$NON-NLS-1$
    			refreshGUI();
  			}
  			else if (propertyName.equals("target")) { //$NON-NLS-1$
    			LibraryResource record = null;
    			if (e.getNewValue() instanceof LibraryTreeNode) {
  	  			LibraryTreeNode node = (LibraryTreeNode)e.getNewValue();
  	  			if (node.record instanceof LibraryCollection) {
  		  			setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
  	  				if (!LibraryComPADRE.loadResources(node)) {
    		  			setCursor(Cursor.getDefaultCursor());
  	  		  		JOptionPane.showMessageDialog(LibraryBrowser.this, 
  	  		  				ToolsRes.getString("LibraryBrowser.Dialog.NoResources.Message"), //$NON-NLS-1$
  	  		  				ToolsRes.getString("LibraryBrowser.Dialog.NoResources.Title"), //$NON-NLS-1$
  	  		  				JOptionPane.PLAIN_MESSAGE); 
  	  		  		return;
  	  				}
  	  				if (node.createChildNodes()) {
  	  					node.setTarget(null);
  	  				}
  	  				node.record.setType(LibraryResource.COLLECTION_TYPE);
  	  				node.record.setDescription(null);
  	  				LibraryTreePanel.htmlPanesByNode.remove(node);
  	  				getSelectedTreePanel().setSelectedNode(node);
  		  			setCursor(Cursor.getDefaultCursor());
  	  				return;
  	  			}
  	  			record = node.record.getClone();
  	  			record.setBasePath(node.getBasePath());
    			}
    			else record = (LibraryResource)e.getNewValue();
    			
    			String target = record.getTarget();
    			if (target!=null && (target.toLowerCase().endsWith(".pdf") //$NON-NLS-1$
    					 || target.toLowerCase().endsWith(".html") //$NON-NLS-1$
    					 || target.toLowerCase().endsWith(".htm"))) { //$NON-NLS-1$
    				target = XML.getResolvedPath(target, record.getBasePath());
    				target = ResourceLoader.getURIPath(target);
		  			setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
    				OSPDesktop.displayURL(target);
		  			setCursor(Cursor.getDefaultCursor());
    			}
    			else {
	    			// forward the event to browser listeners
	    			firePropertyChange("target", e.getOldValue(), record); //$NON-NLS-1$
    			}
  			}
  		}
  	};
    
    // create edit button
    editButton = new JButton();
    editButton.setOpaque(false);
	  editButton.setBorder(buttonBorder);
    editButton.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
      	final LibraryTreePanel treePanel = getSelectedTreePanel();
      	if (!treePanel.isEditing()) {
        	treePanel.setEditing(true);
        	refreshGUI();
      	}
      	else if (!treePanel.isCollectionChanged()) {
        	treePanel.setEditing(false);
        	refreshGUI();
      	}
      	else {
        	JPopupMenu popup = new JPopupMenu(); 
          JMenuItem item = new JMenuItem(ToolsRes.getString("LibraryBrowser.MenuItem.SaveEdits")); //$NON-NLS-1$
          popup.add(item);
          item.addActionListener(new ActionListener() {
      		  public void actionPerformed(ActionEvent e) {
            	save();
            	treePanel.setEditing(false);
            	refreshGUI();
            }
          });
          item = new JMenuItem(ToolsRes.getString("LibraryBrowser.MenuItem.Discard")); //$NON-NLS-1$
          popup.add(item);
          item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
            	treePanel.setEditing(false);
            	treePanel.revert();
            	refreshGUI();
      		  }		
      		});
          popup.show(editButton, 0, editButton.getHeight());      		
      	}
      }
    });
    
    toolbar = new JToolBar();
    toolbar.setFloatable(false);
    Border empty = BorderFactory.createEmptyBorder(1, 2, 1, 2);
    Border etched = BorderFactory.createEtchedBorder();
    toolbar.setBorder(BorderFactory.createCompoundBorder(etched, empty));
    toolbar.add(editButton);
    toolbar.addSeparator();
    toolbar.add(commandLabel);
    toolbar.add(commandField);
    toolbar.add(commandButton);

    add(toolbar, BorderLayout.NORTH);
    
    // menu items
    fileMenu = new JMenu();
    menubar.add(fileMenu);
    newItem = new JMenuItem();
    int mask = Toolkit.getDefaultToolkit().getMenuShortcutKeyMask();
    newItem.setAccelerator(KeyStroke.getKeyStroke('N', mask));
    newItem.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
      	createNewCollection();
      }
    });
    openItem = new JMenuItem();
    openItem.setAccelerator(KeyStroke.getKeyStroke('O', mask));
    openItem.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        open();
      }
    });
    saveItem = new JMenuItem();
    saveItem.setAccelerator(KeyStroke.getKeyStroke('S', mask));
    saveItem.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
      	save();
      }
    });
    saveAsItem = new JMenuItem();
    saveAsItem.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
      	saveAs();
      }
    });
    exportItem = new JMenuItem();
    exportItem.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
      	export(library);
      }
    });
    exitItem = new JMenuItem();
    exitItem.setAccelerator(KeyStroke.getKeyStroke('Q', mask));
    exitItem.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
      	exit();
      }
    });
    
    
    collectionsMenu = new JMenu();
    menubar.add(collectionsMenu);
    libraryManagerItem = new JMenuItem();
    libraryManagerItem.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
      	library.getManager(LibraryBrowser.this).setVisible(true);
      }
    });
    helpMenu = new JMenu();
    menubar.add(helpMenu);
    helpItem = new JMenuItem();
    helpItem.setAccelerator(KeyStroke.getKeyStroke('H', mask));
    helpItem.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        showHelp();
      }
    });
    helpMenu.add(helpItem);
    helpMenu.addSeparator();
    logItem = new JMenuItem();
    logItem.setAccelerator(KeyStroke.getKeyStroke('L', mask));
    logItem.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        Point p0 = new Frame().getLocation();
        JFrame frame = OSPLog.getOSPLog();
        if((frame.getLocation().x==p0.x)&&(frame.getLocation().y==p0.y)) {
          Point p = getLocation();
          frame.setLocation(p.x+28, p.y+28);
        }
        frame.setVisible(true);
      }
    });
    helpMenu.add(logItem);
    helpMenu.addSeparator();
    aboutItem = new JMenuItem();
    aboutItem.setAccelerator(KeyStroke.getKeyStroke('A', mask));
    aboutItem.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        showAboutDialog();
      }
    });
    helpMenu.add(aboutItem);
    
    // create html about-browser pane
    htmlAboutPane = new JTextPane() {
      public void paintComponent(Graphics g) {
        if(OSPRuntime.antiAliasText) {
          Graphics2D g2 = (Graphics2D) g;
          RenderingHints rh = g2.getRenderingHints();
          rh.put(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON);
          rh.put(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        }
        super.paintComponent(g);
      }    	
    };
    htmlAboutPane.setEditable(false);
    htmlAboutPane.setFocusable(false);
    htmlAboutPane.setContentType("text/html"); //$NON-NLS-1$
    htmlAboutPane.addHyperlinkListener(LibraryTreePanel.hyperlinkListener);
    htmlAboutPane.setText("<h2>"+ToolsRes.getString("LibraryBrowser.Title")+"</h2>"); //$NON-NLS-1$ //$NON-NLS-2$ //$NON-NLS-3$
		htmlScroller = new JScrollPane(htmlAboutPane);
		Resource res = ResourceLoader.getResource(htmlAboutPath);
		if (res!=null) {
			try {
				URL url = res.getURL();
				htmlAboutPane.setPage(url);
			} catch (Exception ex) {}
		}

		// check for internet connectivity
  	webConnected = ResourceLoader.isURLAvailable("http://www.opensourcephysics.org"); //$NON-NLS-1$
  	if (!webConnected) {
  		JOptionPane.showMessageDialog(this, 
  				ToolsRes.getString("LibraryBrowser.Dialog.ServerUnavailable.Message"), //$NON-NLS-1$
  				ToolsRes.getString("LibraryBrowser.Dialog.ServerUnavailable.Title"), //$NON-NLS-1$
  				JOptionPane.WARNING_MESSAGE);  		
  	}
  	  	
  	library.load(libraryPath);
  	  	
  	final Runnable runner = new Runnable() {
  		public void run() {
  	    // load libraries and add previously open tabs
  	  	if (library.openTabPaths!=null) {  	  		
  	  		for (String path: library.openTabPaths) {
  	  			addTab(path, false);
  	  		}
  	  		LibraryTreePanel treePanel = getSelectedTreePanel();
  	  		if (treePanel!=null) {
  		  		commandField.setText(treePanel.pathToCollection);
  		  		commandField.setCaretPosition(0);
  	      	commandField.setBackground(Color.white);
  	      	commandField.setForeground(LibraryTreePanel.defaultForeground);
  	  		}
  	  	}
  			
  		}
  	};
    if (externalDialog!=null) {
    	externalDialog.addWindowListener(new java.awt.event.WindowAdapter() {
        public void windowOpened(WindowEvent e) {
    			runner.run();
        }
      });
    }
    else {
    	frame.addWindowListener(new java.awt.event.WindowAdapter() {
        public void windowOpened(WindowEvent e) {
    			runner.run();
        }
      });    	
    }
  	
  }
  
  /**
   * Refreshes the GUI, including locale-dependent resources strings.
   */
  protected void refreshGUI() {
  	if (tabbedPane.getTabCount()==0) {
  		remove(tabbedPane);
      add(htmlScroller, BorderLayout.CENTER);
      validate();
  	}
  	else {
  		remove(htmlScroller);
      add(tabbedPane, BorderLayout.CENTER);
  	}
  	// set text strings
  	setTitle(getTitle());
    fileMenu.setText(ToolsRes.getString("Menu.File")); //$NON-NLS-1$
    newItem.setText(ToolsRes.getString("LibraryBrowser.MenuItem.New")); //$NON-NLS-1$
    openItem.setText(ToolsRes.getString("LibraryBrowser.MenuItem.Open")); //$NON-NLS-1$
    saveItem.setText(ToolsRes.getString("LibraryBrowser.MenuItem.Save")); //$NON-NLS-1$
    saveAsItem.setText(ToolsRes.getString("LibraryBrowser.MenuItem.SaveAs")); //$NON-NLS-1$
    exportItem.setText(ToolsRes.getString("LibraryBrowser.MenuItem.Export")); //$NON-NLS-1$
    exitItem.setText(ToolsRes.getString("MenuItem.Exit")); //$NON-NLS-1$
    collectionsMenu.setText(ToolsRes.getString("LibraryBrowser.Menu.Collections")); //$NON-NLS-1$
    libraryManagerItem.setText(ToolsRes.getString("LibraryBrowser.MenuItem.Organize")+"..."); //$NON-NLS-1$ //$NON-NLS-2$
    helpMenu.setText(ToolsRes.getString("Menu.Help")); //$NON-NLS-1$
    helpItem.setText(ToolsRes.getString("LibraryBrowser.MenuItem.Help"));                    //$NON-NLS-1$
    logItem.setText(ToolsRes.getString("MenuItem.Log"));                               //$NON-NLS-1$
    aboutItem.setText(ToolsRes.getString("MenuItem.About"));                           //$NON-NLS-1$
    commandLabel.setText(ToolsRes.getString("LibraryTreePanel.Label.Target")); //$NON-NLS-1$
  	commandButton.setText(ToolsRes.getString("LibraryTreePanel.Button.Load")); //$NON-NLS-1$
  	commandField.setToolTipText(ToolsRes.getString("LibraryBrowser.Field.Command.Tooltip")); //$NON-NLS-1$
  	// rebuild file menu
    fileMenu.removeAll();
    fileMenu.add(newItem);
    fileMenu.add(openItem);
    fileMenu.addSeparator();
    fileMenu.add(saveItem);
    fileMenu.add(saveAsItem);
//    fileMenu.addSeparator();
//    fileMenu.add(exportItem);
    fileMenu.addSeparator();
    fileMenu.add(exitItem);
    
  	LibraryTreePanel treePanel = getSelectedTreePanel();
    if (treePanel!=null) {
      editButton.setText(!treePanel.isEditing()?
      		ToolsRes.getString("LibraryBrowser.Button.OpenEditor"): //$NON-NLS-1$
      		ToolsRes.getString("LibraryBrowser.Button.CloseEditor")); //$NON-NLS-1$
      editButton.setEnabled(treePanel.isEditable());
      saveItem.setEnabled(treePanel.isCollectionChanged());
      int i = tabbedPane.getSelectedIndex();
    	String title = tabbedPane.getTitleAt(i);
      if (treePanel.isCollectionChanged() && !title.endsWith("*")) { //$NON-NLS-1$
      	tabbedPane.setTitleAt(i, title+"*");  //$NON-NLS-1$
      }
      else if (!treePanel.isCollectionChanged() && title.endsWith("*")) { //$NON-NLS-1$
      	tabbedPane.setTitleAt(i, title.substring(0, title.length()-1)); 
      }
      treePanel.refreshGUI();
    }
    else {
    	editButton.setText(ToolsRes.getString("LibraryBrowser.Button.OpenEditor")); //$NON-NLS-1$
      saveItem.setEnabled(false);
      editButton.setEnabled(false);
      commandField.setText(null);
      commandButton.setEnabled(false);
    }
    repaint();
  }
  
  /**
   * Gets the title of this browser.
   * @return the title
   */
  public String getTitle() {
  	return ToolsRes.getString("LibraryBrowser.Title"); //$NON-NLS-1$
  }
  
  /**
   * Sets the title.
   * @param title the title
   */
  public void setTitle(String title) {
     if (frame!=null) {
    	frame.setTitle(title);
    }
    else if (externalDialog!=null) {
    	externalDialog.setTitle(title);
    }
  }
  
  /**
   * Opens a collection file using a file chooser.
   */
  protected void open() {
    File file = GUIUtils.showOpenDialog(LibraryBrowser.this);
		if (file !=null) {
			String path = file.getAbsolutePath();
			XMLControl control = new XMLControlElement(path);
			if (!control.failedToRead()) {
				if (control.getObjectClass()==LibraryCollection.class) {
					path = XML.forwardSlash(path);
	    		loadCollection(path);
	    		refreshGUI();
				}
//				else if (control.getObjectClass()==Library.class) {
//					Library newLibrary = (Library)control.loadObject(null);
//					if (library.importLibrary(path, newLibrary)) {
//						library.rebuildCollectionsMenu(browser);
//		        JOptionPane.showMessageDialog(LibraryBrowser.this, 
//		        		ToolsRes.getString("LibraryBrowser.Dialog.LibraryAdded.Message") //$NON-NLS-1$
//		        		+":\n\""+newLibrary.getName()+"\""); //$NON-NLS-1$ //$NON-NLS-2$
//					}					
//				}
			}
		}
  }
  
  /**
   * Saves the selected LibraryTreePanel collection.
   */
  protected void save() {
  	LibraryTreePanel treePanel = getSelectedTreePanel();
		treePanel.save();
		refreshGUI();
  }

  /**
   * Exports a library by saving a clone of this browser's library
   * after the user is prompted to give the clone a descriptive name.
   * 
   * @param library the library to export
   */
  protected void export(Library library) {
  	Library lib = library.getCloneForExport();
  	String message = ToolsRes.getString("LibraryBrowser.Dialog.NameLibrary.Message1") //$NON-NLS-1$
    		+ "\n"+ToolsRes.getString("LibraryBrowser.Dialog.NameLibrary.Message2");  //$NON-NLS-1$//$NON-NLS-2$
  	Object input = JOptionPane.showInputDialog(this, message,
    		ToolsRes.getString("LibraryBrowser.Dialog.NameLibrary.Title"),   //$NON-NLS-1$
        JOptionPane.QUESTION_MESSAGE, null, null, lib.getName());
    if(input==null || input.equals("")) {                            //$NON-NLS-1$
      return;
    }
    lib.setName(input.toString());
  	String path = getSavePath();
		if (path!=null) {
	  	// save library
	  	lib.save(path);
		}
  }
  
  public void setVisible(boolean vis) {
  	super.setVisible(vis);
  	if (externalDialog!=null) {
  		externalDialog.setVisible(vis);
  	}
  	else frame.setVisible(vis);
  }

  /**
   * Exits this browser.
   * @return true if exited, false if cancelled by user
   */
  public boolean exit() {
  	for (int i=0; i < tabbedPane.getTabCount(); i++) {
  		LibraryTreePanel treePanel = getTreePanel(i);
      if (!treePanel.saveChanges(getTabTitle(i))) return false; // true unless the user cancels      		
  	}
  	// inform library of open tabs
  	int n = tabbedPane.getTabCount();
  	if (n>0) {
  		library.openTabPaths = new String[n];
  		for (int i=0; i<n; i++) {
  			library.openTabPaths[i] = getTreePanel(i).pathToCollection;
  		}
  	}
  	else {
  		library.openTabPaths = null;
  	}
  	// save library
  	library.save(libraryPath);
  	
  	if (exitOnClose) {
      System.exit(0);
    } else {
    	refreshGUI();
      setVisible(false);
    }
  	return true;
  }

  /**
   * Saves the current collection as a new file.
   */
  protected void saveAs() {
  	String path = getSavePath();
		if (path!=null) {
	  	LibraryTreePanel treePanel = getSelectedTreePanel();
			treePanel.setCollection(treePanel.collection, path, true);
			save();
			treePanel.setEditing(true);
		}
  }
  
  /**
   * Returns a path to which a collection or library file can be saved.
   * This adds the extension ".xml", if none, and checks for duplicates.
   * 
   * @return the save path
   */
  protected String getSavePath() {
		File file = GUIUtils.showSaveDialog(this);
		if (file ==null) return null;
		String path = file.getAbsolutePath();
    String extension = XML.getExtension(path);
    if (extension==null) {
    	path = XML.stripExtension(path)+".xml"; //$NON-NLS-1$
    	file = new File(path);
      if(file.exists()) {
        int response = JOptionPane.showConfirmDialog(this, 
        		DisplayRes.getString("DrawingFrame.ReplaceExisting_message") //$NON-NLS-1$
        		+" "+file.getName()+"?", //$NON-NLS-1$ //$NON-NLS-2$
        		DisplayRes.getString("DrawingFrame.ReplaceFile_option_title"), //$NON-NLS-1$
            JOptionPane.YES_NO_CANCEL_OPTION);
        if(response!=JOptionPane.YES_OPTION) {
          return null;
        }
      }
    }
    return path;  	
  }
  
  /**
   * Adds a collection to this browser's library after prompting the user to 
   * assign it a name.
   * 
   * @param path the path to the collection
   */
  protected void addToCollections(String path) {
  	if (library.containsPath(path, true)) {
  		return;
  	}
    String proposed = getTabTitle(path);
    if (proposed==null) {
    	LibraryCollection collection = loadXML(path);
    	if (collection!=null) proposed = collection.getName();
    }
    if (proposed.equals("")) { //$NON-NLS-1$
    	proposed = XML.getName(path); // filename
    }
    
//    String message = ToolsRes.getString("LibraryBrowser.Dialog.NameCollection.Message"); //$NON-NLS-1$
//    Object input = JOptionPane.showInputDialog(this, 
//    		message+"\n"+path, //$NON-NLS-1$
//    		ToolsRes.getString("LibraryBrowser.Dialog.NameCollection.Title"),   //$NON-NLS-1$
//        JOptionPane.QUESTION_MESSAGE, null, null, proposed);
//    if(input==null || input.equals("")) {                            //$NON-NLS-1$
//      return;
//    }
//    proposed = input.toString().trim();
    
//    String newName = getUniqueName(data, proposed, true);
    library.addCollection(path, proposed);
    rebuildCollectionsMenu();
  	refreshGUI();
  }

  /**
   * Creates a new empty LibraryCollection file.
   */
  protected void createNewCollection() {
  	String path = getSavePath();
		if (path!=null) {
			LibraryCollection collection = new LibraryCollection(null);
			String base = XML.getDirectoryPath(path);
			collection.setBasePath(XML.forwardSlash(base));
			// save new collection
			XMLControl control = new XMLControlElement(collection);
			control.write(path);
			path = XML.forwardSlash(path);
  		int n = addTab(path, false);
    	if (n>-1) {
    		tabbedPane.setSelectedIndex(n);
    		LibraryTreePanel treePanel = getSelectedTreePanel();
    		treePanel.setEditing(true);
    	}
		}
  }
  
  /**
   * Returns a name that is not a duplicate of an existing name.
   * 
   * @param proposed a proposed name
   * @param nameToIgnore a name that is ignored when comparing
   * @return a unique name that is the proposed name plus a possible suffix
   */
  protected String getUniqueName(String proposed, String nameToIgnore) {
  	proposed = proposed.trim();
  	if (isDuplicateName(proposed, nameToIgnore)) {
  		int i = 2;
  		String s = proposed+" ("+i+")"; //$NON-NLS-1$ //$NON-NLS-2$
  		while (isDuplicateName(s, nameToIgnore)) {
  			i++;
  			s = proposed+" ("+i+")"; //$NON-NLS-1$ //$NON-NLS-2$
  		}
  		return s;
  	}
  	return proposed;
  }
  
  /**
   * Determines if a name duplicates an existing name.
   * 
   * @param name the proposed name
   * @param nameToIgnore a name that is ignored when comparing
   * @return true if name is a duplicate
   */
  protected boolean isDuplicateName(String name, String nameToIgnore) {
  	// compare with existing names in library and tabbedPane
  	for (String next: library.getNames()) {
  		if (next.equals(nameToIgnore)) continue;
  		if (name.equals(next)) return true;
  	}
  	for (int i=0; i<tabbedPane.getTabCount(); i++) {
  		String title = tabbedPane.getTitleAt(i);
  		if (title.endsWith("*")) //$NON-NLS-1$
  			title = title.substring(0, title.length()-1);
  		if (title.equals(nameToIgnore)) continue;
  		if (name.equals(title)) return true; 		
  	}
  	return false;
  }
  
  /**
   * Creates a new empty LibraryTreePanel.
   * @return the library tree panel
   */
  protected LibraryTreePanel createLibraryTreePanel() {
  	LibraryTreePanel treePanel = new LibraryTreePanel(this);
    treePanel.addPropertyChangeListener(treePanelListener);
  	return treePanel;
  }
  
  /**
   * Shows the about dialog.
   */
  protected void showAboutDialog() {
    String aboutString = getTitle()+" 1.0,  Aug 2011\n"   //$NON-NLS-1$
                         +"Open Source Physics Project\n" //$NON-NLS-1$
                         +"www.opensourcephysics.org";    //$NON-NLS-1$
    JOptionPane.showMessageDialog(this, aboutString, ToolsRes.getString("Dialog.About.Title")+" "+getTitle(), //$NON-NLS-1$ //$NON-NLS-2$
      JOptionPane.INFORMATION_MESSAGE);
  }

  /**
   * Shows the help frame and displays a help HTML page.
   */
  protected void showHelp() {
    String helpPath = XML.getResolvedPath(helpName, helpBase);
    if(ResourceLoader.getResource(helpPath)==null) {
      String classBase = "/org/opensourcephysics/resources/tools/html/"; //$NON-NLS-1$
      helpPath = XML.getResolvedPath(helpName, classBase);
    }
    if((helpFrame==null)||!helpPath.equals(helpFrame.getTitle())) {
      helpFrame = new TextFrame(helpPath);
      helpFrame.enableHyperlinks();
      helpFrame.setSize(760, 560);
      // center on the screen
      Dimension dim = Toolkit.getDefaultToolkit().getScreenSize();
      int x = (dim.width-helpFrame.getBounds().width)/2;
      int y = (dim.height-helpFrame.getBounds().height)/2;
      helpFrame.setLocation(x, y);
    }
    helpFrame.setVisible(true);
  }
  
//______________________________ inner classes _________________________________
  
  class TabTitle extends JPanel {
  	JLabel titleLabel, iconLabel;
  	Icon normalIcon, boldIcon;
  	Action action;
  	
  	TabTitle(Icon lightIcon, Icon heavyIcon) {
  		super(new BorderLayout());
  		this.setOpaque(false);
  		titleLabel = new JLabel();
  		iconLabel = new JLabel();
  		iconLabel.setBorder(BorderFactory.createEmptyBorder(0, 6, 0, 0));
  		iconLabel.addMouseListener(new MouseAdapter() {
        public void mouseClicked(MouseEvent e) {
        	int i = getTabIndexFromTitle(titleLabel.getText());
        	if (i>-1 && tabbedPane.getSelectedIndex()!=i) tabbedPane.setSelectedIndex(i);
        	action.actionPerformed(null);
        }
        public void mouseEntered(MouseEvent e) {
      		iconLabel.setIcon(boldIcon);
        }
        public void mouseExited(MouseEvent e) {
      		iconLabel.setIcon(normalIcon);
        }
  		});
  		add(titleLabel, BorderLayout.WEST);
  		add(iconLabel, BorderLayout.EAST);
  		setIcons(lightIcon, heavyIcon);
  	}
  	
  	void setTitle(String title) {
  		titleLabel.setText(title);
  	}
  	
  	void setIcons(Icon lightIcon, Icon heavyIcon) {
  		normalIcon = lightIcon;
  		boldIcon = heavyIcon;
  		iconLabel.setIcon(normalIcon);
  	}
  	
  	void setAction(Action action) {
  		this.action = action;
  	}

  }

//______________________________ static methods ____________________________

  /**
   * Entry point when run as an independent application.
   * 
   * @param args String[] ignored
   */
  public static void main(String[] args) {
  	final LibraryBrowser browser = LibraryBrowser.getBrowser();
  	String trackerLibraryPath = "http://www.cabrillo.edu/~dbrown/tracker/library/tracker_library.xml"; //$NON-NLS-1$
  	browser.addOSPLibrary(trackerLibraryPath);
  	browser.addComPADRECollection(LibraryComPADRE.EJS_SERVER_TREE);
  	browser.addComPADRECollection(LibraryComPADRE.TRACKER_SERVER_TREE);
  	browser.rebuildCollectionsMenu();
//  	browser.addPropertyChangeListener("target", new PropertyChangeListener() { //$NON-NLS-1$
//  		public void propertyChange(PropertyChangeEvent e) {
//  			LibraryResource record = (LibraryResource)e.getNewValue();
//				String target = XML.getResolvedPath(record.getTarget(), record.getBasePath());
//					  				
//  			ArrayList<String> extensions = new ArrayList<String>();
//  			for (String ext: VideoIO.getVideoExtensions()) {
//  				extensions.add(ext);
//  			}
//  			extensions.add("trk"); //$NON-NLS-1$
//  			extensions.add("zip"); //$NON-NLS-1$
//  			for (String ext: extensions) {
//  				if (target.endsWith("."+ext)) { //$NON-NLS-1$
//  			    Tracker tracker = Tracker.getTracker();
//  			    final TFrame frame = tracker.getFrame();
//  			    frame.setVisible(true);
//            try {
//        			target = ResourceLoader.getURIPath(target);
//							URL url = new URL(target);
//							TrackerIO.open(url, new TrackerPanel(), frame);
//						} catch (Exception ex) {ex.printStackTrace();}
//     			}
//  			}
//  		}
//  	});
 	
  	browser.exitOnClose = true;
    Dimension dim = Toolkit.getDefaultToolkit().getScreenSize();
    int x = (dim.width - browser.getBounds().width) / 2;
    int y = (dim.height - browser.getBounds().height) / 2;
    browser.setLocation(x, y);
    browser.setVisible(true);
  }
    
}
