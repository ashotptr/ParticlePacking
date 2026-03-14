/*
 * Open Source Physics software is free software as described near the bottom of this code file.
 *
 * For additional information and documentation on Open Source Physics please see:
 * <http://www.opensourcephysics.org/>
 */

package org.opensourcephysics.tools;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusAdapter;
import java.awt.event.FocusEvent;
import java.awt.font.FontRenderContext;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;

import javax.swing.AbstractListModel;
import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JTextField;
import javax.swing.JToolBar;
import javax.swing.ListModel;
import javax.swing.ListSelectionModel;
import javax.swing.SwingConstants;
import javax.swing.WindowConstants;
import javax.swing.border.Border;
import javax.swing.border.TitledBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import org.opensourcephysics.controls.XML;
import org.opensourcephysics.controls.XMLControl;
import org.opensourcephysics.controls.XMLControlElement;

/**
 * A Library for a LibraryBrowser. Maintains its own collections and provides
 * access to collections from comPADRE, OSP and imported libraries.
 *
 * @author Douglas Brown
 */
public class Library {
	
	private String name; // name of the library
	private ArrayList<String> pathList = new ArrayList<String>();
	private HashMap<String, String> pathToNameMap = new HashMap<String, String>();
	private ArrayList<String> comPADREPathList = new ArrayList<String>();
	private HashMap<String, String> comPADREPathToNameMap = new HashMap<String, String>();
	private ArrayList<String> ospPathList = new ArrayList<String>();
	private HashMap<String, Library> ospPathToLibraryMap = new HashMap<String, Library>();
	private ArrayList<String> importedPathList = new ArrayList<String>();
	private HashMap<String, Library> importedPathToLibraryMap = new HashMap<String, Library>();
	private ArrayList<String> subPathList = new ArrayList<String>();
	private HashMap<String, Library> subPathToLibraryMap = new HashMap<String, Library>();
	protected String[] openTabPaths;
	protected Manager manager;
	protected String chooserDir;
	protected LibraryBrowser browser;
	
	/**
	 * Adds an OSP-sponsored library. OSP libraries are not under user control.
	 * 
	 * @param path the library path
	 * @return true if successfully added
	 */
	public boolean addOSPLibrary(String path) {		
  	if (ospPathList.contains(path))
  		return false;
  	XMLControl control = new XMLControlElement(path);
  	if (control.failedToRead() || control.getObjectClass()!=Library.class)
  		return false;
  	Library library = new Library();
  	control.loadObject(library);
  	library.browser = this.browser;
  	ospPathList.add(path);
  	ospPathToLibraryMap.put(path, library);
  	return true;
	}

	/**
	 * Imports a library. Imported libraries are managed by the user.
	 * 
	 * @param path the library path
	 * @return true if successfully imported
	 */
	public boolean importLibrary(String path) {		
  	if (importedPathList.contains(path))
  		return false;
  	XMLControl control = new XMLControlElement(path);
  	if (control.failedToRead() || control.getObjectClass()!=Library.class)
  		return false;
  	Library library = new Library();
  	library.browser = this.browser;
  	control.loadObject(library);
  	return importLibrary(path, library);
	}
	
	/**
	 * Adds a comPADRE collection. ComPADRE collections are not under user control.
	 * 
	 * @param path the comPADRE query
	 * @param name the name of the collection
	 * @return true if successfully added
	 */
	public boolean addComPADRECollection(String path, String name) {
		path = path.trim();
		// don't add duplicate paths
  	if (comPADREPathList.contains(path))
  		return false;
		comPADREPathList.add(path);
		comPADREPathToNameMap.put(path, name.trim());
		return true;
	}
	
	/**
	 * Adds a sublibrary. Sublibraries are shown as submenus in a Library's Collections menu.
	 * Sublibraries are not under user control.
	 * 
	 * @param lib the library
	 * @return true if successfully added
	 */
	public boolean addSubLibrary(String path) {		
  	if (subPathList.contains(path))
  		return false;
  	XMLControl control = new XMLControlElement(path);
  	if (control.failedToRead() || control.getObjectClass()!=Library.class)
  		return false;
  	Library library = new Library();
  	library.browser = this.browser;
  	control.loadObject(library);
  	subPathList.add(path);
  	subPathToLibraryMap.put(path, library);
  	return true;
	}

	/**
	 * Returns a string representation of this library.
	 * 
	 * @return the name of the library
	 */
	@Override
  public String toString() {
  	return getName();
  }
  
//_____________________ protected and private methods _________________________
	
	/**
	 * Sets the name of this library.
	 * 
	 * @param name the name
	 */
	protected void setName(String name) {
		if (name==null) {
      name = System.getProperty("user.home").replace('\\', '/'); //$NON-NLS-1$
      if(name.endsWith("/")) {                                         //$NON-NLS-1$ 
        name = name.substring(0, name.length()-1); 
      }
      name = XML.getName(name)+" "+ToolsRes.getString("Library.Name");  //$NON-NLS-1$//$NON-NLS-2$
		}
		this.name = name;
	}

	/**
	 * Gets the name of this library.
	 * 
	 * @return the name
	 */
	protected String getName() {
		return name;
	}

	/**
	 * Saves this library in an xml file.
	 * 
	 * @param path the path to the saved file
	 */
	protected void save(String path) {
		if (path==null) return;
  	XMLControl control = new XMLControlElement(this);
  	control.write(path);
	}
	
	/**
	 * Loads this library from an xml file.
	 * 
	 * @param path the path to the file
	 */
	protected void load(String path) {
		if (path==null) return;
  	XMLControl control = new XMLControlElement(path);
  	control.loadObject(this);
	}
	
	/**
	 * Gets the names of all collections maintained by this library.
	 * 
	 * @return a collection of names
	 */
	protected Collection<String> getNames() {
		return pathToNameMap.values();
	}
	
	/**
	 * Returns true if this library has no collections.
	 * 
	 * @return true if empty
	 */
	protected boolean isEmpty() {
		return pathList.isEmpty();
	}
	
	/**
	 * Returns true if this library contains a collection path.
	 * 
	 * @param path the collection path
	 * @param allLists true to search in all collection lists
	 * @return true if this contains the path
	 */
	protected boolean containsPath(String path, boolean allLists) {
		path = path.trim();
		int n = path.indexOf(LibraryComPADRE.PRIMARY_ONLY);
		if (n>-1)
			path = path.substring(0, n);
		boolean containsPath = pathList.contains(path);
		if (allLists) {
			containsPath = containsPath 
					|| comPADREPathList.contains(path) 
					|| ospPathList.contains(path);
		}
		return containsPath;
	}
	
	/**
	 * Adds a collection to this library.
	 * 
	 * @param path the path to the collection
	 * @param name the menu item name for the collection
	 */
	protected void addCollection(String path, String name) {
		path = path.trim();
		// don't add duplicate paths
  	if (pathList.contains(path))
  		return;
		pathList.add(path);
		pathToNameMap.put(path, name.trim());
	}
	
	/**
	 * Renames a collection.
	 * 
	 * @param path the path to the collection
	 * @param newName the new name
	 */
	protected void renameCollection(String path, String newName) {
		path = path.trim();
		// change only paths that have already been added
  	if (!pathList.contains(path))
  		return;
		pathToNameMap.put(path, newName.trim());
  }
	
	/**
	 * Gets a clone of this library that is suitable for exporting. The exported
	 * library has no OSP libraries, ComPADRE collections or imported libraries.
	 * 
	 * @return a Library for export
	 */
	protected Library getCloneForExport() {
  	Library lib = new Library();
  	lib.pathList = pathList;
  	lib.pathToNameMap = pathToNameMap;
  	lib.name = name;
  	return lib;
	}

	/**
	 * Gets the collections manager for this library.
	 * 
	 * @param browser a LibraryBrowser
	 * @return the collections manager
	 */
	protected Manager getManager(LibraryBrowser browser) {
		if (manager==null) {
			if (LibraryBrowser.externalDialog!=null)
				manager = new Manager(browser, LibraryBrowser.externalDialog);
			else
				manager = new Manager(browser, LibraryBrowser.frame);
	    // center on screen
	    Dimension dim = Toolkit.getDefaultToolkit().getScreenSize();
	    int x = (dim.width-manager.getBounds().width)/2;
	    int y = (dim.height-manager.getBounds().height)/2;
	    manager.setLocation(x, y);
		}
		if (pathList.size()>0 && manager.collectionList.getSelectedIndex()==-1) {
			manager.collectionList.setSelectedIndex(0);
		}
		if (importedPathList.size()>0 && manager.guestList.getSelectedIndex()==-1) {
			manager.guestList.setSelectedIndex(0);
		}
		manager.refreshGUI();
		return manager;
	}
  
	/**
	 * Rebuilds the Collections menu.
	 * 
	 * @param browser the LibraryBrowser that will display the menu
	 */
  synchronized protected void rebuildCollectionsMenu(LibraryBrowser browser) {
  	JMenu menu = browser.collectionsMenu;
  	menu.removeAll();
		JMenu myLibraryMenu = new JMenu(ToolsRes.getString("Library.Name.Local")); //$NON-NLS-1$
		menu.add(myLibraryMenu);
  	if (!pathList.isEmpty()) {
	  	for (String path: pathList) {
	  		String name = pathToNameMap.get(path);
	      JMenuItem item = new JMenuItem(name);
	      myLibraryMenu.add(item);
	      item.addActionListener(browser.loadCollectionAction);
	      item.setToolTipText(path);
	      item.setActionCommand(path);
	  	}
  	}
  	if (!comPADREPathList.isEmpty()) {
  		JMenu submenu = new JMenu(ToolsRes.getString("Library.Name.ComPADRE")); //$NON-NLS-1$
  		menu.add(submenu);
	  	for (String path: comPADREPathList) {
	  		String name = comPADREPathToNameMap.get(path);      
	  		JMenuItem item = new JMenuItem(name);
	  		submenu.add(item);
	      item.addActionListener(browser.loadCollectionAction);
//	  		if (LibraryComPADRE.primary_only)
//	  			path += LibraryComPADRE.PRIMARY_ONLY;
	      item.setToolTipText(path);
	      item.setActionCommand(path);
	  	}
  	}
  	if (!ospPathList.isEmpty()) {
	  	for (String path: ospPathList) {
	  		Library library = ospPathToLibraryMap.get(path);
	  		JMenu submenu = new JMenu(library.getName());
	  		menu.add(submenu);
	  		populateSubMenu(submenu, library, browser);
	  	}
  	}
  	if (!importedPathList.isEmpty()) {
	  	menu.addSeparator();
	  	for (String path: importedPathList) {
	  		Library library = importedPathToLibraryMap.get(path);
	  		JMenu submenu = new JMenu(library.getName());
	  		menu.add(submenu);
	  		for (String next: library.pathList) {
	    		String name = library.pathToNameMap.get(next);
	        JMenuItem item = new JMenuItem(name);
	        submenu.add(item);
	        item.addActionListener(browser.loadCollectionAction);
	        item.setToolTipText(next);
	        item.setActionCommand(next);	  			
	  		}
	  	}
  	}
  	menu.addSeparator();
  	menu.add(browser.libraryManagerItem);
  }
  
	/**
	 * Rebuilds the Collections menu.
	 * 
	 * @param browser the LibraryBrowser that will display the menu
	 */
  private void populateSubMenu(JMenu menu, Library library, LibraryBrowser browser) {
		for (String next: library.pathList) {
  		String name = library.pathToNameMap.get(next);
      JMenuItem item = new JMenuItem(name);
      menu.add(item);
      item.addActionListener(browser.loadCollectionAction);
      item.setToolTipText(next);
      item.setActionCommand(next);	  			
		}
  	if (!library.subPathList.isEmpty()) {
	  	for (String path: library.subPathList) {
	  		Library lib = library.subPathToLibraryMap.get(path);
	  		JMenu submenu = new JMenu(lib.getName());
	  		menu.add(submenu);
	  		populateSubMenu(submenu, lib, browser);
	  	}
  	}

  }

    
	/**
	 * Imports a Library if not already imported.
	 * 
	 * @param path the path to the library
	 * @param library the library
	 * @return true if imported
	 */
	protected boolean importLibrary(String path, Library library) {		
  	if (importedPathList.contains(path))
  		return false;
  	importedPathList.add(path);
  	importedPathToLibraryMap.put(path, library);
  	return true;
	}
	
	/**
	 * A class that enables the user to manage this library's collections.
	 */
  class Manager extends JDialog {
  	
  	LibraryBrowser browser;
  	JTabbedPane tabbedPane;
  	JPanel collectionsPanel, importedPanel;
  	JList collectionList, guestList;
  	JTextField nameField, pathField;
  	ActionListener nameAction, pathAction;
  	JButton moveUpButton, moveDownButton, addButton, removeButton, renameButton;
  	JToolBar buttonbar;
  	Box namebar, pathbar, menuItemPanel;
  	JLabel nameLabel, pathLabel;
  	TitledBorder menuTitleBorder, importTitleBorder;
  	
  	/**
  	 * Constructor for a frame
  	 * 
  	 * @param browser a LibraryBrowser
  	 * @param frame the frame
  	 */
  	protected Manager(LibraryBrowser browser, JFrame frame) {
  		super(frame, false);
  		this.browser = browser;
      setDefaultCloseOperation(WindowConstants.HIDE_ON_CLOSE);
      createGUI();
      setSize(new Dimension(500, 300));  		
  	}
  	
  	/**
  	 * Constructor for a dialog
  	 * 
  	 * @param browser a LibraryBrowser
  	 * @param dialog the dialog
  	 */
  	protected Manager(LibraryBrowser browser, JDialog dialog) {
  		super(dialog, false);
  		this.browser = browser;
      setDefaultCloseOperation(WindowConstants.HIDE_ON_CLOSE);
      createGUI();
      setSize(new Dimension(500, 300));  		
  	}
  	
  	/**
  	 * Creates the GUI.
  	 */
  	protected void createGUI() {
  		// create collections list
  		ListModel collectionListModel = new AbstractListModel() {
        public int getSize() {
        	return pathList.size();
        }
        public Object getElementAt(int i) { 
        	String path = pathList.get(i);
        	return pathToNameMap.get(path);
        }
      };
  		collectionList = new JList(collectionListModel);
  		collectionList.addListSelectionListener(new ListSelectionListener() {
  			public void valueChanged(ListSelectionEvent e) {
  				refreshGUI();
  			}
  		});
  		collectionList.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);

  		// create import list
  		ListModel importListModel = new AbstractListModel() {
        public int getSize() {
        	return importedPathList.size();
        }
        public Object getElementAt(int i) { 
        	String path = importedPathList.get(i);
        	return importedPathToLibraryMap.get(path).getName();
        }
      };
  		guestList = new JList(importListModel);
  		guestList.addListSelectionListener(new ListSelectionListener() {
  			public void valueChanged(ListSelectionEvent e) {
  				refreshGUI();
  			}
  		});
  		guestList.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);

  		// create name action, field and label
  		nameAction = new ActionListener() {
        public void actionPerformed(ActionEvent e) {
        	String path = pathField.getText();
    			String prev = pathToNameMap.get(path);
        	String input = nameField.getText().trim();
          if(input==null || input.equals("") || input.equals(prev)) { //$NON-NLS-1$
            return;
          }
          renameCollection(path, input);
    			rebuildCollectionsMenu(browser);
    			collectionList.repaint();
    			refreshGUI();
       	}
      };
	  	nameField = new LibraryTreePanel.EntryField();
	  	nameField.addActionListener(nameAction);
	  	nameField.addFocusListener(new FocusAdapter() {
	      public void focusGained(FocusEvent e) {
	      	nameField.selectAll();
	      }
	      public void focusLost(FocusEvent e) {
	      	nameAction.actionPerformed(null);
	      }
	    });
//  		nameField.setEditable(false);
	  	nameField.setBackground(Color.white);
	
	  	nameLabel = new JLabel();
	  	nameLabel.setBorder(BorderFactory.createEmptyBorder(0, 4, 0, 2));
	  	nameLabel.setHorizontalAlignment(SwingConstants.TRAILING);
  	
  		// create path action, field and label
  		pathAction = new ActionListener() {
        public void actionPerformed(ActionEvent e) {
  	    	int i = collectionList.getSelectedIndex();
        	String path = pathList.get(i);
    			String name = pathToNameMap.get(path);
        	String input = pathField.getText().trim();
          if(input==null || input.equals("") || input.equals(path)) { //$NON-NLS-1$
            return;
          }
          pathList.remove(i);
          pathList.add(i, input);
      		pathToNameMap.remove(path);
      		pathToNameMap.put(input, name);

//          renameCollection(path, input);
    			rebuildCollectionsMenu(browser);
    			collectionList.repaint();
    			refreshGUI();
       	}
      };
    	pathField = new LibraryTreePanel.EntryField();
	  	pathField.addActionListener(pathAction);
    	pathField.addFocusListener(new FocusAdapter() {
        public void focusGained(FocusEvent e) {
        	pathField.selectAll();
        }
        public void focusLost(FocusEvent e) {
        	pathAction.actionPerformed(null);
        }
      });
//    	pathField.setEditable(false);
    	pathField.setBackground(Color.white);

    	pathLabel = new JLabel();
      pathLabel.setBorder(BorderFactory.createEmptyBorder(0, 4, 0, 2));
      pathLabel.setHorizontalAlignment(SwingConstants.TRAILING);
    	
      // create buttons
      moveUpButton = new JButton();
      moveUpButton.setOpaque(false);
      moveUpButton.setBorder(LibraryBrowser.buttonBorder);
      moveUpButton.addActionListener(new ActionListener() {
        public void actionPerformed(ActionEvent e) {
        	boolean isImports = tabbedPane.getSelectedComponent()==importedPanel;
        	JList list = isImports? guestList: collectionList;
        	ArrayList<String> paths = isImports? importedPathList: pathList;
  				int i = list.getSelectedIndex();
  				String path = paths.get(i);
  				paths.remove(path);
  				paths.add(i-1, path);
        	list.setSelectedIndex(i-1);
    			rebuildCollectionsMenu(browser);
    			browser.refreshGUI();
       	}
      });
      moveDownButton = new JButton();
      moveDownButton.setOpaque(false);
      moveDownButton.setBorder(LibraryBrowser.buttonBorder);
      moveDownButton.addActionListener(new ActionListener() {
        public void actionPerformed(ActionEvent e) {
        	boolean isImports = tabbedPane.getSelectedComponent()==importedPanel;
        	JList list = isImports? guestList: collectionList;
        	ArrayList<String> paths = isImports? importedPathList: pathList;
  				int i = list.getSelectedIndex();
  				String path = paths.get(i);
  				paths.remove(path);
  				paths.add(i+1, path);
        	list.setSelectedIndex(i+1);
    			rebuildCollectionsMenu(browser);
    			browser.refreshGUI();
       	}
      });
      addButton = new JButton();
      addButton.setOpaque(false);
      addButton.setBorder(LibraryBrowser.buttonBorder);
      addButton.addActionListener(new ActionListener() {
        public void actionPerformed(ActionEvent e) {
        	boolean imported = tabbedPane.getSelectedComponent()==importedPanel;
        	String message = imported? 
        			ToolsRes.getString("LibraryBrowser.Dialog.AddLibrary.Message"): //$NON-NLS-1$
        			ToolsRes.getString("LibraryBrowser.Dialog.AddCollection.Message"); //$NON-NLS-1$
          String title = imported? 
            	ToolsRes.getString("LibraryBrowser.Dialog.AddLibrary.Title"): //$NON-NLS-1$
            	ToolsRes.getString("LibraryBrowser.Dialog.AddCollection.Title"); //$NON-NLS-1$
          Object input = JOptionPane.showInputDialog(browser, 
          		message, title, JOptionPane.QUESTION_MESSAGE, null, null, null);
          if(input==null || input.equals("")) {                            //$NON-NLS-1$
            return;
          }
          String path = input.toString();
          path = XML.forwardSlash(path);
          path = ResourceLoader.getNonURIPath(path);
    			XMLControl control = new XMLControlElement(path);
    			if (!control.failedToRead()) { 
            // see if path points to a collection
    				if (!imported && control.getObjectClass()==LibraryCollection.class) {
    					
  	      		browser.addToCollections(path);
  	      		refreshGUI();
  	      		collectionList.repaint();
  	      		collectionList.setSelectedIndex(pathList.size()-1);
  						rebuildCollectionsMenu(browser);
  						return;
  	  			}
      			// see if path points to a library
    				if (imported && Library.class==control.getObjectClass()) {
    			  	Library newLibrary = new Library();
    			  	newLibrary.browser = Library.this.browser;
    			  	control.loadObject(newLibrary);
    					if (importLibrary(path, newLibrary)) {
    	      		refreshGUI();
    	      		guestList.repaint();
    	        	guestList.setSelectedIndex(importedPathList.size()-1);
    						rebuildCollectionsMenu(browser);
    					}
    					return;
    				}
     			}
  		  	String s = ToolsRes.getString("LibraryBrowser.Dialog.CollectionNotFound.Message"); //$NON-NLS-1$
  		  	JOptionPane.showMessageDialog(manager, 
  		  			s+":\n"+path, //$NON-NLS-1$
  						ToolsRes.getString("LibraryBrowser.Dialog.CollectionNotFound.Title"), //$NON-NLS-1$
  						JOptionPane.WARNING_MESSAGE);  		
          
       	}
      });
      removeButton = new JButton();
      removeButton.setOpaque(false);
      removeButton.setBorder(LibraryBrowser.buttonBorder);
      removeButton.addActionListener(new ActionListener() {
        public void actionPerformed(ActionEvent e) {
        	boolean isImports = tabbedPane.getSelectedComponent()==importedPanel;
        	JList list = isImports? guestList: collectionList;
        	ArrayList<String> paths = isImports? importedPathList: pathList;
  				int i = list.getSelectedIndex();
  				String path = paths.get(i);
  				paths.remove(path);
  				if (isImports)
  					importedPathToLibraryMap.remove(path);
  				else
  					pathToNameMap.remove(path);
  				list.repaint();
        	if (i>=paths.size()) {
        		list.setSelectedIndex(paths.size()-1);
        	}
    			rebuildCollectionsMenu(browser);
    			refreshGUI();
    			browser.refreshGUI();
       	}
      });
      
//      renameButton = new JButton();
//      renameButton.setOpaque(false);
//      renameButton.setBorder(LibraryBrowser.buttonBorder);
//      renameButton.addActionListener(new ActionListener() {
//        public void actionPerformed(ActionEvent e) {
//        	String path = pathField.getText();
//    			String proposed = pathToNameMap.get(path);
//          Object input = JOptionPane.showInputDialog(browser, 
//          		ToolsRes.getString("LibraryManager.Dialog.Rename.Message")+":", //$NON-NLS-1$ //$NON-NLS-2$
//          		ToolsRes.getString("LibraryManager.Dialog.Rename.Title"),   //$NON-NLS-1$
//              JOptionPane.QUESTION_MESSAGE, null, null, proposed);
//          if(input==null || input.equals("") || input.equals(proposed)) { //$NON-NLS-1$
//            return;
//          }
//          renameCollection(path, input.toString());
//    			rebuildCollectionsMenu(browser);
//    			collectionList.repaint();
//    			refreshGUI();
//       	}
//      });

      // assemble components
  		collectionsPanel = new JPanel(new BorderLayout());
  		JScrollPane scroller = new JScrollPane(collectionList);
      Border empty = BorderFactory.createEmptyBorder(2, 2, 2, 2);
      Border etched = BorderFactory.createEtchedBorder();
  		collectionList.setBorder(BorderFactory.createCompoundBorder(etched, empty));
  		menuTitleBorder = BorderFactory.createTitledBorder(""); //$NON-NLS-1$
  		scroller.setBorder(menuTitleBorder);
  		collectionsPanel.add(scroller, BorderLayout.CENTER);
  		
      namebar = Box.createHorizontalBox();
      namebar.setBorder(BorderFactory.createEmptyBorder(0, 0, 2, 4));
      namebar.add(nameLabel);
      namebar.add(nameField);    
      pathbar = Box.createHorizontalBox();
      pathbar.setBorder(BorderFactory.createEmptyBorder(0, 0, 2, 4));
      pathbar.add(pathLabel);
      pathbar.add(pathField); 
      menuItemPanel = Box.createVerticalBox();
      menuItemPanel.add(namebar);
      menuItemPanel.add(pathbar);
      
  		collectionsPanel.add(menuItemPanel, BorderLayout.SOUTH);
  		buttonbar = new JToolBar();
      buttonbar.setFloatable(false);
      empty = BorderFactory.createEmptyBorder(1, 2, 1, 2);
      buttonbar.setBorder(BorderFactory.createCompoundBorder(etched, empty));
      buttonbar.add(moveUpButton);
      buttonbar.add(moveDownButton);
      buttonbar.add(addButton);
      buttonbar.add(removeButton);
//      buttonbar.add(renameButton);
  		collectionsPanel.add(buttonbar, BorderLayout.NORTH);
  		
  		importedPanel = new JPanel(new BorderLayout());  		
  		scroller = new JScrollPane(guestList);
      empty = BorderFactory.createEmptyBorder(2, 2, 2, 2);
  		guestList.setBorder(BorderFactory.createCompoundBorder(etched, empty));
  		importTitleBorder = BorderFactory.createTitledBorder(""); //$NON-NLS-1$
  		scroller.setBorder(importTitleBorder);
  		importedPanel.add(scroller, BorderLayout.CENTER);
      
  		// create tabbedPane
  		tabbedPane = new JTabbedPane();
  		setContentPane(tabbedPane);  		
  		tabbedPane.addTab("", collectionsPanel); //$NON-NLS-1$
//  		tabbedPane.addTab("", importedPanel); //$NON-NLS-1$

      // add change listener last
  		tabbedPane.addChangeListener(new ChangeListener() {
        public void stateChanged(ChangeEvent e) {
        	if (tabbedPane.getSelectedComponent()==collectionsPanel) {
        		collectionsPanel.add(buttonbar, BorderLayout.NORTH);
        		collectionsPanel.add(menuItemPanel, BorderLayout.SOUTH);
        		refreshGUI();
        	}
        	else if (tabbedPane.getSelectedComponent()==importedPanel) {
        		importedPanel.add(buttonbar, BorderLayout.NORTH);
        		importedPanel.add(menuItemPanel, BorderLayout.SOUTH);
        		refreshGUI();
        	}
        }
      });
  	}
  	
  	/**
  	 * Refreshes the GUI including locale-based resource strings.
  	 */
  	protected void refreshGUI() {
  		setTitle(ToolsRes.getString("LibraryManager.Title")); //$NON-NLS-1$
  		addButton.setText(ToolsRes.getString("LibraryManager.Button.Add")); //$NON-NLS-1$
  		removeButton.setText(ToolsRes.getString("LibraryManager.Button.Remove")); //$NON-NLS-1$
//  		renameButton.setText(ToolsRes.getString("LibraryManager.Button.Rename")); //$NON-NLS-1$
    	moveUpButton.setText(ToolsRes.getString("LibraryTreePanel.Button.Up")); //$NON-NLS-1$
    	moveDownButton.setText(ToolsRes.getString("LibraryTreePanel.Button.Down")); //$NON-NLS-1$
    	nameLabel.setText(ToolsRes.getString("LibraryManager.Label.Name")+":"); //$NON-NLS-1$ //$NON-NLS-2$
    	pathLabel.setText(ToolsRes.getString("LibraryManager.Label.Path")+":"); //$NON-NLS-1$ //$NON-NLS-2$
			menuTitleBorder.setTitle(ToolsRes.getString("LibraryManager.Title.MenuItems")+":"); //$NON-NLS-1$ //$NON-NLS-2$
			importTitleBorder.setTitle(ToolsRes.getString("LibraryManager.Title.Import")+":"); //$NON-NLS-1$ //$NON-NLS-2$
    	tabbedPane.setTitleAt(0, ToolsRes.getString("Library.Name.Local")); //$NON-NLS-1$
//    	tabbedPane.setTitleAt(1, ToolsRes.getString("LibraryManager.Tab.Import")); //$NON-NLS-1$
			
  		// adjust size of labels so they right-align
      int w = 0;
      Font font = nameLabel.getFont();
      FontRenderContext frc = new FontRenderContext(null, false, false); 
      Rectangle2D rect = font.getStringBounds(nameLabel.getText()+" ", frc); //$NON-NLS-1$
      w = Math.max(w, (int) rect.getWidth()+4);
      rect = font.getStringBounds(pathLabel.getText()+" ", frc); //$NON-NLS-1$
      w = Math.max(w, (int) rect.getWidth()+4);

      Dimension labelSize = new Dimension(w, 20);
      nameLabel.setPreferredSize(labelSize);
      nameLabel.setMinimumSize(labelSize);
      pathLabel.setPreferredSize(labelSize);
      pathLabel.setMinimumSize(labelSize);

    	pathField.setForeground(LibraryTreePanel.defaultForeground);
    	if (tabbedPane.getSelectedComponent()==collectionsPanel) {
	    	int i = collectionList.getSelectedIndex();
	  		moveDownButton.setEnabled(i<pathList.size()-1);
	  		moveUpButton.setEnabled(i>0);
	  		if (i>-1 && pathList.size()>i) {
		  		removeButton.setEnabled(true);
//		  		renameButton.setEnabled(true);
					String path = pathList.get(i);
					pathField.setText(path);
					pathField.setCaretPosition(0);
					String name = pathToNameMap.get(path);
					nameField.setText(name);
		      Resource res = ResourceLoader.getResource(path);
		      if (res==null) {
		      	pathField.setForeground(LibraryTreePanel.darkRed);
//		      	pathField.setBackground(LibraryTreePanel.lightRed);
		      }
	  		}
	  		else {
		  		removeButton.setEnabled(false);
//		  		renameButton.setEnabled(true);
					pathField.setText(null);  			
	      	pathField.setBackground(Color.white);
	  		}
			}
			else if (tabbedPane.getSelectedComponent()==importedPanel) {
	    	int i = guestList.getSelectedIndex();
	  		moveDownButton.setEnabled(i<importedPathList.size()-1);
	  		moveUpButton.setEnabled(i>0);
//	  		renameButton.setEnabled(false);
	  		if (i>-1 && importedPathList.size()>i) {
		  		removeButton.setEnabled(true);
					String path = importedPathList.get(i);
					pathField.setText(path);
					pathField.setCaretPosition(0);
	  		}
	  		else {
		  		removeButton.setEnabled(false);
					pathField.setText(null);  			
	  		}
			}
			nameField.setBackground(Color.white);
    	pathField.setBackground(Color.white);
  	}
  }

  /**
   * Returns an ObjectLoader to save and load data for this class.
   *
   * @return the object loader
   */
  public static XML.ObjectLoader getLoader() {
    return new Loader();
  }
	
  /**
   * A class to save and load data for this class.
   */
  static class Loader implements XML.ObjectLoader {
  	
    /**
     * Saves an object's data to an XMLControl.
     *
     * @param control the control to save to
     * @param obj the object to save
     */
    public void saveObject(XMLControl control, Object obj) {
    	Library library = (Library)obj;
    	control.setValue("name", library.getName()); //$NON-NLS-1$
    	if (!library.pathList.isEmpty()) {
	    	String[] paths = library.pathList.toArray(new String[0]);
	    	control.setValue("collection_paths", paths); //$NON-NLS-1$
	    	String[] names = new String[paths.length];
	    	for (int i=0; i< paths.length; i++) {
	    		names[i] = library.pathToNameMap.get(paths[i]);
	    	}
	    	control.setValue("collection_names", names); //$NON-NLS-1$
    	}
    	if (!library.subPathList.isEmpty()) {
	    	String[] paths = library.subPathList.toArray(new String[0]);
	    	control.setValue("sublibrary_paths", paths); //$NON-NLS-1$
    	}
    	if (!library.importedPathList.isEmpty()) {
	    	String[] paths = library.importedPathList.toArray(new String[0]);
	    	control.setValue("imported_library_paths", paths); //$NON-NLS-1$
    	}
    	control.setValue("open_tabs", library.openTabPaths); //$NON-NLS-1$
    	control.setValue("chooser_directory", library.chooserDir); //$NON-NLS-1$
    }
    
    /**
     * Creates a new object.
     *
     * @param control the XMLControl with the object data
     * @return the newly created object
     */
    public Object createObject(XMLControl control) {
      return new Library();
    }

    /**
     * Loads an object with data from an XMLControl.
     *
     * @param control the control
     * @param obj the object
     * @return the loaded object
     */
    public Object loadObject(XMLControl control, Object obj) {
    	final Library library = (Library)obj;
    	library.setName(control.getString("name")); //$NON-NLS-1$
    	String[] paths = (String[])control.getObject("collection_paths"); //$NON-NLS-1$
    	if (paths!=null) {
      	String[] names = (String[])control.getObject("collection_names"); //$NON-NLS-1$
      	library.pathList.clear();
      	library.pathToNameMap.clear();
    		for (int i=0; i<paths.length; i++) {
    			if (paths[i]==null || names[i]==null) continue;
    			library.pathList.add(paths[i]);
    			library.pathToNameMap.put(paths[i], names[i]);
    		}
    	}
    	paths = (String[])control.getObject("sublibrary_paths"); //$NON-NLS-1$
    	if (paths!=null) {
    		final String[] subs = paths;
    		Runnable runner = new Runnable() {
    			public void run() {
    				for (String path: subs) {
    					library.addSubLibrary(path);
    				}
    				if (library.browser!=null) {
    					library.browser.rebuildCollectionsMenu();
    				}
    			}
    		};
    		new Thread(runner).start();
    	}
    	paths = (String[])control.getObject("imported_library_paths"); //$NON-NLS-1$
    	if (paths!=null) {
    		final String[] imports = paths;
    		Runnable runner = new Runnable() {
    			public void run() {
    				for (String path: imports) {
    					library.importLibrary(path);
    				}
    				if (library.browser!=null) {
    					library.browser.rebuildCollectionsMenu();
    				}
    			}
    		};
    		new Thread(runner).start();
    	}
    	library.openTabPaths = (String[])control.getObject("open_tabs"); //$NON-NLS-1$
    	library.chooserDir = control.getString("chooser_directory"); //$NON-NLS-1$
    	return obj;
    }
  }  	
}

