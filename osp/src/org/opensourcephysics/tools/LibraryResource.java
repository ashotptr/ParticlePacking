/*
 * Open Source Physics software is free software as described near the bottom of this code file.
 *
 * For additional information and documentation on Open Source Physics please see:
 * <http://www.opensourcephysics.org/>
 */

package org.opensourcephysics.tools;

import java.util.ArrayList;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import javax.swing.Icon;
import javax.swing.ImageIcon;

import org.opensourcephysics.controls.XML;
import org.opensourcephysics.controls.XMLControl;
import org.opensourcephysics.media.core.VideoIO;

/**
 * This represents a library resource.
 *
 * @author Douglas Brown
 * @version 1.0
 */
public class LibraryResource {
	
  // static fields
  protected static final String UNKNOWN_TYPE = "Unknown"; //$NON-NLS-1$
  protected static final String COLLECTION_TYPE = "Collection"; //$NON-NLS-1$
  protected static final String TRACKER_TYPE = "Tracker"; //$NON-NLS-1$
  protected static final String EJS_TYPE = "EJS"; //$NON-NLS-1$
  protected static final String VIDEO_TYPE = "Video"; //$NON-NLS-1$
  protected static final String HTML_TYPE = "HTML"; //$NON-NLS-1$
  protected static final String PDF_TYPE = "PDF"; //$NON-NLS-1$
  protected static String[] resourceTypes 
  		= {TRACKER_TYPE, EJS_TYPE, VIDEO_TYPE, HTML_TYPE, PDF_TYPE, UNKNOWN_TYPE};
  protected static ArrayList<String> allTypes = new ArrayList<String>();
  protected static Icon htmlIcon, videoIcon, trackerIcon, ejsIcon, pdfIcon, unknownIcon;

  static {
  	allTypes.add(COLLECTION_TYPE);
  	for (String next: resourceTypes) allTypes.add(next);
    String imageFile = "/org/opensourcephysics/resources/tools/images/html.gif";        //$NON-NLS-1$
    htmlIcon = new ImageIcon(LibraryResource.class.getResource(imageFile));
    imageFile = "/org/opensourcephysics/resources/tools/images/pdf.gif";        //$NON-NLS-1$
    pdfIcon = new ImageIcon(LibraryResource.class.getResource(imageFile));
    imageFile = "/org/opensourcephysics/resources/tools/images/video.gif";        //$NON-NLS-1$
    videoIcon = new ImageIcon(LibraryResource.class.getResource(imageFile));
    imageFile = "/org/opensourcephysics/resources/tools/images/trackericon.gif"; //$NON-NLS-1$
    trackerIcon = new ImageIcon(LibraryResource.class.getResource(imageFile));
    imageFile = "/org/opensourcephysics/resources/tools/images/ejsicon.gif";        //$NON-NLS-1$
    ejsIcon = new ImageIcon(LibraryResource.class.getResource(imageFile));
    imageFile = "/org/opensourcephysics/resources/tools/images/question_mark.gif";        //$NON-NLS-1$
    unknownIcon = new ImageIcon(LibraryResource.class.getResource(imageFile));
  }
		
  // instance fields
	private String name=""; //$NON-NLS-1$
	private String htmlPath=""; // rel or abs path to html that describes this resource //$NON-NLS-1$
	private String basePath=""; // base path for target and/or info //$NON-NLS-1$
	private String description=""; //$NON-NLS-1$
	private String target=""; // rel or abs path to target //$NON-NLS-1$
	private String type=UNKNOWN_TYPE;
  protected String displayName;
  private Map<String, Object> properties = new TreeMap<String, Object>();
  private String thumbnail;
	
  /**
   * Constructor.
   *
   * @param name the name of the resource
   */
	public LibraryResource(String name) {
		setName(name);
	}
	
  /**
   * Gets the name of this resource (never null).
   *
   * @return the name
   */
	public String getName() {
		return name;
	}
	
  /**
   * Sets the name of this resource.
   * 
   * @param aName the name
   * @return true if changed
   */
	public boolean setName(String aName) {
		aName = aName==null? "": aName.trim(); //$NON-NLS-1$
		if (!aName.equals(name)) {
			name = aName;
			return true;
		}
		return false;
	}
	
  /**
   * Gets the base path.
   *
   * @return the base path
   */
	public String getBasePath() {
		return basePath;
	}
	
  /**
   * Sets the base path of this resource.
   * 
   * @param path the base path
   * @return true if changed
   */
	public boolean setBasePath(String path) {
		path = path==null? "": path.trim(); //$NON-NLS-1$
		if (!path.equals(basePath)) {
			basePath = path;
			return true;
		}
		return false;
	}
	
  /**
   * Gets the target of this resource (file name or comPADRE command).
   *
   * @return the target
   */
	public String getTarget() {
		return "".equals(target)? null: target; //$NON-NLS-1$
	}

  /**
   * Sets the target of this resource.
   * 
   * @param path the target path
   * @return true if changed
   */
	public boolean setTarget(String path) {
		path = path==null? "": path.trim(); //$NON-NLS-1$
		if (!path.equals(target)) {
			thumbnail = null;
			target = path;
			path = path.toUpperCase();
			if (path.endsWith(".TRK")) //$NON-NLS-1$
				setType(LibraryResource.TRACKER_TYPE);
			else if (path.endsWith(".PDF")) //$NON-NLS-1$
				setType(LibraryResource.PDF_TYPE);
			else if (path.indexOf("EJS")>-1) { //$NON-NLS-1$
				setType(LibraryResource.EJS_TYPE);
  		}
  		else if (path.endsWith(".ZIP")) { //$NON-NLS-1$
		    Runnable runner = new Runnable() {
		      public void run() {
	  				String zipPath = XML.getResolvedPath(target, getBasePath());
						Set<String> files = ResourceLoader.getZipContents(zipPath);
						for (String next: files) {
							if (next.toUpperCase().endsWith(".TRK")) { //$NON-NLS-1$
								setType(LibraryResource.TRACKER_TYPE);
								break;
							}
						}
		      }
		    };
		    new Thread(runner).start();
			}
  		else if (path.equals("")) { //$NON-NLS-1$
  			if (getHTMLPath()==null)
  				setType(LibraryResource.UNKNOWN_TYPE);
  			else setType(LibraryResource.HTML_TYPE);
  		}
  		else for (String ext: VideoIO.getVideoExtensions()) {
				if (path.indexOf("."+ext.toUpperCase())>-1) { //$NON-NLS-1$
					setType(LibraryResource.VIDEO_TYPE);
					
	  		}
			}
			return true;
		}
		return false;
	}
	
  /**
   * Gets the path to the html page displayed in the browser.
   *
   * @return the html path
   */
	public String getHTMLPath() {
		return htmlPath;
	}
	  	
  /**
   * Sets the html path of this resource.
   * 
   * @param path the html path
   * @return true if changed
   */
	public boolean setHTMLPath(String path) {
		path = path==null? "": path.trim(); //$NON-NLS-1$
		if (!path.equals(htmlPath)) {
			htmlPath = path;
			if (!(this instanceof LibraryCollection) 
					&& getTarget()==null) {
				if (path.equals("")) {//$NON-NLS-1$
					setType(LibraryResource.UNKNOWN_TYPE);
				}
				else setType(LibraryResource.HTML_TYPE);
			}
			return true;
		}
		return false;
	}
	
  /**
   * Gets the description, which must be in html code.
   *
   * @return the description
   */
	public String getDescription() {
		return description;
	}
  	
  /**
   * Sets the description of this resource.
   * Note: the description must be in html code, since it is displayed
   * in the html pane of the LibraryTreePanel if the html path is empty.
   * 
   * @param desc the description in HTML code
   * @return true if changed
   */
	public boolean setDescription(String desc) {
		desc = desc==null? "": desc.trim(); //$NON-NLS-1$
		if (!desc.equals(description)) {
			description = desc;
			return true;
		}
		return false;
	}
  	
  /**
   * Gets the type of resource.
   *
   * @return the one of the static constant types defined in this class
   */
	public String getType() {
		return type;
	}
  	
  /**
   * Sets the type of this resource.
   * The types are static constants defined in this class.
   * 
   * @param type the type
   * @return true if changed
   */
	public boolean setType(String type) {
		if (this.type.equals(type))
			return false;
		for (String next: allTypes) {
			if (next.equals(type)) {
				this.type = next;
				return true;
			}
		}
		return false;
	}
	
  /**
   * Sets an arbitrary property.
   * 
   * @param name the name of the property
   * @param value the value of the property
   */
	public void setProperty(String name, Object value) {
		properties.put(name, value);
	}
	
  /**
   * Gets a property value. May return null.
   * 
   * @param name the name of the property
   * @return the value of the property
   */
	public Object getProperty(String name) {
		return properties.get(name);
	}
	
  /**
   * Returns the names of all defined properties.
   * @return a set of names
   */
	public Set<String> getPropertyNames() {
		return properties.keySet();
	}
	
  /**
   * Gets the icon for the tree node associated with this resource.
   *
   * @return the icon
   */
	public Icon getIcon() {
		if (type==TRACKER_TYPE) return trackerIcon;
		if (type==EJS_TYPE) return ejsIcon;
		if (type==VIDEO_TYPE) return videoIcon;
		if (type==HTML_TYPE) return htmlIcon;
		if (type==PDF_TYPE) return pdfIcon;
		return null;
	}
	
  /**
   * Gets the thumbnail of this resource, if any.
   *
   * @return the thumbnail
   */
	public String getThumbnail() {
		return thumbnail;
	}

  /**
   * Gets the thumbnail of this resource, if any.
   *
   * @param imagePath the path to a thumbnail image
   */
	public void setThumbnail(String imagePath) {
		thumbnail = imagePath;
	}

  /**
   * Gets a clone of this resource.
   *
   * @return the clone
   */
	public LibraryResource getClone() {
		LibraryResource resource = new LibraryResource(name);
  	resource.setBasePath(basePath);
  	resource.setTarget(target);
  	resource.setHTMLPath(htmlPath);
  	resource.setDescription(description);
  	resource.setType(type);
  	for (String next: getPropertyNames()) {
  		resource.setProperty(next, getProperty(next));
  	}
		return resource;
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
   * The ObjectLoader class to save and load LibraryResource data.
   */
  static class Loader implements XML.ObjectLoader {

    /**
     * Saves an object's data to an XMLControl.
     *
     * @param control the control to save to
     * @param obj the object to save
     */
    public void saveObject(XMLControl control, Object obj) {
    	LibraryResource res = (LibraryResource)obj;
    	control.setValue("name", res.name); //$NON-NLS-1$
    	if (!"".equals(res.description)) //$NON-NLS-1$
    		control.setValue("description", res.description); //$NON-NLS-1$
    	if (!"".equals(res.htmlPath)) //$NON-NLS-1$
    		control.setValue("html_path", res.htmlPath); //$NON-NLS-1$
    	if (!"".equals(res.basePath)) //$NON-NLS-1$
    		control.setValue("base_path", res.basePath); //$NON-NLS-1$
    	if (!"".equals(res.target)) //$NON-NLS-1$
    		control.setValue("target", res.getTarget()); //$NON-NLS-1$
    	if (!UNKNOWN_TYPE.equals(res.type))
    		control.setValue("type", res.type); //$NON-NLS-1$
    }
    
    /**
     * Creates a new object.
     *
     * @param control the XMLControl with the object data
     * @return the newly created object
     */
    public Object createObject(XMLControl control){
    	String name = control.getString("name"); //$NON-NLS-1$
      return new LibraryResource(name);
    }

    /**
     * Loads an object with data from an XMLControl.
     *
     * @param control the control
     * @param obj the object
     * @return the loaded object
     */
    public Object loadObject(XMLControl control, Object obj) {
    	LibraryResource res = (LibraryResource)obj;
    	// name is loaded in createObject() method
    	res.setDescription(control.getString("description")); //$NON-NLS-1$
    	res.setBasePath(control.getString("base_path")); //$NON-NLS-1$
    	String target = control.getString("target"); //$NON-NLS-1$
    	if (target!=null) res.target = target;
    	res.setHTMLPath(control.getString("html_path")); //$NON-NLS-1$
    	res.setType(control.getString("type")); //$NON-NLS-1$
    	return res;
    }
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
