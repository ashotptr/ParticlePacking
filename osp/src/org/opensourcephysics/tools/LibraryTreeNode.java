/*
 * Open Source Physics software is free software as described near the bottom of this code file.
 *
 * For additional information and documentation on Open Source Physics please see:
 * <http://www.opensourcephysics.org/>
 */

package org.opensourcephysics.tools;

import java.awt.Dimension;
import java.awt.Graphics2D;
import java.awt.geom.AffineTransform;
import java.awt.image.BufferedImage;
import java.io.File;
import java.lang.reflect.Method;
import java.net.URL;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Set;

import javax.imageio.ImageIO;
import javax.swing.SwingWorker;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.TreePath;

import org.opensourcephysics.controls.OSPLog;
import org.opensourcephysics.controls.XML;
import org.opensourcephysics.media.core.VideoIO;
import org.opensourcephysics.media.gif.GifDecoder;

/**
 * A DefaultMutableTreeNode for a library collection tree.
 *
 * @author Douglas Brown
 * @version 1.0
 */
public class LibraryTreeNode extends DefaultMutableTreeNode {

  protected static HashMap<String, URL> htmlURLs = new HashMap<String, URL>();
  protected static HashMap<String, URL> targetURLs = new HashMap<String, URL>();
  protected static Dimension thumbnailDimension = new Dimension(160, 120);


  protected LibraryResource record;
  protected boolean editable = true;
  protected LibraryTreePanel treePanel;
  protected ArrayList<LibraryResource> resources = new ArrayList<LibraryResource>();

  /**
   * Constructs a node with a LibraryResource.
   *
   * @param resource the resource
   * @param treePanel the LibraryTreePanel that will use the node
   */
  protected LibraryTreeNode(LibraryResource resource, LibraryTreePanel treePanel) {
  	this.record = resource;
  	this.treePanel = treePanel;
  	if (treePanel.tree!=null)
  		createChildNodes();
    setUserObject(this);
  }
  
  /**
   * Creates the child nodes of this node if this is a collection node.
   *
   * @return true if children were added
   */
  protected boolean createChildNodes() {
  	ArrayList<String> children = new ArrayList<String>();
  	for (int i=0; i< getChildCount(); i++) {
  		children.add(this.getChildAt(i).toString());
  	}
  	boolean changed = false;
    if (record instanceof LibraryCollection) {
    	LibraryCollection collection = (LibraryCollection)record;
    	for (LibraryResource next: collection.getResources()) {
    		if (next!=null && !children.contains(next.getName())) {
    			LibraryTreeNode newNode = new LibraryTreeNode(next, treePanel);
        	if (treePanel.insertChildAt(newNode, this, getChildCount())) {
        		changed = true;
        	}
    		}
    	}
    }
    if (changed) treePanel.collectionChanged();
    return changed;
  }
  
  /**
   * Returns the name of this node's resource.
   * 
   * @return the name
   */
  protected String getName() {
  	return record.getName();
  }
  
  /**
   * Returns the base path of this node's resource.
   * 
   * @return the base path
   */
  protected String getBasePath() {
  	String base = record.getBasePath();
  	if (!base.equals("")) //$NON-NLS-1$
  		return base;
  	LibraryTreeNode parent = (LibraryTreeNode)getParent();
  	if (parent!=null)
  		return parent.getBasePath();
  	return ""; //$NON-NLS-1$
  }

  /**
   * Returns the html path of this node's resource.
   * 
   * @return the html path
   */
  protected String getHTMLPath() {
  	String path = record.getHTMLPath();
		if (path!=null && !path.trim().equals("")) { //$NON-NLS-1$
	  	path = XML.getResolvedPath(path, getBasePath());
			return path;
		}
  	return null;
  }
  
  /**
   * Returns the html URL of this node's resource, or null if html path
   * is empty or invalid.
   * 
   * @return the URL
   */
  protected URL getHTMLURL() {
  	String path = getHTMLPath();
  	if (path==null) return null;
  	// first try to get URL with raw path
  	if (htmlURLs.keySet().contains(path)) { 
  		URL url = htmlURLs.get(path);
  		if (url!=null) return url;
  	}
  	else {
  		Resource res = ResourceLoader.getResource(path);
  		if (res!=null) {
  			URL url = res.getURL();
  			htmlURLs.put(path, url);
  			return url;
  		}
    	// try with URI form of path
    	String uriPath = ResourceLoader.getURIPath(path);
  		res = ResourceLoader.getResource(uriPath);
  		if (res!=null) {
  			URL url = res.getURL();
  			htmlURLs.put(path, url);
  			return url;
  		}
  		htmlURLs.put(path, null);
  	}
  	return null;
  }
  
  /**
   * Returns an HTML string that describes this node's resource.
   * This is displayed if no html URL is available.
   * 
   * @return the html string
   */
  protected String getHTMLString() {
  	if (!record.getDescription().equals("")) { //$NON-NLS-1$
  		return record.getDescription();
  	}
  	StringBuffer buffer = new StringBuffer();
    buffer.append("<h2>"+record.getName()+"</h2>"); //$NON-NLS-1$ //$NON-NLS-2$
    
  	String typeStr = "LibraryResource.Type."+record.getType(); //$NON-NLS-1$
    typeStr = ToolsRes.getString("LibraryResource.Type")+": "+ToolsRes.getString(typeStr); //$NON-NLS-1$ //$NON-NLS-2$
    buffer.append("<blockquote>"+typeStr); //$NON-NLS-1$ 
    
    boolean isThumbnailType = record.getType().equals(LibraryResource.VIDEO_TYPE) && record.getTarget()!=null;
    String thumb = isThumbnailType? record.getThumbnail(): null;
    if (isThumbnailType && thumb==null) {
			final String source = XML.getResolvedPath(getTarget(), getBasePath());
			final String ext = XML.getExtension(source);
			File thumbFile = getThumbnailFile(source);
			final String thumbPath = thumbFile.getAbsolutePath();
			if (thumbFile.exists()) {
				thumb = thumbFile.getAbsolutePath();
				record.setThumbnail(thumb);
			}
			else {
				OSPLog.finer("needs new thumbnail "+thumbFile.getAbsolutePath()); //$NON-NLS-1$
				class ThumbnailLoader extends SwingWorker<File, Object> {
		       @Override
		       public File doInBackground() {
						// create a new thumbnail
						File thumbFile = null; 
						if (ext!=null && "GIF".equals(ext.toUpperCase())) { //$NON-NLS-1$
							GifDecoder decoder = new GifDecoder();
						  int status = decoder.read(source);
						  if(status!=0) { // error
								OSPLog.fine("failed to create thumbnail for GIF "+thumbPath); //$NON-NLS-1$
						  }
						  else {
						  	BufferedImage image = decoder.getImage();
						  	thumbFile = createThumbnailFile(image, thumbPath);
						  }
						}
						else if (ext!=null && "PNG".equals(ext.toUpperCase())) { //$NON-NLS-1$
							try {
								URL url = new URL(ResourceLoader.getURIPath(source));
								BufferedImage image = ImageIO.read(url);
						  	thumbFile = createThumbnailFile(image, thumbPath);
							} catch (Exception e) {
								OSPLog.fine("failed to create thumbnail for PNG "+thumbPath); //$NON-NLS-1$
							}
						}
		 				else {
		 		      String className = "org.opensourcephysics.media.xuggle.XuggleThumbnailTool"; //$NON-NLS-1$
		 		      Class<?>[] types = new Class<?>[] {Dimension.class, String.class, String.class};
		 		      Object[] values = new Object[] {thumbnailDimension, source, thumbPath};
		 			    try {
		 			      Class<?> xuggleClass = Class.forName(className);
		 			      Method method=xuggleClass.getMethod("createThumbnailFile", types); //$NON-NLS-1$
		 			      thumbFile = (File)method.invoke(null, values);
		 					} catch(Exception ex) {
		 						OSPLog.fine("failed to create thumbnail: "+ex.toString()); //$NON-NLS-1$
		 					} catch(Error err) {
		 					}
		 				 }
						return thumbFile;
		       }

		       @Override
		       protected void done() {
	           try {
	          	 File thumbFile = get();
	     	       record.setThumbnail(thumbFile==null || !thumbFile.exists()? null: thumbFile.getAbsolutePath());
	     	       if (treePanel.getSelectedNode()==LibraryTreeNode.this) {
	     	      	 LibraryTreePanel.htmlPanesByNode.remove(LibraryTreeNode.this);      			
	     	      	 treePanel.showInfo(LibraryTreeNode.this);
	     	       }
	           } catch (Exception ignore) {
	           }
		       }
		   }
		 
		   new ThumbnailLoader().execute();
//				// create a new thumbnail
//				if (ext!=null && "GIF".equals(ext.toUpperCase())) { //$NON-NLS-1$
//					GifDecoder decoder = new GifDecoder();
//			    int status = decoder.read(source);
//			    if(status!=0) { // error
//						OSPLog.fine("failed to create thumbnail for GIF "+thumbPath); //$NON-NLS-1$
//			    }
//			    else {
//			    	BufferedImage image = decoder.getImage();
//			    	thumbFile = createThumbnailFile(image, thumbPath);
//			    }
//				}
//				else if (ext!=null && "PNG".equals(ext.toUpperCase())) { //$NON-NLS-1$
//					try {
//						URL url = new URL(ResourceLoader.getURIPath(source));
//						BufferedImage image = ImageIO.read(url);
//			    	thumbFile = createThumbnailFile(image, thumbPath);
//					} catch (Exception e) {
//						OSPLog.fine("failed to create thumbnail for PNG "+thumbPath); //$NON-NLS-1$
//					}
//				}
//				else {
//		      String className = "org.opensourcephysics.media.xuggle.XuggleThumbnailTool"; //$NON-NLS-1$
//		      Class<?>[] types = new Class<?>[] {Dimension.class, String.class, String.class};
//		      Object[] values = new Object[] {thumbnailDimension, source, thumbPath};
//			    try {
//			      Class<?> xuggleClass = Class.forName(className);
//			      Method method=xuggleClass.getMethod("createThumbnailFile", types); //$NON-NLS-1$
//			      thumbFile = (File)method.invoke(null, values);
//					} catch(Exception ex) {
//						OSPLog.fine("failed to create thumbnail: "+ex.toString()); //$NON-NLS-1$
//					} catch(Error err) {
//					}
//				}
//	      record.setThumbnail(thumbFile==null || !thumbFile.exists()? null: thumbFile.getAbsolutePath());
//				thumb = thumbFile==null || !thumbFile.exists()? null: thumbFile.getAbsolutePath();
			}
    }
    if (thumb!=null) {
    	thumb = XML.forwardSlash(thumb);
    	thumb = ResourceLoader.getURIPath(thumb);
      buffer.append ("<p><img src=\""+thumb+"\" alt=\""+record.getName()+"\"></p>"); //$NON-NLS-1$ //$NON-NLS-2$ //$NON-NLS-3$
    }
    
    String collection = " "+ToolsRes.getString("LibraryResource.Type.Collection.Description"); //$NON-NLS-1$ //$NON-NLS-2$
    for (String type: LibraryResource.allTypes) {
    	if (type.equals(LibraryResource.UNKNOWN_TYPE)) continue;
    	if (type.equals(LibraryResource.PDF_TYPE)) continue;
    	String[] types = new String[] {type};    	
    	if (type.equals(LibraryResource.HTML_TYPE)) {
    		type = "Other"; //$NON-NLS-1$
    		types = new String[] {LibraryResource.HTML_TYPE, LibraryResource.PDF_TYPE, LibraryResource.UNKNOWN_TYPE};
    	}
    	ArrayList<LibraryResource> children = getChildResources(types);
      if (!children.isEmpty()) { // node has children
      	String s = "LibraryResource.Type."+type+".List"; //$NON-NLS-1$ //$NON-NLS-2$
        buffer.append("<p>"+ToolsRes.getString(s) //$NON-NLS-1$
        		+" "+toString()+collection+":</p>\n");  //$NON-NLS-1$//$NON-NLS-2$
        buffer.append("<ul>\n"); //$NON-NLS-1$
        for (LibraryResource next: children) {
        	String name = next.getName();
        	if (name.equals("")) //$NON-NLS-1$
        		name = ToolsRes.getString("LibraryResource.Name.Default"); //$NON-NLS-1$
          buffer.append("<li>"+name+"</li>\n"); //$NON-NLS-1$ //$NON-NLS-2$      	
        }
        buffer.append("</ul>\n"); //$NON-NLS-1$
      }   	
    }
    buffer.append("</blockquote>"); //$NON-NLS-1$
    
    return buffer.toString();
  }
  
  /**
   * Returns the target of this node's resource.
   * The target may be absolute or relative to base path.
   * 
   * @return the target
   */
  protected String getTarget() {
		return record.getTarget();
  }

  /**
   * Returns the absolute target path of this node's resource.
   * 
   * @return the absolute target path
   */
  protected String getAbsoluteTarget() {
  	if (getTarget()==null) return null;
  	if (record instanceof LibraryCollection) {
	  	return getBasePath()+getTarget();
  	}
		return XML.getResolvedPath(getTarget(), getBasePath());
  }

  /**
   * Returns the target URL of this node's resource, or null if target path
   * is empty or invalid.
   * 
   * @return the URL
   */
  protected URL getTargetURL() {
  	String path = getAbsoluteTarget();
  	if (path==null) return null;
  	// first try to get URL with raw path
  	if (targetURLs.keySet().contains(path)) { 
  		URL url = targetURLs.get(path);
  		if (url!=null) return url;
  	}
  	else {
  		Resource res = ResourceLoader.getResource(path);
  		if (res!=null) {
  			URL url = res.getURL();
  			targetURLs.put(path, url);
  			return url;
  		}
    	// try with URI form of path
    	String uriPath = ResourceLoader.getURIPath(path);
  		res = ResourceLoader.getResourceZipURLsOK(uriPath);
  		if (res!=null) {
  			URL url = res.getURL();
  			targetURLs.put(path, url);
  			return url;
  		}
  		targetURLs.put(path, null);
  	}
  	return null;
  }
  
  /**
   * Used by the tree node to get the display name.
   *
   * @return the display name of the node
   */
  @Override
  public String toString() {
  	String name = getName();
  	if (name.equals("")) { //$NON-NLS-1$
  		if (record.displayName!=null)
  			return record.displayName;
  		if (record instanceof LibraryCollection)
  			return ToolsRes.getString("LibraryCollection.Name.Default"); //$NON-NLS-1$
			return ToolsRes.getString("LibraryResource.Name.Default"); //$NON-NLS-1$
  	}
  	return name;
  }
  
  /**
   * Determines if this node is editable.
   * Note: returns true only if this and its parent are editable.
   * 
   * @return true of editable
   */
  protected boolean isEditable() {
  	if (isRoot()) return editable;
  	LibraryTreeNode parent = (LibraryTreeNode)getParent();
  	return editable && parent.isEditable();
  }
  
  /**
   * Sets the editable property for this node.
   * 
   * @param edit true to make this node editable
   */
  protected void setEditable(boolean edit) {
  	editable = edit;
  }
  
  /**
   * Sets the name of this node's resource.
   * 
   * @param name the name
   */
  protected void setName(String name) {
		if (record.setName(name)) {
	  	treePanel.tree.getModel().valueForPathChanged(new TreePath(getPath()), name);
			treePanel.showInfo(this);
  		treePanel.collectionChanged();
		}
  }
  
  /**
   * Sets the target of this node's resource. May be absolute or relative path.
   * 
   * @param path the target path
   * @return true if changed
   */
  protected boolean setTarget(String path) {
		if (record.setTarget(path)) {
			// target has changed
			if (path==null) path = ""; //$NON-NLS-1$
			if (path.endsWith(".trk")) //$NON-NLS-1$
				setType(LibraryResource.TRACKER_TYPE);
			else if (path.indexOf("EJS")>-1) { //$NON-NLS-1$
				setType(LibraryResource.EJS_TYPE);
  		}
  		else if (path.endsWith(".zip")) { //$NON-NLS-1$
		    Runnable runner = new Runnable() {
		      public void run() {
	  				String zipPath = getAbsoluteTarget();
						Set<String> files = ResourceLoader.getZipContents(zipPath);
						for (String next: files) {
							if (next.endsWith(".trk")) { //$NON-NLS-1$
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
				if (path.toUpperCase().endsWith("."+ext.toUpperCase())) { //$NON-NLS-1$
					setType(LibraryResource.VIDEO_TYPE);					
	  		}
			}
  		LibraryTreePanel.htmlPanesByNode.remove(this);
  		record.setThumbnail(null);
			treePanel.showInfo(this);
  		treePanel.collectionChanged();
  		return true;
		}
		return false;
  }
  
  /**
   * Sets the html path of this node's resource.
   * 
   * @param path the html path
   */
  protected void setHTMLPath(String path) {
  	if (record.setHTMLPath(path)) {
  		treePanel.showInfo(this);
  		treePanel.collectionChanged();
  	}
  }
  
  /**
   * Sets the base path of this node's resource.
   * 
   * @param path the base path
   */
  protected void setBasePath(String path) {
		if (record.setBasePath(path)) {
  		LibraryTreePanel.htmlPanesByNode.remove(this);      			
  		record.setThumbnail(null);
			treePanel.showInfo(this);
  		treePanel.collectionChanged();
		}
  }
  
  /**
   * Sets the type of this node's resource.
   * The types are static constants defined by LibraryResource.
   * 
   * @param type the type
   */
  protected void setType(String type) {
  	if (record.setType(type)) {
  		LibraryTreePanel.htmlPanesByNode.remove(this);      			
  		treePanel.showInfo(this);
  		treePanel.collectionChanged();
  	}
  }
  
  /**
   * Returns this node's child resources, if any, of a given set of types.
   * The types are static constants defined by LibraryResource.
   * 
   * @param types an array of resource types
   * @return a list of LibraryResources
   */
  protected ArrayList<LibraryResource> getChildResources(String[] types) {
  	resources.clear();
  	for (String type: types) {
	    for (int i=0; i<getChildCount(); i++) {
	    	LibraryTreeNode child = (LibraryTreeNode)getChildAt(i);
	    	if (child.record.getType().equals(type))
	    		resources.add(child.record);
	    }
  	}
  	return resources;
  }
  
  /**
   * Returns a File that points to the cached thumbnail, if any, for a video file.
   * Note: the thumbnail file may not exist--this just determines where it should be.
   * 
   * @param videoURL the path to the video file 
   * @return the thumbnail File
   */
  protected File getThumbnailFile(String videoURL) {
		String name = XML.stripExtension(XML.getName(videoURL))+".png"; //$NON-NLS-1$
		File thumbFile = ResourceLoader.getOSPCacheHostDirectory(videoURL);
		thumbFile = new File(thumbFile, "thumbnails"); //$NON-NLS-1$
		return new File(thumbFile, name);
  }
  
  /**
   * Creates a thumbnail image and writes it to a specified path.
   * @param image the full-size image from which to create the thumbnail 
   * @param path the path for the thumbnail image file
   * @return the thumbnail File, or null if failed
   */
  protected File createThumbnailFile(BufferedImage image, String path) {
  	// determine image resize factor
    double widthFactor = thumbnailDimension.getWidth()/image.getWidth();
    double heightFactor = thumbnailDimension.getHeight()/image.getHeight();
    double factor = Math.min(widthFactor, heightFactor);

    // determine dimensions of thumbnail image
    int w = (int)(image.getWidth()*factor);
    int h = (int)(image.getHeight()*factor);
  	
		// create and draw thumbnail image
    BufferedImage thumbnailImage = new BufferedImage(w, h, BufferedImage.TYPE_3BYTE_BGR);
		Graphics2D g = thumbnailImage.createGraphics();

    AffineTransform transform = AffineTransform.getScaleInstance(factor, factor);
    g.setTransform(transform);
    g.drawImage(image, 0, 0, null);
    
    // write thumbnail image to file
    return VideoIO.writeImageFile(thumbnailImage, path);
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
