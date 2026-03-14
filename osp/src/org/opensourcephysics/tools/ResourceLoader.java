/*
 * Open Source Physics software is free software as described near the bottom of
 * this code file.
 *
 * For additional information and documentation on Open Source Physics please
 * see: <http://www.opensourcephysics.org/>
 */

package org.opensourcephysics.tools;
import java.applet.AudioClip;
import java.awt.Image;
import java.awt.image.BufferedImage;
import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileFilter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.Reader;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLClassLoader;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

import javax.swing.ImageIcon;
import javax.swing.JOptionPane;

import org.opensourcephysics.controls.OSPLog;
import org.opensourcephysics.controls.XML;
import org.opensourcephysics.display.OSPRuntime;

/**
 * This defines static methods for loading resources.
 *
 * @author Douglas Brown
 * @version 1.0
 */
public class ResourceLoader {
	
  @SuppressWarnings("javadoc")
	public static final FileFilter OSP_CACHE_FILTER;
  @SuppressWarnings("javadoc")
	public static final String WIN_XP_DEFAULT_CACHE;
  @SuppressWarnings("javadoc")
	public static final String WINDOWS_DEFAULT_CACHE;
  @SuppressWarnings("javadoc")
	public static final String OSX_DEFAULT_CACHE;
  @SuppressWarnings("javadoc")
	public static final String LINUX_DEFAULT_CACHE;

  protected static ArrayList<String> searchPaths = new ArrayList<String>();                        // search paths
  protected static ArrayList<String> appletSearchPaths = new ArrayList<String>();                  // search paths for apples
  protected static int maxPaths = 20;                                                              // max number of paths in history
  protected static Hashtable<String, Resource> resources = new Hashtable<String, Resource>();      // cached resources
  protected static boolean cacheEnabled=false, canceled=false;
  protected static Map<String, URLClassLoader> zipLoaders = new TreeMap<String, URLClassLoader>(); // maps zip to zipLoader
  protected static URLClassLoader xsetZipLoader; // zipLoader of current xset
  protected static ArrayList<String> extractExtensions = new ArrayList<String>();
  protected static ArrayList<String> pathsNotFound = new ArrayList<String>();
  protected static File ospCache;
  protected static boolean zipURLsOK;
  
  static {
  	OSP_CACHE_FILTER = new FileFilter() {
  		public boolean accept(File file) {
  			return file.isDirectory() && file.getName().startsWith("osp-"); //$NON-NLS-1$
  		}
  	};
  	WIN_XP_DEFAULT_CACHE = "/Local Settings/Application Data/OSP/Cache"; //$NON-NLS-1$
  	WINDOWS_DEFAULT_CACHE = "/AppData/Local/OSP/Cache"; //$NON-NLS-1$
  	OSX_DEFAULT_CACHE = "/Library/Caches/OSP"; //$NON-NLS-1$
  	LINUX_DEFAULT_CACHE = "/Downloads/OSP/Cache"; //$NON-NLS-1$
  }

  /**
   * Private constructor to prevent instantiation.
   */
  private ResourceLoader() {
    /** empty block */
  }

  /**
   * Gets a resource specified by name. If no resource is found using the name
   * alone, the searchPaths are searched.
   *
   * @param name the file or URL name
   * @return the Resource, or null if none found
   */
  public static Resource getResource(String name) {
    return getResource(name, true);
  }

  /**
   * Gets a resource specified by name. If no resource is found using the name
   * alone, the searchPaths are searched. This will find a zip file as a URL
   * resource, unlike the getResource(String) method.
   *
   * @param name the file or URL name
   * @return the Resource, or null if none found
   */
  public static Resource getResourceZipURLsOK(String name) {
  	zipURLsOK = true;
    Resource res = getResource(name, true);
  	zipURLsOK = false;
    return res;
  }

  /**
   * Gets a resource specified by name. If no resource is found using
   * the name alone, the searchPaths are searched.
   * Files are searched only if searchFile is true.
   *
   * @param name the file or URL name
   * @param searchFiles true to search files
   * @return the Resource, or null if none found
   */
  public static Resource getResource(String name, boolean searchFiles) {
    try {
      URL url = getAppletResourceURL(name); // added by W. Christian
      if(url!=null) {
        return new Resource(url);
      }
    } catch(Exception ex) {}
    return getResource(name, Resource.class, searchFiles);
  }

  /**
   * Gets a resource specified by name and Class. If no resource is found using
   * the name alone, the searchPaths are searched.
   *
   * @param name the file or URL name
   * @param type the Class providing default ClassLoader resource loading
   * @return the Resource, or null if none found
   */
  public static Resource getResource(String name, Class<?> type) {
    return getResource(name, type, true);
  }

  /**
   * Gets a resource specified by name and Class. If no resource is found using
   * the name alone, the searchPaths are searched.
   * Files are searched only if searchFile is true.
   *
   * @param name the file or URL name
   * @param type the Class providing default ClassLoader resource loading
   * @param searchFiles true to search files
   * @return the Resource, or null if none found
   */
  public static Resource getResource(String name, Class<?> type, boolean searchFiles) {
    if((name==null)||name.equals("")) { //$NON-NLS-1$
      return null;
    }
    pathsNotFound.clear();
    // Remove leading and trailing inverted commas (added by Paco)
    if(name.startsWith("\"")) { //$NON-NLS-1$   
      name = name.substring(1);
    }
    if(name.endsWith("\"")) { //$NON-NLS-1$   
      name = name.substring(0, name.length()-1);
    }
    if(name.startsWith("./")) { //$NON-NLS-1$   
      name = name.substring(2);
    }
    if(OSPRuntime.isAppletMode()||(OSPRuntime.applet!=null)) { // added by Paco
      Resource appletRes = null;
      // following code added by Doug Brown 2009/11/14
      if(type==OSPRuntime.applet.getClass()) {
        try {
          URL url = type.getResource(name);
          appletRes = createResource(url);
          if(appletRes!=null) {
            return appletRes;
          }
        } catch(Exception ex) {}
      }  // end code added by Doug Brown 2009/11/14
      for(Iterator<String> it = searchPaths.iterator(); it.hasNext(); ) {
        String path = getPath(it.next(), name);
        appletRes = findResourceInClass(path, type, searchFiles);
        if(appletRes!=null) {
          return appletRes;
        }
      }
      appletRes = findResourceInClass(name, type, searchFiles);
      if(appletRes!=null) {
        return appletRes;
      }
    }
    // look for resource with name only
    Resource res = findResource(name, type, searchFiles);
    if(res!=null) {
      return res;
    }
    pathsNotFound.add(name);
    StringBuffer err = new StringBuffer("Not found: "+name); //$NON-NLS-1$
    err.append(" [searched "+name); //$NON-NLS-1$
    // look for resource in searchPaths
    for(String next: searchPaths) {
      String path = getPath(next, name);
    	if (pathsNotFound.contains(path))
    		continue;
      res = findResource(path, type, searchFiles);
      if(res!=null) {
        return res;
      }
      pathsNotFound.add(path);
      err.append(";"+path); //$NON-NLS-1$
    }
    err.append("]"); //$NON-NLS-1$
    OSPLog.fine(err.toString());
    return null;
  }

  /**
   * Gets a resource specified by base path and name. If base path is relative
   * and no resource is found using the base alone, the searchPaths are
   * searched.
   *
   * @param basePath the base path
   * @param name the file or URL name
   * @return the Resource, or null if none found
   */
  public static Resource getResource(String basePath, String name) {
    return getResource(basePath, name, Resource.class);
  }

  /**
   * Gets a resource specified by base path and name. If base path is relative
   * and no resource is found using the base alone, the searchPaths are
   * searched. Files are searched only if searchFile is true.
   *
   * @param basePath the base path
   * @param name the file or URL name
   * @param searchFiles true to search files
   * @return the Resource, or null if none found
   */
  public static Resource getResource(String basePath, String name, boolean searchFiles) {
    return getResource(basePath, name, Resource.class, searchFiles);
  }

  /**
   * Gets a resource specified by base path, name and class. If base path is
   * relative and no resource is found using the base alone, the searchPaths
   * are searched.
   *
   * @param basePath the base path
   * @param name the file or URL name
   * @param type the Class providing ClassLoader resource loading
   * @return the Resource, or null if none found
   */
  public static Resource getResource(String basePath, String name, Class<Resource> type) {
    return getResource(basePath, name, type, true);
  }

  /**
 * Gets a resource specified by base path, name and class. If base path is
 * relative and no resource is found using the base alone, the searchPaths
 * are searched. Files are searched only if searchFile is true.
 *
 * @param basePath the base path
 * @param name the file or URL name
 * @param type the Class providing ClassLoader resource loading
 * @param searchFiles true to search files
 * @return the Resource, or null if none found
 */
  public static Resource getResource(String basePath, String name, Class<Resource> type, boolean searchFiles) {
  	if(basePath==null) {
      return getResource(name, type);
    }
    if(name.startsWith("./")) { //$NON-NLS-1$   
      name = name.substring(2);
    }
    // look for resource with basePath and name
    pathsNotFound.clear();
    String path = getPath(basePath, name);
    Resource res = findResource(path, type, searchFiles);
    if(res!=null) {
      return res;
    }
    // keep looking only if base path is relative
    if(basePath.startsWith("/")||(basePath.indexOf(":/")>-1)) { //$NON-NLS-1$ //$NON-NLS-2$
      return null;
    }
    pathsNotFound.add(path);
    StringBuffer err = new StringBuffer("Not found: "+path); //$NON-NLS-1$
    err.append(" [searched "+path); //$NON-NLS-1$
    if(OSPRuntime.applet!=null) {                  // applet mode
      String docBase = OSPRuntime.applet.getDocumentBase().toExternalForm();
      docBase = XML.getDirectoryPath(docBase)+"/"; //$NON-NLS-1$
      path = getPath(getPath(docBase, basePath), name);
    	if (!pathsNotFound.contains(path)) {
	      res = findResource(path, type, searchFiles);
	      if(res!=null) {
	        return res;
	      }
	      pathsNotFound.add(path);
	      err.append(";"+path);                        //$NON-NLS-1$
    	}
      String codeBase = OSPRuntime.applet.getCodeBase().toExternalForm();
      if(!codeBase.equals(docBase)) {
        path = getPath(getPath(codeBase, basePath), name);
      	if (!pathsNotFound.contains(path)) {
	        res = findResource(path, type, searchFiles);
	        if(res!=null) {
	          return res;
	        }
	        pathsNotFound.add(path);
	        err.append(";"+path);                      //$NON-NLS-1$
      	}
      }
    }
    // look for resource in searchPaths
    for(Iterator<String> it = searchPaths.iterator(); it.hasNext(); ) {
      path = getPath(getPath(it.next(), basePath), name);
    	if (pathsNotFound.contains(path))
    		continue;
      res = findResource(path, type, searchFiles);
      if(res!=null) {
        return res;
      }
      pathsNotFound.add(path);
      err.append(";"+path); //$NON-NLS-1$
    }
    err.append("]"); //$NON-NLS-1$
    OSPLog.fine(err.toString());
    return null;
  }

  /**
   * Adds a path at the beginning of the searchPaths list.
   *
   * @param base the base path to add
   */
  public static void addSearchPath(String base) {
    if((base==null)||base.equals("")||(maxPaths<1)) { //$NON-NLS-1$
      return;
    }
    synchronized(searchPaths) {
      if(searchPaths.contains(base)) {
        searchPaths.remove(base);
      } else {
        OSPLog.fine("Added path: "+base);   //$NON-NLS-1$
      }
      searchPaths.add(0, base);
      while(searchPaths.size()>Math.max(maxPaths, 0)) {
        base = searchPaths.get(searchPaths.size()-1);
        OSPLog.fine("Removed path: "+base); //$NON-NLS-1$
        searchPaths.remove(base);
      }
    }
  }

  /**
   * Removes a path from the searchPaths list.
   *
   * @param base the base path to remove
   */
  public static void removeSearchPath(String base) {
    if((base==null)||base.equals("")) { //$NON-NLS-1$
      return;
    }
    synchronized(searchPaths) {
      if(searchPaths.contains(base)) {
        OSPLog.fine("Removed path: "+base); //$NON-NLS-1$
        searchPaths.remove(base);
      }
    }
  }

  /**
   * Adds a search path at the beginning of the applet's search path list.
   * Added by Wolfgang Christian.
   *
   * @param base the base path to add
   */
  public static void addAppletSearchPath(String base) {
    if((base==null)||(maxPaths<1)) {
      return;
    }
	base=base.trim();
	if(!base.endsWith("/"))base=base+"/";  //$NON-NLS-1$//$NON-NLS-2$
    synchronized(appletSearchPaths) {
      if(appletSearchPaths.contains(base)) {
        appletSearchPaths.remove(base);                 // search path will be added to top of list later
      } else {
        OSPLog.fine("Applet search path added: "+base); //$NON-NLS-1$
      }
      appletSearchPaths.add(0, base);
      while(appletSearchPaths.size()>Math.max(maxPaths, 0)) {
        base = appletSearchPaths.get(appletSearchPaths.size()-1);
        OSPLog.fine("Removed path: "+base);             //$NON-NLS-1$
        appletSearchPaths.remove(base);
      }
    }
  }

  /**
   * Removes a path from the applet search path list.
   * Added by Wolfgang Christian.
   *
   * @param base the base path to remove
   */
  public static void removeAppletSearchPath(String base) {
    if((base==null)||base.equals("")) { //$NON-NLS-1$
      return;
    }
    synchronized(appletSearchPaths) {
      if(appletSearchPaths.contains(base)) {
        OSPLog.fine("Applet search path removed: "+base); //$NON-NLS-1$
        appletSearchPaths.remove(base);
      }
    }
  }

  /**
   * Sets the cacheEnabled property.
   *
   * @param enabled true to enable the cache
   */
  public static void setCacheEnabled(boolean enabled) {
    cacheEnabled = enabled;
  }

  /**
   * Gets the cacheEnabled property.
   *
   * @return true if the cache is enabled
   */
  public static boolean isCacheEnabled() {
    return cacheEnabled;
  }

  /**
   * Adds an extension to the end of the extractExtensions list.
   * Files with this extension found inside jars are extracted before loading.
   *
   * @param extension the extension to add
   */
  public static void addExtractExtension(String extension) {
    if((extension==null)||extension.equals("")) { //$NON-NLS-1$
      return;
    }
    if(!extension.startsWith(".")) { //$NON-NLS-1$
      extension = "."+extension;     //$NON-NLS-1$
    }
    OSPLog.finest("Added extension: "+extension); //$NON-NLS-1$
    synchronized(extractExtensions) {
      extractExtensions.add(extension);
    }
  }

  /**
   * Cancels the current operation when true.
   *
   * @param cancel true to cancel
   */
  public static void setCanceled(boolean cancel) {
  	canceled = cancel;
  }
  
  /**
   * Determines if the current operation is canceled.
   *
   * @return true if canceled
   */
  public static boolean isCanceled() {
  	return canceled;
  }
  
  // ___________________________ convenience methods _________________________
  
  /**
   * Opens and returns an input stream. May return null.
   * 
   * @param path the path
   * @return the input stream
   */
  public static InputStream openInputStream(String path) {
    Resource res = getResource(path);
    return(res==null) ? null : res.openInputStream();
  }

  /**
   * Opens and returns a reader. May return null.
   * 
   * @param path the path
   * @return the reader
   */
  public static Reader openReader(String path) {
    Resource res = getResource(path);
    return(res==null) ? null : res.openReader();
  }

  /**
   * Gets a string. May return null.
   * 
   * @param path the path
   * @return the string
   */
  public static String getString(String path) {
    Resource res = getResource(path);
    return(res==null) ? null : res.getString();
  }

  /**
   * Gets an image icon. May return null.
   * 
   * @param path the path
   * @return the icon
   */
  public static ImageIcon getIcon(String path) {
    URL url = getAppletResourceURL(path); // added by W. Christian
    if(url!=null) {
      return new ImageIcon(url);
    }
    Resource res = getResource(path);
    return(res==null) ? null : res.getIcon();
  }

  /**
   * Gets an image. May return null.
   * 
   * @param path the path
   * @return the image
   */
  public static Image getImage(String path) {
    URL url = getAppletResourceURL(path); // added by W. Christian
    if(url!=null) {
      return new ImageIcon(url).getImage();
    }
    Resource res = getResource(path);
    return(res==null) ? null : res.getImage();
  }

  /**
   * Gets a buffered image. May return null.
   * 
   * @param path the path
   * @return the image
   */
  public static BufferedImage getBufferedImage(String path) {
    Resource res = getResource(path);
    return(res==null) ? null : res.getBufferedImage();
  }

  /**
  * Gets a buffered image. May return null.
  * 
  * @param path the path
   * @param bufferedImageType one of the types defined by the BufferedImage class
  * @return the image
  */
 public static BufferedImage getBufferedImage(String path, int bufferedImageType) {
   Resource res = getResource(path);
   return(res==null) ? null : res.getBufferedImage(bufferedImageType);
 }

  /**
   * Gets an audio clip. May return null.
   * 
   * @param path the path
   * @return the audio clip
   */
  public static AudioClip getAudioClip(String path) {
    Resource res = getResource(path);
    return(res==null) ? null : res.getAudioClip();
  }
  
  /**
   * Sets the directory for cached files.
   * 
   * @param cache the desired cache directory
   */
  public static void setOSPCache(File cache) {
    if (!cache.exists() || !cache.isDirectory()) {
    	cache.mkdirs();
    }
  	ospCache = cache;
  }

  /**
   * Gets the directory for cached files.
   * 
   * @return the OSP cache
   */
  public static File getOSPCache() {
  	return ospCache;
  }
  
  /**
   * Gets a subdirectory for cached files associated with a particular URL path.
   * 
   * @param urlPath the path to the file
   * @return the subdirectory file, or the temp directory if there is no cache
   */
  public static File getOSPCacheHostDirectory(String urlPath) { 
  	File cache = getOSPCache();
		String host = ""; //$NON-NLS-1$
		if (cache!=null) {
			try {
				URL url = new URL(urlPath);
				host = url.getHost().replace('.', '_');
			} catch (MalformedURLException e) {				
			}
			if (host.equals("")) //$NON-NLS-1$
				host = "local_machine"; //$NON-NLS-1$
			cache = new File(cache, "osp-"+host); //$NON-NLS-1$
		}
		else {
			cache = new File(System.getProperty("java.io.tmpdir")); //$NON-NLS-1$			
 		}
		return cache;
  }

  /**
   * Downloads a file from the web to the OSP Cache.
   * 
   * @param urlPath the path to the file
   * @param fileName the name to assign the downloaded file
   * @param subdirectory a subdirectory path within the cache (may be null)
   * @param alwaysOverwrite true to overwrite an existing file, if any
   * @return the downloaded file, or null if failed to download
   */
  public static File downloadToOSPCache(String urlPath, String fileName, 
  		String subdirectory, boolean alwaysOverwrite) {
		if (fileName==null) return null;
  	File cache = getOSPCacheHostDirectory(urlPath);
		if (subdirectory!=null) {
			cache = new File(cache, subdirectory);
		}
		File file = ResourceLoader.download(urlPath, fileName, cache, alwaysOverwrite);
		if (file==null) {
			JOptionPane.showMessageDialog(null, 
					ToolsRes.getString("ResourceLoader.Dialog.FailedToDownload.Message1") //$NON-NLS-1$ 
					+"\n"+ToolsRes.getString("ResourceLoader.Dialog.FailedToDownload.Message2") //$NON-NLS-1$ //$NON-NLS-2$ 
					+"\n"+ToolsRes.getString("ResourceLoader.Dialog.FailedToDownload.Message3"), //$NON-NLS-1$ //$NON-NLS-2$ 
					ToolsRes.getString("ResourceLoader.Dialog.FailedToDownload.Title"), //$NON-NLS-1$ 
					JOptionPane.ERROR_MESSAGE);
		}
		return file;
  }
  
  /**
   * Clears an OSP cache directory.
   * 
   * @param cache the cache to clear
   * @return true if successfully cleared
   */
  public static boolean clearOSPCache(File cache) {
  	if (cache==null || !cache.canWrite()) return true;
  	boolean success = true;
  	File[] files = cache.listFiles(OSP_CACHE_FILTER);
  	if (files==null) return true;
    for(File next: files) {
    	success = success && deleteFile(next);
    }
    if (!success) {
			JOptionPane.showMessageDialog(null, 
					ToolsRes.getString("ResourceLoader.Dialog.UnableToClearCache.Message1") //$NON-NLS-1$ 
					+"\n"+ToolsRes.getString("ResourceLoader.Dialog.UnableToClearCache.Message2"), //$NON-NLS-1$ //$NON-NLS-2$ 
					ToolsRes.getString("ResourceLoader.Dialog.UnableToClearCache.Title"), //$NON-NLS-1$ 
					JOptionPane.WARNING_MESSAGE);
    }
    return success;
  }

  /**
   * Deletes a file or folder. In case of a folder, deletes all contents
   * and the folder itself.
   * 
   * @param file the file to delete
   */
  private static boolean deleteFile(File file) {
    if(file.isDirectory()) {
      File[] files = file.listFiles();
      for(File next: files) {
        deleteFile(next);
      }
    }
    return file.delete();
  }
  
  /**
   * Gets the contents of a zip file.
   * 
   * @param zipPath the path to the zip file
   * @return a set of file names in alphabetical order
   */
  public static Set<String> getZipContents(String zipPath) {
    Set<String> fileNames = new TreeSet<String>();
    try {
    	URL url = new URL(getURIPath(zipPath));    	
    	OSPLog.finest("zip url: "+url.toExternalForm()); //$NON-NLS-1$
      BufferedInputStream bufIn = new BufferedInputStream(url.openStream());
      ZipInputStream input = new ZipInputStream(bufIn);
      ZipEntry zipEntry=null;
      while ((zipEntry=input.getNextEntry())!=null) {
      	OSPLog.finest("zip entry: "+zipEntry); //$NON-NLS-1$
        if (zipEntry.isDirectory()) continue;
        String fileName = zipEntry.getName();
        fileNames.add(fileName);
      }
      input.close();
    }
    catch (Exception ex) {}    
    return fileNames;
  }
  
  /**
   * Unzips a ZIP file into the given directory. ZIP file may be on a server.
   * Can be canceled using the static setCanceled(boolean) method.
   * Note this does not warn of possible overwrites.
   * 
   * @param zipPath the (url) path to the zip file
   * @param targetDir target directory to save the extracted files
   * @param alwaysOverwrite true to overwrite existing files, if any
   * @return the Set of extracted files
   */
  public static Set<File> unzip(String zipPath, File targetDir, boolean alwaysOverwrite) {
  	if (targetDir==null)
			targetDir = new File(System.getProperty("java.io.tmpdir")); //$NON-NLS-1$
  	OSPLog.finer("unzipping "+zipPath+" to "+targetDir); //$NON-NLS-1$ //$NON-NLS-2$
  	try {
    	URL url = new URL(getURIPath(zipPath));    	
      BufferedInputStream bufIn = new BufferedInputStream(url.openStream());
      ZipInputStream input = new ZipInputStream(bufIn);
      ZipEntry zipEntry=null;
      Set<File> fileSet = new HashSet<File>();
      byte[] buffer = new byte[1024];
      setCanceled(false);
      while ((zipEntry=input.getNextEntry()) != null) {
        if (zipEntry.isDirectory()) continue;
        if (isCanceled()) {
          input.close();
          return null;
        }
        String filename = zipEntry.getName();
        File file = new File(targetDir, filename);
        if (!alwaysOverwrite && file.exists()) {
          fileSet.add(file);
        	continue;
        }
        file.getParentFile().mkdirs();
        int bytesRead;
        FileOutputStream output = new FileOutputStream(file);
        while ((bytesRead=input.read(buffer)) != -1) 
        	output.write(buffer, 0, bytesRead);
        output.close();
        input.closeEntry();
        fileSet.add(file);
      }
      input.close();
      return fileSet;
    }
    catch (Exception ex) { 
      ex.printStackTrace();
      return null;
    }    
  }  
  
  /**
   * Downloads a file from the web.
   * 
   * @param urlPath the path to the file
   * @param fileName the name to assign the downloaded file
   * @param targetDir the directory in which the downloaded file will be written
   * @param alwaysOverwrite true to overwrite an existing file, if any
   * @return the downloaded file
   */
  public static File download(String urlPath, String fileName, File targetDir, boolean alwaysOverwrite) {
  	urlPath = getURIPath(urlPath);
		if (!targetDir.exists())
			targetDir.mkdirs();
  	String localPath = XML.getResolvedPath(fileName, targetDir.getAbsolutePath());
		localPath = XML.forwardSlash(localPath);
    File localFile = new File(localPath);
    if (alwaysOverwrite || !localFile.exists()) {
	  	OSPLog.finer("downloading "+urlPath+" to "+targetDir+" as "+fileName); //$NON-NLS-1$ //$NON-NLS-2$ //$NON-NLS-3$
  		try {
				URL url = new URL(urlPath);
				InputStream reader = url.openStream();
        FileOutputStream writer = new FileOutputStream(localPath);
        byte[] buffer = new byte[65536]; // 2^14 = 64K
        int bytesRead = 0;
        while ((bytesRead = reader.read(buffer)) > 0) { 
        	writer.write(buffer, 0, bytesRead);
        }  
        writer.close();
        reader.close();
			} catch (Exception ex) {
				ex.printStackTrace();
			} 
    }
    if (localFile.exists()) return localFile; 
    return null;
  }
  
  /**
   * Determines if a url path is available (ie both valid and connected).
   *
   * @param urlPath the path in URI form
   * @return true if available
   */
  public static boolean isURLAvailable(String urlPath) {
	  try {
      // make a URL, open a connection, get content
      URL url = new URL(urlPath);
      HttpURLConnection urlConnect = (HttpURLConnection)url.openConnection();
      urlConnect.getContent();	
	  } catch (Exception ex) {
		  return false;
	  }
	  return true;
  }

  /**
   * Removes protocol and "%20" from URI paths.
   *
   * @param uriPath the path in URI form
   * @return the path
   */
  public static String getNonURIPath(String uriPath) {
  	if (uriPath==null) return null;
  	String path = uriPath;
//		String path = XML.forwardSlash(uriPath.trim());
  	// remove file protocol, if any
    if (path.startsWith("file:")) {                     //$NON-NLS-1$
      path = path.substring(5);
    }
    // remove all but one leading slash
    while (path.startsWith("//")) { //$NON-NLS-1$
      path = path.substring(1);
    }
    // remove last leading slash if drive is specified
    if (path.startsWith("/") && path.indexOf(":")>-1) { //$NON-NLS-1$ //$NON-NLS-2$
      path = path.substring(1);
    }
    // replace "%20" with space
    int j = path.indexOf("%20");                       //$NON-NLS-1$
    while(j>-1) {
      String s = path.substring(0, j);
      path = s+" "+path.substring(j+3);                //$NON-NLS-1$
      j = path.indexOf("%20");                         //$NON-NLS-1$
    }
//    // replace "%26" with "&"
//    j = path.indexOf("%26");                           //$NON-NLS-1$
//    while(j>-1) {
//      String s = path.substring(0, j);
//      path = s+"&"+path.substring(j+3);                //$NON-NLS-1$
//      j = path.indexOf("%26");                         //$NON-NLS-1$
//    }
  	return path;
  }

  /**
   * Converts a path to URI form (spaces replaced by "%20", etc).
   *
   * @param path the path
   * @return the path in URI form
   */
  public static String getURIPath(String path) {
  	if (path==null) return null;
		// trim and change backslashes to forward slashes
		path = XML.forwardSlash(path.trim());
		// add forward slash at end if needed
		if (!path.equals("")  //$NON-NLS-1$
				&& XML.getExtension(path)==null
				&& !path.endsWith("/")) //$NON-NLS-1$
			path += "/"; //$NON-NLS-1$
    // replace spaces with "%20"
    int i = path.indexOf(" ");                       //$NON-NLS-1$
    while(i>-1) {
      String s = path.substring(0, i);
      path = s+"%20"+path.substring(i+1);            //$NON-NLS-1$
      i = path.indexOf(" ");                         //$NON-NLS-1$
    }
//    // replace "&" with "%26"
//    i = path.indexOf("&");                           //$NON-NLS-1$
//    while(i>-1) {
//      String s = path.substring(0, i);
//      path = s+"%26"+path.substring(i+1);            //$NON-NLS-1$
//      i = path.indexOf("&");                         //$NON-NLS-1$
//    }
    // add file protocol if path to local file
		if (!path.equals("") && !path.startsWith("http:") && !path.startsWith("file:/")) { //$NON-NLS-1$ //$NON-NLS-2$ //$NON-NLS-3$
			String protocol = OSPRuntime.isWindows()? "file:/": "file://"; //$NON-NLS-1$ //$NON-NLS-2$
			path = protocol+path;
		}
  	return path;
  }

  // ______________________________ private methods ___________________________

  /**
   * Gets the resource URL using the applet's class loader.
   * Added by Wolfgang Christian.
   *
   * @param name of the resource
   * @return URL of the Resource, or null if none found
   */
  private static URL getAppletResourceURL(String name) {
    if((OSPRuntime.applet==null)||(name==null)||name.trim().equals("")) { //$NON-NLS-1$
      return null;
    }
    if(name.startsWith("http:")||name.startsWith("https:")){ //$NON-NLS-1$  //$NON-NLS-2$ // open a direct connection for http and https resources
    	try {
			return new java.net.URL(name);
		} catch (MalformedURLException e) {
			//e.printStackTrace();
		} 
    }
    name=name.trim();  // remove whitespace
    if(!name.startsWith("/")) { //$NON-NLS-1$ // try applet search paths for relative paths
      for(Iterator<String> it = appletSearchPaths.iterator(); it.hasNext(); ) {
        String path = it.next();
        String tempName=name;                   // tempName may change
    	if(tempName.startsWith("../")) {        //$NON-NLS-1$   
    		  tempName = tempName.substring(3);     //remove prefix
              path=path.substring(0, path.length()-1); // drop trailing slash
              int last=path.lastIndexOf("/"); //$NON-NLS-1$ // find last directory slash
              path=(last>0)?path.substring(0, last):"/";   //$NON-NLS-1$ // drop last directory if it exists 
        }else if(tempName.startsWith("./")) {         //$NON-NLS-1$   
        	tempName = tempName.substring(2);     //remove reference to current directory
        } 
        URL url = OSPRuntime.applet.getClass().getResource(path+tempName);
        if(url!=null) {
          return url;
        }
      }
    }
    return OSPRuntime.applet.getClass().getResource(name); // url not found in applet search paths
  }

  /**
   * Creates a Resource from a file.
   *
   * @param path the file path
   * @return the resource, if any
   */
  static private Resource createFileResource(String path) {
      // don't create file resources when in applet mode
    if(OSPRuntime.applet!=null) {
      return null;
    }
    // ignore paths that refer to zip or jar files
    if((path.indexOf(".zip")>-1)||(path.indexOf(".jar")>-1)) { //$NON-NLS-1$ //$NON-NLS-2$
      return null;
    }
    File file = new File(path);
    try {
      if(file.exists()&&file.canRead()) {
        Resource res = new Resource(file);
        if(path.endsWith("xset")) {                                    //$NON-NLS-1$
          xsetZipLoader = null;
        }
        OSPLog.finer("File: "+XML.forwardSlash(res.getAbsolutePath())); //$NON-NLS-1$
        return res;
      }
    } catch(Exception ex) {
      /** empty block */
    }
    return null;
  }

  /**
   * Creates a Resource from a URL.
   *
   * @param path the url path
   * @return the resource, if any
   */
  static private Resource createURLResource(String path) {
    // ignore paths that refer to zip or jar files unless explicitly OK
    if(!zipURLsOK && 
    		(path.indexOf(".zip")>-1 || path.indexOf(".jar")>-1)) { //$NON-NLS-1$ //$NON-NLS-2$
      return null;
    }
    Resource res = null;
    // following added by Doug Brown 2009/11/14
    if(OSPRuntime.applet!=null) {
      try { // let applet class try to get it first
        //URL url = OSPRuntime.applet.getClass().getResource(path);
        URL url = getAppletResourceURL(path);
        res = createResource(url);
      } catch(Exception ex) {
        /** empty block */
      }
    } // end code added by Doug Brown 2009/11/14
    if(res==null) {
      // if path includes protocol, use it directly
      if(path.indexOf(":/")>-1) {  //$NON-NLS-1$
        try {
//          URL url = new URL(path); // changed to use URI path 2011/09/11 DB
          URL url = new URL(getURIPath(path));
          res = createResource(url);
        } catch(Exception ex) {
          /** empty block */
        }
      }
      // else if applet mode and relative path, search document and code base
      else {
        if((OSPRuntime.applet!=null)&&!path.startsWith("/")) {             //$NON-NLS-1$
          // first check document base
          URL docBase = OSPRuntime.applet.getDocumentBase();
          try {
            // following added by Doug Brown 2009/11/14
            String basePath = docBase.toString();
            // strip query, if any, from document base
            int n = basePath.indexOf("?");                                 //$NON-NLS-1$
            if(n>-1) {
              docBase = new URL(basePath.substring(0, n));
            }
            // end code added by Doug Brown 2009/11/14
            URL url = new URL(docBase, path);
            res = createResource(url);
          } catch(Exception ex) {
            /** empty block */
          }
          if(res==null) {
            URL codeBase = OSPRuntime.applet.getCodeBase();
            String s = XML.getDirectoryPath(docBase.toExternalForm())+"/"; //$NON-NLS-1$
            if(!codeBase.toExternalForm().equals(s)) {
              try {
                URL url = new URL(codeBase, path);
                res = createResource(url);
              } catch(Exception ex) {
                /** empty block */
              }
            }
          }
        }
      }
    }
    if(res!=null) {
      if(path.endsWith(".xset")) {                                  //$NON-NLS-1$
        xsetZipLoader = null;
      }
      OSPLog.finer("URL: "+XML.forwardSlash(res.getAbsolutePath())); //$NON-NLS-1$
    }
    return res;
  }

  /**
   * Creates a Resource from within a zip or jar file.
   *
   * @param path the file path
   * @return the resource, if any
   */
  @SuppressWarnings("resource")
static private Resource createZipResource(String path) {
  	// convert to non-URI form
  	path = getNonURIPath(path);

    // get separate zip base and relative file name
    String base = null;
    String fileName = path;
    // look for zip or jar base path
    int i = path.indexOf("zip!/"); //$NON-NLS-1$
    if(i==-1) {
      i = path.indexOf("jar!/"); //$NON-NLS-1$
    }
    if(i==-1) {
      i = path.indexOf("exe!/"); //$NON-NLS-1$
    }
    if(i>-1) {
      base = path.substring(0, i+3);
      fileName = path.substring(i+5);
    }
    if(base==null) {
      if(path.endsWith(".zip")                           //$NON-NLS-1$
        ||path.endsWith(".jar")                          //$NON-NLS-1$
        ||path.endsWith(".exe")) {                       //$NON-NLS-1$
        String name = XML.stripExtension(XML.getName(path));
        base = path;
        fileName = name+".xset";                         //$NON-NLS-1$
      } else if(path.endsWith(".xset")) {                //$NON-NLS-1$
        base = path.substring(0, path.length()-4)+"zip"; //$NON-NLS-1$
      }
    }
    
    // if loading from a web file, download to OSP cache
  	boolean isZip = base!=null && 
  			(base.endsWith(".zip")||base.endsWith(".jar")); //$NON-NLS-1$ //$NON-NLS-2$
  	boolean deleteOnExit = ospCache==null;
  	if (isZip && path.startsWith("http:")) {             //$NON-NLS-1$
  		String zipFileName = XML.getName(base);
      File zipFile = downloadToOSPCache(base, zipFileName, null, false);
      if (zipFile!=null) {
      	if (deleteOnExit)
      		zipFile.deleteOnExit();
      	base = zipFile.getAbsolutePath();
      	path = base+"!/"+fileName; //$NON-NLS-1$
      }
  		
  		
//  		
//  		
//  		File zipDir = null;
//	  	String zipName = XML.getName(base);
//  		if (ospCache!=null) {
//				try {
//					URL url = new URL(getURIPath(path));
//					String host = url.getHost().replace('.', '_')+"/"; //$NON-NLS-1$
//					host += zipName.replace('.', '_') + "/"; //$NON-NLS-1$
//					zipDir = new File(ospCache, host);
//				} catch (MalformedURLException e) {}			
//  		}
//  		if (zipDir==null) {
//  			zipDir = new File(System.getProperty("java.io.tmpdir")); //$NON-NLS-1$			
//      	deleteOnExit = true;
//  		}
//      File zipFile = download(base, zipName, zipDir, false);
//      if (zipFile!=null) {
//      	if (deleteOnExit)
//      		zipFile.deleteOnExit();
//      	base = zipFile.getAbsolutePath();
//      	path = base+"!/"+fileName; //$NON-NLS-1$
//      }
    }
  	 	
    URLClassLoader zipLoader = null;
    URL url = null;
    if(base!=null) {
      // use existing zip loader, if any
      zipLoader = zipLoaders.get(base);
      if(zipLoader!=null) {
        url = zipLoader.findResource(fileName);
      } else {
        try {
          // create new zip loader
          URL[] urls = new URL[] {new URL("file", null, base)};  //$NON-NLS-1$
          zipLoader = new URLClassLoader(urls);
          url = zipLoader.findResource(fileName);
          if(url==null) {                                        // workaround works in IE?
            URL classURL = Resource.class.getResource("/"+base); //$NON-NLS-1$
            if(classURL!=null) {
              urls = new URL[] {classURL};
              zipLoader = new URLClassLoader(urls);
              url = zipLoader.findResource(fileName);
            }
          }
          if(url!=null) {
            zipLoaders.put(base, zipLoader);
          }
        } catch(Exception ex) {
          /** empty block */
        }
      }
    }
    // if not found, use xset zip loader, if any
    if((url==null)&&(xsetZipLoader!=null)) {
      url = xsetZipLoader.findResource(fileName);
      if(url!=null) {
        Iterator<String> it = zipLoaders.keySet().iterator();
        while(it.hasNext()) {
          Object key = it.next();
          if(zipLoaders.get(key)==xsetZipLoader) {
            base = (String) key;
            break;
          }
        }
      }
    }
    String launchJarPath = OSPRuntime.getLaunchJarPath();
    // if still not found, use launch jar loader, if any
    if((url==null)&&(launchJarPath!=null)) {
      zipLoader = zipLoaders.get(launchJarPath);
      if(zipLoader!=null) {
        url = zipLoader.findResource(fileName);
      } else {
        try {
          // create new zip loader
          URL[] urls = new URL[] {new URL("file", null, launchJarPath)};  //$NON-NLS-1$
          zipLoader = new URLClassLoader(urls);
          url = zipLoader.findResource(fileName);
          if(url==null) {                                                 // workaround works in IE?
            URL classURL = Resource.class.getResource("/"+launchJarPath); //$NON-NLS-1$
            if(classURL!=null) {
              urls = new URL[] {classURL};
              zipLoader = new URLClassLoader(urls);
              url = zipLoader.findResource(fileName);
            }
          }
          if(url!=null) {
            zipLoaders.put(launchJarPath, zipLoader);
          }
        } catch(Exception ex) {
          /** empty block */
        }
      }
      if(url!=null) {
        base = launchJarPath;
      }
    }
    if(url!=null) {     // successfully found url
      // extract file if extension is flagged for extraction
      Iterator<String> it = extractExtensions.iterator();
      while(it.hasNext()) {
        String ext = it.next();
        if(url.getFile().endsWith(ext)) {
          File zipFile = new File(base);
        	String targetPath = fileName;
        	String parent = zipFile.getParent();
        	// if target path is relative, resolve wrt parent folder of zip file
        	if(parent!=null && !targetPath.startsWith("/") //$NON-NLS-1$
        			&& fileName.indexOf(":/")==-1) { //$NON-NLS-1$
        		targetPath = XML.getResolvedPath(fileName, parent);
        	}
          File target = new File(targetPath);
          if(!target.exists()) {
            target = JarTool.extract(zipFile, fileName, targetPath);
            if (deleteOnExit)
            	target.deleteOnExit();
          }
          return createFileResource(target.getAbsolutePath());
        }
      }
      try {
        Resource res = createResource(url);
        if((res==null)||(res.getAbsolutePath().indexOf(path)==-1)) {
          return null;
        }
        if(fileName.endsWith("xset")) {                               //$NON-NLS-1$
          xsetZipLoader = zipLoader;
        }
        OSPLog.finer("Zip: "+XML.forwardSlash(res.getAbsolutePath())); //$NON-NLS-1$
        return res;
      } catch(IOException ex) {
        /** empty block */
      }
    }
    return null;
  }

  /**
   * Creates a Resource from a class resource, typically in a jar file.
   *
   * @param name the resource name
   * @param type the class providing the classloader
   * @return the resource, if any
   */
  static private Resource createClassResource(String name, Class<?> type) {
    // ignore any name that has a protocol
    if(name.indexOf(":/")!=-1) { //$NON-NLS-1$
      return null;
    }
    String fullName = name;
    int i = name.indexOf("jar!/"); //$NON-NLS-1$
    if(i==-1) {
      i = name.indexOf("exe!/"); //$NON-NLS-1$
    }
    if(i!=-1) {
      name = name.substring(i+5);
    }
    Resource res = null;
    try {                                   // check relative to root of jarfile containing specified class
      URL url = type.getResource("/"+name); //$NON-NLS-1$
      res = createResource(url);
    } catch(Exception ex) {
      /** empty block */
    }
    if(res==null) {
      try { // check relative to specified class
        URL url = type.getResource(name);
        res = createResource(url);
      } catch(Exception ex) {
        /** empty block */
      }
    }
    // if resource is found, log and set launchJarName if not yet set
    if(res!=null) {
      String path = XML.forwardSlash(res.getAbsolutePath());
      // don't return resources from Java runtime system jars
      if((path.indexOf("/jre")>-1)&&(path.indexOf("/lib")>-1)) { //$NON-NLS-1$ //$NON-NLS-2$
        return null;
      }
      // don't return resources that don't contain original name
      if(path.indexOf(fullName)==-1) {
        return null;
      }
      if(name.endsWith("xset")) {                                //$NON-NLS-1$
        xsetZipLoader = null;
      }
      OSPLog.finer("Class resource: "+path);                      //$NON-NLS-1$
      OSPRuntime.setLaunchJarPath(path);
    }
    return res; // may be null
  }

  /**
   * Creates a Resource.
   *
   * @param url the URL
   * @return the resource, if any
   * @throws IOException
   */
  static private Resource createResource(URL url) throws IOException {
    if(url==null) {
      return null;
    }
    // check that url is accessible
    InputStream stream = url.openStream();
    if(stream.read()==-1) {
      return null;
    }
    stream.close();
    return new Resource(url);
  }

  /**
   * Finds the resource using only the class resource loader
   */
  private static Resource findResourceInClass(String path, Class<?> type, boolean searchFiles) { // added by Paco
    path = path.replaceAll("/\\./", "/"); // This eliminates any embedded /./ //$NON-NLS-1$ //$NON-NLS-2$
    if(type==null) {
      type = Resource.class;
    }
    Resource res = null;
    // look for cached resource
    if(cacheEnabled) {
      res = resources.get(path);
      if((res!=null)&&(searchFiles||(res.getFile()==null))) {
        OSPLog.finest("Found in cache: "+path); //$NON-NLS-1$
        return res;
      }
    }
    if((res = createClassResource(path, type))!=null) {
      if(cacheEnabled) {
        resources.put(path, res);
      }
      return res;
    }
    return null;
  }

  private static Resource findResource(String path, Class<?> type, boolean searchFiles) {
  	path = path.replaceAll("/\\./", "/"); // This eliminates any embedded /./ //$NON-NLS-1$ //$NON-NLS-2$
    if(type==null) {
      type = Resource.class;
    }
    Resource res = null;
    // look for cached resource
    if(cacheEnabled) {
      res = resources.get(path);
      if((res!=null)&&(searchFiles||(res.getFile()==null))) {
        OSPLog.finest("Found in cache: "+path); //$NON-NLS-1$
        return res;
      }
    }
    // try to load resource in file/url/zip/class order
    // search files only if flagged
    if((searchFiles&&(res = createFileResource(path))!=null)
    		||(res = createURLResource(path))!=null
    		||(res = createZipResource(path))!=null
    		||(res = createClassResource(path, type))!=null) {
      if(cacheEnabled) {
        resources.put(path, res);
      }
      return res;
    }
    return null;
  }

  /**
   * Gets a path from a base path and file name.
   *
   * @param base the base path
   * @param name the file name
   * @return the path
   */
  private static String getPath(String base, String name) {
    if(base==null) {
      base = ""; //$NON-NLS-1$
    }
    if(base.endsWith(".jar")||base.endsWith(".zip")) { //$NON-NLS-1$ //$NON-NLS-2$
      base += "!";                                     //$NON-NLS-1$
    }
    String path = XML.getResolvedPath(name, base);
    // correct the path so that it works with Mac
    if(OSPRuntime.isMac()&&path.startsWith("file:/")&&!path.startsWith("file:///")) { //$NON-NLS-1$ //$NON-NLS-2$
      path = path.substring(6);
      while(path.startsWith("/")) {                                                   //$NON-NLS-1$
        path = path.substring(1);
      }
      path = "file:///"+path;                                                         //$NON-NLS-1$
    }
    return path;
  }

}

/*
 * Open Source Physics software is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License (GPL) as
 * published by the Free Software Foundation; either version 2 of the License,
 * or(at your option) any later version.
 *
 * Code that uses any portion of the code in the org.opensourcephysics package
 * or any subpackage (subdirectory) of this package must must also be be
 * released under the GNU GPL license.
 *
 * This software is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this; if not, write to the Free Software Foundation, Inc., 59 Temple Place,
 * Suite 330, Boston MA 02111-1307 USA or view the license online at
 * http://www.gnu.org/copyleft/gpl.html
 *
 * Copyright (c) 2007 The Open Source Physics project
 * http://www.opensourcephysics.org
 */
