package org.opensourcephysics.tools;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

import org.opensourcephysics.controls.XML;
import org.opensourcephysics.display.OSPRuntime;

	/**
	 * ExtensionsManager manages Java extensions directories.
	 * Its primary use is to copy Xuggle jars into appropriate ext directories.
	 * 
	 * @author Douglas Brown
	 * @version 1.0
	 */
public class ExtensionsManager {

	static String xuggleHome;
	static ExtensionsFilter extFilter;
	
	static {
		xuggleHome = System.getenv("XUGGLE_HOME"); //$NON-NLS-1$
		extFilter = new ExtensionsFilter();		
	}
	
	/**
	 * Main method when used as a stand-alone application.
	 * @param args ignored
	 */
	public static void main(String[] args) {
		// "print" extensions to stdout for use by Tracker and Xuggle installers
		ExtensionsManager manager = new ExtensionsManager();
		manager.printExtensionDirectoriesForXuggle(" "); //$NON-NLS-1$		
	}
	
  /**
   * Finds extension directories for Xuggle and prints a delimited list to System.out. 
   * A single space delimiter is parsable by Bitrock installers.
   * 
   * @param separator the delimiter between items in the list
   */
	public void printExtensionDirectoriesForXuggle(String separator) {
		Collection<File> extDirs = findJavaExtensionDirectoriesForXuggle();
		StringBuffer buf = new StringBuffer();
		for (File next: extDirs) {
			String fileName = XML.forwardSlash(next.getPath());
			buf.append("\""+fileName+"\""+separator); //$NON-NLS-1$ //$NON-NLS-2$
		}
		String s = buf.toString();
		// remove last separator
		if (s.length()>=separator.length()) {
			s = s.substring(0, s.length()-separator.length());
		}
		System.out.print(s);
	}

	
  /**
   * Copies Xuggle jar files to a target directory. Does not overwrite newer files
   * unless overwrite flag is true.
   *
   * @param dir the directory
   * @param overwrite true to always overwrite
   * @return true if copied or existing files are up to date
   */
	public boolean copyXuggleJarsTo(File dir, boolean overwrite) {
		if (xuggleHome==null || dir==null) {
			return false;
		}
	  if (!new File(xuggleHome+"/share/java/jars/xuggle-xuggler.jar").exists()) { //$NON-NLS-1$
	  	return false;
	  }
    File xuggleJarDir = new File(xuggleHome+"/share/java/jars"); //$NON-NLS-1$
    String[] jarNames = DiagnosticsForXuggle.xuggleJarNames;
	  File xuggleFile = new File(xuggleJarDir, jarNames[0]);
    long modified = xuggleFile.lastModified();
    File extFile = new File(dir, jarNames[0]);
    // copy xuggle jars
    if (overwrite || !extFile.exists() || extFile.lastModified()!=modified) {
//    	if (extFile.exists()) extFile.delete();
	    for (String next: jarNames) {
	      xuggleFile = new File(xuggleJarDir, next);
	      extFile = new File(dir, next);
	      if (!copyFile(xuggleFile, extFile)) {
	      	return false;
	      }
	    }
	    
      return true;
    }
		return true;
	}
	
  /**
   * Removes Xuggle files from a target directory.
   *
   * @param extensionDir the directory
   * @return true if removed
   */
	public boolean removeXuggleFrom(File extensionDir) {
    String[] jarNames = {"xuggle-xuggler.jar","slf4j-api.jar", //$NON-NLS-1$ //$NON-NLS-2$
    		"logback-core.jar","logback-classic.jar"}; //$NON-NLS-1$ //$NON-NLS-2$
    for (String next: jarNames) {
    	File extFile = new File(extensionDir, next);
      if (extFile.exists() && !extFile.delete()) {
      	return false;
      }
    }
    return true;
	}
	
  /**
   * Copies QTJava.zip to a target directory. Does not overwrite newer files
   * unless overwrite flag is true.
   *
   * @param extensionDir the directory
   * @param overwrite true to always overwrite
   * @return true if copied
   */
	public boolean copyQTJavaTo(File extensionDir, boolean overwrite) {
  	String qtJarName = "QTJava.zip"; //$NON-NLS-1$
    String[] folderNames = {
    		"C:/Program Files/QuickTime/QTSystem/", //$NON-NLS-1$
    		"C:/windows/system32/", //$NON-NLS-1$
    		"C:/windows/system/", //$NON-NLS-1$
    		"C:/winNT/system32/", //$NON-NLS-1$
    		"system/library/java/extensions/"}; //$NON-NLS-1$
    // look for most recent QTJava.zip in system folders
    long modified = 0;
    File qtSource = null; // file to be copied
    for (String next: folderNames) {
      File qtFile = new File(next+qtJarName);
    	if (!qtFile.exists()) continue;
    	long date = qtFile.lastModified();
    	if (date>modified) {
    		modified = Math.max(modified, date);
    		qtSource = qtFile;
    	}     
    }
	  if (qtSource==null)	return false;
    File extFile = new File(extensionDir, qtJarName);
    // copy xuggle jars
    if (overwrite || !extFile.exists() || extFile.lastModified()<modified) {
	    if (!copyFile(qtSource, extFile))	return false;
      return true;
    }
		return false;
	}
	
  /**
   * Removes QTJava from a target directory.
   *
   * @param extensionDir the directory
   * @return true if removed
   */
	public boolean removeQTJavaFrom(File extensionDir) {
		File extFile = new File(extensionDir, "QTJava.zip"); //$NON-NLS-1$
    if (extFile.exists() && !extFile.delete())
      return false;
    return true;
	}
	
  /**
   * Finds all java extension directories on the current machine.
   *
   * @return a collection of java extension directory files
   */
	public Collection<File> findAllJavaExtensionDirectories() {
    // set of all Java extensions directories on this machine
    Set<File> extDirs = new HashSet<File>();
    // set of extension directories used by the running Java VM 
    Set<String> vmExtDirs = new HashSet<String>();
    // set of "Java level" directories to search
    Set<File> searchPaths = new HashSet<File>();

		// get and parse system extension directories property into vmExtDirs
    String paths = XML.forwardSlash(System.getProperty("java.ext.dirs")); //$NON-NLS-1$   
    String separator = System.getProperty("path.separator"); //$NON-NLS-1$
    int n = paths.indexOf(separator);
    while (n>-1) {
    	vmExtDirs.add(paths.substring(0, n));
    	paths = paths.substring(n+1);
      n = paths.indexOf(separator);
    }
    if (!"".equals(paths)) {//$NON-NLS-1$
    	vmExtDirs.add(paths);
    }
       
    for (String next: vmExtDirs) {
    	File dir = new File(next);
    	if (!dir.exists()) continue;
  		extDirs.add(dir);
    	
      if (OSPRuntime.isMac()) {
  			// search path: /JavaVirtualMachines
  			while (dir.getPath().indexOf("/JavaVirtualMachines")>-1) { //$NON-NLS-1$
  				if (dir.getName().equals("JavaVirtualMachines")) { //$NON-NLS-1$
  					searchPaths.add(dir);
  					break;
  				}
  				dir = dir.getParentFile();
  			}
      }
      else if (OSPRuntime.isLinux()) {
    	  // search path: /jvm
  			while (dir.getPath().indexOf("/jvm")>-1) { //$NON-NLS-1$
  				if (dir.getName().equals("jvm")) { //$NON-NLS-1$
  					searchPaths.add(dir);
  					break;
  				}
  				dir = dir.getParentFile();
  			}
      }
    }
    
  	if (OSPRuntime.isWindows()) {
	    String progfiles = System.getenv("ProgramFiles"); //$NON-NLS-1$
			String w6432 = System.getenv("ProgramW6432"); //$NON-NLS-1$
			String x86 = System.getenv("ProgramFiles(x86)"); //$NON-NLS-1$
			// add Program Files (may or may not be x86) Java directory to search path if it exists
			if (progfiles!=null) {
				File file = new File(progfiles, "Java"); //$NON-NLS-1$
    		if (file.exists()) searchPaths.add(file);   				
			}    		
			// add "Program Files" Java directory to search path if it exists
			if (w6432!=null) { // 64-bit Windows
				File file = new File(w6432, "Java"); //$NON-NLS-1$
    		if (file.exists()) searchPaths.add(file);   				
			}    		
			// add "Program Files (x86)" Java directory to search path if it exists
			if (x86!=null) { // 64-bit Windows
				// add x86 Java directory to search path if it exists
				File file = new File(x86, "Java"); //$NON-NLS-1$
    		if (file.exists()) searchPaths.add(file);   				
			}
  	}
    
    // search all searchPaths and add all extensions directories found
    for (File next: searchPaths) {
			findJavaExtensionDirectories(next, extDirs);
    }    
    return extDirs;
	}
	
  /**
   * Finds java extension directories for Xuggle jars. Limits results to Java 1.6 or greater.
   *
   * @return a collection of java extension directory files
   */
	public Collection<File> findJavaExtensionDirectoriesForXuggle() {
    // set of all Java extension directories on this machine
    Collection<File> extDirs = findAllJavaExtensionDirectories();
    // set of extension directories to remove 
    Set<File> toExclude = new HashSet<File>();
    for (File next: extDirs) {
    	File javaFile = next.getParentFile().getParentFile();
    	int n = javaFile.getPath().indexOf("jdk"); //$NON-NLS-1$
    	n = Math.max(n, javaFile.getPath().indexOf("jre")); //$NON-NLS-1$   	
    	int oldVersion = javaFile.getPath().indexOf("1.5."); //$NON-NLS-1$
    	oldVersion = Math.max(oldVersion, javaFile.getPath().indexOf("-5-")); //$NON-NLS-1$   	
    	oldVersion = Math.max(oldVersion, javaFile.getPath().indexOf("1.4.")); //$NON-NLS-1$   	
    	oldVersion = Math.max(oldVersion, javaFile.getPath().indexOf("1.3.")); //$NON-NLS-1$   	
    	oldVersion = Math.max(oldVersion, javaFile.getPath().indexOf("1.2.")); //$NON-NLS-1$   	
    	oldVersion = Math.max(oldVersion, javaFile.getPath().indexOf("1.3.")); //$NON-NLS-1$   	
    	oldVersion = Math.max(oldVersion, javaFile.getPath().indexOf("1.2.")); //$NON-NLS-1$   	
  		if (n>-1 && oldVersion>-1) {
  			toExclude.add(next);
  		}
    }
    extDirs.removeAll(toExclude);
    return extDirs;
	}
	
  /**
   * Finds all java extension directory files (recursively) in a directory.
   * Extension directories are added to a collection of Files and returned.
   * Does not search symbolic links.
   *
   * @param dir the directory
   * @param extDirs the collection to add to
   * @return the collection
   */
	public Collection<File> findJavaExtensionDirectories(File dir, Set<File> extDirs) {
		try { // don't search symlinks
			if (!dir.getCanonicalPath().equals(dir.getAbsolutePath()))
				return extDirs;
		} catch (IOException e) {
		}
			
		// search all children contained in the directory
		String[] fileNames = dir.list();
		if (fileNames!=null && fileNames.length>0) {
			for (String next: fileNames) {
				File subDir = new File(dir, next);
				// if subdirectory is an extensions folder, add it
				if (extFilter.accept(subDir, subDir.getName())) {
					extDirs.add(subDir);
				}
				// else search the next level down
				else {
					findJavaExtensionDirectories(subDir, extDirs);
				}
			}
		}
		return extDirs;
	}
	
  /**
   * Finds public jres in which Xuggle should work.
   * @param requires32Bit true if Xuggle requires a 32-bit Java VM (eg, Xuggle 3.4 on Windows)
   * @return a collection of java jre files
   */
	public Collection<File> findPublicJREsForXuggle(boolean requires32Bit) {
    Collection<File> extDirs = findJavaExtensionDirectoriesForXuggle();
    Set<File> jreDirs = new HashSet<File>();
		String x86 = System.getenv("ProgramFiles(x86)"); //$NON-NLS-1$
    // iterate through extDirs to fill jreDirs
    for (File next: extDirs) {
    	// move up two levels from lib/ext
    	File javaFile = next.getParentFile().getParentFile();
    	if (OSPRuntime.isWindows()) {
    		// if 32-bit required, eliminate 64-bit
    		if (requires32Bit && x86!=null && next.getPath().indexOf(x86)==-1) {
    			continue;
    		}

    		// eliminate non-public jre
    		if (next.getPath().indexOf("jdk")>-1) continue; //$NON-NLS-1$
    		javaFile = new File(javaFile, "bin/java.exe"); //$NON-NLS-1$
    		if (!javaFile.exists()) continue;
    		String jrePath = OSPRuntime.getJREPath(javaFile);
    		jreDirs.add(new File(jrePath));
    	}
      else {
    		javaFile = new File(javaFile, "bin/java"); //$NON-NLS-1$
    		if (!javaFile.exists()) continue;
    		String jrePath = OSPRuntime.getJREPath(javaFile);
    		jreDirs.add(new File(jrePath));
      }
    }
    return jreDirs;
	}
	
  /**
   * Copies a source file to a target file.
   *
   * @param inFile the source
   * @param outFile the target
   * @return true if successfully copied
   */
  public boolean copyFile(File inFile, File outFile) {
  	byte[] buffer = new byte[100000];
    try {
    	InputStream in = new FileInputStream(inFile);
    	OutputStream out = new FileOutputStream(outFile);
			while (true) {
				synchronized (buffer) {
					int amountRead = in.read(buffer);
					if (amountRead == -1) {
						break;
					}
					out.write(buffer, 0, amountRead);
				}
			}
			in.close();
			out.close();
			// following line sometimes fails on Windows 7??
			outFile.setLastModified(inFile.lastModified());
		}                   
    catch (IOException ex) {
    	return false;
    }
  	return true;
  }

	/**
	 * ExtensionsFilter identifies Java extensions directories.
	 * 
	 * Windows:
	 * 		typical jdk: Program Files\Java\jdkX.X.X_XX\jre\lib\ext
   *    typical jre: Program Files\Java\jreX.X.X_XX\lib\ext
   *             and Program Files\Java\jreX\lib\ext
   *		on 64-bit Windows, 32-bit VMs are in: Program Files (x86)\Java\...
   *		non-jre: Sun\Java\lib\ext.
   *		exclude: \Program Files (x86)\Java\Java3D\x.x.x\lib\ext
   *		jre search in: \Java
	 *
	 * OS X:
	 * 		typical: /System/Library/Java/JavaVirtualMachines/X.X.X.jdk/Contents/Home/lib/ext
	 *    non-jre: /Library/Java/Extensions and/or /System/Library/Java/Extensions
	 * 		jre search in: /JavaVirtualMachines
	 * 
	 * Linux:
	 * 		typical: /usr/lib/jvm/java-X-openjdk/jre/lib/ext
	 * 		      or /usr/lib/jvm/java-X.X.X-openjdk/jre/lib/ext
	 * 		      or /usr/lib/jvm/java-X-sun-X.X.X.XX/jre/lib/ext
	 * 		      or /usr/lib/jvm/java-X.X.X-sun/jre/lib/ext
	 * 		non-jre: /usr/java/packages/lib/ext
	 * 		jre search in: /jvm
	 */
  static class ExtensionsFilter implements FilenameFilter {  	
    public boolean accept(File dir, String name) {
    	if (!dir.isDirectory()) return false;
    	// standardize paths to forward slash
      String path = XML.forwardSlash(dir.getPath());
      // accept jre extensions directories on all platforms
      if (path.endsWith("/lib/ext")) { //$NON-NLS-1$
        if (path.endsWith("/jre/lib/ext")) { //$NON-NLS-1$
        	return true;
        }
        String jre = XML.getName(path.substring(0, path.length()-8));
        if (jre.indexOf("jre")>-1) //$NON-NLS-1$
        	return true;
        if (path.indexOf("jdk")>-1 && path.indexOf("/Java/")>-1) //$NON-NLS-1$ //$NON-NLS-2$
        	return true;
        if (path.indexOf("jre")>-1 && path.indexOf("/Java/")>-1) //$NON-NLS-1$ //$NON-NLS-2$
        	return true;
      }
      // accept non-jre extensions directories
      // Linux
      if (path.endsWith("java/packages/lib/ext")) return true; //$NON-NLS-1$
      // Windows
      if (path.endsWith("Java/lib/ext")) return true; //$NON-NLS-1$
      // Mac OSX
      if (path.endsWith("Java/Extensions")) return true; //$NON-NLS-1$
      return false;
    }
  }

}
