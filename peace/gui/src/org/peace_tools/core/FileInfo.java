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
// Authors:   Dhananjai M. Rao              raodm@muohio.edu
//
//---------------------------------------------------------------------

package org.peace_tools.core;

/**
 * This is a simple class that is used to encapsulate information
 * about a given file. This class was introduced to enable 
 * handing back consistent information about files both on local
 * and remote servers.
 */
public class FileInfo {
	/**
	 * Bit value indicating if this entry is a directory. If
	 * this bit value is set in the file attributes then this
	 * entry is a directory.
	 */
	public static final int DIR_ATTRIB = 0x1;
	
	/**
	 * Bit value indicating if this entry is a regular file. If
	 * this bit value is set in the file attributes then this
	 * entry is a regular file.
	 */
	public static final int FILE_ATTRIB = 0x2;
	
	/**
	 * Bit value indicating if this entry is readable. If
	 * this bit value is set in the file attributes then this
	 * entry can be read by the user.
	 */
	public static final int READ_ATTRIB = 0x4;
	
	/**
	 * Bit value indicating if this entry is writable. If
	 * this bit value is set in the file attributes then this
	 * entry can be written to by the user.
	 */
	public static final int WRITE_ATTRIB = 0x8;
	
	/**
	 * Bit value indicating if this entry is executable. If
	 * this bit value is set in the file attributes then this
	 * entry can be executed by the user.
	 */
	public static final int EXEC_ATTRIB = 0x10;
	
	/**
	 * Create a complete file information class with all the 
	 * information.
	 * 
	 * This constructor must be used to create a complete object
	 * that contains all the associated information regarding
	 * a given file entry.
	 * 
	 * @param absolutePath The absolute path to the given file on
	 * a given host.
	 * 
	 * @param lastModified A long value representing the time the 
	 * file was last modified, measured in milliseconds since the 
	 * epoch (00:00:00 GMT, January 1, 1970), or 0L if the file does
	 * not exist or if an I/O error occurs.
	 *  
	 * @param size The size of the file in bytes. This value may
	 * be meaningful only for regular files.
	 * 
	 * @param attributes The various attributes defining the file. 
	 * These attributes are built by bitwise OR-ing together the
	 * various attributes defined in this class (ex: <code>
	 * FileInfo.FILE_ATTRIB | FileInfo.READ_ATTRIB</code>). If the
	 * file does not exist then the attributes must be 0. 
	 */
	public FileInfo(String absolutePath, long lastModified, 
			long size, int attributes) {
		this.path         = absolutePath;
		this.lastModified = lastModified;
		this.size         = size;
		this.attributes   = attributes;
	}

	/**
	 * The path to the file represented by this object.
	 * 
	 * @return The absolute path to the file.
	 */
	public String getPath() { return path; }

	/**
	 * The size of the file in bytes.
	 * 
	 * @return The size of the file in bytes. If the file does not
	 * exist then this method returns -1.
	 */
	public long getSize() { return size; }

	/** 
	 * Determine if this entry exists.
	 * 
	 * @return This method returns true if the file corresponding to
	 * this abstract entry exists.
	 */
	public boolean exists() { return attributes != 0; }
	
	/** Determine if this entry is a directory.
	 * 
	 * @return This method returns true if the entry is a directory.
	 * The return value is valid only if the entry exists.
	 */
	public boolean isDirectory() { return (attributes & DIR_ATTRIB) != 0; }

	/** Determine if this entry is a regular file.
	 * 
	 * @return This method returns true if the entry is a regular file.
	 * The return value is valid only if the entry exists.
	 */
	public boolean isFile() { return (attributes & FILE_ATTRIB) != 0; }

	/** Determine if this entry is readable.
	 * 
	 * @return This method returns true if the entry is readable by the
	 * user. The return value is valid only if the entry exists.
	 */
	public boolean canRead() { return (attributes & READ_ATTRIB) != 0; }

	/** Determine if this entry is writable.
	 * 
	 * @return This method returns true if the entry is writable by the
	 * user. The return value is valid only if the entry exists.
	 */
	public boolean canWrite() { return (attributes & WRITE_ATTRIB) != 0; }

	/** Determine if this entry is executable.
	 * 
	 * @return This method returns true if the entry is executable by the
	 * user. The return value is valid only if the entry exists.
	 */
	public boolean canExecute() { return (attributes & EXEC_ATTRIB) != 0; }

	/**
	 *  Returns the time that the file denoted by this 
	 *  abstract file information was last modified. 
	 *  
	 * @return A long value representing the time the 
	 * file was last modified, measured in milliseconds since the 
	 * epoch (00:00:00 GMT, January 1, 1970), or 0L if the file does
	 * not exist.
	 */
	public long getLastModified() { return lastModified; }

	@Override
	public String toString() {
		return "FileInfo for " + path + 
			": last modified: " + lastModified +
			", size: " + size + 
			", attributes 0x" + 
			Integer.toHexString(attributes);
	}
	
	/**
	 * The absolute path to the file referenced by this file attribute.
	 * This value is set in the constructor. 
	 */
	private final String path;
	
	/** A long value representing the time the 
	 * file was last modified, measured in milliseconds since the 
	 * epoch (00:00:00 GMT, January 1, 1970), or 0L if the file does
	 * not exist or if an I/O error occurs	 
	 */
	private final long lastModified;
	
	/**
	 * The size of the file in bytes. This value is -1 if the file
	 * does not exist.
	 */
	private long size;
	
	/**
	 * The various attributes defining the file. 
	 * These attributes are built by bitwise OR-ing together the
	 * various attributes defined in this class (ex: <code>
	 * FileInfo.FILE_ATTRIB | FileInfo.READ_ATTRIB</code>). If the
	 * file does not exist then the attributes must be 0. 
	 */
	private int attributes;
}
