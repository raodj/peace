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

package org.peace_tools.workspace;

import java.util.EventObject;

public class WorkspaceEvent extends EventObject {
	/**
	 * Enumeration for the different types of Workpace entries 
	 * that can be reported via a WorkspaceEvent object. 
	 *
	 */
	public enum EntryType {
		/** Enumeration used to report information regarding a data set
		 * (that consists of a EST file, a set of MSTData objects, and a
		 * MSTClusterData file). The source object in the event is a 
		 * DataSet object.
		 */
		DATA_SET,
		/** Enumeration used to report information regarding a MST
		 * data file. The source object in the event is a MSTData object.
		 */
		MST_DATA,
		/** Enumeration used to report information regarding a cluster
		 * data file. The source object in the event is a MSTClusterData object.
		 */
		MST_CLUSTER_DATA,
		/** Enumeration used to report information regarding a job or job-list entry.
		 * The source object in the event is a JobOrListWrapper object.
		 */
		JOB,
		/** Enumeration used to report information regarding a Server entry.
		 * The source object in the event is a Server object.
		 */
		SERVER,
		/**
		 * Enumeration used to report that the list of DB classifiers in this
		 * work space have all changed. Currently, the classifiers are all
		 * changed in one fail swoop. 
		 */
		CLASSIFIER_LIST,
		/**
		 * Enumeration used to broadcast/report the list of file entries
		 * in this work space have all changed.
		 */
		GENERATED_FILE_LIST,
		/**
		 * Enumeration used to report that a single file entry with a 
		 * generated file list (within a data set) has changed.
		 */
		FILE_ENTRY,
		/**
		 * Enumeration used to report that a sublist or an entry
		 * in the sublist has changed. The object in this case
		 * is a JobList (for top-level entry) or a 
		 * JobOrListWrapper (for sublist).
		 */
		JOB_LIST
	};
	
	/**
	 * Enumeration for the different types of operations on a Workspace 
	 * entity that can be reported via a WorkspaceEvent object. 
	 *
	 */
	public enum Operation {
		/** Enumeration used to report that a new entry of a given type
		 * (a data set, a Job, or Server) has been added to the Workspace.
		 */
		INSERT,
		/** Enumeration used to report that the status (or information) 
		 * associated with an entry of a given type (a mst, a cluster, a Job, 
		 * or Server) has changed.
		 */
		UPDATE,
		/** Enumeration used to report that an existing entry of a given type
		 * (a data set, a Job, or Server) has been removed from the Workspace.
		 */
		DELETE 
	};
	
	/**
	 * Constructor to create an event that can be used to report change in the
	 * status of a complete data set.
	 * 
	 * @param ds The data set that has been inserted, deleted, or updated.
	 *  
	 * @param operation The type of operation (insert, delete, or update) that
	 * has already occurred to the data set.
	 */
	public WorkspaceEvent(DataSet ds, WorkspaceEvent.Operation operation) {
		super(ds);
		this.entryType    = EntryType.DATA_SET;
		this.operation    = operation;
		this.indexPos     = -1;
		this.deletedEntry = null;
	}

	/**
	 * Constructor to create an event that can be used to report change in the
	 * status of a job entry.
	 * 
	 * @param job The job entry that has been inserted, deleted, or 
	 * updated.
	 *  
	 * @param operation The type of operation (insert, delete, or update) that
	 * has already occurred to the job.
	 * 
	 * @param container The top-level job list or JobOrListWrapper for sublist
	 * that contains the job entry that was modified.
	 */
	public WorkspaceEvent(Job job, WorkspaceEvent.Operation operation) {
		super(job);
		this.entryType    = EntryType.JOB;
		this.operation    = operation;
		this.indexPos     = -1;
		this.deletedEntry = null;
	}

	/**
	 * Constructor to create an event that can be used to report change in the
	 * entries contained by a job-list.
	 * 
	 * @param jobList The job list entry whose contents have been modified.
	 *  
	 * @param operation The type of operation (insert, delete, or update) that
	 * has already occurred to the job list.
	 * 
	 * @param indexPos The index position within the job list where the 
	 * change occurred.
	 * 
	 * @param entry The object at the given index position that is being
	 * inserted, updated, or removed.
	 */
	public WorkspaceEvent(JobList jList, WorkspaceEvent.Operation operation, 
			int indexPos, Object entry) {
		super(jList);
		this.entryType    = EntryType.JOB_LIST;
		this.operation    = operation;
		this.indexPos     = indexPos;
		this.deletedEntry = entry;
	}

	/**
	 * Constructor to create an event that can be used to report change in the
	 * status of a server entry.
	 * 
	 * @param server The server entry that has been inserted, deleted, or 
	 * updated.
	 *  
	 * @param operation The type of operation (insert, delete, or update) that
	 * has already occurred to the server entry.
	 */
	public WorkspaceEvent(Server server, WorkspaceEvent.Operation operation) {
		super(server);
		this.entryType    = EntryType.SERVER;
		this.operation    = operation;
		this.indexPos     = -1;
		this.deletedEntry = null;
	}
	
	/**
	 * Constructor to create an event that can be used to report change to the
	 * list of DB classifiers. This method always sets the operation to
	 * an UPDATE as DB classifiers are changed in one swoop and individual
	 * entries are not updated (due to the nature of the classifiers).
	 * 
	 * @param list The classifier list that has been changed.
	 */
	public WorkspaceEvent(ClassifierList list) {
		super(list);
		this.entryType    = EntryType.CLASSIFIER_LIST;
		this.operation    = Operation.UPDATE;
		this.indexPos     = -1;
		this.deletedEntry = null;
	}

	/**
	 * Constructor to create an event that can be used to report change in the
	 * status of a generated file list.
	 * 
	 * @param gfl The generated file list entry that has been inserted, deleted, or 
	 * updated.
	 *  
	 * @param operation The type of operation (insert, delete, or update) that
	 * has already occurred to the entry.
	 */
	public WorkspaceEvent(GeneratedFileList gfl, WorkspaceEvent.Operation operation) {
		super(gfl);
		this.entryType    = EntryType.GENERATED_FILE_LIST;
		this.operation    = operation;
		this.indexPos     = -1;
		this.deletedEntry = null;
	}
	
	/**
	 * Constructor to create an event that can be used to report change in the
	 * status of a generated file entry.
	 * 
	 * @param gfl The generated file entry that has been inserted, deleted, or 
	 * updated.
	 *  
	 * @param operation The type of operation (insert, delete, or update) that
	 * has already occurred to the entry.
	 */
	public WorkspaceEvent(FileEntry fe, WorkspaceEvent.Operation operation) {
		super(fe);
		this.entryType    = EntryType.FILE_ENTRY;
		this.operation    = operation;
		this.indexPos     = -1;
		this.deletedEntry = null;
	}
	
	/**
	 * Obtain the index position within a job list where a change occurred.
	 * 
	 * The return value from this method is meaningful only when the
	 * {@link #getEntryType()} method returns {@link EntryType#JOB_LIST}.
	 * 
	 * @return The index position within the job list where a change 
	 * occurred.
	 */
	public int getIndexPos() { return indexPos; }
	
	/**
	 * Obtain the entry that has been removed from a list.
	 * 
	 * This method can be used to obtain the entry that has been removed
	 * from a list. A valid non-null object is returned by this method 
	 * when the {@link #getEntryType()} method returns {@link EntryType#JOB_LIST}.
	 * In other cases this method returns null.
	 * 
	 * @return This method returns the object that has been removed from a
	 * list. The index position from where this object was removed can be
	 * obtained from {@link #getIndexPos()} method.
	 */
	public Object getDeletedEntry() { return deletedEntry; }
	
	/**
	 * Determine the type of entry regarding which a status change is being
	 * reported.
	 * 
	 * <p><b>Note:</b>  The entry type determines the data type of the source object included
	 * in the event as indicated in the enumeration definition(s).</p>
	 * 
	 * @return One of the predefined enumerations identifying the type of entity
	 * that this event pertains to.
	 */
	public WorkspaceEvent.EntryType getEntryType() { return entryType; }

	/**
	 * Determine the type of operation reported by this event.
	 * 
	 * @return The resulting type of operation that is being reported by this
	 * event.
	 */
	public WorkspaceEvent.Operation getOperation() { return operation; }
	
	/**
	 * This instance variable is used to contain the type of the entry 
	 * regarding which a status change is being reported.
	 */
	private final EntryType entryType;
	
	/**
	 * This instance variable is used to contain the type of status change 
	 * being reported via this event.
	 */
	private final Operation operation;

	/**
	 * The index position within the {@link #getSource()} source where an entry
	 * was added, removed, or updated. If the index position is not
	 * applicable, then it is set to -1.
	 */
	private final int indexPos;
	
	/**
	 * This instance variable is used to hold the object that has been
	 * removed from a list. The list can be obtained via call to 
	 * {@link #getSource()} method. This object is currently used
	 * only when {@link #entryType} is set to {@link EntryType#JOB_LIST}.
	 * In other cases it is set to NULL. 
	 */
	private final Object deletedEntry;
	
	/**
	 * A generated serialization UID (required as the base class is defined 
	 * to be serializable).
	 */
	private static final long serialVersionUID = -2063315046609995997L;
}
