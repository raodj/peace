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
		/** Enumeration used to report information regarding a job entry.
		 * The source object in the event is a Job object.
		 */
		JOB,
		/** Enumeration used to report information regarding a Server entry.
		 * The source object in the event is a Server object.
		 */
		SERVER,
		/**
		 * Enumeration used to report that the list of DB classifiers in this
		 * work space have all changed. Currently, the classifiers are all
		 * changed in one flail swoop. 
		 */
		CLASSIFIER_LIST
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
	 * has already occured to the data set.
	 */
	public WorkspaceEvent(DataSet ds, WorkspaceEvent.Operation operation) {
		super(ds);
		this.entryType = EntryType.DATA_SET;
		this.operation = operation;
	}

	/**
	 * Constructor to create an event that can be used to report change in the
	 * status of a mst data file.
	 * 
	 * @param mstFile The MST data file that has been inserted, deleted, or updated.
	 *  
	 * @param operation The type of operation (insert, delete, or update) that
	 * has already occured to the MST data.
	 */
	public WorkspaceEvent(MSTData mst, WorkspaceEvent.Operation operation) {
		super(mst);
		this.entryType = EntryType.MST_DATA;
		this.operation = operation;
	}

	/**
	 * Constructor to create an event that can be used to report change in the
	 * status of a mst cluster data file.
	 * 
	 * @param mstFile The MST cluster data file that has been inserted, 
	 * deleted, or updated.
	 *  
	 * @param operation The type of operation (insert, delete, or update) that
	 * has already occured to the MST cluster data.
	 */
	public WorkspaceEvent(MSTClusterData clusters, WorkspaceEvent.Operation operation) {
		super(clusters);
		this.entryType = EntryType.MST_CLUSTER_DATA;
		this.operation = operation;
	}

	/**
	 * Constructor to create an event that can be used to report change in the
	 * status of a job entry.
	 * 
	 * @param mstFile The job entry that has been inserted, deleted, or 
	 * updated.
	 *  
	 * @param operation The type of operation (insert, delete, or update) that
	 * has already occurred to the job.
	 */
	public WorkspaceEvent(Job job, WorkspaceEvent.Operation operation) {
		super(job);
		this.entryType = EntryType.JOB;
		this.operation = operation;
	}

	/**
	 * Constructor to create an event that can be used to report change in the
	 * status of a server entry.
	 * 
	 * @param mstFile The server entry that has been inserted, deleted, or 
	 * updated.
	 *  
	 * @param operation The type of operation (insert, delete, or update) that
	 * has already occurred to the server entry.
	 */
	public WorkspaceEvent(Server server, WorkspaceEvent.Operation operation) {
		super(server);
		this.entryType = EntryType.SERVER;
		this.operation = operation;
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
		this.entryType = EntryType.CLASSIFIER_LIST;
		this.operation = Operation.UPDATE;
	}

	/**
	 * Determine the type of entry regarding which a status change is being
	 * reported.
	 * 
	 * @note The entry type determines the data type of the source object included
	 * in the event as indicated in the enmeration definition(s).
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
	 * A generated serialization UID (required as the base class is defined 
	 * to be serializable).
	 */
	private static final long serialVersionUID = -2063315046609995997L;
}
