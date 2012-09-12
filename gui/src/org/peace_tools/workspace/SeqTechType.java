package org.peace_tools.workspace;

/**
 * Sequencing technology enumerations to match enumerations in
 * the PEACE schema file.
 * 
 * The following enumeration type is used internally within 
 * PEACE and DECAGON (particularly, the Synthetic Dataset 
 * Generator (SDG) wizard) to refer to the
 * three different types of sequencing technologies that are 
 * currently supported by SDG. The strings associated with
 * the enumerations refer to the names of the parameters in
 * the parameter set XML files used for defining the set of
 * parameters used by SDG.
 */
public enum SeqTechType {
	/**
	 * Enumeration used to identify Sanger sequencing technology.
	 */
	SANGER("Sanger"),
	/**
	 * Enumeration used to identify Roche 454 sequencing technology.
	 */
	R454("454"), 
	/**
	 * Enumeration used to identify Illumina sequencing technology.
	 */
	ILLUMINA("Illumina"),
	/**
	 * Enumeration used to identify unknown/unspecified sequencing technology.
	 */
	UNKNOWN("Unknown");
	
	
	/**
	 * The parameter prefix string that is used to in the Parameter Set
	 * XML description file associated with this technology. This string
	 * is case sensitive.
	 */
	private final String name;
	
	/**
	 * The constructor used to create enumeration constants.
	 * 
	 * @param name The parameter prefix string that is used to in the 
	 * Parameter Set XML description file associated with this technology.
	 * This string is case sensitive. 
	 */
	private SeqTechType(final String name) {
		this.name = name;
	}
	
	@Override
	public String toString() {
		return name;
	}
}