//
// This file was generated by the JavaTM Architecture for XML Binding(JAXB) Reference Implementation, vJAXB 2.1.10 
// See <a href="http://java.sun.com/xml/jaxb">http://java.sun.com/xml/jaxb</a> 
// Any modifications to this file will be lost upon recompilation of the source schema. 
// Generated on: 2012.06.18 at 06:29:58 PM EDT 
//


package org.peace_tools.decagon.jaxb;

import java.util.ArrayList;
import java.util.List;
import javax.xml.bind.annotation.XmlAccessType;
import javax.xml.bind.annotation.XmlAccessorType;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlType;


/**
 * 
 * 	This type defines the basic constructs for conditional parameters in
 * 	a DECAGON configuration or assembler parameter definition XML file.
 * 	Conditions can have nested sub-conditions and in such cases they are
 * 	called compound-conditions. Conditions that do not have nested
 * 	sub-conditions are called simple conditions. When a list of conditions
 * 	are specified connected by condition connectors (namely: "and" or "or"),
 * 	the connectors are applied in the order in which they are specified. This
 * 	is different than traditional Boolean expressions where "and" has higher
 * 	precedence over "or". An example of a conditional fragment in a
 * 	parameter is shown below:
 * 	
 * 	<blockquote><pre>
 * 	<decagon:parameter name="AsmErrorFlag">
 *         <decagon:condition var="OS" comp="eq" value="Linux" connector="and"/>
 *         <decagon:condition var="SANGER_SUB_ERR" comp="gt" value="0" connector="and"/>
 * 
 * 	    <decagon:kind>BOOLEAN</decagon:kind>
 * 	    <decagon:cmdLine>--has-subst-errs</decagon:cmdLine>
 * 	    <decagon:cmdLineKind>IMPLICIT_VALUE</decagon:cmdLineKind>
 * 	    <decagon:value>true</decagon:value>
 * 	    <decagon:mutable>false</decagon:mutable>
 * 	    <decagon:hidden>true</decagon:hidden>
 * 	    <decagon:description>Flag passed to assembler's command-line when 
 * 	    sanger sequences have substitution errors in them</decagon:description>
 *     </decagon:parameter>
 *     </pre></blockquote>
 * 			
 * 
 * 
 * 
 */
@XmlAccessorType(XmlAccessType.FIELD)
@XmlType(name = "Condition", namespace = "http://www.peace-tools.org/decagon/", propOrder = {
    "condition"
})
public class Condition {

    @XmlElement(namespace = "http://www.peace-tools.org/decagon/")
    protected List<Condition> condition;
    @XmlAttribute(required = true)
    protected String var;
    @XmlAttribute(required = true)
    protected ComparatorKind comp;
    @XmlAttribute
    protected Double numValue;
    @XmlAttribute
    protected String strValue;
    @XmlAttribute
    protected ConnectorKind connector;

    /**
     * 
     * 	An optional list of sub-conditions associated with this condition. The
     * 	sub-conditions are evaluated first prior to applying connectors with the
     * 	main condition. 
     * 					Gets the value of the condition property.
     * 
     * <p>
     * This accessor method returns a reference to the live list,
     * not a snapshot. Therefore any modification you make to the
     * returned list will be present inside the JAXB object.
     * This is why there is not a <CODE>set</CODE> method for the condition property.
     * 
     * <p>
     * For example, to add a new item, do as follows:
     * <pre>
     *    getCondition().add(newItem);
     * </pre>
     * 
     * 
     * <p>
     * Objects of the following type(s) are allowed in the list
     * {@link Condition }
     * 
     * 
     */
    public List<Condition> getCondition() {
        if (condition == null) {
            condition = new ArrayList<Condition>();
        }
        return this.condition;
    }

    /**
     * 
     * 	Each condition is based on the presence (is the variable defined
     * 	or not) or the value of a given DECAGON variable. This attribute
     * 	refers to the variable to be used to evaluate a condition.
     * 				
     * 
     * @return
     *     possible object is
     *     {@link String }
     *     
     */
    public String getVar() {
        return var;
    }

    /**
     * Sets the value of the var property.
     * 
     * @param value
     *     allowed object is
     *     {@link String }
     *     
     */
    public void setVar(String value) {
        this.var = value;
    }

    /**
     * 
     * 	Each condition requires a comparator that defines the type of 
     * 	check (or comparison) to be performed. This attribute specifies the
     * 	type of check to be performed. DECAGON supports a variety of 
     * 	predefined conditions that can be used as the value for this
     * 	attribute.
     * 				
     * 
     * @return
     *     possible object is
     *     {@link ComparatorKind }
     *     
     */
    public ComparatorKind getComp() {
        return comp;
    }

    /**
     * Sets the value of the comp property.
     * 
     * @param value
     *     allowed object is
     *     {@link ComparatorKind }
     *     
     */
    public void setComp(ComparatorKind value) {
        this.comp = value;
    }

    /**
     * 
     * 	Comparators other than "def" and "ndef" require a value to compare a
     * 	variable against. This attribute must be used to specify the constant
     * 	numeric value to be used in the comparsion operation. This value is
     * 	used only if the variable used in the condition also has a numeric
     * 	value. Otherwise this value is ignored and the strValue (if specified)
     * 	is used for comparison.
     * 				
     * 
     * @return
     *     possible object is
     *     {@link Double }
     *     
     */
    public Double getNumValue() {
        return numValue;
    }

    /**
     * Sets the value of the numValue property.
     * 
     * @param value
     *     allowed object is
     *     {@link Double }
     *     
     */
    public void setNumValue(Double value) {
        this.numValue = value;
    }

    /**
     * 
     * 	Comparators other than "def" and "ndef" require a value to compare a
     * 	variable against. This attribute must be used to specify the constant
     * 	string to be used in the comparsion operation. If intValue attribute
     * 	is specified then it is given preference over the strValue.
     * 				
     * 
     * @return
     *     possible object is
     *     {@link String }
     *     
     */
    public String getStrValue() {
        return strValue;
    }

    /**
     * Sets the value of the strValue property.
     * 
     * @param value
     *     allowed object is
     *     {@link String }
     *     
     */
    public void setStrValue(String value) {
        this.strValue = value;
    }

    /**
     * 
     * 	Connectors must be used to connector sequence of Conditions together to 
     * 	form a single compound boolean expression. The default connector is "and".
     * 				
     * 
     * @return
     *     possible object is
     *     {@link ConnectorKind }
     *     
     */
    public ConnectorKind getConnector() {
        if (connector == null) {
            return ConnectorKind.AND;
        } else {
            return connector;
        }
    }

    /**
     * Sets the value of the connector property.
     * 
     * @param value
     *     allowed object is
     *     {@link ConnectorKind }
     *     
     */
    public void setConnector(ConnectorKind value) {
        this.connector = value;
    }

}
