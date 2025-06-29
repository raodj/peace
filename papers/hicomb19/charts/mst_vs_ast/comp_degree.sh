#!/bin/bash

# Convenience script to print the average node degree from an MST
if [ $# -ne 1 ]; then
    echo "Specify the MST/AST file to be used"
    exit 2
fi

# Setup a variable to make the script readable
mst="$1"

# First get the nodes that have child-nodes associated with them
parentNodes=( `grep "^[0-9]" "${mst}" | cut -d"," -f1 | sort -nu` )
echo "Number of parent nodes: ${#parentNodes[@]}"

# Count the degree of each of the parent nodes to determine average
# degree.
degreeSum=0
for id in ${parentNodes[@]}
do
    nodeDegree=`grep -c "^${id}," "${mst}"`
    # echo "${id}, ${nodeDegree}"
    degreeSum=$(( degreeSum + nodeDegree ))
done

# Compute average degree
avgDegree=`echo "${degreeSum} / ${#parentNodes[@]}" | bc -l`
printf "Average node degree: %.2f\n" ${avgDegree}

# End of script
