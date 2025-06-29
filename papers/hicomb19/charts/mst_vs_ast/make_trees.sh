#!/bin/bash

# Script to generate plots of MST
~/miami/research/peace/src/peace --fastaFiles HA_500.fa --output-mst-file HA_500.mst --output-cls-file HA_500.cls --clusterMaker mst --flat-print 2> /dev/null

# Create a DOT output using mst and cls files
echo "digraph HA_500_mst {" > "ha_500_digraph.dot"
echo 'node[label="",width=0.6,shape=circle];' >>  "ha_500_digraph.dot"
grep "^[0-9]" HA_500.mst | cut -d"," -f1,2 | sed 's/,/ -> /g' | sed 's/$/ [dir=none];/g' >> "ha_500_digraph.dot"
# Add coloring information to the dot file.
cut -d"," -f1,3 HA_500.cls | tr "," "=" | cut -d"=" -f2,4 | tr "=" " " | awk '{ print $2 " [style=filled,colorscheme=", ($1 < 13) ? "set312,color=" $1 : "rdylbu11,color=" $1-12, "];"}' >> "ha_500_digraph.dot"
echo "}" >> "ha_500_digraph.dot"

# Generate PDF
neato -Tpdf -Gsize=6,6\! ha_500_digraph.dot -o ha_500_digraph.pdf

#-----------------------------------------------------------
# Script to generate plots of AST
~/miami/research/peace/src/peace --fastaFiles HA_500.fa --output-mst-file HA_500_ast.mst --output-cls-file HA_500_ast.cls --clusterMaker mst --flat-print --maxUse 1.0 2> /dev/null

# Create a DOT output using mst and cls files
echo "digraph HA_500_ast {" > "ha_500_ast_digraph.dot"
echo 'node[label="",width=0.3,shape=circle];' >>  "ha_500_ast_digraph.dot"
grep "^[0-9]" HA_500_ast.mst | cut -d"," -f1,2 | sed 's/,/ -> /g' | sed 's/$/ [dir=none];/g' >> "ha_500_ast_digraph.dot"
# Add coloring information to the dot file.
cut -d"," -f1,3 HA_500_ast.cls | tr "," "=" | cut -d"=" -f2,4 | tr "=" " " | awk '{ print $2 " [style=filled,colorscheme=", ($1 < 13) ? "set312,color=" $1 : "rdylbu11,color=" $1-12, "];"}' >> "ha_500_ast_digraph.dot"
echo "}" >> "ha_500_ast_digraph.dot"

# Generate PDF
neato -Tpdf -Gsize=6,6\! ha_500_ast_digraph.dot -o ha_500_ast_digraph.pdf



# end of script
