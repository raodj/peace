#
# This file requires gnuplot 5.0
#

# Set output file formats
set terminal pdfcairo enhanced color size 6in, 3in font ",18"
# set terminal pdfcairo enhanced color size 5in, 6in font ",18"
set output "maxuse_quality.pdf"

# Set input file formats
dataFile = "maxuse_anal.csv"
set datafile separator ','

set key left center
# set key top right

set xlabel "AST-threshold value" font " Bold"
set ylabel "NMI/Purity" font " Bold"
set y2label "Number of clusters" font " Bold" tc lt 3

set yrange [0:1]
set ytics nomirror

# set y2range [30:*]
set y2tics nomirror tc lt 3

plot dataFile using 1:6 with linespoints lw 2 title "NMI",\
     '' using 1:7 with linespoints lw 2 title "Purity",\
     '' using 1:5 axes x1y2 with linespoints lw 2 title "#Clusters"

# Plot chart on runtime characteristics

set output "maxuse_runtime.pdf"

set yrange [0:*]
set y2range [*:*]
set y2tics nomirror tc lt 6

set ylabel "Runtime (seconds)" font " Bold"
set y2label "Peak memory (MB)" font " Bold" tc lt 6

plot dataFile using 1:3 with linespoints lc 4 lw 2 title "Runtime (sec)",\
     '' using 1:($4/1000) axes x1y2 with linespoints lc 6 lw 2 title "Memory (MB)"
