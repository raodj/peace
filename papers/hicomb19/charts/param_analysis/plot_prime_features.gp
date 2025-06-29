#
# This file requires gnuplot 5.0
#

# Set output file formats
set terminal pdfcairo enhanced color size 6in, 3in font ",18"
# set terminal pdfcairo enhanced color size 5in, 6in font ",18"
set output "primes_features_quality.pdf"

# Set input file formats
dataFile = "prime_features_anal.csv"
set datafile separator ','

set key left bottom

set xlabel "Number of primes features" font " Bold"
set ylabel "NMI/Purity" font " Bold"
set y2label "Number of clusters" font " Bold" tc lt 3

set yrange [0:1]
set ytics nomirror

set y2range [30:*]
set y2tics nomirror tc lt 3

# Setup log scale for x-axis
# set logscale x
set xrange [1:32]
red = "#d95319"
# set label "(note: log-scale)" at graph 0.5, 0.05 center font " Bold" tc rgb red

plot dataFile using 2:6 with linespoints lw 2 title "NMI",\
     '' using 2:7 with linespoints lw 2 title "Purity",\
     '' using 2:5 axes x1y2 with linespoints lw 2 title "#Clusters"

# Plot chart on runtime characteristics

set output "primes_features_runtime.pdf"
set key left center

set yrange [0:*]
set y2range [*:*]
set y2tics nomirror tc lt 6

set ylabel "Runtime (seconds)" font " Bold"
set y2label "Peak memory (MB)" font " Bold" tc lt 6

plot dataFile using 2:3 with linespoints lc 4 lw 2 title "Runtime",\
     '' using 2:($4/1000) axes x1y2 with linespoints lc 6 lw 2 title "Memory (MB)"
