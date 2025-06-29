#
# This file requires gnuplot 5.0
#

# Set  and output file formats
set terminal pdfcairo enhanced color size 6in, 3.5in
set output "maxUse_convergence.pdf"

set datafile separator ','

# Enable drawing plot inside a plot
set multiplot

set key maxrows 3
set yrange [0:65070]

set ylabel "Reads to process (n)" font " Bold"
set xlabel "Iteration in AST construction" font " Bold"

getData(maxUse) = sprintf("< grep -n '^@' maxUse_%s.o* | tr ':' ',' | cut -d',' -f1,5", maxUse)


plot getData("1.0") using ($1-5):(65069-$2) with linespoints lw 2 ps 0.5 title "AST=1.0",\
     getData("0.95") using ($1-5):(65069-$2) with linespoints lw 2 ps 0.3 title "AST=0.95",\
     getData("0.90") using ($1-5):(65069-$2) with linespoints lw 2 ps 0.3 title "AST=0.90",\
     getData("0.80") using ($1-5):(65069-$2) with linespoints lw 2 ps 0.3 title "AST=0.80",\
     getData("0.50") using ($1-5):(65069-$2) with linespoints lw 2 ps 0.3 title "AST=0.50",\
     getData("0.30") using ($1-5):(65069-$2) with linespoints lw 2 ps 0.3 title "AST=0.30",\
     getData("0.10") using ($1-5):(65069-$2) with linespoints lw 2 ps 0.3 title "AST=0.10",\


# Plot runtimes as bar charts inside the main plot
set origin 0.52, 0.35
set size 0.45, 0.5

set yrange [0:*]
set xrange [*:*]

set xtics rotate right

set ylabel "Runtime (sec)" offset 2,0
set xlabel "AST threshold"

set style data histograms
set boxwidth 0.9 absolute
set style fill solid 0.75

# Generate sequential color numbers
col = 0
getNextColor(c)=(col=col+1,col)

reg(x) =  (440.3217426 * x * x) - (734.9822713 * x) + 320.4111546

plot "runtimes.csv" using ($0):2:(getNextColor(0)):xticlabels(1) with boxes lc variable notitle

# End multiplot
unset multiplot
