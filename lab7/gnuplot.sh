#!/usr/bin/gnuplot 
set term png enhanced size 600,300 

set size ratio -1

set o "u-1000.png"
set xrange [0:2]
set yrange [0:0.9]
plot "q-1000.txt" u 1:2:5 with image notitle

set o "v-1000.png"
set xrange [0:2]
set yrange [0:0.9]
plot "q-1000.txt" u 1:2:6 with image notitle

set o "u-4000.png"
set xrange [0:2]
set yrange [0:0.9]
plot "q-4000.txt" u 1:2:5 with image notitle

set o "v-4000.png"
set xrange [0:2]
set yrange [0:0.9]
plot "q-4000.txt" u 1:2:6 with image notitle


set o "u+4000.png"
set xrange [0:2]
set yrange [0:0.9]
plot "q+4000.txt" u 1:2:5 with image notitle

set o "v+4000.png"
set xrange [0:2]
set yrange [0:0.9]
plot "q+4000.txt" u 1:2:6 with image notitle







set o "fi-1000.png"
set contour
set cbr [-55:-50]
set cntrparam levels increment -55,0.2,-50
unset surface
set view map
unset key
sp "q-1000.txt" u 1:2:3:3 w l lw -1 palette  t '' 


set o "zeta-1000.png"
set contour
set cbr [-200:400]
set cntrparam levels increment -200,10,400
unset surface
set view map
unset key
sp "q-1000.txt" u 1:2:4:4 w l lw -1 palette  t '' 


set o "fi-4000.png"
set contour
set cbr [-218:-202]
set cntrparam levels increment -218,1,-202
unset surface
set view map
unset key
sp "q-4000.txt" u 1:2:3:3 w l lw -1 palette  t '' 


set o "zeta-4000.png"
set contour
set cbr [-800:1200]
set cntrparam levels increment -800,10,1200
unset surface
set view map
unset key
sp "q-4000.txt" u 1:2:4:4 w l lw -1 palette  t ''



set o "fi+4000.png"
set contour
set cbr [202:218]
set cntrparam levels increment 202,1,218
unset surface
set view map
unset key
sp "q+4000.txt" u 1:2:3:3 w l lw -1 palette  t '' 


