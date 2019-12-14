set term png 

set out "wyn1.png"
set title "nx=ny=4, e1 = e2 = 1"
plot "v1.txt" u 1:2:3 with image

set out "wyn2.png"
set title "nx=ny=50, e1 = e2 = 1"
set xrange [0:5]
set yrange [0:5]
plot "v2.txt" u 1:2:3 with image

set out "wyn3.png"
set title "nx=ny=100, e1 = e2 = 1"
set xrange [0:10]
set yrange [0:10]
plot "v3.txt" u 1:2:3 with image

set out "wyn4.png"
set title "nx=ny=200, e1 = e2 = 1"
set xrange [0:20]
set yrange [0:20]
plot "v4.txt" u 1:2:3 with image

set out "wyn5.png"
set title "nx=ny=100, e1 = e2 = 1"
set xrange [0:10]
set yrange [0:10]
plot "v5.txt" u 1:2:3 with image

set out "wyn6.png"
set title "nx=ny=100, e1 =1, e2 = 2"
set xrange [0:10]
set yrange [0:10]
plot "v6.txt" u 1:2:3 with image

set out "wyn7.png"
set title "nx=ny=100, e1 = 1, e2 = 10"
set xrange [0:10]
set yrange [0:10]
plot "v7.txt" u 1:2:3 with image
