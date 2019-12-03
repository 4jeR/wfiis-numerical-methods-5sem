set term png 

set xl "x"
set yl "y"

set xl "it"
set yl "S(it)"
set logscale x


set title "wynik-global"
set out "global.png"
p "global-06.txt" u 1:2 w p t "S(it)", \
  "global-10.txt" u 1:2 w p t "S2(it)"
 
set title "wynik-local"
set out "local.png"
p "local-10.txt" u 1:2 w l t "S(it)", \
  "local-14.txt" u 1:2 w l t "S2(it)", \
  "local-18.txt" u 1:2 w l t "S3(it)", \
  "local-19.txt" u 1:2 w l t "S4(it)"
 



set out "map_globalv.png"
set title "Relaxed V"
plot "map-global-06.txt" u 1:2:3 t "w = 0.6" with image,\
     "map-global-10.txt" u 1:2:3 t "w = 1.0" with image

set out "map_global_error.png"
set title "Error"
plot "map-global-06.txt" u 1:2:4 t "w = 0.6" with image,\
     "map-global-10.txt" u 1:2:4 t "w = 1.0" with image
