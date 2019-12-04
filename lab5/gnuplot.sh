set term png 
set xl "iter"
set yl "S(iter)"
set out "result_s_it.png"
set title "S(it)"
p "k16.txt" u 1:2 w l t "S(it), k = 16",\
  "k8.txt" u 1:2 w l t "S(it), k = 8",\
  "k4.txt" u 1:2 w l t "S(it), k = 4",\
  "k2.txt" u 1:2 w l t "S(it), k = 2",\
  "k1.txt" u 1:2 w l t "S(it), k = 1"


set xl "x"
set yl "y"

set out "v_map_k16.png"
set title "V-map-k16"
plot "map16.txt" u 1:2:3 t "k = 16" with image

set out "v_map_k8.png"
set title "V-map-k8"
plot "map8.txt" u 1:2:3 t "k = 8" with image

set out "v_map_k4.png"
set title "V-map-k4"
plot "map4.txt" u 1:2:3 t "k = 4" with image

set out "v_map_k2.png"
set title "V-map-k2"
plot "map2.txt" u 1:2:3 t "k = 2" with image

set out "v_map_k1.png"
set title "V-map-k1"
plot "map1.txt" u 1:2:3 t "k = 1" with image