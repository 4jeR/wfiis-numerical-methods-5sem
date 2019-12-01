set term png 



set xl "x"
set yl "y"


set out "global-map06.png"
set title "map-wynik-global06"
p "map-global-06.txt" u 1:3 w p pt 2 t "d(x)", \
  "map-global-06.txt" u 2:3 w p pt 2 t "S(y)"


set out "global-map10.png"
set title "map-wynik-global10"
p "map-global-10.txt" u 1:3 w p pt 2 t "d(x)", \
  "map-global-10.txt" u 2:3 w p pt 2 t "S(y)"

set xl "it"
set yl "S(it)"
set logscale x


set title "wynik-global"
set out "global.png"
p "global-06.txt" u 1:2 w p pt 4 t "S(it)", \
  "global-10.txt" u 1:2 w p pt 4 t "S2(it)"
 
set title "wynik-local"
set out "local.png"
p "local-10.txt" u 1:2 w p pt 4 t "S(it)", \
  "local-14.txt" u 1:2 w p pt 4 t "S2(it)", \
  "local-18.txt" u 1:2 w p pt 4 t "S3(it)", \
  "local-19.txt" u 1:2 w p pt 4 t "S4(it)"
 


