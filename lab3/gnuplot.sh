# task
set term png 
set title "RK_2"
set xl "t"

set out "RK2_xtt.png"
p "RK2_a.txt" u 1:3 w p pt 7 t "x(t) TOL = 10 ^ -2", \
  "RK2_b.txt" u 1:3 w p pt 7 t "x(t) TOL = 10 ^ -5", 
  
set out "RK2_dtt.png"
p "RK2_a.txt" u 1:2 w p pt 7 t "dt(t) TOL = 10 ^ -2", \
  "RK2_b.txt" u 1:2 w p pt 7 t "dt(t) TOL = 10 ^ -5", 

set out "RK2_vt.png"
p "RK2_a.txt" u 1:4 w p pt 7 t "v(t) TOL = 10 ^ -2", \
  "RK2_b.txt" u 1:4 w p pt 7 t "v(t) TOL = 10 ^ -5", 

set out "RK2_vx.png"
p "RK2_a.txt" u 3:4 w p pt 7 t "v(x) TOL = 10 ^ -2", \
  "RK2_b.txt" u 3:4 w p pt 7 t "v(x) TOL = 10 ^ -5"

set title "Trapezy"
set out "Trapez_xtt.png"
p "Trapez_a.txt" u 1:3 w p pt 7 t "x(t) TOL = 10 ^ -2", \
  "Trapez_b.txt" u 1:3 w p pt 7 t "x(t) TOL = 10 ^ -5", 
  
set out "Trapez_dtt.png"
p "Trapez_a.txt" u 1:2 w p pt 7 t "dt(t) TOL = 10 ^ -2", \
  "Trapez_b.txt" u 1:2 w p pt 7 t "dt(t) TOL = 10 ^ -5", 

set out "Trapez_vt.png"
p "Trapez_a.txt" u 1:4 w p pt 7 t "v(t) TOL = 10 ^ -2", \
  "Trapez_b.txt" u 1:4 w p pt 7 t "v(t) TOL = 10 ^ -5", 

set out "Trapez_vx.png"
p "Trapez_a.txt" u 3:4 w p pt 7 t "v(x) TOL = 10 ^ -2", \
  "Trapez_b.txt" u 3:4 w p pt 7 t "v(x) TOL = 10 ^ -5"
