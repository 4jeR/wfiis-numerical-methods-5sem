# task 1 picord
set term png 
set out "picord.png"
set xl "t"
set yl "y(t)"
set xrange [0:100]
set title "Picord"
p "picard.txt" u 1:2 w p pt 2 t "u(t)", \
  "picard.txt" u 1:3 w p pt 2 t "z(t) = N - u(t)"

# task 2 newton
set term png 
set out "newton.png"
set xl "t"
set yl "y(t)"
set xrange [0:100]
set title "Newton"
p "newton.txt" u 1:2 w p pt 2 t "u(t)", \
  "newton.txt" u 1:3 w p pt 2 t "z(t) = N - u(t)"

# task 3 RK2
set term png 
set out "rk2.png"
set xl "t"
set yl "y(t)"
set xrange [0:100]
set title "Niejawna metoda RK2"
p "rk2.txt" u 1:2 w p pt 2 t "u(t)", \
  "rk2.txt" u 1:3 w p pt 2 t "z(t) = N - u(t)"