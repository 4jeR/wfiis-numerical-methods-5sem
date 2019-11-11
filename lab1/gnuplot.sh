# task 1 euler
set term png 
set out "zad1.png"
set xl "t"
set yl "y(t)"
set title "Metoda Eulera"
p "wynik-1a.txt" u 1:2 w p pt 2 t "dt=0.01", \
  "wynik-1b.txt" u 1:2 w p pt 3 t "dt=0.1", \
  "wynik-1c.txt" u 1:2 w p pt 4 t "dt=1.0", \
  "wynik-1-analitics.txt" u 1:2 w p pt 6 t "analityczne"


# task 1 error
set term png 
set out "zad1-errors.png"
set xl "t"
set yl "error = y_num(t) - y_dok(t)"
set title "Blad globalny"
p "error-1a.txt" u 1:2 w p pt 2 t "dt=0.01", \
  "error-1b.txt" u 1:2 w p pt 3 t "dt=0.1", \
  "error-1c.txt" u 1:2 w p pt 4 t "dt=1.0", 


# task 2 rk2 
set term png 
set out "zad2.png"
set xl "t"
set yl "y(t)"
set title "Metoda Eulera"
p "wynik-2a.txt" u 1:2 w p pt 2 t "dt=0.01", \
  "wynik-2b.txt" u 1:2 w p pt 3 t "dt=0.1", \
  "wynik-2c.txt" u 1:2 w p pt 4 t "dt=1.0", \
  "wynik-2-analitics.txt" u 1:2 w p pt 6 t "analityczne"


# task 2 error
set term png 
set out "zad2-errors.png"
set xl "t"
set yl "error = y_num(t) - y_dok(t)"
set title "Blad globalny"
p "error-2a.txt" u 1:2 w p pt 2 t "dt=0.01", \
  "error-2b.txt" u 1:2 w p pt 3 t "dt=0.1", \
  "error-2c.txt" u 1:2 w p pt 4 t "dt=1.0", 


# task 3 rk4
set term png 
set out "zad3.png"
set xl "t"
set yl "y(t)"
set title "Metoda Eulera"
p "wynik-3a.txt" u 1:2 w p pt 2 t "dt=0.01", \
  "wynik-3b.txt" u 1:2 w p pt 3 t "dt=0.1", \
  "wynik-3c.txt" u 1:2 w p pt 4 t "dt=1.0", \
  "wynik-3-analitics.txt" u 1:2 w p pt 6 t "analityczne"


# task 3 error
set term png 
set out "zad3-errors.png"
set xl "t"
set yl "error = y_num(t) - y_dok(t)"
set title "Blad globalny"
p "error-3a.txt" u 1:2 w p pt 2 t "dt=0.01", \
  "error-3b.txt" u 1:2 w p pt 3 t "dt=0.1", \
  "error-3c.txt" u 1:2 w p pt 4 t "dt=1.0", 



# task 4
# Q(t)
set term png 
set out "zad4_qt.png"
set xl "t"
set yl "Q(t)"
set title "Metoda RK4, Q(t)"
p "zad4_0.5wo.txt" u 1:2 w p pt 2 t "0.5 omega_0", \
  "zad4_0.8wo.txt" u 1:2 w p pt 3 t "0.8 omega_0", \
  "zad4_1.0wo.txt" u 1:2 w p pt 4 t "1.0 omega_0", \
  "zad4_1.2wo.txt" u 1:2 w p pt 6 t "1.2 omega_0"


# I(t)
set term png 
set out "zad4_it.png"
set xl "t"
set yl "I(t)"
set title "Metoda RK4, I(t)"
p "zad4_0.5wo.txt" u 1:3 w p pt 2 t "0.5 omega_0", \
  "zad4_0.8wo.txt" u 1:3 w p pt 3 t "0.8 omega_0", \
  "zad4_1.0wo.txt" u 1:3 w p pt 4 t "1.0 omega_0", \
  "zad4_1.2wo.txt" u 1:3 w p pt 6 t "1.2 omega_0"