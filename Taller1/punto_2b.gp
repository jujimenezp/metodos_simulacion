set grid

set terminal pdf enhanced
set out 'data/Taller1/punto_2b.pdf'

set xl 'λ'
set yl 'f(λ)'

plot 'data/Taller1/punto_2b.dat' w l lc rgb 'blue' notitle
