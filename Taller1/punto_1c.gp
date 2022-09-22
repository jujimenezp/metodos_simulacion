set grid
set style line 1 lc rgb 'black' lt 1 lw 1 pt 7 ps 0.4

set xl 'Ritmo reproductivo básico R_0'
set yl 'Susceptibles finales'

set term pdf

set o 'data/Taller1/punto_1c.pdf'

plot 'data/Taller1/punto_1c.dat' u 1:2 w lp ls 1 notitle
