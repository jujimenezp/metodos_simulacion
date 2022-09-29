set g
set xtics 10
set xl 'Tiempo'
set yl 'Coordenada Y'

set term pdf
set o 'data/Taller1/punto_5b.pdf'
set arrow from 60,25 to 60,80 nohead lc rgb "blue"
set label "t_{eq} = 60" at graph 0.32,0.8

plot 'data/Taller1/punto_5b.dat' u 1:2 w l title 'Y promedio'
