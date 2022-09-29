set g
set xtics 10
set xl 'Tiempo'
set yl 'Coordenada Y'

set term pdf
set o 'data/Taller1/punto_5b.pdf'

plot 'data/Taller1/punto_5b.dat' u 1:2 w l title 'Y promedio'
