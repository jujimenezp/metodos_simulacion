set grid
set xl 'Tiempo'
set yl 'Coordenada y'
set title 'Y promedio'

set term pdf
set o 'data/Taller1/punto_5h.pdf'
plot 'data/Taller1/punto_5h.dat' w l notitle
