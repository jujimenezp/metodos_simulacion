set g
set xl 'Temperatura k_B T'
set yl 'Presión'
set term pdf
set o 'data/Taller1/punto_5i_presion.pdf'

plot 'data/Taller1/punto_5i_presion.dat' w l notitle
