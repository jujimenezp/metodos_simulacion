set g
set xl 'Tiempo (s)'
set yl 'Torque (Dina*cm)'
set xrange [0.174:0.178]
set ytics 1e7

set term pdf font "Arial,10"
set o 'data/Taller1/punto_4a.pdf'

plot 'data/Taller1/punto_4a_0.1e10.dat' u 1:2 w l title 'K=0.1e10',\
     'data/Taller1/punto_4a_0.2e10.dat' u 1:2 w l title 'K=0.2e10',\
     'data/Taller1/punto_4a_0.5e10.dat' u 1:2 w l title 'K=0.5e10',\
     'data/Taller1/punto_4a_1e10.dat' u 1:2 w l title 'K=1e10',\
     'data/Taller1/punto_4a_2.5e10.dat' u 1:2 w l title 'K=2.5e10',\
     'data/Taller1/punto_4a_10e10.dat' u 1:2 w l title 'K=10e10'
