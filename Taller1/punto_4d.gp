set g
set xl 'Tiempo (s)'
set yl 'Torque (Dina*cm)'
set xrange [-1:15]

set term pdf font "Arial,10"
set o 'data/Taller1/punto_4d.pdf'

plot 'data/Taller1/punto_4d_1e09.dat' u 1:2 w l lw 4.5 title 'K=0.1e10',\
     'data/Taller1/punto_4d_2e09.dat' u 1:2 w l lw 4 title 'K=0.2e10',\
     'data/Taller1/punto_4d_5e09.dat' u 1:2 w l lw 3.5 title 'K=0.5e10',\
     'data/Taller1/punto_4d_1e10.dat' u 1:2 w l lw 3 title 'K=1e10',\
     'data/Taller1/punto_4d_2.5e10.dat' u 1:2 w l lw 2.5 title 'K=2.5e10',\
     'data/Taller1/punto_4d_1e11.dat' u 1:2 w l lw 2 title 'K=10e10'
