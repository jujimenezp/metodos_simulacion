set grid
set title 'Radios de giro'
set term jpeg
set o 'data/3.jpg'

set xl 'Tiempo'
set yl 'R_{giro}'

plot 'data/3_0.05.dat' u 1:2 w l title 'k_BT=0.05',\
     'data/3_0.5.dat' u 1:2 w l title 'k_BT=0.5',\
     'data/3_10.dat' u 1:2 w l title 'k_BT=10'
