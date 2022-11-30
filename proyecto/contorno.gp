reset
set term pdf
set out "data/proyecto/antena.pdf"
set contour
unset surface
set cntrparam levels incr -2,0.05,2

DATA="data/proyecto/antena_50.dat"

set xrange[5:95]
set yrange[5:95]

set view map

set dgrid3d 100,100,4

set table "contour.txt"
splot DATA
unset table

unset contour
set surface
set table "dgrid.txt"
splot DATA
unset table

reset
set pm3d map
unset key
set palette defined (0 '#352a87', 1 '#0363e1',2 '#1485d4', 3 '#06a7c6', 4 '#38b99e', 5 '#92bf73', 6 '#d9ba56', 7 '#fcce2e', 8 '#f9fb0e')
set autoscale fix
set grid

splot 'dgrid.txt' w pm3d, 'contour.txt' w l lc rgb "black"
