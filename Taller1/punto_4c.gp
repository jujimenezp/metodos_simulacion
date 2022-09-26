set grid
set term pdf
set o 'data/Taller1/punto_4c.pdf'
set xl 'K'

set logscale
T(x)=A*x**a
fit T(x) 'data/Taller1/punto_4b.dat' u 1:2 via A,a

t(x)=B*x**b
B=1
b=-0.01
fit t(x) 'data/Taller1/punto_4b.dat' u 1:3 via B,b

Result = sprintf("T(k)=%.2f*k^{%.2f} \nt(k)=%.2f*k^{%.2f}",A,a,B,b)
set obj 1 rect from graph 0, 1 to graph 0.24, 0.86 fc rgb "white" front
set lab 1 Result at graph 0.01, 0.96 front

plot 'data/Taller1/punto_4b.dat' u 1:2 w lp pt 7 ps 0.2 title 'T medido',\
     '' u 1:3 w lp pt 7 ps 0.4 title 't medido',\
     T(x) w l title "T ajustado",\
     t(x) w l title "t ajustado"
