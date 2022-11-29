set grid
set term jpeg
set o 'data/4.jpg'

set xl 'Factor de atenuación D'
set yl 'Coeficiente de absorción C'

f(x) = A*x+B
fit f(x) 'data/4.dat' via A,B

Result = sprintf("C(D)=A*D+B \nA=%.2f \nB=%.2f",A,B)
set obj 1 rect from graph 0, 1 to graph 0.2, 0.86 fc rgb "white" front
set lab 1 Result at graph 0.01, 0.96 front

plot 'data/4.dat' w p, f(x)
