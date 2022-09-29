set term pdf
set o 'data/Taller1/punto_5f_histogram.pdf'
set style fill solid border -1
width=0.1
set boxwidth width
set xr [-11:11]
set xtics 1
sigma=sqrt(10)
miu=0
Gauss(x)=4e5*exp(-((x-miu)/sigma)**2/2)/(sigma*sqrt(2*pi))

plot 'data/Taller1/punto_5f_histogram.dat' using ($1+width/2):3 ti col smooth frequency with boxes,\
     Gauss(x) w l
