set term pdf
set o 'data/Taller1/punto_5f_histogram.pdf'
set style fill solid border -1
width=0.1
set boxwidth width
set xr [-11:11]
set xtics 1
sigma=sqrt(10)
miu=0
Gauss(x)=A*sigma*exp(-((x-miu)/sigma)**2/2)/(sigma*sqrt(2*pi))

fit Gauss(x) 'data/Taller1/punto_5f_histogram.dat' using 1:3 via A

set label "G(x) = A(σ²2π)^{-1/2}exp[(x-μ)²/2σ²]" at graph 0.05,0.9
plot 'data/Taller1/punto_5f_histogram.dat' using ($1+width/2):3 ti col smooth frequency with boxes,\
     Gauss(x) w l
