#!/bin/bash

echo "Recompilando Taller1/punto_4a.cpp"
g++ -std=c++17 Taller1/punto_4a.cpp -o Taller1/punto_4a.x

for k in 0.1e10 0.2e10 0.5e10 1e10 2.5e10 10e10
do
    echo "Corriendo Taller1/punto_4a.x con k=$k"
    ./Taller1/punto_4a.x $k > data/Taller1/punto_4a_$k.dat
done

gnuplot Taller1/punto_4a.gp
