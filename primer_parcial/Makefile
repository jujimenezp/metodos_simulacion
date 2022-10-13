#Compiling parameters
PDF_VIEWER=okular
FLAGS = -std=c++17
SOURCES = $(shell find . -type f -name "*.cpp" | sed -e 's/\.cpp//g' -e 's/\.\///g')
BCyan=\033[1;36m
NC=\033[0m

.PHONY: punto1
.PRECIOUS: %.x

all:
	@echo -e "${BCyan}PROGRAMAS:${NC}"
	@echo $(SOURCES) | tr " " "\n"
	@echo -e "${BCyan}GRAPH:${NC}"
	@echo "make graph gp=gnuplot_script corre el script gp.gp y abre el PDF resultado gp.pdf con PDF_VIEWER (okular)"
	@echo -e "${BCyan}CLEAN:${NC}"
	@echo "make clean borra los .x y .dat en todas las subcarpetas"

%: %.cpp
	@echo 'Compiling $@.cpp'
	@g++ $(FLAGS) $< -o $@.x
	@echo -e 'Running $@.x'
	@./$@.x > data/$@.dat
	@echo -e "Data stored in data/$@.dat"

%.x: %.cpp
	@g++ $(FLAGS) $< -o $@

punto1:	1
	@echo -e "Graficando trayectoria desde 1.gp"
	@gnuplot 1.gp
	@echo -e "Animando desde data/1.dat"
	@gnuplot data/1.dat

punto2: 2
	@echo -e "Animando desde data/2.dat"
	@gnuplot data/2.dat

punto3: 3
	@echo -e "Graficando el radio de giro desde data/3.dat"
	@gnuplot 3.gp

clean:
	@find . -type f -name "*.x" -delete
	@find . -type f -name "*.dat" -delete
	@find . -type f -name "*.log" -delete
