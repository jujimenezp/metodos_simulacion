#Compiling parameters
PDF_VIEWER=okular
FLAGS = -Wall -fsanitize=address -fsanitize=leak -std=c++17 -lgsl -lgslcblas -lm -O3
SOURCES = $(shell find . -type f -name "*.cpp" | sed -e 's/\.cpp//g' -e 's/\.\///g')
BCyan=\033[1;36m
NC=\033[0m

.PHONY: graph
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
	@echo 'Data stored in data/$@.dat'

%.x: %.cpp
	@g++ $(FLAGS) $< -o $@

graph: $(gp)
	@gnuplot $(gp).gp
	@okular data/$(gp).pdf &

clean:
	@find . -type f -name "*.x" -delete
	@find . -type f -name "*.dat" -delete
	@find . -type f -name "*.log" -delete
	@find . -type f -name "*.csv" -delete
