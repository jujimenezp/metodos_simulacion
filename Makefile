#Compiling parameters
FLAGS = -std=c++17
SOURCES = $(shell find . -type f -name "*.cpp" | sed -e 's/\.cpp//g' -e 's/\.\///g')
BCyan=\033[1;36m
NC=\033[0m

.PHONY: graph
.PRECIOUS: %.x

all:
	@echo -e "${BCyan}PROGRAMAS:${NC}"
	@echo $(SOURCES) | tr " " "\n"

%: %.cpp
	@echo 'Compiling $@.x'
	@g++ $(FLAGS) $< -o $@.x
	@echo -e 'Running $< \nData stored in data/$@.dat'
	@./$@.x > data/$@.dat

%.x: %.cpp
	@g++ $(FLAGS) $< -o $@

clean:
	@find . -type f -name "*.x" -delete
	@find . -type f -name "*.dat" -delete
