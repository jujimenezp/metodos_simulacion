#Compiling parameters
FLAGS = -std=c++11
SOURCES = $(wildcard ./*.cpp)

.PHONY: all $(SOURCES)
.PRECIOUS: %.x

all:
	@bash -c "echo $(SOURCES) | sed 's/cpp/a/'"
	@echo $(SOURCES)

%: %.cpp
	@echo 'Compiling $@.x'
	@g++ $(FLAGS) $< -o $@.x
	@echo -e 'Running $< \nData stored in $@'
	@./$@.x > data/$@.dat

%.x: %.cpp
	@g++ $(FLAGS) $< -o $@

clean:
	@rm -rf *.x
	@rm -rf data/*.dat
