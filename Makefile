#Compiling parameters
FLAGS = -std=c++11

.PRECIOUS: %.x

all: euler

data/%.dat: %.x
	@echo -e 'Running $< \nData stored in $@'
	@./$< > $@

%: %.cpp
	@echo 'Compiling $@.x'
	@g++ $(FLAGS) $< -o $@.x
	@echo 'Running $@.x'
	@./$@.x

%.x: %.cpp
	@g++ $(FLAGS) $< -o $@

clean:
	@rm -rf *.x
	@rm -rf data/*.dat
