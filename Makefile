CC := gcc
FLAGS := -O3 -Wall
BINARIES := bin
SOURCE := src
BODIES = 500
ITERATIONS = 50000#110000 för del 1
SIZE = 1000
FAR =  1000
THREADS = 8
SEED = 100
DEFINES = -DVISUALIZE

benchmark:
	@mkdir -p $(BINARIES)
	@rm -rf outfiles
	@rm -rf benchmark
	$(CC) $(DEFINES) $(FLAGS) -o $(BINARIES)/1 $(SOURCE)/1.c
	$(CC) $(DEFINES) $(FLAGS) -o $(BINARIES)/2 $(SOURCE)/2.c
	$(CC) $(DEFINES) $(FLAGS) -o $(BINARIES)/3 $(SOURCE)/3.c
	$(CC) $(DEFINES) $(FLAGS) -o $(BINARIES)/4 $(SOURCE)/4.c
	./$(SOURCE)/benchmark $(ITERATIONS) $(SIZE) $(THREADS) $(FAR) $(SEED)
	gnuplot -p $(SOURCE)/1.plt
	gnuplot -p $(SOURCE)/2.plt
	gnuplot -p $(SOURCE)/3.plt
	gnuplot -p $(SOURCE)/4.plt

benchmark_1:
	@mkdir -p $(BINARIES)
	@rm -rf outfiles
	@rm -f benchmark/1.dat
	$(CC) $(DEFINES) $(FLAGS) -o $(BINARIES)/1 $(SOURCE)/1.c
	./$(SOURCE)/benchmark_1 $(ITERATIONS) $(SIZE) $(SEED)
	gnuplot -p $(SOURCE)/1.plt

benchmark_2:
	@mkdir -p $(BINARIES)
	@rm -rf outfiles
	@rm -f benchmark/2_*.dat
	$(CC) $(DEFINES) $(FLAGS) -o $(BINARIES)/2 $(SOURCE)/2.c
	./$(SOURCE)/benchmark_2 $(ITERATIONS) $(SIZE) $(THREADS) $(SEED)
	gnuplot -p $(SOURCE)/2.plt

benchmark_3:
	@mkdir -p $(BINARIES)
	@rm -rf outfiles
	@rm -f benchmark/3.dat
	$(CC) $(DEFINES) $(FLAGS) -o $(BINARIES)/3 $(SOURCE)/3.c
	./$(SOURCE)/benchmark_3 $(ITERATIONS) $(SIZE) $(FAR) $(SEED)
	gnuplot -p $(SOURCE)/3.plt

benchmark_4:
	@mkdir -p $(BINARIES)
	@rm -rf outfiles
	@rm -f benchmark/4_*.dat
	$(CC) $(DEFINES) $(FLAGS) -o $(BINARIES)/4 $(SOURCE)/4.c
	./$(SOURCE)/benchmark_4 $(ITERATIONS) $(SIZE) $(THREADS) $(FAR) $(SEED)
	gnuplot -p $(SOURCE)/4.plt

.PHONY: benchmark benchmark_1 benchmark_2 benchmark_3 benchmark_4

4: $(SOURCE)/4.c
	@mkdir -p $(BINARIES)
	@rm -rf outfiles
	@make -C visualize/
	$(CC) $(DEFINES) $(FLAGS) -o $(BINARIES)/$@ $^
	./$(BINARIES)/4 $(BODIES) $(ITERATIONS) $(THREADS) $(FAR) $(SIZE) $(SEED)
	./visualize/bin/Visualize ./outfiles/4.out

3: $(SOURCE)/3.c
	@mkdir -p $(BINARIES)
	@rm -rf outfiles
	@make -C visualize/
	$(CC) $(DEFINES) $(FLAGS) -o $(BINARIES)/$@ $^
	./$(BINARIES)/3 $(BODIES) $(ITERATIONS) $(FAR) $(SIZE) $(SEED)
	./visualize/bin/Visualize ./outfiles/3.out

2: $(SOURCE)/2.c
	@mkdir -p $(BINARIES)
	@rm -rf outfiles
	@make -C visualize/
	$(CC) $(DEFINES) $(FLAGS) -o $(BINARIES)/$@ $^
	./$(BINARIES)/2 $(BODIES) $(ITERATIONS) $(THREADS) $(SIZE) $(SEED)
	./visualize/bin/Visualize ./outfiles/2.out

1: $(SOURCE)/1.c
	@mkdir -p $(BINARIES)
	@rm -rf outfiles
	@make -C visualize/
	$(CC) $(DEFINES) $(FLAGS) -o $(BINARIES)/$@ $^
	./$(BINARIES)/1 $(BODIES) $(ITERATIONS) $(SIZE) $(SEED)
	./visualize/bin/Visualize ./outfiles/1.out