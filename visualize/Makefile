SOURCE = src
BINARIES = bin

Visualize: $(SOURCE)/Visualize.cpp
	mkdir -p $(BINARIES)
	g++ $^ -o $(BINARIES)/$@ -lGLFW -framework OpenGL