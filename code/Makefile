# Default target
all: image

# Compile image.cpp into executable named image
image: image.cpp
	g++ $< -o image

# Generic rule for compiling any .cpp file into the image executable
%: %.cpp
	g++ $< -o image

format: 
	find . -name '*.cpp' -o -name '*.h' | xargs clang-format -i

clean: 
	rm -f image