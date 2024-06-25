# Default target
all: image

# Compile image.cpp into executable named image
image: image.cpp
	g++ $< -o image

# Generic rule for compiling any .cpp file into the image executable
%: %.cpp
	g++ $< -o image

clean: 
	rm -f image

format: 
	clang-format -i *.cpp *.h