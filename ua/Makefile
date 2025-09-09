# Makefile for UnitedAtom in SurfCo
# Based on NPCoronaPredict's UnitedAtom

CXX = g++
CXXFLAGS = -O3 -std=c++11 -Wall -I./src
LDFLAGS = -lm -lboost_filesystem -lboost_system

# Source files
SOURCES = $(wildcard src/*.cpp)
OBJECTS = $(SOURCES:src/%.cpp=build/%.o)
TARGET = build/UnitedAtom

# Default target
all: directories $(TARGET)

# Create necessary directories
directories:
	@mkdir -p build

# Build the executable
$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJECTS) $(LDFLAGS)
	@echo "UnitedAtom built successfully!"
	@echo "Creating symbolic link in SurfCo root..."
	@ln -sf ua/build/UnitedAtom ../UnitedAtom

# Build object files
build/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean build files
clean:
	rm -rf build/*.o $(TARGET)

# Install to SurfCo root
install: $(TARGET)
	@ln -sf ua/build/UnitedAtom ../UnitedAtom
	@echo "UnitedAtom installed to SurfCo root"

# Test the build
test: $(TARGET)
	$(TARGET) --help 2>&1 || echo "UA compiled successfully (no --help option)"

.PHONY: all clean install test directories