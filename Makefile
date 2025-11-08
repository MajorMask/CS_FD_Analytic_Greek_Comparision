# Makefile for Black-Scholes Greeks Validation
# Usage:
#   make          - Compile the program
#   make run      - Compile and run
#   make clean    - Remove generated files
#   make analyze  - Run Python analysis script

CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall
LDFLAGS = -lm

TARGET = bs_greeks_validation
SOURCE = bs_greeks_validation.cpp
HEADERS = bs_call_price.h InverseCumulativeNormal.h

# CSV output files
CSV_FILES = bs_fd_vs_complex_scenario1.csv bs_fd_vs_complex_scenario2.csv

# Default target
all: $(TARGET)

# Compile the program
$(TARGET): $(SOURCE) $(HEADERS)
	@echo "Compiling $(TARGET)..."
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SOURCE) $(LDFLAGS)
	@echo "✓ Compilation successful!"

# Run the validation
run: $(TARGET)
	@echo "Running validation..."
	./$(TARGET)

# Generate plots and analysis
analyze: $(CSV_FILES)
	@echo "Generating plots and statistical analysis..."
	python3 analyze_results.py

# Clean generated files
clean:
	@echo "Cleaning generated files..."
	rm -f $(TARGET) $(CSV_FILES) greeks_error_analysis.png
	@echo "✓ Clean complete"

# Help
help:
	@echo "Available targets:"
	@echo "  make          - Compile the program"
	@echo "  make run      - Compile and run validation"
	@echo "  make analyze  - Generate plots (requires Python)"
	@echo "  make clean    - Remove generated files"
	@echo "  make help     - Show this help"

.PHONY: all run analyze clean help
