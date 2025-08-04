CXX = g++
CXXFLAGS = -O3 -std=c++20 -I. -Ireal_algebraic_numbers
LDFLAGS = -lboost_system -lboost_thread -lntl
SRC = Rational.cpp Polynomial.cpp RealAlgebraicNumber.cpp MyTimer.cpp real_algebraic_numbers/main.cpp
OBJ = $(SRC:.cpp=.o)
TARGET = main.exe

all: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	del /Q *.o *.exe 2>nul || rm -f *.o *.exe

.PHONY: all clean
