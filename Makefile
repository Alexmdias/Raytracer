CC := g++
CFLAGS := -Wall -Wextra -W -std=c++17
CPPFLAGS := 


BIN := raytrace

SRC := $(shell find . -type f -name "*.cpp" ! -name "main.cpp")
OBJ := $(SRC:.cpp=.o)


$(BIN): main.cpp $(OBJ)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $^


.PHONY: test
test: $(BIN)
	./$(BIN) scenes/basic.ray


.PHONY: clean
clean:
	rm *.o $(BIN)


%.o: %.cpp
	$(CC) -c $^ $(CFLAGS)
