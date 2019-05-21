CC = cc
SRC = $(wildcard src/*.c)
OBJ = $(SRC:.c=.o)
BIN = ./bin
TARGET = project.out

LDFLAGS = -lm 
CFLAGS = -I./include -g -Wall -O3

all: dir $(BIN)/$(TARGET)

dir: ${BIN}

${BIN}:
	mkdir -p $(BIN)

%.o: %.c
	$(CC) -fopenmp  $(CFLAGS) -c -o $@ $<

$(BIN)/$(TARGET): $(OBJ)
	$(CC) -fopenmp  $(LDFLAGS) -o $@ $^ -lm

.PHONY: clean
clean:
	rm -f $(OBJ) $(BIN)/$(TARGET)
	
