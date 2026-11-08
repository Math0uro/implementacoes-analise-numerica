CC = gcc
CFLAGS = -std=c11 -O2 -Wall -Wextra -Wpedantic

all: programa

programa: main.c
	$(CC) $(CFLAGS) main.c -lm -o programa

interface: programa
	python interface.py

clean:
	rm -f programa

.PHONY: all interface clean
