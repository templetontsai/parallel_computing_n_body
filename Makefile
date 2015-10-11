CC=gcc
CFLAGS=-std=c99
LDFLAGS=-lm
OBJ=nbody_serial.o
EXEC=nbody_serial
nbody_serial: $(OBJ)
	$(CC) $(CFLAGS) -o $(EXEC) nbody_serial.c $(LDFLAGS)
all: nbody_serial
.PHONY: clean
clean:
	rm -rf *.o nbody_serial
