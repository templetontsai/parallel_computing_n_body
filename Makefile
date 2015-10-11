CC=gcc
CFLAGS=-Wall
LDFLAGS=-lm -lrt
OBJ=nbody_serial.o
EXEC=nbody_serial
nbody_serial: $(OBJ)
	$(CC) $(CFLAGS) -o $(EXEC) nbody_serial.c $(LDFLAGS)
all: nbody_serial
.PHONY: clean
clean:
	rm -rf *.o nbody_serial
