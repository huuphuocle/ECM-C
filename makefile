# the compiler: gcc for C program
CC = gcc

# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
CFLAGS  = -g -Wall

# library flags:
#	-lgmp	gmp.h
#	-lm 	math.h
LIB = -lgmp -lm

# C files
FILES = main.c ecm.c factor.c ladder.c

# headers
HEADERS = ecm.h

# the build target executable:
TARGET = myecm

all: $(FILES) $(HEADERS)
	$(CC) $(CFLAGS) -o $(TARGET) $(FILES) $(LIB)

clean: 
	rm -f $(TARGET)

run: $(TARGET)
	./$(TARGET)