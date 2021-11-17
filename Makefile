###############################################################################
# Makefile for Malaria Project, Spring 2021
###############################################################################
CC = mpicc
CFLAGS = -std=c99 -g -O3
LIBS = -lm

BIN = project prop malaria
files = project.c prop.c

all: $(files)
	$(CC) $(files) -o malaria $(CFLAGS) $(LIBS)

clean:
	$(RM) $(BIN)

	

