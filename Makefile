CC = gcc
CFLAGS = -Wall -Wconversion -O3

LIBS = -lm
EXEC = solver
TEST = test
OBJS = linalg.o mtxio.o 
INCS = linalg.h mtxio.h

all: solver test

$(EXEC): main.o $(OBJS)
	$(CC) -o $@ $^ $(LIBS)

$(TEST): unittests.o $(OBJS)
	$(CC) -o $@ $^ $(LIBS)

%.o: %.c $(INCS)
	$(CC) -c -o $@ $< $(CFLAGS)

.PHONY: clean
clean:
	rm -f main.o unittests.o $(OBJS) $(EXEC) $(TEST) *.~ *.mtx
