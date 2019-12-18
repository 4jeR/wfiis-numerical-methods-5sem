CC=g++
CXXFLAGS= -g -Wall -pedantic
OBJS = lab7.o
OUT = exe

$(OUT): $(OBJS)
	$(CC)  -o $@  $^ $(CXXFLAGS)

.c.o:
	$(CC) -c -o $@ $< $(CXXFLAGS) 

.PHONY: clean run

clean:
	rm -f *o
	rm -f $(OUT)
	rm -f *.txt
	rm -f *.png
	clear

run:
	./$(OUT