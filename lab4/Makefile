CC=g++
CXXFLAGS= -Ofast -Wall -pedantic
OBJS = lab4.o
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