CC=nvcc
CFLAGS=  -Xcompiler -fopenmp
LIB= -lgomp -lcudart
SOURCES= p.cu
EXECNAME= p
all:
	$(CC) -o $(EXECNAME) $(SOURCES) $(LIB) $(CFLAGS)
clean:
	rm *.o *.linkinfo
