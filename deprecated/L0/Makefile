exec:=gena
prefix='./'
arch=anaconda

CFLAGS = -O3 --std=c++17 -I$(INCLUDE)
LIB = -lpython3.12

ifeq ($(arch),miniconda)
	CC = g++
    INCLUDE=/opt/miniconda3/include/python3.12
    LIB_PATH=/opt/miniconda3/lib/
else
	ifeq ($(arch),anaconda)
	CC = g++
    INCLUDE=$(HOME)/anaconda3/include/python3.12
    LIB_PATH=$(HOME)/anaconda3/lib/
endif
endif
${exec}: main.o filter.o io.o packets.o structures.hpp
	$(CC) $(CFLAGS) -L$(LIB_PATH) -o $(exec)  main.o filter.o io.o packets.o  $(LIB) -Wl,-rpath $(LIB_PATH)

io.o:%.o:%.cpp ${io.o} structures.hpp
	$(CC) $(CFLAGS) -c $<

packets.o:%.o:%.cpp ${packets.o} structures.hpp
	$(CC) $(CFLAGS) -c $<

filter.o:%.o:%.cpp ${filter.o} structures.hpp
	$(CC) $(CFLAGS) -c $<

main.o:%.o:%.cpp ${main.o} structures.hpp
	$(CC) $(CFLAGS) -c $< 

clean :
	$(RM) *.o

