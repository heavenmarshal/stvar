CC = g++
CFLAG = -O3 -fPIC -Wall
DEBUG = -g
INPATH = -I /usr/local/include/eigen3/ -I /usr/share/R/include -I /home/ruzhang/Dropbox/personal/projects/gpmodel
LIB = -L /usr/lib/R/lib -lR
stvar.so: stvarwrapper.o stvar.o stvarAr1.o stvarBase.o correlation.o
	$(CC) -shared -o $@ $^ $(LIB)
stvarBase.o: stvarBase.cpp stvarBase.hpp kernel.hpp
	$(CC) -o $@ $(CFLAG) $(DEBUG) $(INPATH) -c $<
stvar.o: stvar.cpp stvarBase.hpp stvar.hpp
	$(CC) -o $@ $(CFLAG) $(DEBUG) $(INPATH) -c $<
stvarAr1.o: stvarAr1.cpp stvarBase.hpp stvarAr1.hpp
	$(CC) -o $@ $(CFLAG) $(DEBUG) $(INPATH) -c $<
stvarwrapper.o: stvarwrapper.cpp stvarBase.hpp stvar.hpp stvarAr1.hpp
	$(CC) -o $@ $(CFLAG) $(DEBUG) $(INPATH) -c $<
all: stvar.so
.PHONY: clean
clean:
	rm -f *.o
