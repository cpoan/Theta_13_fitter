all: unified nGd all
unified: unified_main.cc
	g++ -g -o unified unified_main.cc -lEG `root-config --cflags` `root-config --libs` -lTreePlayer -lMathMore -lMinuit 
nGd: nGd_main.cc
	g++ -g -o nGd nGd_main.cc -lEG `root-config --cflags` `root-config --libs` -lTreePlayer -lMathMore -lMinuit 
all: main.cc
	g++ -g -o all main.cc -lEG `root-config --cflags` `root-config --libs` -lTreePlayer -lMathMore -lMinuit 

clean:
	rm ./*~ ./*.o ./main
