install: ../bin/gapdis ../bin/getpairs ../bin/joinOrder ../bin/joinGap ../bin/perfectjoin ../bin/collect2 ../bin/iLSLS ../bin/collect ../bin/scfwtmp ../bin/scfwtmp2 ../bin/scf ../bin/combine ../bin/scaffold mt2
../bin/gapdis:gapdis.cpp
	g++ -o ../bin/gapdis gapdis.cpp -std=c++0x -mcmodel=medium -g 
../bin/joinGap:joinGap.cpp
	g++ -o ../bin/joinGap joinGap.cpp -std=c++0x -mcmodel=medium -g 
../bin/joinOrder:joinOrder.cpp
	g++ -o ../bin/joinOrder joinOrder.cpp -std=c++0x -mcmodel=medium -g 
../bin/perfectjoin:perfectjoin.cpp
	g++ -o ../bin/perfectjoin perfectjoin.cpp -std=c++0x -mcmodel=medium -g 
../bin/getpairs:getpairs.cpp
	g++ -o ../bin/getpairs getpairs.cpp -std=c++0x -mcmodel=medium -g 
../bin/collect2:collect2.cpp
	g++ -o ../bin/collect2 collect2.cpp -std=c++0x -mcmodel=medium -g 
../bin/collect:collect.cpp
	g++ -o ../bin/collect collect.cpp -std=c++0x -mcmodel=medium -g 
../bin/iLSLS: iLSLS.cpp iLSLS.h
	g++ -o ../bin/iLSLS iLSLS.cpp -std=c++0x -mcmodel=medium -g 
../bin/combine:combine.cpp
	g++ -o ../bin/combine combine.cpp -std=c++0x -mcmodel=medium -g 
../bin/scaffold:scaffold.cpp
	g++ -o ../bin/scaffold scaffold.cpp -std=c++0x -mcmodel=medium -g 
mt2:
	sh compile.sh
