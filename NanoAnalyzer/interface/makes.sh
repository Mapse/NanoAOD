g++ -o VertexFitter.o -c vxlite/VertexFitter.cc libCLHEP.a libtLite.a
g++ -o DAFVertexFinder.o -c vxlite/DAFVertexFinder.cc libCLHEP.a libtLite.a
ar cr libFitterFinder.a DAFVertexFinder.o VertexFitter.o
g++ -o texMain.o -c texMain.C
g++ -o binary texMain.o -L. -lFitterFinder -ltLite -lCLHEP
