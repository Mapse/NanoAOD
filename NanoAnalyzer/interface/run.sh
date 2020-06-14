g++ -o readTupleexecutable.o maintest.C libCLHEP.a libtLite.a `root-config --glibs --cflags` -g
./readTupleexecutable.o
