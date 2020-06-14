g++ -fPIC -o DAFVertexFinder.so -c vxlite/DAFVertexFinder.cc
g++ -shared -Wl,-soname, -o libDAFVertexFinder.so DAFVertexFinder.so libtLite.a

g++ -fPIC -o VertexFitter.so -c vxlite/VertexFitter.cc
g++ -shared -Wl,-soname, -o libVertexFitter.so VertexFitter.so libtLite.a

cp libDAFVertexFinder.so ../../../../lib/slc5_amd64_gcc434/
cp libVertexFitter.so ../../../../lib/slc5_amd64_gcc434/
rm *.so
