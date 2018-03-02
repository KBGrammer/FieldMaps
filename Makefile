CCC         = g++
OPT         = -O0
CFLAGS=-I/usr/include/gsl
LDFLAGS=-lgsl -lgslcblas

all: fieldmap

MagFieldMap.o: MagFieldMap.cpp MagFieldMap.h
	$(CCC) $(OPT) $(CFLAGS) -c MagFieldMap.cpp

rng.o: rng.cpp rng.h
	$(CCC) $(OPT) $(CFLAGS) -c rng.cpp

mymath.o: mymath.cpp globals.h mymath.h
	$(CCC) $(OPT) $(CFLAGS) -c mymath.cpp

Vector3.o: Vector3.cpp globals.h Vector3.h mymath.o
	$(CCC) $(OPT) $(CFLAGS) -c Vector3.cpp mymath.o

Matrix.o: Matrix.cpp globals.h Matrix.h  Vector3.h mymath.o
	$(CCC) $(OPT) $(CFLAGS) -c Matrix.cpp mymath.o

EFieldPoint.o: EFieldPoint.cpp globals.h EFieldPoint.h
	$(CCC) $(OPT) $(CFLAGS) -c EFieldPoint.cpp

ChargedParticle.o: ChargedParticle.cpp Matrix.h globals.h ChargedParticle.h
	$(CCC) $(OPT) $(CFLAGS) -c ChargedParticle.cpp

EFieldMap.o: EFieldMap.cpp EFieldMap.h globals.h nanoflann.hpp EFieldPoint.o
	$(CCC) $(OPT) $(CFLAGS) -c EFieldMap.cpp EFieldPoint.o

fieldmap: main.cpp Vector3.o EFieldPoint.o EFieldMap.o mymath.o MagFieldMap.o ChargedParticle.o rng.o Matrix.o
	$(CCC) $(OPT) $(CFLAGS) -o fieldmap main.cpp Vector3.o EFieldMap.o EFieldPoint.o mymath.o MagFieldMap.o ChargedParticle.o Matrix.o rng.o $(LDFLAGS)

nanoflanntest: pointcloud_example.cpp nanoflann.hpp
	$(CCC) $(OPT) $(CFLAGS) -o nanoflanntest pointcloud_example.cpp

clean:
	rm -f *.o *~ */*.e4* */*.pe4* */*.po4*