INCLUDES := -I/home/wpornpat/Codes/C/chemC-ATLAS/include -I/usr/include
LIBFLAGS := -L/usr/lib/atlas-base -latlas -lcblas -llapack_atlas -lm
SRC := /home/wpornpat/Codes/C/chemC-ATLAS/src

all: icf_calc

icf_calc: $(SRC)/icf_main.o $(SRC)/utils.o $(SRC)/pbc.o $(SRC)/matrix.o $(SRC)/icf2d.o $(SRC)/distance.o $(SRC)/coordnum.o $(SRC)/chemmatrixaux.o
	gcc -o icf_calc $(SRC)/icf_main.o $(SRC)/utils.o $(SRC)/pbc.o $(SRC)/matrix.o $(SRC)/icf2d.o $(SRC)/distance.o $(SRC)/coordnum.o \
	   $(SRC)/chemmatrixaux.o $(INCLUDES) $(LIBFLAGS)

$(SRC)/icf_main.o: $(SRC)/icf_main.c
	gcc -c $(SRC)/icf_main.c -o $(SRC)/icf_main.o $(INCLUDES) $(LIBFLAGS)

$(SRC)/utils.o: $(SRC)/utils.c
	gcc -c $(SRC)/utils.c -o $(SRC)/utils.o $(INCLUDES)

$(SRC)/pbc.o: $(SRC)/pbc.c
	gcc -c $(SRC)/pbc.c -o $(SRC)/pbc.o $(INCLUDES) $(LIBFLAGS)

$(SRC)/matrix.o: $(SRC)/matrix.c
	gcc -c $(SRC)/matrix.c -o $(SRC)/matrix.o $(INCLUDES) $(LIBFLAGS)

$(SRC)/icf2d.o: $(SRC)/icf2d.c
	gcc -c $(SRC)/icf2d.c -o $(SRC)/icf2d.o $(INCLUDES) $(LIBFLAGS)

$(SRC)/distance.o: $(SRC)/distance.c
	gcc -c $(SRC)/distance.c -o $(SRC)/distance.o $(INCLUDES) $(LIBFLAGS)

$(SRC)/coordnum.o: $(SRC)/coordnum.c
	gcc -c $(SRC)/coordnum.c -o $(SRC)/coordnum.o $(INCLUDES) $(LIBFLAGS)

$(SRC)/chemmatrixaux.o: $(SRC)/chemmatrixaux.c
	gcc -c $(SRC)/chemmatrixaux.c -o $(SRC)/chemmatrixaux.o $(INCLUDES) $(LIBFLAGS)
