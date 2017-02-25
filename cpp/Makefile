IDIR =./
CC=g++
ODIR=obj
CFLAGS=-std=c++11 -O3
INCLUDES=-I/Users/ashivnishekhawat/Software/Suitesparse/SuiteSparse/include
LDIRS=-L/Users/ashivnishekhawat/Software/Suitesparse/SuiteSparse/lib -L/Users/ashivnishekhawat/Software/Suitesparse/SuiteSparse/CXSParse/Lib

LIBS=/Users/ashivnishekhawat/Software/Suitesparse/SuiteSparse/CSparse/Lib/libcsparse.a -lcxsparse -lm -lumfpack

_DEPS = hierarchicalMaterials.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = hierarchicalMaterials.o bond.o node.o hexagon.o coordinates_2d.o hexagonal_unit_cell.o\
       hexagonal_grid.o hierarchical_grid.o linear_system.o st_to_cc.o CPPSparse.o util_funcs.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(INCLUDES)

hierarchicalMaterials: $(OBJ)
	$(CC) $(CFLAGS) $(INCLUDES) $(LDIRS) -o $@ $^  $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core
