IDIR =./
CC=g++
ODIR=obj
CFLAGS=-std=c++11 -O3

_DEPS = hierarchicalMaterials.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = hierarchicalMaterials.o bond.o node.o hexagon.o coordinates_2d.o hexagonal_unit_cell.o hexagonal_grid.o hierarchical_grid.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

hierarchicalMaterials: $(OBJ)
	$(CC) -o  $@ $^ $(CFLAGS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core
