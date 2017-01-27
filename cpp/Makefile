IDIR =./
CC=g++
ODIR=obj

_DEPS = hierarchicalMaterials.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = hierarchicalMaterials.o bond.o node.o hexagon.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $<

hierarchicalMaterials: $(OBJ)
	$(CC) -o  $@ $^

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core