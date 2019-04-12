## Compiler, tools and options
CC         = gcc
LINKER     = gfortran
CCFLAGS    = -O1 -mavx #-Q --help=optimizers
LDFLAGS    = 

LAPACKHOME = /home/jascho/lapack-3.8.0
INCPATH    =  -I$(LAPACKHOME)/CBLAS/include
LIBS       =  -L$(LAPACKHOME) -lcblas -lrefblas 


## Files
OBJECTS = c_timer.o lab4.o 
TARGET  = lab4


## Implicit rules
.SUFFIXES: .c
.c.o:
	$(CC) -c  $(CCFLAGS) $(INCPATH) $<


## Build rules
all: $(TARGET)


$(TARGET): $(OBJECTS)
	$(LINKER) -o $@  $(OBJECTS)  $(CCFLAGS) $(LDFLAGS) $(LIBS)

clean:
	rm -f $(OBJECTS) $(TARGET)
	rm -f *~ core
