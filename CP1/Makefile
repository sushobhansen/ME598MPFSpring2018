#Makefile for CP1

PROGRAM = cp1
SOURCE = cp1.c initialize.c ppecoeffs.c ppesource.c momentum.c bc.c sor.c momentumcorr.c writedata.c bilinterp.c particleintegrate.c
OBJECTS = cp1.o initialize.o ppecoeffs.o ppesource.o momentum.o bc.o sor.o momentumcorr.o writedata.o bilinterp.o particleintegrate.o
ARCHIVE = cp1.tar

CC = gcc
CFLAGS = -lm -I.
DEPS = cp1headers.h
OPTIONS = 
#OPTIONS = -g -debug extended -traceback -O0

%.o: %.c $(DEPS)
	$(CC) $(OPTIONS) -c -o $@ $< $(CFLAGS)

$(PROGRAM):	$(OBJECTS)
	$(CC) $(OPTIONS) -o $(PROGRAM).out $(OBJECTS) $(CFLAGS)
	
clean: 
	rm -f $(OBJECTS) $(PROGRAM).out *.tar *.plt
	
rebuild:
	clean $(PROGRAM)
	
archive:
	@echo Creating code archive tar-file $(ARCHIVE) ...
	tar cf $(ARCHIVE) $(SOURCE) Makefile
	@ls -l $(ARCHIVE)

help:
	@echo Try:
	@echo make $(PROGRAM) .... to build the program named $(PROGRAM)
	@echo make clean .... to clean up, removing object files and program $(PROGRAM)
	@echo make archive .... to make an archive tar file you can transfer or submit
