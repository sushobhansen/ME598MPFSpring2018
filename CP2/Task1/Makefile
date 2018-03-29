#Makefile for CP2Task1

PROGRAM = cp2task1
SOURCE = cp2task1.c writedata.c initialize.c bc.c advect.c reinitialize.c
OBJECTS = cp2task1.o writedata.o initialize.o bc.o advect.o reinitialize.o
ARCHIVE = cp2task1.tar

CC = gcc
CFLAGS = -lm -I.
DEPS = cp2task1headers.h
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
