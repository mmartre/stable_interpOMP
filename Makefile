# Compiling and linking options
CC      = gcc
CFLAGS  = -Wall -fPIC -std=c99 -O3 -march=native -ffast-math -fopenmp
LDLIBS  = -lgsl -lgslcblas -lm -L/usr/local/lib/


SOURCES = $(shell ls -1 *.c* | xargs)

EXEC_SOURCES = $(shell grep -l "^int main" $(SOURCES) | xargs)

EXECS = talk
#$(shell echo $(EXEC_SOURCES) | sed -e 's:\.c[p]*::g')
DEPS = $(shell echo $(SOURCES) | sed -e 's:\.c[p]*:\.d:g')

all:	main clean

main: main.o stable_common.o mcculloch.o stable_dist.o stable_fit.o stable_interp.o

%.d : %.c
	@set -e; $(CC) -MM $(CFLAGS) $< \
	| sed 's/\($*\)\.o[ :]*/\1.o $@ : /g' > $@; \
	[ -s $@ ] || rm -f $@

-include $(DEPS)

# Executable
%.o :	%.c
	@echo -n compilando objeto \'$<\'...
	@$(CC) $(CFLAGS) $< -c
	@echo [OK]


% :	%.o
	@echo -n compilando ejecutable \'$@\'...
	@$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)
	@echo [OK]

# Clean
clean:	
	@rm -f $(wildcard *.o *.d core* *.P) $(EXECS)

