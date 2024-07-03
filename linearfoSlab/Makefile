PROG	= sfem

C_SRCS	= fem.c inputData.c FEMmatrix.c eigenModeCalculation.c \
	  Initialfield.c calcNeff.c propagation.c analysis.c \
	  calcDispersion.c material.c utility.c power.c intpol.c
C_OBJS	= fem.o inputData.o FEMmatrix.o eigenModeCalculation.o \
	  Initialfield.o calcNeff.o propagation.o analysis.o \
	  calcDispersion.o material.o utility.o power.o intpol.o

F_SRCS	= bhlu1c.c eigv4p.c
F_OBJS	= bhlu1c.o eigv4p.o

FLAG	= -O -g
CFLAGS = $(FLAG)
FFLAGS = $(FLAG)

LIBS	= -L/home/fujisawa/utility/lib -lMFC -lm
CC = g++

all : $(PROG)

$(PROG) : $(C_OBJS) $(F_OBJS)
	g++ $(FLAG) $(C_OBJS) $(F_OBJS) -o $@ $(LIBS)
#$(C_OBJS) $(F_OBJS) :
#	g++ -c -I/mnt8/denpa/tsuji/include $<

clean :
	/bin/rm -f $(PROG) $(C_OBJS) $(F_OBJS)
