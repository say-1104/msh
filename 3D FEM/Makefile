PROG = optfem3d

C_SRCS = optfem3d.cpp inputData.cpp analysis.cpp elementMatrix.cpp \
	solveMatrixEquation.cpp intpol.cpp elementMatrix2D.cpp \
	output.cpp utility.cpp bicgstab.cpp pwsmp_pre.cpp \
	intpol3D.cpp matrix-calc.cpp pzssmp_run.cpp \
	mozaic.cpp analysis2.cpp port.cpp wfm.cpp calcFitness.cpp random.c

C_OBJS = optfem3d.o inputData.o analysis.o elementMatrix.o \
	solveMatrixEquation.o intpol.o elementMatrix2D.o \
	output.o utility.o bicgstab.o pwsmp_pre.o \
	intpol3D.o matrix-calc.o pzssmp_run.o \
	mozaic.o analysis2.o port.o wfm.o calcFitness.o random.o


CXEROOT = /home/muratsubaki/intel/oneapi
CC = $(CXEROOT)/mpi/latest/bin/mpiicpc
MKLROOT = $(CXEROOT)/mkl/latest
MPIROOT = $(CXEROOT)/mpi/latest


INCLUDE_PATH = -I$(CXEROOT)/compiler/include -I$(MPIROOT)/include



FLAGS = -gxx-name=g++ -O3 -shared-intel -m64 -parallel -qopenmp -ip -ipo -DMKL_ILP64 -mcmodel=large


LIBS = -Wl,--start-group \
           $(MKLROOT)/lib/intel64/libmkl_intel_ilp64.so \
           $(MKLROOT)/lib/intel64/libmkl_intel_thread.so \
           $(MKLROOT)/lib/intel64/libmkl_core.so  -Wl,--end-group
LIBS_PATH = -L$(MKLROOT)/lib/intel64 -L/usr/lib64 #-L$(CXEROOT)/compiler/lib/intel64
INCLUDE_PATH = -I$(MKLROOT)/../compiler/include/intel64/ \
               -I$(MKLROOT)/include/ \
               -I$(MKLROOT)/include/intel64 \
               -I$(MKLROOT)/include/intel64/ilp64


$(PROG) : $(C_OBJS) $(F_OBJS)
	$(CC) $(FLAGS) $(C_OBJS) -o $@ $(LIBS_PATH) $(LIBS) $(INCLUDE_PATH)
$(C_OBJS) : $(C_SRCS)
	$(CC) $(FLAGS) -c $(C_SRCS) $(INCLUDE_PATH)


clean :
	/bin/rm -f $(C_OBJS) $(F_OBJS) $(PROG)
    
