all: dgeqrf_test.out

PATH_ROOT=$(PEODE)
PATH_TEST=$(PATH_ROOT)/test/dgeqrf

COMPILER=gfortran
STATIC_LIB_COMPILER=$(COMPILER) -c

dgeqrf_test.out: dgeqrf_test.o
	$(COMPILER) dgeqrf_test.o \
	  -llapack \
	  -o dgeqrf_test.out

dgeqrf_test.o: $(PATH_TEST)/dgeqrf_test.f90
	$(STATIC_LIB_COMPILER) $(PATH_TEST)/dgeqrf_test.f90 \
	  -o dgeqrf_test.o

clear:
	rm *.o *.out
