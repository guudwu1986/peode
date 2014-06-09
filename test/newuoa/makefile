all: newuoa_test.out

PATH_ROOT=$(PEODE)
PATH_NEWUOA=$(PATH_ROOT)/solver/newuoa_dp
PATH_TEST=$(PATH_ROOT)/test/newuoa

COMPILER=gfortran -fdefault-real-8
STATIC_LIB_COMPILER=$(COMPILER) -c

newuoa_test.out: bigden.o biglag.o newuoa.o newuob.o trsapp.o update.o \
  calfun.o newuoa_test.o
	$(COMPILER) bigden.o biglag.o newuoa.o newuob.o \
	  trsapp.o update.o \
	  calfun.o newuoa_test.o \
	  -o newuoa_test.out

bigden.o: $(PATH_NEWUOA)/bigden.f
	$(STATIC_LIB_COMPILER) $(PATH_NEWUOA)/bigden.f \
	  -o bigden.o
biglag.o: $(PATH_NEWUOA)/biglag.f
	$(STATIC_LIB_COMPILER) $(PATH_NEWUOA)/biglag.f \
	  -o biglag.o
newuoa.o: $(PATH_NEWUOA)/newuoa.f
	$(STATIC_LIB_COMPILER) $(PATH_NEWUOA)/newuoa.f \
	  -o newuoa.o
newuob.o: $(PATH_NEWUOA)/newuob.f
	$(STATIC_LIB_COMPILER) $(PATH_NEWUOA)/newuob.f \
	  -o newuob.o
trsapp.o: $(PATH_NEWUOA)/trsapp.f
	$(STATIC_LIB_COMPILER) $(PATH_NEWUOA)/trsapp.f \
	  -o trsapp.o
update.o: $(PATH_NEWUOA)/update.f
	$(STATIC_LIB_COMPILER) $(PATH_NEWUOA)/update.f \
	  -o update.o
calfun.o: $(PATH_TEST)/calfun.f90
	$(STATIC_LIB_COMPILER) $(PATH_TEST)/calfun.f90 \
	  -o calfun.o
newuoa_test.o: $(PATH_TEST)/newuoa_test.f90
	$(STATIC_LIB_COMPILER) $(PATH_TEST)/newuoa_test.f90 \
	  -o newuoa_test.o

clear:
	rm *.o *.out *.mod
