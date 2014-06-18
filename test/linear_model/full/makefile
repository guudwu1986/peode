all: linear_model_full_test.out

PATH_ROOT=$(PEODE)
PATH_RKF45=$(PATH_ROOT)/solver/rkf45
PATH_LINEAR_MODEL=$(PATH_ROOT)/model/linear
PATH_TEST=$(PATH_ROOT)/test/linear_model/full

COMPILER=gfortran -fdefault-real-8
STATIC_LIB_COMPILER=$(COMPILER) -c

linear_model_full_test.out: rkf45.o machine.o linear_ode_mod.o \
  linear_model_full_test.o
	$(COMPILER) rkf45.o machine.o linear_ode_mod.o \
	  linear_model_full_test.o \
	  -o linear_model_full_test.out

rkf45.o: $(PATH_RKF45)/rkf45.f
	$(STATIC_LIB_COMPILER) $(PATH_RKF45)/rkf45.f \
	  -o rkf45.o
machine.o: $(PATH_RKF45)/machine.f
	$(STATIC_LIB_COMPILER) $(PATH_RKF45)/machine.f \
	  -o machine.o
linear_ode_mod.o: $(PATH_LINEAR_MODEL)/linear_ode_mod.f90
	$(STATIC_LIB_COMPILER) $(PATH_LINEAR_MODEL)/linear_ode_mod.f90 \
	  -o linear_ode_mod.o
linear_model_full_test.o: $(PATH_TEST)/linear_model_full_test.f90
	$(STATIC_LIB_COMPILER) $(PATH_TEST)/linear_model_full_test.f90 \
	  -o linear_model_full_test.o

clear:
	rm *.o *.out *.mod
