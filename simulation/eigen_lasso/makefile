all: eigen_lasso_sim.out

PATH_ROOT=$(PEODE)
PATH_NEWUOA=$(PATH_ROOT)/solver/newuoa_dp
PATH_RKF45=$(PATH_ROOT)/solver/rkf45
PATH_LINEAR_MODEL=$(PATH_ROOT)/model/linear
PATH_LINEAR_INVERSE=$(PATH_ROOT)/solver/linear_inverse
PATH_TEST=$(PATH_ROOT)/simulation/eigen_lasso

COMPILER=gfortran
STATIC_LIB_COMPILER=$(COMPILER) -c

eigen_lasso_sim.out: \
  bigden.o biglag.o newuoa.o newuob.o trsapp.o update.o \
  rkf45.o machine.o \
  linear_ode_mod.o \
  linear_inverse_eigen_mod.o \
  linear_inverse_mod.o \
  eigen_lasso_sim.o
	$(COMPILER) \
	  bigden.o biglag.o newuoa.o newuob.o \
	  trsapp.o update.o \
	  rkf45.o machine.o \
	  linear_ode_mod.o \
	  linear_inverse_eigen_mod.o \
	  linear_inverse_mod.o \
	  eigen_lasso_sim.o \
	  -llapack \
	  -o eigen_lasso_sim.out

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
rkf45.o: $(PATH_RKF45)/rkf45.f
	$(STATIC_LIB_COMPILER) $(PATH_RKF45)/rkf45.f \
	  -o rkf45.o
machine.o: $(PATH_RKF45)/machine.f
	$(STATIC_LIB_COMPILER) $(PATH_RKF45)/machine.f \
	  -o machine.o
linear_ode_mod.o: $(PATH_LINEAR_MODEL)/linear_ode_mod.f90
	$(STATIC_LIB_COMPILER) $(PATH_LINEAR_MODEL)/linear_ode_mod.f90 \
	  -o linear_ode_mod.o
linear_inverse_eigen_mod.o: \
  $(PATH_LINEAR_INVERSE)/linear_inverse_eigen_mod.f90
	$(STATIC_LIB_COMPILER) \
	  $(PATH_LINEAR_INVERSE)/linear_inverse_eigen_mod.f90 \
	  -o linear_inverse_eigen_mod.o
linear_inverse_mod.o: \
  $(PATH_LINEAR_INVERSE)/linear_inverse_mod.f90
	$(STATIC_LIB_COMPILER) \
	  $(PATH_LINEAR_INVERSE)/linear_inverse_mod.f90 \
	  -o linear_inverse_mod.o
eigen_lasso_sim.o: $(PATH_TEST)/eigen_lasso_sim.f90
	$(STATIC_LIB_COMPILER) \
	  $(PATH_TEST)/eigen_lasso_sim.f90 \
	  -o eigen_lasso_sim.o

clear:
	rm *.o *.out *.mod
