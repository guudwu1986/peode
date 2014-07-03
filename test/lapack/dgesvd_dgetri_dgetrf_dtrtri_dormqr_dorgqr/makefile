all: dgesvd_dgetri_dgetrf_dtrtri_dormqr_dorgqr_test.out

PATH_ROOT=$(PEODE)
PATH_TEST=$(PATH_ROOT)/test/lapack/dgesvd_dgetri_dgetrf_dtrtri_dormqr_dorgqr

COMPILER=gfortran
STATIC_LIB_COMPILER=$(COMPILER) -c

dgesvd_dgetri_dgetrf_dtrtri_dormqr_dorgqr_test.out: \
  dgesvd_dgetri_dgetrf_dtrtri_dormqr_dorgqr_test.o
	$(COMPILER) dgesvd_dgetri_dgetrf_dtrtri_dormqr_dorgqr_test.o \
	  -llapack \
	  -o dgesvd_dgetri_dgetrf_dtrtri_dormqr_dorgqr_test.out

dgesvd_dgetri_dgetrf_dtrtri_dormqr_dorgqr_test.o: \
  $(PATH_TEST)/dgesvd_dgetri_dgetrf_dtrtri_dormqr_dorgqr_test.f90
	$(STATIC_LIB_COMPILER) \
	  $(PATH_TEST)/dgesvd_dgetri_dgetrf_dtrtri_dormqr_dorgqr_test.f90 \
	  -o dgesvd_dgetri_dgetrf_dtrtri_dormqr_dorgqr_test.o

clear:
	rm *.o *.out
