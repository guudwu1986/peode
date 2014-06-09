program main

  implicit none

  integer , parameter &
    :: dim_fun = 5

  double precision , dimension(dim_fun) &
    :: x

  integer , parameter &
    :: num_point = 2*dim_fun+1
  double precision &
    :: rhobeg = 2
  double precision &
    :: rhoend = 1e-8
  integer &
    :: iprint = 0
  integer &
    :: maxfun = 1000
  double precision , &
    dimension((num_point+13)*(num_point+dim_fun)+3*dim_fun*(dim_fun+3)/2) &
    :: w

  interface
    SUBROUTINE NEWUOA (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(*),W(*)
    end SUBROUTINE NEWUOA
  end interface

  x = 0

  call newuoa &
    ( &
      dim_fun &
      , num_point &
      , x &
      , rhobeg &
      , rhoend &
      , iprint &
      , maxfun &
      , w &
    )

  write(*,*) x

  stop

end program main
