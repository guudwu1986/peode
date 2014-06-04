program main

  use Linear_Ode_mod
  implicit none

  interface
    subroutine rkf45(f,neqn,y,t,tout,relerr,abserr,iflag,work,iwork)
      integer neqn,iflag,iwork(5)
      double precision y(neqn),t,tout,relerr,abserr,work(1)
      external f
    end subroutine rkf45
  end interface

  integer &
    :: ind

  integer &
    :: dim_ode = 3
  integer &
    :: dim_time = 5
  double precision , dimension(:) , allocatable &
    :: time_point
  double precision &
    :: tol = 1e-8
  double precision , dimension(:,:) , allocatable &
    :: linear
  double precision , dimension(:) , allocatable &
    :: constant
  double precision , dimension(:) , allocatable &
    :: initial
  double precision , dimension(:,:) , allocatable &
    :: ode_result

  allocate ( time_point(dim_time) )
  time_point = (/ ( ind , ind=0,(dim_time-1) ) /)
  allocate ( linear(dim_ode,dim_ode) )
  linear = 0
  allocate ( constant(dim_ode) )
  constant = 0
  allocate ( initial(dim_ode) )
  initial = (/ ( ind , ind=1,dim_ode ) /)
  allocate ( ode_result(dim_ode,dim_time) )

  do ind = 1,dim_ode
    linear(ind,ind) = ind
    constant(ind) = ind
  end do

  call OdeSolve &
    ( &
      ode_result &
      , initial &
      , linear &
      , constant &
      , time_point &
      , tol &
    )

  write(*,*) ode_result

  stop

end program main
