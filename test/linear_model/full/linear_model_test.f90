! Testing program for "OdeSolve" in module Linear_Ode_mod

! Read data from "data", compute solution, then compare to truth.

! Order of file "data":
!   Dimension of ODE
!   Dimension of time
!   Time points
!   True solution
!   Initial condition
!   Linear term
!   Constant term

program main

  use Linear_Ode_mod
  implicit none

  integer &
    :: ind
  integer &
    :: NUM_UNIT

  integer &
    :: dim_ode
  integer &
    :: dim_time
  double precision , dimension(:) , allocatable &
    :: time_point
  double precision &
    :: tol = 1e-8
  double precision , dimension(:) , allocatable , target &
    :: linear
  double precision , dimension(:) , allocatable , target &
    :: constant
  double precision , dimension(:) , allocatable &
    :: initial
  double precision , dimension(:) , allocatable &
    :: ode_result
  double precision , dimension(:) , allocatable &
    :: truth
  integer , dimension(:) , allocatable &
    :: iter
  integer , dimension(:) , allocatable &
    :: info

! Initialization!{{{

  open ( newunit = NUM_UNIT , file = 'data' ,&
    action = 'read' , status = 'old' )

  read ( NUM_UNIT , * ) dim_ode
  read ( NUM_UNIT , * ) dim_time

  allocate ( time_point(dim_time) )
  read ( NUM_UNIT , * ) time_point

  allocate ( truth(dim_ode*dim_time) )
  read ( NUM_UNIT , * ) truth

  allocate ( initial(dim_ode) )
  read ( NUM_UNIT , * ) initial

  allocate ( linear(dim_ode*dim_ode) )
  read ( NUM_UNIT , * ) linear

  allocate ( constant(dim_ode) )
  read ( NUM_UNIT , * ) constant

  close ( unit = NUM_UNIT )

  allocate ( ode_result(dim_ode*dim_time) )
  allocate ( iter(dim_time) )
  allocate ( info(dim_time) )
!}}}

  call OdeSolve &
    ( &
      dim_ode &
      , dim_time &
      , ode_result &
      , initial &
      , linear &
      , constant &
      , time_point &
      , tol &
      , iter &
      , info &
    )

  ode_result = ode_result - truth
  write(*,*) sum(ode_result**2)
  write(*,*) iter
  write(*,*) info

  stop

end program main
