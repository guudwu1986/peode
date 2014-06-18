! Testing program for "OdeSolveSparse" in module Linear_Ode_mod

! Read data from "data", compute solution, then compare to truth.

! Order of file "data":
!   Dimension of ODE
!   Dimension of time
!   Number of non-zeroes in linear term
!   Number of non-zeroes in constant term
!   Time points
!   True solution
!   Initial condition
!   Index of non-zeroes in linear term
!   Non-zero entries in linear term
!   Index of non-zeroes in constant term
!   Non-zero entries in constant term

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
  integer &
    :: dim_linear
  integer &
    :: dim_constant
  double precision , dimension(:) , allocatable &
    :: timepoint
  double precision &
    :: tol = 1e-8
  integer , dimension(:) , allocatable &
    :: nonzero_linear
  double precision , dimension(:) , allocatable &
    :: linear
  integer , dimension(:) , allocatable &
    :: nonzero_constant
  double precision , dimension(:) , allocatable &
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
  read ( NUM_UNIT , * ) dim_linear
  read ( NUM_UNIT , * ) dim_constant

  allocate ( timepoint(dim_time) )
  read ( NUM_UNIT , * ) timepoint

  allocate ( truth(dim_ode*dim_time) )
  read ( NUM_UNIT , * ) truth

  allocate ( initial(dim_ode) )
  read ( NUM_UNIT , * ) initial

  allocate ( nonzero_linear(dim_linear) )
  read ( NUM_UNIT , * ) nonzero_linear

  allocate ( linear(dim_linear) )
  read ( NUM_UNIT , * ) linear

  allocate ( nonzero_constant(dim_constant) )
  read ( NUM_UNIT , * ) nonzero_constant

  allocate ( constant(dim_constant) )
  read ( NUM_UNIT , * ) constant

  close ( unit = NUM_UNIT )

  allocate ( ode_result(dim_ode*dim_time) )
  allocate ( iter(dim_time) )
  allocate ( info(dim_time) )
!}}}

  call OdeSolveSparse &
    ( &
      dim_ode &
      , dim_time &
      , dim_linear &
      , dim_constant &
      , ode_result &
      , initial &
      , nonzero_linear &
      , linear &
      , nonzero_constant &
      , constant &
      , timepoint &
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
