! Testing program for subroutine "InverseNewuoaSparse"
! in module "Linear_Ode_Inverse_mod".

! Read data from "data", compute parameters from observation data,
! then compare with true values.
! Also compare estimated system solution with observation data.

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

  use Linear_Ode_Inverse_mod
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
    :: rhobeg = 10
  double precision &
    :: rhoend = 1e-3
  integer &
    :: maxfun = 1e4
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
    :: par
  double precision , dimension(:) , allocatable &
    :: observation
  double precision , dimension(:) , allocatable &
    :: estimate_result

  integer , dimension(:) , allocatable &
    :: info
  integer , dimension(:) , allocatable &
    :: iter

! Initialization!{{{

  open ( newunit = NUM_UNIT , file = 'data' ,&
    action = 'read' , status = 'old' )

  read ( NUM_UNIT , * ) dim_ode
  read ( NUM_UNIT , * ) dim_time
  read ( NUM_UNIT , * ) dim_linear
  read ( NUM_UNIT , * ) dim_constant

  allocate ( timepoint(dim_time) )
  read ( NUM_UNIT , * ) timepoint

  allocate ( observation(dim_ode*dim_time) )
  read ( NUM_UNIT , * ) observation

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
!}}}

  allocate ( par(dim_ode+dim_linear+dim_constant) )
  par = 0
  par ( 1:dim_ode ) = observation ( 1:dim_ode )
!  par ( (dim_ode+1):(dim_ode+dim_linear) ) = linear
!  par ( (dim_ode+dim_linear+1):(dim_ode+dim_linear+dim_constant) ) &
!    = constant

  call InverseNewuoaSparse &
    ( &
      dim_ode &
      , dim_time &
      , dim_linear &
      , dim_constant &
      , nonzero_linear &
      , nonzero_constant &
      , par &
      , timepoint &
      , observation &
      , rhobeg &
      , rhoend &
      , maxfun &
      , tol &
    )

  write(*,*) 'Initial:'
  write(*,*) initial - par ( 1:dim_ode )
  write(*,*) 'Linear:'
  write(*,*) linear - par ( (dim_ode+1):(dim_ode+dim_linear) )
  write(*,*) 'Constant:'
  write(*,*) &
    constant &
    - par ( (dim_ode+dim_linear+1):(dim_ode+dim_linear+dim_constant) )

  allocate ( estimate_result(dim_ode*dim_time) )
  allocate ( iter(dim_time) )
  allocate ( info(dim_time) )

  call OdeSolveSparse &
    ( &
      dim_ode &
      , dim_time &
      , dim_linear &
      , dim_constant &
      , estimate_result &
      , par ( 1 : dim_ode ) &
      , nonzero_linear &
      , par ( (dim_ode+1) : (dim_ode+dim_linear) ) &
      , nonzero_constant &
      , par ( &
          (dim_ode+dim_linear+1) &
          : (dim_ode+dim_linear+dim_constant) &
        ) &
      , timepoint &
      , tol &
      , iter &
      , info &
    )
  write(*,*) 'Observation:'
  write(*,*) observation - estimate_result

  stop

end program main
