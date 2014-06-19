! Testing program for subroutine "InverseNewuoa"
! in module "Linear_Ode_Inverse_mod".

! Read data from "data", compute parameters from observation data,
! then compare with true values.
! Also compare estimated system solution with observation data.

! Order of file "data":
!   Dimension of ODE
!   Dimension of time
!   Time points
!   True solution
!   Initial condition
!   Linear term
!   Constant term

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
  double precision , dimension(:) , allocatable &
    :: linear
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

! Read!{{{

  open ( newunit = NUM_UNIT , file = 'data' ,&
    action = 'read' , status = 'old' )

  read ( NUM_UNIT , * ) dim_ode
  read ( NUM_UNIT , * ) dim_time

  allocate ( timepoint(dim_time) )
  read ( NUM_UNIT , * ) timepoint

  allocate ( observation(dim_ode*dim_time) )
  read ( NUM_UNIT , * ) observation

  allocate ( initial(dim_ode) )
  read ( NUM_UNIT , * ) initial

  allocate ( linear(dim_ode*dim_ode) )
  read ( NUM_UNIT , * ) linear

  allocate ( constant(dim_ode) )
  read ( NUM_UNIT , * ) constant

  close ( unit = NUM_UNIT )
!}}}

  allocate ( par(dim_ode**2+2*dim_ode) )
  par = 0
  par ( 1:dim_ode ) = observation ( 1:dim_ode )
!  par ( (dim_ode+1):(dim_ode**2+dim_ode) ) = linear
!  par ( (dim_ode**2+dim_ode+1):(dim_ode**2+2*dim_ode) ) = constant

  call InverseNewuoa &
    ( &
      dim_ode &
      , dim_time &
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
  write(*,*) linear - par ( (dim_ode+1):(dim_ode**2+dim_ode) )
  write(*,*) 'Constant:'
  write(*,*) &
    constant - par ( (dim_ode**2+dim_ode+1):(dim_ode**2+2*dim_ode) )

  allocate ( estimate_result(dim_ode*dim_time) )
  allocate ( iter(dim_time) )
  allocate ( info(dim_time) )

  call OdeSolve &
    ( &
      dim_ode &
      , dim_time &
      , estimate_result &
      , par ( 1:dim_ode ) &
      , par ( (dim_ode+1):(dim_ode**2+dim_ode) ) &
      , par ( (dim_ode**2+dim_ode+1):(dim_ode**2+2*dim_ode) ) &
      , timepoint &
      , tol &
      , iter &
      , info &
    )

  write(*,*) 'Observation:'
  write(*,*) observation - estimate_result

  stop

end program main
