! Simulation for eigenvalue-updating algorithm.

program main

  use Linear_Ode_Inverse_mod

  implicit none

  integer &
    :: NUM_UNIT

  integer &
    :: dim_ode
  integer &
    :: dim_time
  double precision , dimension(:) , allocatable &
    :: timepoint
  double precision , dimension(:) , allocatable &
    :: eigen
  double precision , dimension(:) , allocatable &
    :: scaling
  double precision , dimension(:) , allocatable &
    :: observation
  double precision &
    :: ridge_parameter = 1d-10
  double precision , dimension(:) , allocatable &
    :: lasso_parameter
  double precision &
    :: rhobeg = 1d-0
  double precision &
    :: rhoend = 1d-3
  integer &
    :: maxfun = 100000

  double precision , dimension(:) , allocatable &
    :: linear_e
  double precision , dimension(:) , allocatable &
    :: initial_e
  double precision , dimension(:) , allocatable &
    :: curve_e

  integer &
    :: ind

! Read!{{{

  open ( newunit = NUM_UNIT , file = 'data' ,&
    action = 'read' , status = 'old' )

  read ( NUM_UNIT , * ) dim_ode
  read ( NUM_UNIT , * ) dim_time

  allocate ( timepoint(dim_time) )
  read ( NUM_UNIT , * ) timepoint

  allocate ( observation(dim_ode*dim_time) )
  read ( NUM_UNIT , * ) observation

  allocate ( eigen(dim_ode) )
  eigen = 0
  allocate ( scaling(dim_ode) )
  scaling ( 1 : (dim_ode/2) ) = -0.7
  scaling ( (dim_ode/2+1) : dim_ode ) &
    = (/(ind,ind=1,dim_ode/2)/) * 2 * 3.14
  scaling = scaling / ( timepoint(dim_time) - timepoint(1) )

  close ( unit = NUM_UNIT )
!}}}

  allocate ( lasso_parameter ( dim_ode*dim_ode ) )
  lasso_parameter = 0d0

  allocate ( linear_e ( dim_ode*dim_ode ) )
  allocate ( initial_e ( dim_ode ) )
  allocate ( curve_e ( dim_ode*dim_time ) )

  call InverseEigenLasso &
    ( &
      dim_ode &
      , dim_time &
      , eigen &
      , scaling &
      , timepoint &
      , observation &
      , ridge_parameter &
      , lasso_parameter &
      , rhobeg &
      , rhoend &
      , maxfun &
      , linear_e &
      , initial_e &
      , curve_e &
    )

  write(*,*) eigen
  write(*,*) linear_e

  stop

end program main
