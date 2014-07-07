! Test program for "InverseEigen" in module
! "Linear_Ode_Inverse_Eigen_mod"

program main

  use Linear_Ode_Inverse_Eigen_mod

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
    :: linear
  double precision , dimension(:) , allocatable &
    :: eigen
  double precision , dimension(:) , allocatable &
    :: scaling
  double precision , dimension(:) , allocatable &
    :: initial
  double precision , dimension(:) , allocatable &
    :: observation
  double precision &
    :: ridge_parameter = 1e-10
  double precision &
    :: rhobeg = 1
  double precision &
    :: rhoend = 1e-5
  integer &
    :: maxfun = 1000

  double precision , dimension(:) , allocatable &
    :: linear_e
  double precision , dimension(:) , allocatable &
    :: initial_e
  double precision , dimension(:) , allocatable &
    :: curve_e

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

  allocate ( eigen(dim_ode) )
  eigen = 1
  allocate ( scaling(dim_ode) )
  read ( NUM_UNIT , * ) scaling ( 1 : (dim_ode/2) )
  read ( NUM_UNIT , * ) scaling ( (dim_ode/2+1) : dim_ode )

  close ( unit = NUM_UNIT )
!}}}

  allocate ( linear_e ( dim_ode*dim_ode ) )
  allocate ( initial_e ( dim_ode ) )
  allocate ( curve_e ( dim_ode*dim_time ) )

  call InverseEigen &
    ( &
      dim_ode &
      , dim_time &
      , eigen &
      , scaling &
      , timepoint &
      , observation &
      , ridge_parameter &
      , rhobeg &
      , rhoend &
      , maxfun &
      , linear_e &
      , initial_e &
      , curve_e &
    )

  write(*,*) maxval ( abs ( linear_e - linear ) )
  write(*,*) maxval ( abs ( initial_e - initial ) )
  write(*,*) maxval ( abs ( curve_e - observation ) )

  stop

end program main
