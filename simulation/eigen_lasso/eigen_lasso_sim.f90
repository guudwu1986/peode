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
    :: ridge_parameter = 1d-10
  double precision , dimension(:) , allocatable &
    :: lasso_parameter
  double precision &
    :: rhobeg = 1d-1
  double precision &
    :: rhoend = 1d-3
  integer &
    :: maxfun = 100000
  double precision &
    :: tp_weight = 1
  double precision &
    :: tp
  double precision &
    :: tn
  double precision &
    :: tp_trial
  double precision &
    :: tn_trial

  double precision , dimension(:) , allocatable &
    :: linear_e
  double precision , dimension(:) , allocatable &
    :: initial_e
  double precision , dimension(:) , allocatable &
    :: curve_e

  double precision , dimension(:) , allocatable &
    :: zero
  double precision , dimension(:) , allocatable &
    :: nonzero
  double precision , dimension(:) , allocatable &
    :: zero_index
  double precision , dimension(:) , allocatable &
    :: nonzero_index
  double precision &
    :: max_zero
  double precision &
    :: min_nonzero

  integer &
    :: ind
  integer &
    :: ind1
  integer &
    :: ind2

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

  allocate ( lasso_parameter ( dim_ode*dim_ode ) )
  where ( linear .eq. 0 )
    lasso_parameter = 1d3
  elsewhere
    lasso_parameter = 0
  end where

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

!  write(*,*) maxval ( abs ( linear_e - linear ) )
!  write(*,*) maxval ( abs ( initial_e - initial ) )
!  write(*,*) maxval ( abs ( curve_e - observation ) )

  allocate ( nonzero ( 2*dim_ode ) )
  allocate ( zero ( dim_ode**2 - 2*dim_ode ) )
  allocate ( nonzero_index ( 2*dim_ode ) )
  allocate ( zero_index ( dim_ode**2 - 2*dim_ode ) )

  do ind = 2 , dim_ode-4 , 2
    ind1 = ind*(dim_ode+1)+1
    nonzero ( 2*ind+1 : 2*ind+2 ) &
      = abs ( linear_e ( ind1 : ind1+1 ) )
    nonzero ( 2*ind+3 : 2*ind+4 ) &
      = abs ( linear_e ( ind1+dim_ode : ind1+dim_ode+1 ) )
    ind1 = ind*dim_ode
    ind2 = ind*(dim_ode-2)
    zero ( ind2+1 : ind2+ind ) &
      = abs ( linear_e ( ind1+1 : ind1+ind ) )
    zero ( ind2+ind+1 : ind2+dim_ode-2 ) &
      = abs ( linear_e ( ind1+ind+3 : ind1+dim_ode ) )
    ind1 = ind1+dim_ode
    ind2 = ind2+(dim_ode-2)
    zero ( ind2+1 : ind2+ind ) &
      = abs ( linear_e ( ind1+1 : ind1+ind ) )
    zero ( ind2+ind+1 : ind2+dim_ode-2 ) &
      = abs ( linear_e ( ind1+ind+3 : ind1+dim_ode ) )
  end do

  zero ( 1 : dim_ode-2 ) &
    = abs ( linear_e ( 3 : dim_ode ) )
  zero ( (dim_ode-1) : (dim_ode-2)*2 ) &
    = abs ( linear_e ( dim_ode+3 : 2*dim_ode ) )
  ind1 = (dim_ode-2)*dim_ode
  ind2 = (dim_ode-2)*(dim_ode-2)
  zero ( ind2+1 : ind2+ind ) &
    = abs ( linear_e ( ind1+1 : ind1+ind ) )
  ind1 = ind1+dim_ode
  ind2 = ind2+(dim_ode-2)
  zero ( ind2+1 : ind2+ind ) &
    = abs ( linear_e ( ind1+1 : ind1+ind ) )

  nonzero ( 1 : 2 ) = abs ( linear_e ( 1 : 2 ) )
  nonzero ( 3 : 4 ) = abs ( linear_e ( dim_ode+1 : dim_ode+2 ) )
  nonzero ( 2*dim_ode-3 : 2*dim_ode-2 ) &
    = abs ( linear_e ( dim_ode**2-dim_ode-1 : dim_ode**2-dim_ode ) )
  nonzero ( 2*dim_ode-1 : 2*dim_ode ) &
    = abs ( linear_e ( dim_ode**2-1 : dim_ode**2 ) )

  max_zero = maxval ( zero )
  min_nonzero = minval ( nonzero )

  tp = 0
  tn = 0

  if ( min_nonzero .gt. max_zero ) then
    tp = 1
    tn = 1
  else
    do ind = 1 , dim_ode*2
      where ( nonzero .ge. nonzero(ind) )
        nonzero_index = 1
      elsewhere
        nonzero_index = 0
      end where
      tp_trial = sum(nonzero_index) / size(nonzero_index)
      if ( nonzero(ind) .gt. max_zero ) then
        tn_trial = 1
      else
        where ( zero .lt. nonzero(ind) )
          zero_index = 1
        elsewhere
          zero_index = 0
        end where
        tn_trial = sum(zero_index) / size(zero_index)
      end if
!      write(*,*) 'test:'
!      write(*,*) tp_trial
!      write(*,*) tn_trial
      if ( tp_weight*tp_trial+tn_trial .gt. tp_weight*tp+tn ) then
        tp = tp_trial
        tn = tn_trial
      end if
    end do
  end if

  write(*,*) tp
  write(*,*) tn

  stop

end program main
