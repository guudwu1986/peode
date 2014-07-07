! Test program for "ConstructODE" in module
! "Linear_Ode_Inverse_Eigen_mod"

module test
  contains
subroutine ConstructODE &!{{{
  ( &
    Dim_Ode &
    , Dim_Time &
    , Eigen &
    , Timepoint &
    , Observation &
    , Ridge_Parameter &
    , Linear &
    , Initial &
    , Info &
  )

  implicit none

  integer , intent(in) &
    :: Dim_Ode
  integer , intent(in) &
    :: Dim_Time
  double precision , dimension(Dim_Ode) , intent(inout) &
    :: Eigen
  double precision , dimension(Dim_Time) , intent(in) , target &
    :: Timepoint
  double precision , dimension(Dim_Ode*Dim_Time) , intent(in) &
    , target &
    :: Observation
  double precision , intent(in) &
    :: Ridge_Parameter
  double precision , dimension(Dim_Ode*Dim_Ode) , intent(out) &
    , target &
    :: Linear
  double precision , dimension(Dim_Ode) , intent(out) &
    :: Initial
  integer , intent(out) &
    :: Info
! 0: Normal
! 1: LU-factorize Q^T Y
! 2: Inverse Q^T Y
! 3: Inverse R
! 255: QR-factorize basis
! 254: Reconstruct Q
! 253: LU-factorize Q^T Y
! 252: Inverse Q^T Y
! 251: Inverse R

  double precision , dimension(Dim_Ode+Dim_Time,Dim_Ode) &
    :: basis
  double precision , dimension(Dim_Ode+Dim_Time,Dim_Ode) &
    :: q
  double precision , dimension(Dim_Ode,Dim_Ode) &
    :: r

  double precision , dimension(:) , allocatable &
    :: work
  integer , dimension(:) , allocatable &
    :: iwork
  integer &
    :: ind
  double precision , dimension(Dim_Ode) &
    :: tau
  integer &
    :: info_i

  double precision , dimension(:,:) , pointer &
    :: p_linear
  double precision , dimension(:,:) , pointer &
    :: p_observation

  double precision , dimension(Dim_Ode,Dim_Ode) &
    :: xtx
  double precision , dimension(Dim_Ode,Dim_Ode) &
    :: xty

! Main target:
!   Compute similarity transformation operator S and inverse
!   by S = R^{-1} Q^T Y.
!   Compute Linear output by: S Sigma S^{-1}

  p_linear ( 1:Dim_Ode , 1:Dim_Ode ) => Linear
  p_observation ( 1:Dim_Time , 1:Dim_Ode ) => Observation
  Info = 0

! Generate basis:
!   Each column is a basis.
!   A diagonal square matrix is binded to the bottom,
!   used for ridge regression.
!   Generate Sigma, stored in Linear.
!   Generate Initial.

  basis = 0
  p_linear = 0
  Initial = 0

  allocate ( work(Dim_Time) )
  do ind = 1 , Dim_Ode/2
    work = dexp ( Eigen(ind) * Timepoint )
    basis ( 1:Dim_Time , 2*ind-1 ) &
      = work * dsin ( Eigen(ind+Dim_Ode/2) * Timepoint )
    basis ( 1:Dim_Time , 2*ind ) &
      = work * dcos ( Eigen(ind+Dim_Ode/2) * Timepoint )
    p_linear ( 2*ind-1 , 2*ind-1 ) = Eigen(ind)
    p_linear ( 2*ind , 2*ind ) = Eigen(ind)
    p_linear ( 2*ind-1 , 2*ind ) = Eigen(ind+Dim_Ode/2)
    p_linear ( 2*ind , 2*ind-1 ) = -Eigen(ind+Dim_Ode/2)
    Initial ( 2*ind ) = 1
  end do
  deallocate ( work )

  do ind = 1 , Dim_Ode
    basis ( ind+Dim_Time , ind ) = Ridge_Parameter
  end do

  xtx = matmul ( transpose(basis) , basis )!{{{
  xty = matmul ( transpose(basis(1:Dim_Time,:)) , p_observation )

  allocate ( iwork(Dim_Ode) )
  call dgetrf &
    ( &
      Dim_Ode &
      , Dim_Ode &
      , xtx &
      , Dim_Ode &
      , iwork &
      , info_i &
    )
  allocate ( work(1) )
  ind = -1
  call dgetri &
    ( &
      Dim_Ode &
      , xtx &
      , Dim_Ode &
      , iwork &
      , work &
      , ind &
      , info_i &
    )
  ind = work(1)
  deallocate ( work )
  allocate ( work(ind) )
  call dgetri &
    ( &
      Dim_Ode &
      , xtx &
      , Dim_Ode &
      , iwork &
      , work &
      , ind &
      , info_i &
    )
  deallocate ( work )
  deallocate ( iwork )

  xtx = matmul ( xtx , xty )
  xtx = transpose(xtx)

  p_linear = matmul ( xtx , p_linear )
  Initial = matmul ( xtx , Initial )

  allocate ( iwork(Dim_Ode) )
  call dgetrf &
    ( &
      Dim_Ode &
      , Dim_Ode &
      , xtx &
      , Dim_Ode &
      , iwork &
      , info_i &
    )
  allocate ( work(1) )
  ind = -1
  call dgetri &
    ( &
      Dim_Ode &
      , xtx &
      , Dim_Ode &
      , iwork &
      , work &
      , ind &
      , info_i &
    )
  ind = work(1)
  deallocate ( work )
  allocate ( work(ind) )
  call dgetri &
    ( &
      Dim_Ode &
      , xtx &
      , Dim_Ode &
      , iwork &
      , work &
      , ind &
      , info_i &
    )
  deallocate ( work )
  deallocate ( iwork )
  p_linear = matmul ( p_linear , xtx )!}}}

!! QR-decompose basis!{{{
!
!  allocate ( work(1) )
!  ind = -1
!  call dgeqrf ( &
!    Dim_Ode+Dim_Time &
!    , Dim_Ode &
!    , basis &
!    , Dim_Ode+Dim_Time &
!    , tau &
!    , work &
!    , ind &
!    , info_i &
!  )
!  ind = work(1)
!  deallocate ( work )
!  allocate ( work(ind) )
!  call dgeqrf ( &
!    Dim_Ode+Dim_Time &
!    , Dim_Ode &
!    , basis &
!    , Dim_Ode+Dim_Time &
!    , tau &
!    , work &
!    , ind &
!    , info_i &
!  )
!  deallocate ( work )
!
!  if ( info_i .lt. 0 ) then
!    Info = 255
!    stop
!  end if
!
!! Reconstruct Q
!
!  q = basis
!
!  allocate ( work(1) )
!  ind = -1
!  call dorgqr &
!    ( &
!      Dim_Ode+Dim_Time &
!      , Dim_Ode &
!      , Dim_Ode &
!      , q &
!      , Dim_Ode+Dim_Time &
!      , tau &
!      , work &
!      , ind &
!      , info_i &
!    )
!  ind = work(1)
!  deallocate ( work )
!  allocate ( work(ind) )
!  call dorgqr &
!    ( &
!      Dim_Ode+Dim_Time &
!      , Dim_Ode &
!      , Dim_Ode &
!      , q &
!      , Dim_Ode+Dim_Time &
!      , tau &
!      , work &
!      , ind &
!      , info_i &
!    )
!  deallocate ( work )
!
!  if ( info_i .lt. 0 ) then
!    Info = 254
!    stop
!  end if
!
!! Compute Q^T Y, saved in "r"
!
!  r = matmul ( transpose(q(1:Dim_Time,:)) , p_observation )
!
!! Compute (Q^T Y) Sigma (Q^T Y)^{-1}, (Q^T Y)Initial
!
!  Initial = matmul ( r , Initial )
!
!  p_linear = matmul ( r , p_linear )
!
!  allocate ( iwork(Dim_Ode) )
!  call dgetrf &
!    ( &
!      Dim_Ode &
!      , Dim_Ode &
!      , r &
!      , Dim_Ode &
!      , iwork &
!      , info_i &
!    )
!
!  if ( info_i .lt. 0 ) then
!    Info = 253
!    stop
!  end if
!
!  if ( info_i .gt. 0 ) then
!    Info = 1
!  end if
!
!  allocate ( work(1) )
!  ind = -1
!  call dgetri &
!    ( &
!      Dim_Ode &
!      , r &
!      , Dim_Ode &
!      , iwork &
!      , work &
!      , ind &
!      , info_i &
!    )
!  ind = work(1)
!  deallocate ( work )
!  allocate ( work(ind) )
!  call dgetri &
!    ( &
!      Dim_Ode &
!      , r &
!      , Dim_Ode &
!      , iwork &
!      , work &
!      , ind &
!      , info_i &
!    )
!  deallocate ( work )
!  deallocate ( iwork )
!
!  if ( info_i .lt. 0 ) then
!    Info = 252
!    stop
!  end if
!
!  if ( info_i .gt. 0 ) then
!    Info = 2
!  end if
!
!  p_linear = matmul ( p_linear , r )
!
!! Compute inverse of R,
!! then compute S Sigma S^{-1}, S Initial by
!! R^{-1} Q^T Y Sigma (Q^T Y)^{-1} R, R^{-1} Q^T Y Initial.
!
!  r = 0
!
!  do ind = 1 , Dim_Ode
!    r ( 1 : ind , ind ) = basis ( 1 : ind , ind )
!  end do
!
!  p_linear = matmul ( p_linear , r )
!
!  call dtrtri &
!    ( &
!      'U' &
!      , 'N' &
!      , Dim_Ode &
!      , r &
!      , Dim_Ode &
!      , info_i &
!    )
!
!  if ( info_i .lt. 0 ) then
!    Info = 251
!    stop
!  end if
!
!  if ( info_i .gt. 0 ) then
!    Info = 3
!  end if
!
!  p_linear = matmul ( r , p_linear )
!
!  Initial = matmul ( r , Initial )!}}}

  return

end subroutine ConstructODE!}}}
end module test

program main

!use test
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
  :: initial
double precision , dimension(:) , allocatable &
  :: observation
double precision &
  :: ridge_parameter = 1e-10
double precision , dimension(:) , allocatable &
  :: linear_e
double precision , dimension(:) , allocatable &
  :: initial_e
integer &
  :: info

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
read ( NUM_UNIT , * ) eigen ( 1 : (dim_ode/2) )
read ( NUM_UNIT , * ) eigen ( (dim_ode/2+1) : dim_ode )

close ( unit = NUM_UNIT )
!}}}

allocate ( linear_e(dim_ode*dim_ode) )
allocate ( initial_e(dim_ode) )

call ConstructODE &
( &
  dim_ode &
  , dim_time &
  , eigen &
  , timepoint &
  , observation &
  , ridge_parameter &
  , linear_e &
  , initial_e &
  , info &
)

write(*,*) info
write(*,*) maxval ( abs ( linear_e - linear ) )
write(*,*) maxval ( abs ( initial_e - initial ) )

stop

end program main
