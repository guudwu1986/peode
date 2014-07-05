! A module to solve a linear-ODE inverse problem
! via eigenvalue-updating scheme.

! Used by module Linear_Ode_Inverse_mod

! ResidualSumOfSquares:
!   A private subroutine which computes the residue between observation
!   and system curves when eigenvalues of the system is given.
! ResidueNewuoaObjective:
!   A private subroutine which is a wrapper for "ResidualSumOfSquares",
!   to be called by NEWUOA-solver.
! TestResidueNewuoaObjective:
!   A wrapper for "ResidueNewuoaObjective" testing.
!   Should be commented.
! InverseEigen:
!   A public subroutine to compute an ODE inverse problem
!   via eigenvalue-updating.

module Linear_Ode_Inverse_Eigen_mod

  implicit none
  private

  public &
    :: ResidualSumOfSquares
!  public &
!    :: ResidueNewuoaObjective
!  public &
!    :: InverseEigen

!  public &
!    :: TestResidueNewuoaObjective

  integer &
    :: m_dim_time
  double precision &
    :: m_ridge_parameter

  double precision , dimension(:) , pointer &
    :: mp_timepoint
  double precision , dimension(:) , pointer &
    :: mp_observation
  double precision , dimension(:) , pointer &
    :: mp_scaling

  contains

subroutine ResidualSumOfSquares &!{{{
  ( &
    Dim_Ode &
    , Dim_Time &
    , Eigen &
    , Timepoint &
    , Observation &
    , Ridge_Parameter &
    , Rss &
  )

  implicit none

  integer , intent(in) &
    :: Dim_Ode
  integer , intent(in) &
    :: Dim_Time
  double precision , dimension(Dim_Ode) , intent(in) &
    :: Eigen
  double precision , dimension(Dim_Time) , intent(in) &
    :: Timepoint
  double precision , dimension(Dim_Ode*Dim_Time) , intent(in) &
    , target &
    :: Observation
  double precision , intent(in) &
    :: Ridge_Parameter
  double precision , intent(out) &
    :: Rss

  integer &
    :: ind
  double precision , dimension(Dim_Ode+Dim_Time,Dim_Ode) &
    :: basis
  double precision , &
    dimension ( (Dim_Ode+Dim_Time) , Dim_Ode ) &
    :: q
!  double precision &
!    :: zero_rss

  double precision , dimension(:) , allocatable &
    :: work
  integer &
    :: lwork
  double precision , dimension(Dim_Ode) &
    :: tau
  integer &
    :: info

  double precision , dimension(:,:) , pointer &
    :: p_observation

  p_observation ( 1:Dim_Time , 1:Dim_Ode ) => Observation

! Generate basis:
!   Each column is a basis.
!   A diagonal square matrix is binded to the bottom,
!   used for ridge regression.

  basis = 0

  allocate ( work(Dim_Time) )
  do ind = 1 , Dim_Ode/2
    work = dexp ( Eigen(ind) * Timepoint )
    basis ( : , 2*ind-1 ) &
      = work * dsin ( Eigen(ind+Dim_Ode/2) * Timepoint )
    basis ( : , 2*ind ) &
      = work * dcos ( Eigen(ind+Dim_Ode/2) * Timepoint )
  end do
  deallocate ( work )

  do ind = 1 , Dim_Ode
    basis ( ind+Dim_Time , ind ) = Ridge_Parameter
  end do

! QR-decompose basis

  Rss = 0

  allocate ( work(1) )
  lwork = -1
  call dgeqrf ( &
    Dim_Ode+Dim_Time &
    , Dim_Ode &
    , basis &
    , Dim_Ode+Dim_Time &
    , tau &
    , work &
    , lwork &
    , info &
  )
  lwork = work(1)
  deallocate ( work )
  allocate ( work(lwork) )
  call dgeqrf ( &
    Dim_Ode+Dim_Time &
    , Dim_Ode &
    , basis &
    , Dim_Ode+Dim_Time &
    , tau &
    , work &
    , lwork &
    , info &
  )
  deallocate ( work )

! Reconstruct Q

  q = basis

  allocate ( work(1) )
  lwork = -1
  call dorgqr &
    ( &
      Dim_Ode+Dim_Time &
      , Dim_Ode &
      , Dim_Ode &
      , q &
      , Dim_Ode+Dim_Time &
      , tau &
      , work &
      , lwork &
      , info &
    )
  lwork = work(1)
  deallocate ( work )
  allocate ( work(lwork) )
  call dorgqr &
    ( &
      Dim_Ode+Dim_Time &
      , Dim_Ode &
      , Dim_Ode &
      , q &
      , Dim_Ode+Dim_Time &
      , tau &
      , work &
      , lwork &
      , info &
    )
  deallocate ( work )

! Compute RSS:
!   Residue is ||Y||^2 - ||Q^T Y||^2.

  Rss = sum ( p_observation**2 ) &
    - sum ( matmul ( transpose(q(1:Dim_Time,:)) , p_observation ) **2 )

! Uncomment following block to check
!   cases with 0-eigenvalue with 2 time-points.

!  zero_rss = 0
!  do ind = 1 , Dim_Ode
!    zero_rss = zero_rss &
!      + sum ( &
!          ( p_observation(:,ind) &
!            - sum(p_observation(:,ind))/Dim_Time &
!          ) **2 &
!        )
!  end do
!  write(*,*) zero_rss

  return

end subroutine ResidualSumOfSquares!}}}

subroutine ResidueNewuoaObjective &!{{{
  ( &
    N , X , F &
  )

  implicit none

  integer , intent(in) :: N
  double precision , intent(in) , dimension(N) :: X
  double precision , intent(out) :: F

  call ResidualSumOfSquares &
    ( &
      N &
      , m_dim_time &
      , X * mp_scaling &
      , mp_timepoint &
      , mp_observation &
      , m_ridge_parameter &
      , F &
    )

  return

end subroutine ResidueNewuoaObjective!}}}

!  subroutine TestResidueNewuoaObjective &!{{{
!    ( &
!      Dim_Ode &
!      , Dim_Time &
!      , Eigen &
!      , Scaling &
!      , Timepoint &
!      , Observation &
!      , Ridge_Parameter &
!    )
!
!    implicit none
!
!    integer , intent(in) &
!      :: Dim_Ode
!    integer , intent(in) &
!      :: Dim_Time
!    double precision , dimension(Dim_Ode) , intent(inout) &
!      :: Eigen
!    double precision , dimension(Dim_Ode) , intent(in) , target &
!      :: Scaling
!    double precision , dimension(Dim_Time) , intent(in) , target &
!      :: Timepoint
!    double precision , dimension(Dim_Ode*Dim_Time) , intent(in) &
!      , target &
!      :: Observation
!    double precision , intent(in) &
!      :: Ridge_Parameter
!
!    double precision &
!      :: f
!
!    m_dim_time = Dim_Time
!    m_ridge_parameter = Ridge_Parameter
!
!    mp_scaling => Scaling
!    mp_timepoint => Timepoint
!    mp_observation => Observation
!
!    call ResidueNewuoaObjective ( Dim_Ode , Eigen , f )
!
!    write(*,*) f
!
!    return
!
!  end subroutine TestResidueNewuoaObjective!}}}

!subroutine ConstructODE &!{{{
!  ( &
!    Dim_Ode &
!    , Dim_Time &
!    , Eigen &
!    , Scaling &
!    , Timepoint &
!    , Observation &
!    , Ridge_Parameter &
!  )
!
!  implicit none
!
!  return
!
!end subroutine ConstructODE!}}}

!subroutine InverseEigen &!{{{
!  ( &
!    Dim_Ode &
!    , Dim_Time &
!    , Eigen &
!    , Scaling &
!    , Timepoint &
!    , Observation &
!    , Ridge_Parameter &
!    , Rhobeg &
!    , Rhoend &
!    , Maxfun &
!  )
!
!  implicit none
!
!  integer , intent(in) &
!    :: Dim_Ode
!  integer , intent(in) &
!    :: Dim_Time
!  double precision , dimension(Dim_Ode) , intent(inout) &
!    :: Eigen
!  double precision , dimension(Dim_Ode) , intent(in) , target &
!    :: Scaling
!  double precision , dimension(Dim_Time) , intent(in) , target &
!    :: Timepoint
!  double precision , dimension(Dim_Ode*Dim_Time) , intent(in) &
!    , target &
!    :: Observation
!  double precision , intent(in) &
!    :: Ridge_Parameter
!  double precision , intent(in) &
!    :: Rhobeg
!  double precision , intent(in) &
!    :: Rhoend
!  integer , intent(in) &
!    :: Maxfun
!
!  integer &
!    :: iprint = 3
!
!  double precision , &
!    dimension ( (15*Dim_Ode+97)*Dim_Ode/2+14 ) &
!    :: work
!
!  interface
!    SUBROUTINE NEWUOA (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W,CALFUN)
!      IMPLICIT double precision (A-H,O-Z)
!      DIMENSION X(*),W(*)
!      external CALFUN
!    end SUBROUTINE NEWUOA
!  end interface
!
!  m_dim_time = Dim_Time
!  m_ridge_parameter = Ridge_Parameter
!
!  mp_scaling => Scaling
!  mp_timepoint => Timepoint
!  mp_observation => Observation
!
!  call newuoa &
!    ( &
!      Dim_Ode &
!      , 2*Dim_Ode+1 &
!      , Eigen &
!      , Rhobeg &
!      , Rhoend &
!      , iprint &
!      , Maxfun &
!      , work &
!      , ResidueNewuoaObjective &
!    )
!
!  return
!
!end subroutine InverseEigen!}}}

end module Linear_Ode_Inverse_Eigen_mod
