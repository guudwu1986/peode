! A module to solve a linear-ODE inverse problem
! via eigenvalue-updating scheme.

! Used by module Linear_Ode_Inverse_mod

! ResidualSumOfSquares:
!   A private subroutine which computes the residue between observation
!   and system curves when eigenvalues of the system is given.
! ResidueNewuoaObjective:
!   A private subroutine which is a wrapper for "ResidualSumOfSquares",
!   to be called by NEWUOA-solver.

module Linear_Ode_Inverse_Eigen_mod

  implicit none
  private

!  public &
!    :: ResidualSumOfSquares
!  public &
!    :: ResidueNewuoaObjective
!  public &
!    :: InverseEigen

  public &
    :: TestResidueNewuoaObjective

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
    integer &
      :: ind1
    double precision , dimension(Dim_Ode+Dim_Time,Dim_Ode) &
      :: basis
    double precision &
      :: temp
    double precision , dimension(:) , allocatable &
      :: atemp
    double precision , dimension(Dim_Ode+Dim_Time) &
      :: y

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

    do ind = 1 , Dim_Ode/2
      do ind1 = 1 , Dim_Time
        temp = exp ( Eigen(ind) * Timepoint(ind1) )
        basis ( ind1 , 2*ind-1 ) &
          = temp * sin ( Eigen(ind+Dim_Ode/2) * Timepoint(ind1) )
        basis ( ind1 , 2*ind ) &
          = temp * cos ( Eigen(ind+Dim_Ode/2) * Timepoint(ind1) )
      end do
    end do

    do ind = 1 , Dim_Ode
      basis ( ind+Dim_Time , ind ) = Ridge_Parameter
    end do

! Compute RSS:
!   QR-decompose basis.
!   For each column(curve) Y of observation,
!   Residue is ||Y||^2 - ||QY||^2.

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

    allocate ( atemp(Dim_Ode+Dim_Time) )
    do ind = 1 , Dim_Ode
      y = 0
      y ( 1 : Dim_Time ) = p_observation ( : , ind )
      Rss = Rss + sum(y**2)
      atemp = 0
      do ind1 = Dim_Ode , 1 , -1
        atemp ( (ind1+1) : (Dim_Ode+Dim_Time) ) &
          = basis ( (ind1+1) : (Dim_Ode+Dim_Time) , ind1 )
        atemp ( ind1 ) = 1
        y = y &
          - tau(ind1) * sum(atemp*y) * atemp
      end do
      Rss = Rss - sum(y**2)
    end do
    deallocate ( atemp )

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

  subroutine TestResidueNewuoaObjective &!{{{
    ( &
      Dim_Ode &
      , Dim_Time &
      , Eigen &
      , Scaling &
      , Timepoint &
      , Observation &
      , Ridge_Parameter &
    )

    implicit none

    integer , intent(in) &
      :: Dim_Ode
    integer , intent(in) &
      :: Dim_Time
    double precision , dimension(Dim_Ode) , intent(inout) &
      :: Eigen
    double precision , dimension(Dim_Ode) , intent(in) , target &
      :: Scaling
    double precision , dimension(Dim_Time) , intent(in) , target &
      :: Timepoint
    double precision , dimension(Dim_Ode*Dim_Time) , intent(in) &
      , target &
      :: Observation
    double precision , intent(in) &
      :: Ridge_Parameter

    double precision &
      :: f

    m_dim_time = Dim_Time
    m_ridge_parameter = Ridge_Parameter

    mp_scaling => Scaling
    mp_timepoint => Timepoint
    mp_observation => Observation

    call ResidueNewuoaObjective ( Dim_Ode , Eigen , f )

    write(*,*) f

    return

  end subroutine TestResidueNewuoaObjective!}}}

  subroutine InverseEigen &!{{{
    ( &
      Dim_Ode &
      , Dim_Time &
      , Eigen &
      , Scaling &
      , Timepoint &
      , Observation &
      , Ridge_Parameter &
    )

    implicit none

    integer , intent(in) &
      :: Dim_Ode
    integer , intent(in) &
      :: Dim_Time
    double precision , dimension(Dim_Ode) , intent(inout) &
      :: Eigen
    double precision , dimension(Dim_Ode) , intent(in) , target &
      :: Scaling
    double precision , dimension(Dim_Time) , intent(in) , target &
      :: Timepoint
    double precision , dimension(Dim_Ode*Dim_Time) , intent(in) &
      , target &
      :: Observation
    double precision , intent(in) &
      :: Ridge_Parameter

    m_dim_time = Dim_Time
    m_ridge_parameter = Ridge_Parameter

    mp_scaling => Scaling
    mp_timepoint => Timepoint
    mp_observation => Observation

    return

  end subroutine InverseEigen!}}}

end module Linear_Ode_Inverse_Eigen_mod
