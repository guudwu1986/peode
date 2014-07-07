! A module to solve a linear-ODE inverse problem

! Public subroutine:
!   InverseNewuoa:
!     Solve a linear-ODE system via NEWUOA-solver.
!   InverseNewuoaSparse:
!     Solve a sparse linear-ODE system via NEWUOA-solver.
!   InverseEigen:
!     Solve a linear-ODE system via eigenvalue-updating.

module Linear_Ode_Inverse_mod

  use Linear_Ode_mod
  use Linear_Ode_Inverse_Eigen_mod

  implicit none

  private

  public :: InverseNewuoa
  public :: InverseNewuoaSparse
  public :: InverseEigen

  integer &
    :: m_dim_ode
  integer &
    :: m_dim_time
  double precision , dimension(:) , pointer &
    :: mp_timepoint
  double precision , dimension(:) , pointer &
    :: mp_observation
  double precision &
    :: m_tol

  integer , dimension(:) , pointer &
    :: mp_nonzero_linear
  integer , dimension(:) , pointer &
    :: mp_nonzero_constant

  contains

  subroutine InverseNewuoaObjective &!{{{
    ( &
      N , X , F &
    )

    implicit none

    integer , intent(in) :: N
    double precision , intent(in) , dimension(N) , target :: X
    double precision , intent(out) :: F

    double precision , dimension(:,:) , pointer &
      :: p_linear
    double precision , dimension(:) , pointer &
      :: p_constant
    double precision , dimension(:) , pointer &
      :: p_initial

    double precision , dimension(:) , allocatable &
      :: ode_result
    integer , dimension(:) , allocatable &
      :: iter
    integer , dimension(:) , allocatable &
      :: info

    p_initial => X ( 1:m_dim_ode )
    p_linear ( 1:m_dim_ode , 1:m_dim_ode ) &
      => X ( (m_dim_ode+1) : (m_dim_ode**2+m_dim_ode) )
    p_constant => X ( (m_dim_ode**2+m_dim_ode+1) : (m_dim_ode**2+2*m_dim_ode) )

    allocate ( ode_result(m_dim_ode*m_dim_time) )
    allocate ( iter(m_dim_time) )
    allocate ( info(m_dim_time) )

    call OdeSolve &
      ( &
        m_dim_ode &
        , m_dim_time &
        , ode_result &
        , p_initial &
        , p_linear &
        , p_constant &
        , mp_timepoint &
        , m_tol &
        , iter &
        , info &
      )

    ode_result = ode_result - mp_observation

    F = sum ( ode_result**2 )

  end subroutine InverseNewuoaObjective!}}}

  subroutine InverseNewuoa &!{{{
    ( &
      Dim_Ode &
      , Dim_Time &
      , Par &
      , Timepoint &
      , Observation &
      , Rhobeg &
      , Rhoend &
      , Maxfun &
      , Tol &
    )

    implicit none

    integer , intent(in) &
      :: Dim_Ode
    integer , intent(in) &
      :: Dim_Time
    double precision , intent(inout) &
      , dimension(Dim_Ode**2+2*Dim_Ode) &
      :: Par
    double precision , intent(in) , dimension(Dim_Time) , target &
      :: Timepoint
    double precision , intent(in) , dimension(Dim_Ode*Dim_Time) &
      , target &
      :: Observation
    double precision , intent(in) &
      :: Rhobeg
    double precision , intent(in) &
      :: Rhoend
    integer , intent(in) &
      :: Maxfun
    double precision , intent(in) &
      :: Tol

    double precision , dimension(:) , allocatable &
      :: work

    integer &
      :: iprint = 0

    interface
      SUBROUTINE NEWUOA (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W,CALFUN)
        IMPLICIT double precision (A-H,O-Z)
        DIMENSION X(*),W(*)
        external CALFUN
      end SUBROUTINE NEWUOA
    end interface

    m_dim_ode = Dim_Ode
    m_dim_time = Dim_Time
    mp_observation => Observation
    mp_timepoint => Timepoint
    m_tol = Tol

    allocate ( work ( size(Par)*(size(Par)+7)/2*13+14 ) )

    call newuoa &
      ( &
        size(Par) &
        , 2*size(Par)+1 &
        , par &
        , Rhobeg &
        , Rhoend &
        , iprint &
        , Maxfun &
        , work &
        , InverseNewuoaObjective &
      )

  end subroutine InverseNewuoa!}}}

  subroutine InverseNewuoaSparseObjective &!{{{
    ( &
      N , X , F &
    )

    implicit none

    integer , intent(in) :: N
    double precision , intent(in) , dimension(N) , target :: X
    double precision , intent(out) :: F

    integer &
      :: dim_linear
    integer &
      :: dim_constant

    double precision , dimension(:) , pointer &
      :: p_linear
    double precision , dimension(:) , pointer &
      :: p_constant
    double precision , dimension(:) , pointer &
      :: p_initial

    double precision , dimension(:) , allocatable &
      :: ode_result
    integer , dimension(:) , allocatable &
      :: iter
    integer , dimension(:) , allocatable &
      :: info

    dim_linear = size(mp_nonzero_linear)
    dim_constant = size(mp_nonzero_constant)

    p_initial => X ( 1:m_dim_ode )
    p_linear => X ( (m_dim_ode+1) : (m_dim_ode+dim_linear) )
    p_constant &
      => X ( &
        (m_dim_ode+dim_linear+1) &
        : (m_dim_ode+dim_linear+dim_constant) &
      )

    allocate ( ode_result(m_dim_ode*m_dim_time) )
    allocate ( iter(m_dim_time) )
    allocate ( info(m_dim_time) )

    call OdeSolveSparse &
      ( &
        m_dim_ode &
        , m_dim_time &
        , dim_linear &
        , dim_constant &
        , ode_result &
        , p_initial &
        , mp_nonzero_linear &
        , p_linear &
        , mp_nonzero_constant &
        , p_constant &
        , mp_timepoint &
        , m_tol &
        , iter &
        , info &
      )

    ode_result = ode_result - mp_observation

    F = sum ( ode_result**2 )

  end subroutine InverseNewuoaSparseObjective!}}}

  subroutine InverseNewuoaSparse &!{{{
    ( &
      Dim_Ode &
      , Dim_Time &
      , Dim_Linear &
      , Dim_Constant &
      , Nonzero_Linear &
      , Nonzero_Constant &
      , Par &
      , Timepoint &
      , Observation &
      , Rhobeg &
      , Rhoend &
      , Maxfun &
      , Tol &
    )

    implicit none

    integer , intent(in) &
      :: Dim_Ode
    integer , intent(in) &
      :: Dim_Time
    integer , intent(in) &
      :: Dim_Linear
    integer , intent(in) &
      :: Dim_Constant
    integer , intent(in) , dimension(Dim_Linear) , target &
      :: Nonzero_Linear
    integer , intent(in) , dimension(Dim_Constant) , target &
      :: Nonzero_Constant
    double precision , intent(inout) &
      , dimension(Dim_Ode+Dim_Linear+Dim_Constant) &
      :: Par
    double precision , intent(in) , dimension(Dim_Time) , target &
      :: Timepoint
    double precision , intent(in) , dimension(Dim_Ode*Dim_Time) &
      , target &
      :: Observation
    double precision , intent(in) &
      :: Rhobeg
    double precision , intent(in) &
      :: Rhoend
    integer , intent(in) &
      :: Maxfun
    double precision , intent(in) &
      :: Tol

    double precision , dimension(:) , allocatable &
      :: work

    integer &
      :: iprint = 0

    interface
      SUBROUTINE NEWUOA (N,NPT,X,RHOBEG,RHOEND,IPRINT,MAXFUN,W,CALFUN)
        IMPLICIT double precision (A-H,O-Z)
        DIMENSION X(*),W(*)
        external CALFUN
      end SUBROUTINE NEWUOA
    end interface

    m_dim_ode = Dim_Ode
    m_dim_time = Dim_Time
    mp_observation => Observation
    mp_timepoint => Timepoint
    m_tol = Tol

    mp_nonzero_linear => Nonzero_Linear
    mp_nonzero_constant => Nonzero_Constant

    allocate ( work ( size(Par)*(size(Par)+7)/2*13+14 ) )

    call newuoa &
      ( &
        size(Par) &
        , 2*size(Par)+1 &
        , par &
        , Rhobeg &
        , Rhoend &
        , iprint &
        , Maxfun &
        , work &
        , InverseNewuoaSparseObjective &
      )

  end subroutine InverseNewuoaSparse!}}}

end module Linear_Ode_Inverse_mod
