module Linear_Ode_Inverse_mod

  use Linear_Ode_mod

  implicit none

  private

  public :: InverseNewuoa

  integer &
    :: m_dim_ode
  integer &
    :: m_dim_time
  double precision , dimension(:) , pointer &
    :: m_timepoint
  double precision , dimension(:,:) , pointer &
    :: m_observation
  double precision &
    :: m_tol

  contains

  subroutine InverseNewuoaObjective &!{{{
    ( &
      N , X , F &
    )

    implicit none

    integer , intent(in) :: N
    double precision , intent(in) , dimension(N) , target :: X
    double precision , intent(out) :: F

    integer &
      :: dim_ode
    integer &
      :: dim_time

    double precision , dimension(:,:) , pointer &
      :: p_linear
    double precision , dimension(:) , pointer &
      :: p_constant
    double precision , dimension(:) , pointer &
      :: p_initial

    double precision , dimension(:,:) , allocatable &
      :: ode_result
    integer , dimension(:) , allocatable &
      :: iter
    integer , dimension(:) , allocatable &
      :: info

    dim_ode = size(m_observation,1)
    dim_time = size(m_observation,2)

    p_initial => X ( 1:dim_ode )
    p_linear ( 1:dim_ode , 1:dim_ode ) &
      => X ( (dim_ode+1) : (dim_ode**2+dim_ode) )
    p_constant => X ( (dim_ode**2+dim_ode+1) : (dim_ode**2+2*dim_ode) )

    allocate ( ode_result(dim_ode,dim_time) )
    allocate ( iter(dim_time) )
    allocate ( info(dim_time) )

    call OdeSolve &
      ( &
        dim_ode &
        , dim_time &
        , ode_result &
        , p_initial &
        , p_linear &
        , p_constant &
        , m_timepoint &
        , m_tol &
        , iter &
        , info &
      )

    ode_result = ode_result - m_observation

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

    m_observation ( 1:Dim_Ode , 1:Dim_Time ) => Observation
    m_timepoint => Timepoint
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

end module Linear_Ode_Inverse_mod
