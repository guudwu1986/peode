! Module to solve a linear ODE model

! OdeSolve:
!   Solve a full system.
! OdeSolveSparse:
!   Solve a sparse system.
! Rkf45Derivative, Rkf45DerivativeSparse:
!   Derivative function called by Rkf45-solver.

module Linear_Ode_mod

  implicit none
  private

  public &
    :: OdeSolve
  public &
    :: OdeSolveSparse

! Member variables!{{{
  double precision , dimension(:,:) , pointer &
    :: m_linear
  double precision , dimension(:) , pointer &
    :: m_constant

  integer &
    :: m_dim_ode
  integer , dimension(:) , pointer &
    :: m_nonzero_linear
  integer , dimension(:) , pointer &
    :: m_nonzero_constant
  double precision , dimension(:) , pointer &
    :: m_linear_sparse
  double precision , dimension(:) , pointer &
    :: m_constant_sparse
!}}}

  contains

  subroutine Rkf45Derivative &!{{{
    ( &
      Time &
      , State &
      , Derivative &
    )

    implicit none

    double precision , intent(in) &
      :: Time
    double precision , intent(in) , dimension(:) &
      :: State
    double precision , intent(out) , dimension(:) &
      :: Derivative

    Derivative = matmul ( m_linear , State )

    if ( size(m_constant) .gt. 0 ) then
      Derivative = Derivative + m_constant
    end if

    return

  end subroutine Rkf45Derivative!}}}

  subroutine OdeSolve &!{{{
    ( &
      Dim_Ode &
      , Dim_Time &
      , Ode_Result &
      , Initial &
      , Linear &
      , Constant &
      , Timepoint &
      , Tol &
      , Iter &
      , Info &
    )

    integer , intent(in) &
      :: Dim_Ode
    integer , intent(in) &
      :: Dim_Time
    double precision , intent(out) , dimension(Dim_Ode*Dim_Time) &
      , target &
      :: Ode_Result
    double precision , intent(in) , dimension(Dim_Ode) &
      :: Initial
    double precision , intent(in) , dimension(Dim_Ode*Dim_Ode) , target &
      :: Linear
    double precision , intent(in) , dimension(Dim_Ode) , target &
      :: Constant
    double precision , intent(in) , dimension(Dim_Time) &
      :: Timepoint
    double precision , intent(in) &
      :: Tol
    integer , intent(out) , dimension(Dim_Time) &
      :: Iter
    integer , intent(out) , dimension(Dim_Time) &
      :: Info

    integer &
      :: ind
    double precision , dimension(Dim_Ode) &
      :: state
    double precision &
      :: time
    integer &
      :: iflag
    double precision , dimension(3+6*Dim_Ode) &
      :: work
    integer , dimension(5) &
      :: iwork

    double precision , dimension(:,:) , pointer &
      :: p_ode_result

    interface
      subroutine rkf45 &
        (f,neqn,y,t,tout,relerr,abserr,iflag,work,iwork)
        integer neqn,iflag,iwork(5)
        double precision y(neqn),t,tout,relerr,abserr,work(1)
        external f
      end subroutine rkf45
    end interface

    m_linear ( 1:Dim_Ode , 1:Dim_Ode ) => Linear
    m_constant => Constant

    state = Initial
    time = Timepoint(1)

    iflag = 1
    iwork(1) = 0

    p_ode_result ( 1:Dim_Ode , 1:Dim_Time ) => Ode_Result

    do ind = 1 , size(Timepoint)!{{{

      call rkf45 ( &
        Rkf45Derivative &
        , size(state) &
        , state &
        , time &
        , Timepoint(ind) &
        , Tol &
        , Tol &
        , iflag &
        , work &
        , iwork &
      )

      p_ode_result(:,ind) = state
      time = Timepoint(ind)
      Iter(ind) = iwork(1)
      iwork(1) = 0
      Info(ind) = iflag

    end do!}}}

    return

  end subroutine OdeSolve!}}}

  subroutine Rkf45DerivativeSparse &!{{{
    ( &
      Time &
      , State &
      , Derivative &
    )

    implicit none

    double precision , intent(in) &
      :: Time
    double precision , intent(in) , dimension(:) &
      :: State
    double precision , intent(out) , dimension(:) &
      :: Derivative

    integer &
      :: ind
    integer &
      :: row
    integer &
      :: col

    Derivative = 0
    do ind = 1 , size(m_nonzero_linear)
      row = mod ( m_nonzero_linear(ind)-1 , m_dim_ode ) + 1
      col = ( m_nonzero_linear(ind) - row ) / m_dim_ode + 1
      Derivative(row) &
        = Derivative(row) + m_linear_sparse(ind) * State(col)
    end do

    do ind = 1 , size(m_nonzero_constant)
      Derivative(m_nonzero_constant(ind)) &
        = Derivative(m_nonzero_constant(ind)) + m_constant_sparse(ind)
    end do

    return

  end subroutine Rkf45DerivativeSparse!}}}

  subroutine OdeSolveSparse &!{{{
    ( &
      Dim_Ode &
      , Dim_Time &
      , Dim_Linear &
      , Dim_Constant &
      , Ode_Result &
      , Initial &
      , Nonzero_Linear &
      , Linear &
      , Nonzero_Constant &
      , Constant &
      , Timepoint &
      , Tol &
      , Iter &
      , Info &
    )

    integer , intent(in) &
      :: Dim_Ode
    integer , intent(in) &
      :: Dim_Time
    integer , intent(in) &
      :: Dim_Linear
    integer , intent(in) &
      :: Dim_Constant
    double precision , intent(out) , dimension(Dim_Ode*Dim_Time) &
      , target &
      :: Ode_Result
    double precision , intent(in) , dimension(Dim_Ode) &
      :: Initial
    integer , intent(in) , dimension(Dim_Linear) , target &
      :: Nonzero_Linear
    double precision , intent(in) , dimension(Dim_Linear) , target &
      :: Linear
    integer , intent(in) , dimension(Dim_Constant) , target &
      :: Nonzero_Constant
    double precision , intent(in) , dimension(Dim_Constant) , target &
      :: Constant
    double precision , intent(in) , dimension(Dim_Time) &
      :: Timepoint
    double precision , intent(in) &
      :: Tol
    integer , intent(out) , dimension(Dim_Time) &
      :: Iter
    integer , intent(out) , dimension(Dim_Time) &
      :: Info

    integer &
      :: ind
    double precision , dimension(Dim_Ode) &
      :: state
    double precision &
      :: time
    integer &
      :: iflag
    double precision , dimension(3+6*Dim_Ode) &
      :: work
    integer , dimension(5) &
      :: iwork

    double precision , dimension(:,:) , pointer &
      :: p_ode_result

    interface
      subroutine rkf45 &
        (f,neqn,y,t,tout,relerr,abserr,iflag,work,iwork)
        integer neqn,iflag,iwork(5)
        double precision y(neqn),t,tout,relerr,abserr,work(1)
        external f
      end subroutine rkf45
    end interface

    m_dim_ode = Dim_Ode
    m_nonzero_linear => Nonzero_Linear
    m_nonzero_constant => Nonzero_Constant
    m_linear_sparse => Linear
    m_constant_sparse => Constant

    state = Initial
    time = Timepoint(1)

    iflag = 1
    iwork(1) = 0

    p_ode_result ( 1:Dim_Ode , 1:Dim_Time ) => Ode_Result

    do ind = 1 , size(Timepoint)!{{{

      call rkf45 ( &
        Rkf45DerivativeSparse &
        , size(state) &
        , state &
        , time &
        , Timepoint(ind) &
        , tol &
        , tol &
        , iflag &
        , work &
        , iwork &
      )

      p_ode_result(:,ind) = state
      time = Timepoint(ind)
      Iter(ind) = iwork(1)
      iwork(1) = 0
      Info(ind) = iflag

    end do!}}}

    return

  end subroutine OdeSolveSparse!}}}

end module Linear_Ode_mod
