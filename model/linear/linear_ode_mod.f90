! Module to define a linear ODE model

! OdeSolve:
!   Solve a system.
! Rkf45Derivative:
!   Derivative function called by Rkf45-solver.
! SetParameter:
!   Not used currently.
!   Set member variables, which are used by "Rkf45Derivative".

module Linear_Ode_mod

  implicit none
  private

!  public &
!    :: SetParameter
!  public &
!    :: Rkf45Derivative
  public &
    :: OdeSolve

  ! Member variables!{{{
  double precision , dimension(:,:) , pointer &
    :: m_linear
  double precision , dimension(:) , pointer &
    :: m_constant
!}}}

  contains

!  subroutine SetParameter &!{{{
!    ( &
!      Linear &
!      , Constant &
!    )
!
!    implicit none
!
!    double precision , intent(in) , dimension(:,:) &
!      :: Linear
!    double precision , intent(in) , dimension(:) &
!      :: Constant
!
!    if ( allocated(m_linear) ) then
!      deallocate(m_linear)
!    end if
!    if ( allocated(m_constant) ) then
!      deallocate(m_constant)
!    end if
!
!    allocate ( m_linear(size(Linear,1),size(Linear,2)) )
!    allocate ( m_constant(size(Constant)) )
!    m_linear = Linear
!    m_constant = Constant
!
!    return
!
!  end subroutine SetParameter!}}}

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
      Ode_Result &
      , Initial &
      , Linear &
      , Constant &
      , Time_Point &
      , Tol &
      , Iter &
      , Info &
    )

    double precision , intent(out) , dimension(:,:) &
      :: Ode_Result
    double precision , intent(in) , dimension(:) &
      :: Initial
    double precision , intent(in) , dimension(:,:) , target &
      :: Linear
    double precision , intent(in) , dimension(:) , target &
      :: Constant
    double precision , intent(in) , dimension(:) &
      :: Time_Point
    double precision , intent(in) &
      :: Tol
    integer , intent(out) , dimension(:) &
      :: Iter
    integer , intent(out) , dimension(:) &
      :: Info

    integer &
      :: ind
    double precision , dimension(:) , allocatable &
      :: state
    double precision &
      :: time
    integer &
      :: iflag
    double precision , dimension(:) , allocatable &
      :: work
    integer , dimension(:) , allocatable &
      :: iwork

    interface
      subroutine rkf45 &
        (f,neqn,y,t,tout,relerr,abserr,iflag,work,iwork)
        integer neqn,iflag,iwork(5)
        double precision y(neqn),t,tout,relerr,abserr,work(1)
        external f
      end subroutine rkf45
    end interface

    m_linear => Linear
    m_constant => Constant

    allocate ( state(size(Initial)) )
    state = Initial
    time = time_point(1)

    allocate ( work(3+6*size(Initial)) )
    allocate ( iwork(5) )

    iflag = 1
    iwork(1) = 0

    do ind = 1 , size(Time_Point)!{{{

      call rkf45 ( &
        Rkf45Derivative &
        , size(state) &
        , state &
        , time &
        , time_point(ind) &
        , tol &
        , tol &
        , iflag &
        , work &
        , iwork &
      )

      Ode_Result(:,ind) = state
      time = Time_Point(ind)
      Iter(ind) = iwork(1)
      iwork(1) = 0
      Info(ind) = iflag

    end do!}}}

    return

  end subroutine OdeSolve!}}}

end module Linear_Ode_mod
