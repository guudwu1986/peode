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
      Dim_Ode &
      , Dim_Time &
      , Ode_Result &
      , Initial &
      , Linear &
      , Constant &
      , Time_Point &
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
      :: Time_Point
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
    time = time_point(1)

    iflag = 1
    iwork(1) = 0

    p_ode_result ( 1:Dim_Ode , 1:Dim_Time ) => Ode_Result

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

      p_ode_result(:,ind) = state
      time = Time_Point(ind)
      Iter(ind) = iwork(1)
      iwork(1) = 0
      Info(ind) = iflag

    end do!}}}

    return

  end subroutine OdeSolve!}}}

end module Linear_Ode_mod
