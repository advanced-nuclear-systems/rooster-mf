! ----------------------------------------------------------------------
! Programmer(s): Yutong CHEN @XJTU, MER Dr. Konstantin MIKITYUK @PSI
! ----------------------------------------------------------------------
!
! A module for defining Reactor module, storing global variables
! 
! ----------------------------------------------------------------------
module Reactor_mod
  !======= Inclusions =========
  use Fluid_mod
  use Solid_mod
  use Core_mod
  use Signal_mod
  
  !======= Declarations =========
  implicit none
  
  type :: Reactor
    !======= Member Variables =========
    ! time variables
    real(c_double) :: tstart, tend, dtout
    ! total ode node number
    integer(c_int) :: neq_tot = 0
    ! total signal number
    integer(c_int) :: nsgnl = 0
    ! fluid field object
    type(Fluid)  :: fluid
    ! solid field object
    type(Solid)  :: solid
    ! core object
    type(Core)   :: core
    ! signal object array
    type(Signal), allocatable :: signals(:)
  contains
    !======= Member Procedures ========
    procedure :: write_to_y
    procedure :: read_from_y
    
  end type Reactor
  ! objects
  ! Reactor_glo is the Reactor-type object storing global data
  type(Reactor) :: Reactor_glo

contains

! ----------------------------------------------------------------------
! write_to_y: write unknowns to yval-vector for initialization
! ----------------------------------------------------------------------
  subroutine write_to_y(self, yval)
    !======= Declarations =========
    implicit none
    !== I/O variables ==
    class(Reactor), intent(inout) :: self
    real(c_double), intent(out)   :: yval(:)
    !== local variables ==
    integer(c_int) :: i, j, k

    !======= Internals ============
    k = 1
    ! fluid solver
    if(solve_fluid) then
      ! mdot to be solved by odes
      do i = 1, self%fluid%njunc
        if(self%fluid%juncs(i)%cate=='jun-head') then
          yval(k) = 0.d0 ! initialization
          k = k + 1
        end if
      end do
      ! flen to be solved by odes
      do i = 1, self%fluid%npipe
        if(self%fluid%pipes(i)%cate=='freelevel') then
          yval(k) = self%fluid%pipes(i)%len
          k = k + 1
        end if
      end do
      ! enth to be solved by odes
      do i = 1, self%fluid%npipe
        if(self%fluid%pipes(i)%cate/='break'.and. &
           self%fluid%pipes(i)%cate/='fill') then
          do j = 1, self%fluid%pipes(i)%nz
            yval(k) = self%fluid%pipes(i)%enth(j)
            k = k + 1
          end do
        end if
      end do
    end if
    ! fuelrod solver
    if(solve_fuelrod) then
      do i = 1, self%solid%nfrod
        do j = 1, self%solid%frods(i)%nf
          yval(k) = self%solid%frods(i)%fuel%temp(j)
          k = k + 1
        end do
        do j = 1, self%solid%frods(i)%nc
          yval(k) = self%solid%frods(i)%clad%temp(j)
          k = k + 1
        end do
      end do
    end if
    ! heat structure solver
    if(solve_htstr) then
      do i = 1, self%solid%nhtstr
        do j = 1, self%solid%htstrs(i)%nr
          yval(k) = self%solid%htstrs(i)%temp(j)
          k = k + 1
        end do
      end do
    end if

  end subroutine write_to_y

! ----------------------------------------------------------------------
! read_from_y: read unknowns from y-vector
! ----------------------------------------------------------------------
  subroutine read_from_y(self, yval)
    !======= Declarations =========
    implicit none
    !== I/O variables ==
    class(Reactor), intent(inout) :: self
    real(c_double), intent(in)    :: yval(:)
    !== local variables ==
    integer(c_int) :: i, j, k

    !======= Internals ============
    k = 1
    ! fluid solver
    if(solve_fluid) then
      ! mdot to be solved by odes
      do i = 1, self%fluid%njunc
        if(self%fluid%juncs(i)%cate=='jun-head') then
          self%fluid%mdots(i)= yval(k) 
          k = k + 1
        end if
      end do
      ! flen to be solved by odes
      do i = 1, self%fluid%npipe
        if(self%fluid%pipes(i)%cate=='freelevel') then
          self%fluid%pipes(i)%len = yval(k)
          k = k + 1
        end if
      end do
      ! enth to be solved by odes
      do i = 1, self%fluid%npipe
        if(self%fluid%pipes(i)%cate/='break'.and. &
           self%fluid%pipes(i)%cate/='fill') then
          do j = 1, self%fluid%pipes(i)%nz
            self%fluid%pipes(i)%enth(j) = yval(k)
            k = k + 1
          end do
        end if
      end do
    end if
    ! fuelrod solver
    if(solve_fuelrod) then
      do i = 1, self%solid%nfrod
        do j = 1, self%solid%frods(i)%nf
          self%solid%frods(i)%fuel%temp(j) = yval(k)
          k = k + 1
        end do
        do j = 1, self%solid%frods(i)%nc
          self%solid%frods(i)%clad%temp(j) = yval(k)
          
          k = k + 1
        end do
      end do
    end if
    ! heat structure solver
    if(solve_htstr) then
      do i = 1, self%solid%nhtstr
        do j = 1, self%solid%htstrs(i)%nr
          self%solid%htstrs(i)%temp(j) = yval(k)
          k = k + 1
        end do
      end do
    end if
    
  end subroutine read_from_y


end module Reactor_mod