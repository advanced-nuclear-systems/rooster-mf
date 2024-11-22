! ----------------------------------------------------------------------
! Programmer(s): Yutong CHEN @XJTU, MER Dr. Konstantin MIKITYUK @PSI
! ----------------------------------------------------------------------
!
! A module for defining type(Thermbc)
! 
! ----------------------------------------------------------------------
module Thermbc_mod
  !======= Inclusions ===========
  use Correlate_mod

  !======= Declarations =========
  implicit none
  
  type :: Thermbc
    !======= Member Variables =========
    ! thermbc name, thermbc type
    character(len=20) :: bcname, bctype
    ! thermbc values array, allocatable
    real(c_double), allocatable :: bcval(:)
    ! coupling pipe id, valid only for coupled boundary
    character(len=20) :: pipeid
    ! coupling pipenode, valid only for coupled boundary
    integer(c_int) :: pipenode
  contains
    !======= Member Procedures ========
    procedure :: init_bc      ! constructor
    procedure :: setbcval     ! set the boundary values
  end type Thermbc

contains
  subroutine init_bc(self, bcname, bctype, bcval, pipeid, pipenode)
    !======= Declarations =========
    implicit none
    !== I/O variables ==
    class(Thermbc), intent(inout) :: self
    character(len=20), intent(in) :: bcname, bctype
    real(c_double), intent(in)    :: bcval(:)
    != Optional input variables =
    character(len=20), intent(in), optional :: pipeid
    integer(c_int), intent(in),    optional :: pipenode
    !======= Internals ============
    self%bcname = bcname
    self%bctype = bctype
    select case (self%bctype)
      case('temperature')  ! 1st (Dirichlet) boundary
        ! array storing temperature
        allocate(self%bcval(1))
        self%bcval(1) = bcval(1)
      case('heatflux')     ! 2rd (Neumann) boundary
        ! array storing heatflux
        allocate(self%bcval(1))
        self%bcval(1) = bcval(1)
      case('convective')   ! 3rd (Robin) boundary
        ! array storing htc, ambient temperature
        allocate(self%bcval(2))
        self%bcval = bcval
      case('coupling')     ! fluid-structure coupling boundary
        ! coupling pipe inforamtion
        self%pipeid   = pipeid
        self%pipenode = pipenode
      case default
        print*, 'no match for bctype = ', self%bctype
        stop 'program stopped at Thermbc%init'
    end select
  end subroutine init_bc
  
  subroutine setbcval(self, bcval)
    !======= Declarations =========
    implicit none
    !== I/O variables ==
    class(Thermbc), intent(inout) :: self
    real(c_double), intent(in)    :: bcval(:)
    !======= Internals ============
    select case (self%bctype)
      case('temperature')  ! 1st (Dirichlet) boundary
        ! array storing temperature
        self%bcval(1) = bcval(1)
      case('heatflux')     ! 2rd (Neumann) boundary
        ! array storing heatflux
        self%bcval(1) = bcval(1)
      case('convective')   ! 3rd (Robin) boundary
        ! array storing htc, ambient temperature
        self%bcval = bcval
      case('coupling')     ! fluid-structure coupling boundary
        ! coupling pipe inforamtion
        stop 'warning, coupling thermbc can not be set'
      case default
        print*, 'no match for bctype = ', self%bctype
        stop 'program stopped at Thermbc%setbcval'
    end select
  end subroutine setbcval
  
end module Thermbc_mod