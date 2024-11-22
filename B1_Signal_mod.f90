! ----------------------------------------------------------------------
! Programmer(s): Yutong CHEN @XJTU, MER Dr. Konstantin MIKITYUK @PSI
! ----------------------------------------------------------------------
!
! A module for defining ...
! 
! ----------------------------------------------------------------------
module Signal_mod
  !======= Inclusions =========
  use Correlate_mod
  
  !======= Declarations =========
  implicit none
  
  type :: Signal
    character(len=20) :: id, cate                ! signal id string
    real(c_double), allocatable :: sgntab(:,:)   ! lookup table
  contains
    procedure :: init_sgnl
  end type Signal
contains
  subroutine init_sgnl(self, id, cate, sgntab)
    !======= Declarations =========
    implicit none
    !== I/O variables ==
    class(Signal), intent(inout) :: self
    character(len=20), intent(in):: id, cate
    real(c_double),     optional :: sgntab(:,:)
    !== Local variables ==
    integer(c_int) :: lentab
    !======= Internals ============
    self%id  = id
    self%cate= cate
    select case(self%cate)
      case('constant')
        continue ! do not need any further action
      case('lookup')
        if(present(sgntab)) then
          lentab = size(sgntab(1,:))
          allocate(self%sgntab(2,lentab))
          self%sgntab = sgntab
        else
          print*, 'laking sgnval/sgntab for signal id = ', self%id
          stop 'program stopped at Signal%init_sgnl'
        end if
      case default
        print*, 'no match for signal type = ', self%cate
        stop 'program stopped at Signal%init_sgnl'
    end select
 
  end subroutine init_sgnl
end module Signal_mod