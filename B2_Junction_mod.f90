! ----------------------------------------------------------------------
! Programmer(s): Yutong CHEN @XJTU, MER Dr. Konstantin MIKITYUK @PSI
! ----------------------------------------------------------------------
!
! A module for defining ...
! 
! ----------------------------------------------------------------------
module Junction_mod
  !======= Inclusions =========
  use Correlate_mod
  
  !======= Declarations =========
  implicit none
  
  type :: Junc
    !======= Member Variables =========
    ! category
    character(len=20) cate
    ! from/to_pipeid
    character(len=20) pid_f, pid_t
    ! signal_id (optional)
    character(len=20) sigid
    ! junction values
    real(c_double) :: phead = 0.d0
    real(c_double) :: kfac  = 0.d0
  contains
    !======= Member Procedures =========
    procedure :: init_junc ! initialize junction parameters
  end type Junc

contains

  subroutine init_junc(self, cate, pid_f, pid_t, sigid)
    !======= Declarations =========
    implicit none
    !== I/O variables ==
    class(Junc), intent(inout) :: self
    character(len=20), intent(in)  :: cate, pid_f, pid_t
    character(len=20), optional, intent(in) :: sigid
    !======= Internals ============
    ! required parameters
    self%cate  = cate
    self%pid_f = pid_f
    self%pid_t = pid_t
    ! optional parameters
    if(present(sigid)) then
      self%sigid = sigid
    end if
  end subroutine init_junc


end module Junction_mod