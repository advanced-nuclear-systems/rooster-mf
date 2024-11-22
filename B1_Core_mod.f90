! ----------------------------------------------------------------------
! Programmer(s): Yutong CHEN @XJTU, MER Dr. Konstantin MIKITYUK @PSI
! ----------------------------------------------------------------------
!
! A module for defining ...
! 
! ----------------------------------------------------------------------
module Core_mod
  !======= Inclusions =========
  use Correlate_mod
  
  !======= Declarations =========
  implicit none
  
  type :: Core
    !======= Member Variables =========
    real(c_double) :: power0    ! initial core power
    real(c_double) :: lambda    ! decay constant (1/s)
    real(c_double) :: beta      ! delayed neutron fraction
    real(c_double) :: life      ! prompt neutron lifetime
    real(c_double) :: power     ! real time core power
  contains
    !======= Member Procedures ========
    ! constructor when core power is governed by pointkinetics
    procedure :: init_core_pkcore
  end type Core

contains

  subroutine init_core_pkcore(self, lambda, beta, life)
    !======= Declarations =========
    implicit none
    !== I/O variables ==
    class(Core), intent(inout) :: self
    real(c_double) :: lambda, beta, life
    !======= Internals ============
    self%lambda = lambda
    self%beta   = beta
    self%life   = life
  end subroutine init_core_pkcore
end module Core_mod