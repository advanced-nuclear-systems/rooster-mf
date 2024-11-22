! ----------------------------------------------------------------------
! Programmer(s): Yutong CHEN @XJTU, MER Dr. Konstantin MIKITYUK @PSI
! ----------------------------------------------------------------------
!
! A module for defining type(Gasgap)
! 
! ----------------------------------------------------------------------
module Gasgap_mod
  !======= Inclusions ===========
  use Correlate_mod

  !======= Declarations =========
  implicit none
  
  type :: Gasgap
    real(c_double) :: hgap
  contains
    procedure :: init_gg
  end type Gasgap

contains
  subroutine init_gg(self, hgap)
    !======= Declarations =========
    implicit none
    !== I/O variables ==
    class(Gasgap), intent(inout) :: self
    real(c_double) :: hgap
    !======= Internals ============
    self%hgap = hgap
  end subroutine init_gg

end module Gasgap_mod