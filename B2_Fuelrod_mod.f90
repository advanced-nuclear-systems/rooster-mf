! ----------------------------------------------------------------------
! Programmer(s): Yutong CHEN @XJTU, MER Dr. Konstantin MIKITYUK @PSI
! ----------------------------------------------------------------------
!
! A module for defining Fuelrod 
! 
! ----------------------------------------------------------------------
module Fuelrod_mod
  !======= Inclusions =========
  use Fuel_mod
  use Gasgap_mod
  use Clad_mod

  !======= Declarations =========
  implicit none
  
  type :: Fuelrod
    !======= Member Variables =========
    ! id_string, pipe_id
    character(len=20) :: idstr, pipeid
    ! pipenode, multiplier
    integer(c_int) :: pipenode, mltpl
    ! node number of fuel_pellet, clad, and total
    integer(c_int) :: nf, nc, ntot
    ! heat_flux at clad outer surface, volumetric pellet heat source
    real(c_double) :: qflux_co, qv
    ! radial and axial volumetric heat source factor, axial length
    real(c_double) :: kr, kz, dz
    ! here maybe the member objects can use type names
    ! fuel pellet
    type(Fuel)   :: fuel
    ! gas gap
    type(Gasgap) :: gasgap
    ! cladding
    type(Clad)   :: clad
    
  contains
    !======= Member Procedures ========
    ! here the procedures cannot use procedure names in Fuel/Gasgap/Clad modules
    procedure :: init_frd        ! constructor
    procedure :: comp_rhs_frd    ! calculate the rhs array
  end type Fuelrod

contains
  subroutine init_frd(self, idstr, fuel_obj, gasgap_obj, clad_obj, kr, kz, &
                       pipeid, pipenode, mltpl, dz)
    !======= Declarations =========
    implicit none
    !== I/O variables ==
    class(Fuelrod), intent(inout) :: self
    ! id string, id of corresponding pipe
    character(len=20), intent(in) :: idstr, pipeid
    ! index of pipe node, multiplier
    integer(c_int), intent(in) :: pipenode, mltpl
    ! radial/axial volumetric power density factor, axial length
    real(c_double) :: kr, kz, dz
    ! fuel pellet object
    class(Fuel), intent(in) :: fuel_obj
    ! gas gap object
    class(Gasgap), intent(in) :: gasgap_obj
    ! cladding object
    class(Clad), intent(in)   :: clad_obj
    !======= Internals ============
    self%idstr   = idstr
    self%fuel    = fuel_obj
    self%gasgap  = gasgap_obj
    self%clad    = clad_obj
    self%kr      = kr
    self%kz      = kz
    self%pipeid  = pipeid
    self%pipenode= pipenode
    self%mltpl   = mltpl
    self%dz      = dz
    self%nf      = self%fuel%nr
    self%nc      = self%clad%nr
    self%ntot    = self%fuel%nr + self%clad%nr
  end subroutine init_frd
  
  ! Fuelrod%qflux_co, Fuelrod%qv must be assigned in advance
  subroutine comp_rhs_frd(self, rhs)
    !======= Declarations =========
    implicit none
    !== I/O variables ==
    class(Fuelrod), intent(inout) :: self
    ! rhs array
    real(c_double), intent(out)   :: rhs(:)
    !== local variables ==
    ! fuel->clad heat flux at gas gap
    real(c_double) :: qflux_fc
    ! rhs arrays for fuel and clad
    real(c_double), allocatable :: rhs_f(:), rhs_c(:)
    integer(c_int) :: i, j, k
    !======= Internals ============
    ! check if input is correct
    if(size(rhs)/=self%ntot) then
      stop 'Array size mismatch in Fuelrod%comp_rhs'
    end if
    ! allocate the rhs arrays
    allocate(rhs_f(self%nf))
    allocate(rhs_c(self%nc))
    ! calculate heat flux at gas gap
    qflux_fc = self%gasgap%hgap*(self%fuel%temp(self%nf)-self%clad%temp(1))
    ! patch the volumetric heat sources
    if(self%fuel%nlyr==1) then
      if(self%fuel%isfu(1)) then
        self%fuel%qv = self%qv
      end if
    else
      do i = 1, self%fuel%nlyr
        if(self%fuel%isfu(i)) then
          if(i==1) then
            j = 1
          else
            j = sum(self%fuel%nrs(1:i-1)) + 1
          end if
          k = self%fuel%nrs(i)
          self%fuel%qv(j:j+k-1) = self%qv
        end if
      end do
    end if
    
    ! patch the heat fluxes
    self%fuel%qflux_o = qflux_fc
    self%clad%qflux_i = qflux_fc
    self%clad%qflux_o = self%qflux_co
    ! calculate the rhs arrays
    call self%fuel%calc_rhs_fp(rhs_f)
    call self%clad%calc_rhs_cd(rhs_c)
    ! compose the rhs arrays
    rhs = (/rhs_f, rhs_c/)
    ! deallocate the rhs arrays
    deallocate(rhs_f)
    deallocate(rhs_c)
  end subroutine comp_rhs_frd
end module Fuelrod_mod