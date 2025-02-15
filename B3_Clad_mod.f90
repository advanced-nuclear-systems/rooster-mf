! ----------------------------------------------------------------------
! Programmer(s): Yutong CHEN @XJTU, MER Dr. Konstantin MIKITYUK @PSI
! ----------------------------------------------------------------------
!
! A module for defining ...
! 
! ----------------------------------------------------------------------
module Clad_mod
  !======= Inclusions =========
  use Correlate_mod

  !======= Declarations =========
  implicit none
  
  type :: Clad
    !======= Member Variables =========
    ! id_string and material_id
    character(len=20) :: idstr, matid
    ! inner_radius, outer_radius, node_size
    real(c_double) :: ri, ro, dr
    ! node_quantity
    integer(c_int) :: nr
    ! inner_heatflux, outer_heatflux
    real(c_double) :: qflux_i, qflux_o
    ! node_volume, node_center_positions, inner_boundary_positions
    real(c_double), allocatable :: vol(:), rc(:), rb(:)
   ! density, heat_capacity, heat_conductivity, temperature 
    real(c_double), allocatable :: rho(:), cp(:), kpa(:), temp(:)

  contains
    !======= Member Procedures ========
    procedure :: init_cd       ! constructor
    procedure :: calc_prop     ! calculate material properties
    procedure :: calc_rhs_cd   ! calculate right-hand-side arrays

  end type Clad

contains
  subroutine init_cd(self, idstr, matid, ri, ro, nr, tmp0)
    !======= Declarations =========
    implicit none
    !== I/O variables ==
    class(Clad), intent(inout) :: self
    character(len=20), intent(in) :: idstr, matid
    real(c_double), intent(in) :: ri, ro
    integer(c_int), intent(in) :: nr
    real(c_double), intent(in) :: tmp0   ! Initial temperature
    !== local variables ==
    integer(c_int) :: i
    !======= Internals ============
    self%idstr = idstr
    self%matid = matid
    self%ri = ri
    self%ro = ro
    self%nr = nr
    self%dr = (ro-ri)/real(nr-1, kind=c_double)
    ! allocate arrays
    allocate(self%vol(nr) )
    allocate(self%rc(nr)  )
    allocate(self%rb(nr-1))
    do i = 1, nr
      self%rc(i)  = self%ri + (i-1)*self%dr
      if(i<nr) then
        self%rb(i)  = self%rc(i) + self%dr/2.d0
      endif
      if(i==1) then
        self%vol(i) = pi*(self%rb(i)**(2.d0)-self%ri**(2.d0))
      else if(i<nr) then
        self%vol(i) = pi*(self%rb(i)**(2.d0)-self%rb(i-1)**(2.d0))
      else
        self%vol(i) = pi*(self%ro**(2.d0)-self%rb(i-1)**(2.d0))
      endif
    end do
    ! allocate arrays
    allocate(self%rho(nr) )
    allocate(self%cp(nr)  )
    allocate(self%kpa(nr) )
    allocate(self%temp(nr))
    ! initialize node temperature array
    self%temp = tmp0
    ! initialize thermal physical properties
    call self%calc_prop()
  end subroutine init_cd
  
  subroutine calc_prop(self)
    !======= Declarations =========
    implicit none
    !== I/O variables ==
    class(Clad), intent(inout) :: self
    !======= Internals ============
    call matpro_s(self%matid, self%temp, self%rho, self%cp, self%kpa)
  end subroutine calc_prop

  ! Clad%qflux_i, Clad%qflux_o must be assigned in advance
  subroutine calc_rhs_cd(self, rhs)
    !======= Declarations =========
    implicit none
    !== I/O variables ==
    class(Clad), intent(inout) :: self
    real(c_double), intent(out) :: rhs(:)
    !== local variables ==
    real(c_double) :: qflux(self%nr-1)    ! heatflux at interior cell boundaries
    integer(c_int) :: i
    !======= Internals ============
    ! update the thermal properties
    call self%calc_prop()
    ! heatflux at interior cell boundaries
    do i = 1, self%nr-1
      qflux(i) = (self%temp(i)-self%temp(i+1))/(0.5d0*self%dr*(1.d0/self%kpa(i)+1.d0/self%kpa(i+1)))
    end do
    ! calculate rhs array
    do i = 1, self%nr
      ! heat conduction effect
      if(i==1) then
        rhs(i) = self%qflux_i*2.d0*pi*self%ri - qflux(i)*2.d0*pi*self%rb(i)
      else if(i<self%nr) then
        rhs(i) = qflux(i-1)*2.d0*pi*self%rb(i-1) - qflux(i)*2.d0*pi*self%rb(i)
      else
        rhs(i) = qflux(i-1)*2.d0*pi*self%rb(i-1) - self%qflux_o*2.d0*pi*self%ro
      endif
    end do
    ! divide by heat capacity
    rhs = rhs/self%rho/self%cp/self%vol
    ! tips: direct array manipulation could be more efficient than 'do cycle'
  end subroutine calc_rhs_cd

end module Clad_mod