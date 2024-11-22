! ----------------------------------------------------------------------
! Programmer(s): Yutong CHEN @XJTU, MER Dr. Konstantin MIKITYUK @PSI
! ----------------------------------------------------------------------
!
! A module for defining type(Htstr)
! 
! ----------------------------------------------------------------------
module Htstr_mod
  !======= Inclusions =========
  use Thermbc_mod
  
  !======= Declarations =========
  implicit none
  
  type :: Htstr
    !======= Member Variables =========
    ! id_string, material_id
    character(len=20) :: idstr, matid
    ! inner_radius, outer_radius, node_size
    real(c_double) :: ri, ro, dr
    ! node_quantity
    integer(c_int) :: nr
    ! multiplier
    integer(c_int) :: mltpl
    ! heat_flux on inner/outer boundaries
    real(c_double) :: qflux_i, qflux_o
    ! node_volume, node_center_positions, inner_boundary_positions
    real(c_double), allocatable :: vol(:), rc(:), rb(:)
    ! density, heat_capacity, heat_conductivity, temperature 
    real(c_double), allocatable :: rho(:), cp(:), kpa(:), temp(:)
    ! left_thermbc, right_thermbc
    type(Thermbc)   :: bcleft, bcright

  contains
    !======= Member Procedures =========
    procedure :: init_hs     ! constructor
    procedure :: calc_prop   ! calculate material properties
    procedure :: calc_rhs_hs ! calculate right-hand-side arrays
  
  end type Htstr

contains
  subroutine init_hs(self, idstr, matid, ri, ro, nr, tmp0, bcleft, bcright, mltpl)
    !======= Declarations =========
    implicit none
    !== I/O variables ==
    class(Htstr), intent(inout) :: self
    character(len=20), intent(in) :: idstr, matid
    real(c_double), intent(in) :: ri, ro
    integer(c_int), intent(in) :: nr
    real(c_double), intent(in) :: tmp0   ! Initial temperature
    type(Thermbc)   :: bcleft, bcright
    integer(c_int), intent(in) :: mltpl
    !== local variables ==
    integer(c_int) :: i
    !======= Internals ============
    self%idstr = idstr
    self%matid = matid
    self%ri = ri
    self%ro = ro
    self%nr = nr
    self%dr = (ro-ri)/real(nr-1, kind=c_double)
    ! allocate the arrays
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
    ! allocate the arrays
    allocate(self%rho(nr) )
    allocate(self%cp(nr)  )
    allocate(self%kpa(nr) )
    allocate(self%temp(nr))
    ! initialize node temperature array
    self%temp = tmp0
    ! initialize thermal physical properties
    call self%calc_prop()
    ! initialize the thermbc 
    self%bcleft = bcleft
    self%bcright= bcright
    ! initialize the mltpl
    self%mltpl = mltpl
  
  end subroutine init_hs
  
  subroutine calc_prop(self)
    !======= Declarations =========
    implicit none
    !== I/O variables ==
    class(Htstr), intent(inout) :: self
     
    !======= Internals ============
    call matpro_s(self%matid, self%temp, self%rho, self%cp, self%kpa)

  end subroutine calc_prop
  
  subroutine calc_rhs_hs(self, rhs)
    !======= Declarations =========
    implicit none
    !== I/O variables ==
    class(Htstr), intent(inout) :: self
    real(c_double), intent(out) :: rhs(:)
    !== local variables ==
    real(c_double) :: qflux(self%nr-1)    ! heatflux at inner cell boundaries
    integer(c_int) :: i

    !======= Internals ============
    ! update boundary temperatures if needed
    if(self%bcleft%bctype == 'temperature') then
      self%temp(1) = self%bcleft%bcval(1)
    end if
    if(self%bcright%bctype == 'temperature') then
      self%temp(self%nr) = self%bcright%bcval(1)
    end if
    ! update the thermal properties
    call self%calc_prop()
    ! calculate the heatflux at interior cell boundaries
    do i = 1, self%nr-1
      qflux(i) = (self%temp(i)-self%temp(i+1))/(0.5d0*self%dr*(1.d0/self%kpa(i)+1.d0/self%kpa(i+1)))
    end do
    ! calculate the interior nodes of rhs array
    do i = 2, self%nr-1
      rhs(i) = qflux(i-1)*2.d0*pi*self%rb(i-1) - qflux(i)*2.d0*pi*self%rb(i)
    end do
    ! treat the left thermbc and the left-most node
    select case (self%bcleft%bctype)
      case('temperature')  ! temperature boundary
        rhs(1)  = 0.d0
      case('heatflux')     ! heatflux boundary
        self%qflux_i = self%bcleft%bcval(1)
        rhs(1) = self%qflux_i*2.d0*pi*self%ri - qflux(1)*2.d0*pi*self%rb(1)
      case('convective')   ! convective boundary
        self%qflux_i = self%bcleft%bcval(1)*(self%bcleft%bcval(2)-self%temp(1))
        rhs(1) = self%qflux_i*2.d0*pi*self%ri - qflux(1)*2.d0*pi*self%rb(1)
      case('coupling')     ! coupling boundary
        ! the self%qflux_i should have already been patched
        rhs(1) = self%qflux_i*2.d0*pi*self%ri - qflux(1)*2.d0*pi*self%rb(1)
      case default
        print*, 'no match for left thermbc type = ', self%bcleft%bctype
        stop 'program stopped at Htstr%calc_rhs'
    end select
    ! treat the right thermbc and the right-most node
    select case (self%bcright%bctype)
      case('temperature')  ! temperature boundary
        rhs(self%nr)  = 0.d0
      case('heatflux')     ! heatflux boundary
        self%qflux_o = self%bcright%bcval(1)
        rhs(self%nr) = qflux(self%nr-1)*2.d0*pi*self%rb(self%nr-1)-self%qflux_o*2.d0*pi*self%ro
      case('convective')   ! convective boundary
        self%qflux_o = self%bcright%bcval(1)*(self%temp(self%nr)-self%bcright%bcval(2))
        rhs(self%nr) = qflux(self%nr-1)*2.d0*pi*self%rb(self%nr-1)-self%qflux_o*2.d0*pi*self%ro
      case('coupling')     ! coupling boundary
        ! the self%qflux_o should have already been patched
        rhs(self%nr) = qflux(self%nr-1)*2.d0*pi*self%rb(self%nr-1)-self%qflux_o*2.d0*pi*self%ro
      case default
        print*, 'no match for right thermbc type = ', self%bcright%bctype
        stop 'program stopped at Htstr%calc_rhs'
    end select
    ! divide by heat capacity
    rhs = rhs/self%rho/self%cp/self%vol
    ! tips: direct array manipulation could be more efficient than 'do cycle'
  end subroutine calc_rhs_hs
  
end module Htstr_mod