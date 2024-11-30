! ----------------------------------------------------------------------
! Programmer(s): Yutong CHEN @XJTU, MER Dr. Konstantin MIKITYUK @PSI
! ----------------------------------------------------------------------
!
! A module for defining type(Fuel)
! 
! ----------------------------------------------------------------------
module Fuel_mod
  !======= Inclusions =========
  use Correlate_mod

  !======= Declarations =========
  implicit none
  
  type :: Fuel
    !======= Member Variables =========
    ! id_string, material_id
    character(len=20) :: idstr, matid
    ! material_ids only for multi-layer model
    character(len=20), allocatable :: mats(:)
    ! inner_radius, outer_radius, node_size
    real(c_double) :: ri, ro, dr
    ! interfaces positions, node_sizes of layers
    ! only for multi-layer model
    real(c_double), allocatable :: ris(:), drs(:)
    integer(c_int) :: nr
    ! number of layers
    integer(c_int) :: nlyr = 1
    ! node_quantities of layers, only for multi-layer model
    integer(c_int), allocatable :: nrs(:)
    ! whether the layer is heating material
    logical, allocatable :: isfu(:)
    ! inner_heatflux, outer_heatflux
    real(c_double) :: qflux_i, qflux_o
    ! node_volume, node_center_positions, inner_boundary_positions
    real(c_double), allocatable :: vol(:), rc(:), rb(:)
    ! density, heat_capacity, heat_conductivity, temperature 
    real(c_double), allocatable :: rho(:), cp(:), kpa(:), temp(:)
    ! volumetric_heatsource
    real(c_double), allocatable :: qv(:)
    
  contains
    !======= Member Procedures ========
    procedure :: init_fp1      ! constructor
    procedure :: init_fp2      ! constructor
    procedure :: calc_prop     ! calculate material properties
    procedure :: calc_rhs_fp   ! calculate right-hand-side arrays
    
  end type Fuel

contains
  subroutine init_fp1(self, idstr, matid, ri, ro, nr, isfu, tmp0)
    !======= Declarations =========
    implicit none
    !== I/O variables ==
    class(Fuel), intent(inout) :: self
    character(len=20), intent(in) :: idstr, matid
    real(c_double), intent(in) :: ri, ro
    integer(c_int), intent(in) :: nr
    logical :: isfu
    real(c_double), intent(in) :: tmp0   ! Initial temperature
    !== local variables ==
    integer(c_int) :: i
    !======= Internals ============
    allocate(self%isfu(1))
    self%idstr = idstr
    self%matid = matid
    self%ri = ri
    self%ro = ro
    self%nr = nr
    self%dr = (ro-ri)/real(nr-1, kind=c_double)
    self%isfu(1) = isfu
    ! central hole is assumed to be adiabatic
    self%qflux_i = 0.d0
    ! allocate arrays
    allocate(self%vol(nr) )
    allocate(self%rc(nr)  )
    allocate(self%rb(nr-1))
    do i = 1, nr
      self%rc(i)  = self%ri + (i-1)*self%dr
      if(i<nr) then
        self%rb(i)  = self%rc(i) + self%dr/2.d0
      end if
      if(i==1) then
        self%vol(i) = pi*(self%rb(i)**(2.d0)-self%ri**(2.d0))
      else if(i<nr) then
        self%vol(i) = pi*(self%rb(i)**(2.d0)-self%rb(i-1)**(2.d0))
      else
        self%vol(i) = pi*(self%ro**(2.d0)-self%rb(i-1)**(2.d0))
      end if
    end do
    ! allocate arrays
    allocate(self%rho(nr) )
    allocate(self%cp(nr)  )
    allocate(self%kpa(nr) )
    allocate(self%temp(nr))
    allocate(self%qv(nr)  )
    ! initialize node temperature and qv array
    self%temp = tmp0
    self%qv   = 0.d0
    ! initialize thermal physical properties
    call self%calc_prop()
  end subroutine init_fp1
  
  subroutine init_fp2(self, idstr, mats, ris, nrs, isfu, nlyr, tmp0)
    !======= Declarations =========
    implicit none
    !== I/O variables ==
    class(Fuel), intent(inout) :: self
    character(len=20), intent(in) :: idstr
    character(len=20), intent(in) :: mats(:)
    real(c_double), intent(in) :: ris(:)
    integer(c_int), intent(in) :: nrs(:)
    logical :: isfu(:)
    integer(c_int), intent(in) :: nlyr
    real(c_double), intent(in) :: tmp0   ! Initial temperature
    !== local variables ==
    integer(c_int) :: i, j
    !======= Internals ============
    self%idstr = idstr
    allocate(self%mats(nlyr))
    allocate(self%ris(nlyr+1))
    allocate(self%nrs(nlyr))
    allocate(self%isfu(nlyr))
    self%mats = mats
    self%ris  = ris
    self%ri   = ris(1)
    self%ro   = ris(nlyr+1)
    self%nrs  = nrs
    self%isfu = isfu
    self%nr   = sum(nrs)
    self%nlyr = nlyr
    allocate(self%drs(self%nr))
    j = 1
    do i = 1, nlyr
      if(i==1.or.i==nlyr) then
        self%drs(j:j+nrs(i)-1) = (ris(i+1)-ris(i))/(real(nrs(i),kind=c_double)-0.5d0)
      else
        self%drs(j:j+nrs(i)-1) = (ris(i+1)-ris(i))/(real(nrs(i),kind=c_double))
      end if
      j = j + nrs(i)
    end do
    ! allocate the arrays
    allocate(self%vol(self%nr) )
    allocate(self%rb(self%nr-1))
    self%rb(1) = self%ri + self%drs(1)*0.5d0
    ! inner node boundary position
    do i = 2, self%nr-1
      self%rb(i) = self%rb(i-1) + self%drs(i)
    end do
    ! node volume
    do i = 1, self%nr
      if(i==1) then
        self%vol(i) = pi*(self%rb(i)**(2.d0)-self%ri**(2.d0))
      else if(i==self%nr) then
        self%vol(i) = pi*(self%ro**(2.d0)-self%rb(i-1)**(2.d0))
      else
        self%vol(i) = pi*(self%rb(i)**(2.d0)-self%rb(i-1)**(2.d0))
      end if
    end do
    ! allocate the arrays for properties 
    allocate(self%rho(self%nr) )
    allocate(self%cp(self%nr)  )
    allocate(self%kpa(self%nr) )
    allocate(self%temp(self%nr))
    allocate(self%qv(self%nr)  )
    ! initialize node temperature and qv array
    self%temp = tmp0
    self%qv   = 0.d0
    ! initialize thermal physical properties
    call self%calc_prop()
  end subroutine init_fp2

  subroutine calc_prop(self)
    !======= Declarations =========
    implicit none
    !== I/O variables ==
    class(Fuel), intent(inout) :: self
    !== local variables ==
    integer(c_int) :: i, j, k
    !======= Internals ============
    if(self%nlyr==1) then
      call matpro_s(self%matid, self%temp, self%rho, self%cp, self%kpa)
    else
      j = 1
      do i = 1, self%nlyr
        k = j+self%nrs(i)-1
        call matpro_s(self%mats(i), self%temp(j:k), self%rho(j:k), self%cp(j:k), self%kpa(j:k))
        j = k+1
      end do
    end if
  end subroutine calc_prop

  ! Fuel%qflux_o, Fuel%qv must be assigned in advance
  subroutine calc_rhs_fp(self, rhs)
    !======= Declarations =========
    implicit none
    !== I/O variables ==
    class(Fuel), intent(inout) :: self
    real(c_double), intent(out) :: rhs(:)
    !== local variables ==
    real(c_double) :: qflux(self%nr-1)    ! heatflux at interior cell boundaries
    integer(c_int) :: i
    !======= Internals ============
    ! update the thermal properties
    call self%calc_prop()
    ! update heatflux at interior cell boundaries
    if(self%nlyr==1) then
      ! single-layer model
      do i = 1, self%nr-1
        qflux(i) = (self%temp(i)-self%temp(i+1))/(0.5d0*self%dr*(1.d0/self%kpa(i)+1.d0/self%kpa(i+1)))
      end do
    else
      ! multi-layer model
      do i = 1, self%nr-1
        qflux(i) = (self%temp(i)-self%temp(i+1))/(0.5d0*(self%drs(i)/self%kpa(i)+self%drs(i+1)/self%kpa(i+1)))
      end do
    end if
    ! calculate rhs array
    ! volumetric heat source effect
    rhs = self%qv*self%vol
    do i = 1, self%nr
      ! heat conduction effect
      if(i==1) then
        rhs(i) = rhs(i) + self%qflux_i*2.d0*pi*self%ri - qflux(i)*2.d0*pi*self%rb(i)
      else if(i<self%nr) then
        rhs(i) = rhs(i) + qflux(i-1)*2.d0*pi*self%rb(i-1) - qflux(i)*2.d0*pi*self%rb(i)
      else
        rhs(i) = rhs(i) + qflux(i-1)*2.d0*pi*self%rb(i-1) - self%qflux_o*2.d0*pi*self%ro
      endif
    end do
    ! divide by heat capacity
    rhs = rhs/self%rho/self%cp/self%vol
    ! tips: direct array manipulation could be more efficient than 'do cycle'
  end subroutine calc_rhs_fp

end module Fuel_mod