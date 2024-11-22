! ----------------------------------------------------------------------
! Programmer(s): Yutong CHEN @XJTU, MER Dr. Konstantin MIKITYUK @PSI
! ----------------------------------------------------------------------
!
! A module for defining ...
! 
! ----------------------------------------------------------------------
module Pipe_mod
  !======= Inclusions =========
  use Correlate_mod
  
  !======= Declarations =========
  implicit none
  
  type :: Pipe
    !======= Member Variables =========
    ! required input parameters-----------------------------------------
    ! id_string, material_id, category
    character(len=20) idstr, matid, cate, sigid
    ! optional input parameters-----------------------------------------
    ! hydraulic_dameter, pipe_length, direction, flow_cross_section_area
    real(c_double) :: dhyd, len, dir, areaz
    ! pitch-to-diameter_ratio, bundle_surface-to-total_surface_ratio
    real(c_double) :: p2d, sb2st
    ! number of pipe nodes
    integer(c_int) :: nz
    ! ------------------------------------------------------------------
    ! cell_volume, cell axial length
    real(c_double) :: vol, dz
    ! inlet/outlet enthalpy---------------------------------------------
    real(c_double) :: mh_i, mh_o
    ! pipe mass-flowrate
    real(c_double) :: mdot = 0.d0
    ! calculated array parameters---------------------------------------
    ! cell pressure
    real(c_double), allocatable :: pval(:)
    ! wall-to-fluid net heating power in pipe cell
    real(c_double), allocatable :: qpow(:)
    ! enthalpy, temperature, gas quality
    real(c_double), allocatable :: enth(:), temp(:), xe(:)
    ! density, heat_capacity/conductivity, dynamic/kinetic_viscosity
    real(c_double), allocatable :: rho(:), cp(:), kpa(:), miu(:), niu(:)
    ! dimensionless parameters
    real(c_double), allocatable :: re(:), pr(:), pe(:)
    ! heat transfer mode 
    integer(c_int), allocatable :: htcmod(:)
  contains
    !======= Member Procedures =========
    procedure :: init_pipe     ! calculate material properties
    procedure :: calc_prop     ! calculate material properties
    procedure :: calc_rhs_pipe ! calculate right-hand-side arrays
    
  end type Pipe

contains

  subroutine init_pipe(self, id, mat, cate, nz, inp)
                                          ! dhyd, len, areaz, dir, p2d, sb2st, temp, pout, nz,
    !======= Declarations =========
    implicit none
    !== I/O variables ==
    class(Pipe), intent(inout) :: self
    character(len=20), intent(in) :: id, mat, cate
    integer(c_int), intent(in) :: nz
    real(c_double), intent(in) :: inp(:)    
    !== local variables ==

    !======= Internals ============
    ! required parameters
    self%idstr = id
    self%matid = mat
    self%cate  = cate
    self%nz    = nz
    ! allocatable arrays
    allocate(self%pval(self%nz))
    allocate(self%qpow(self%nz))
    allocate(self%enth(self%nz))
    allocate(self%temp(self%nz))
    allocate(self%xe(self%nz))
    allocate(self%rho(self%nz))
    allocate(self%cp(self%nz))
    allocate(self%kpa(self%nz))
    allocate(self%miu(self%nz))
    allocate(self%niu(self%nz))
    allocate(self%re(self%nz))
    allocate(self%pr(self%nz))
    allocate(self%pe(self%nz))
    allocate(self%htcmod(self%nz))
    self%htcmod = 0
    ! assign optional parameters
    select case (cate)
      case('round')
        ! inp = (/dhyd, len, areaz, dir, temp, pval/)
        self%dhyd = inp(1)
        self%len  = inp(2)
        self%areaz= inp(3)
        self%dir  = inp(4)
        self%temp = inp(5) ! patch temperature array
        self%pval = inp(6) ! patch pressure array
        ! calculate pipe parameters
        self%dz = self%len/self%nz
        self%vol = self%dz*self%areaz
      case('brods')
        ! inp = (/dhyd, len, areaz, dir, p2d, sb2st, temp, pval/)
        self%dhyd = inp(1)
        self%len  = inp(2)
        self%areaz= inp(3)
        self%dir  = inp(4)
        self%p2d  = inp(5)
        self%sb2st= inp(6)
        self%temp = inp(7)
        self%pval = inp(8)
        ! calculate pipe parameters
        self%dz = self%len/self%nz
        self%vol = self%dz*self%areaz
      case('wrods')
        ! inp = (/dhyd, len, areaz, dir, p2d, sb2st, temp, pval/)
        ! temporarily use Rehme model, have the same inp as barerod
        self%dhyd = inp(1)
        self%len  = inp(2)
        self%areaz= inp(3)
        self%dir  = inp(4)
        self%p2d  = inp(5)
        self%sb2st= inp(6)
        self%temp = inp(7)
        self%pval = inp(8)
        ! calculate pipe parameters
        self%dz = self%len/self%nz
        self%vol = self%dz*self%areaz
      case('freelevel')
        ! inp = (/dhyd, len, areaz, temp, pval/)
        self%dhyd = inp(1)
        self%len  = inp(2)
        self%areaz= inp(3)
        self%temp = inp(4)
        self%pval = inp(5)
        ! calculate pipe parameters
        self%dz = self%len/self%nz
        self%vol = self%dz*self%areaz
      case('fill')
        ! inp = (/temp, pval/)
        ! inlet boundary pipe with only 1 cell
        ! intended for the 'from' pipe of 'jun-flow' junction
        ! only used for thermal-hydraulic state calculation
        self%temp = inp(1)
        self%pval = inp(2)
      case('break')
        ! inp = (/temp, pval/)
        ! outlet boundary pipe with only 1 cell
        ! intended for the pressure boundary pipe
        self%temp = inp(1) ! effective only when reverseflow occurs
        self%pval = inp(2) ! boundary pressure
      case default
        print*, 'no match for pipe category = ', cate
        stop 'program stopped at Construct_inputs'
    end select
    ! initialize the cell wall-to-fluid heating power
    self%qpow = 0.d0
    ! 
    self%re = 0.d0
    self%pe = 0.d0
  end subroutine init_pipe

  subroutine calc_prop(self)
    !======= Declarations =========
    implicit none
    !== I/O variables ==
    class(Pipe), intent(inout) :: self
    !== local variables ==

    !======= Internals ============
    self%xe = 0.d0
    select case (self%matid)
      case('LBE')
        call fh_lbe(self%enth, self%temp, self%rho, self%miu, self%niu, &
                    self%cp, self%kpa)
        self%pr = self%cp*self%miu/self%kpa
      case('Na1')
        call fh_na1(self%enth, self%temp, self%rho, self%miu, self%niu, &
                    self%cp, self%kpa)
        self%pr = self%cp*self%miu/self%kpa
      case('H2O')
        call fph_h2o(self%pval, self%enth, self%rho, self%niu, self%miu,&
                    self%kpa, self%cp, self%temp, self%xe)
        self%pr = self%cp*self%miu/self%kpa
      case default
        print*, 'no match for fluid material name = ', self%matid
        stop 'program stopped at pipe%calc_pro'
    end select
  end subroutine calc_prop

  subroutine calc_rhs_pipe(self, rhs)
    !======= Declarations =========
    implicit none
    !== I/O variables ==
    class(Pipe), intent(inout) :: self
    ! rhs array
    real(c_double), intent(out) :: rhs(:)
    !== local variables ==
    real(c_double), allocatable :: rho_vol(:)
    integer(c_int) :: i, nn
    !======= Internals ============
    ! initialization
    rhs = 0.d0
    ! convective effect
    if(self%cate=='freelevel') then
      rhs(1) = rhs(1) + self%mh_i - self%mh_o
      self%vol = self%areaz*self%len
    else ! self%nz>1, not free level
      ! the first cell
      rhs(1) = rhs(1) + self%mh_i 
      if(self%mdot>0.d0) then
        rhs(1) = rhs(1) - self%mdot*self%enth(1) 
      else
        rhs(1) = rhs(1) - self%mdot*self%enth(2)
      end if
      ! the last cell
      nn = self%nz
      if(self%mdot>0.d0) then
        rhs(nn) = rhs(nn) + self%mdot*self%enth(nn-1)
      else
        rhs(nn) = rhs(nn) + self%mdot*self%enth(nn)
      end if
      rhs(nn) = rhs(nn) - self%mh_o
      if(nn>2) then
        do i = 2, nn-1
          if(self%mdot>0.d0) then
            rhs(i) = rhs(i) + self%mdot*self%enth(i-1) - self%mdot*self%enth(i)
          else
            rhs(i) = rhs(i) + self%mdot*self%enth(i) - self%mdot*self%enth(i+1)
          end if
        end do
      end if
      allocate(rho_vol(self%nz))
      rhs = rhs + self%qpow
    end if
    rho_vol = self%rho*self%vol
    rhs = rhs/rho_vol
    deallocate(rho_vol)
  end subroutine calc_rhs_pipe
end module Pipe_mod