! ----------------------------------------------------------------------
! Programmer(s): Yutong CHEN @XJTU, MER Dr. Konstantin MIKITYUK @PSI
! ----------------------------------------------------------------------
!
! A module for defining type solid
! 
! ----------------------------------------------------------------------
module Solid_mod
  !======= Inclusions =========
  use Fuelrod_mod
  use Htstr_mod

  !======= Declarations =========
  implicit none
  
  type :: Solid
    !======= Member Variables =========
    ! total fuel volume
    real(c_double) :: fvol_tot = 0.d0
    ! number of fuelrods, htstr
    integer(c_int) :: nfrod, nhtstr
    ! array storing start/end index in rhs array
    integer(c_int), allocatable :: idx1_f(:), idx2_f(:)
    integer(c_int), allocatable :: idx1_h(:), idx2_h(:)
    ! array storing fuelrod objects
    type(Fuelrod), allocatable :: frods(:)
    ! array storing heat structure objects
    type(Htstr), allocatable   :: htstrs(:)
 
  contains
    !======= Member Procedures ========
    procedure :: init_sld_f      ! constructor for setting frods array
    procedure :: init_sld_h      ! constructor for setting htstrs array
    procedure :: comp_rhs_sld    ! compose the rhs array
  end type Solid
 
contains

  subroutine init_sld_f(self, frod_array)
    !======= Declarations =========
    implicit none
    !== I/O variables ==
    class(Solid), intent(inout) :: self
    class(Fuelrod), intent(in)  :: frod_array(:)
    !== local variables ==
    integer(c_int) :: i, j, k, n
    !======= Internals ============
    self%nfrod = size(frod_array)
    ! allocate and assign the arrays storing structure objects
    allocate(self%frods(self%nfrod))
    self%frods = frod_array
    ! calculate the total number of fuelrod nodes
    allocate(self%idx1_f(self%nfrod))
    allocate(self%idx2_f(self%nfrod))
    do i = 1, self%nfrod
      neq_frod = neq_frod + self%frods(i)%ntot
      if(i==1) then
        self%idx1_f(i) = 1
      else
        self%idx1_f(i) = self%idx2_f(i-1)+1
      end if
      self%idx2_f(i) = self%idx1_f(i) + self%frods(i)%ntot-1
      if(self%frods(i)%fuel%nlyr==1) then
        if(self%frods(i)%fuel%isfu(1)) then
          self%fvol_tot  = self%fvol_tot + sum(self%frods(i)%fuel%vol)*self%frods(i)%dz*self%frods(i)%mltpl
        end if
      else
        do j = 1, self%frods(i)%fuel%nlyr
          if(self%frods(i)%fuel%isfu(j)) then
            if(j==1) then
              k = 1
            else
              k = sum(self%frods(i)%fuel%nrs(1:j-1)) + 1
            end if
            n = self%frods(i)%fuel%nrs(j)
            self%fvol_tot  = self%fvol_tot + sum(self%frods(i)%fuel%vol(k:k+n-1))*self%frods(i)%dz*self%frods(i)%mltpl
          end if
        end do
      end if
    end do
  end subroutine init_sld_f
  
  subroutine init_sld_h(self, htstr_array)
    !======= Declarations =========
    implicit none
    !== I/O variables ==
    class(Solid), intent(inout) :: self
    class(Htstr), intent(in)    :: htstr_array(:)
    !== local variables ==
    integer(c_int) :: i
    !======= Internals ============
    self%nhtstr= size(htstr_array)
    ! allocate and assign the arrays storing structure objects
    allocate(self%htstrs(self%nhtstr))
    self%htstrs= htstr_array
    ! calculate the total number of htstr nodes
    allocate(self%idx1_h(self%nhtstr))
    allocate(self%idx2_h(self%nhtstr))
    do i = 1, self%nhtstr
      neq_hstr = neq_hstr + self%htstrs(i)%nr
      if(i==1) then
        self%idx1_h(i) = 1
      else
        self%idx1_h(i) = self%idx2_h(i-1)+1
      end if
      self%idx2_h(i) = self%idx1_h(i) + self%htstrs(i)%nr-1
    end do
  end subroutine init_sld_h

  subroutine comp_rhs_sld(self, rhs)
    !======= Declarations =========
    implicit none
    !== I/O variables ==
    class(Solid), intent(inout) :: self
    ! rhs array
    real(c_double), intent(out) :: rhs(:)
    !== local variables ==
    integer(c_int) :: i
    real(c_double), allocatable :: rhs_frod(:)
    real(c_double), allocatable :: rhs_htstr(:)
    !======= Internals ============
    ! fuelrods exist in solid field
    if(solve_fuelrod) then
      allocate(rhs_frod(neq_frod))
      do i = 1, self%nfrod
        call self%frods(i)%comp_rhs_frd(rhs_frod(self%idx1_f(i):self%idx2_f(i)))
      end do
    end if
    ! htstrs exist in solid field
    if(solve_htstr) then
      allocate(rhs_htstr(neq_hstr))
      !$OMP PARALLEL NUM_THREADS(num_of_threads)
      !$OMP DO
      do i = 1, self%nhtstr
        call self%htstrs(i)%calc_rhs_hs(rhs_htstr(self%idx1_h(i):self%idx2_h(i)))
      end do
      !$OMP END DO
      !$OMP END PARALLEL
    end if
    ! compose the rhs array
    if(allocated(rhs_frod)) then
      if(allocated(rhs_htstr)) then
        rhs = (/rhs_frod,rhs_htstr/)
        deallocate(rhs_frod)
        deallocate(rhs_htstr)
      else
        rhs = rhs_frod
        deallocate(rhs_frod)
      end if
    else
      if(allocated(rhs_htstr)) then
        rhs = rhs_htstr
        deallocate(rhs_htstr)
      else
        print*, 'no structure object allocated'
        stop 'program stopped at Solid%comp_rhs_sld'
      end if
    end if
  end subroutine comp_rhs_sld

end module Solid_mod