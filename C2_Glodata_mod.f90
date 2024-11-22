! ----------------------------------------------------------------------
! Programmer(s): Yutong CHEN @XJTU, MER Dr. Konstantin MIKITYUK @PSI
! ----------------------------------------------------------------------
!
! A module for defining Global data
! 
! ----------------------------------------------------------------------
module Glodata_mod
  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  !======= Declarations =========
  implicit none
  !== Constant Parameters ==
  ! scientific constants
  real(c_double), parameter :: pi     = 3.141592653589793d0
  real(c_double), parameter :: grav   = 9.8d0
  real(c_double), parameter :: small  = 1.0d-10
  real(c_double), parameter :: minor  = 1.0d-20
  real(c_double), parameter :: large  = 1.0d+10
  real(c_double), parameter :: huge   = 1.0d+20
  real(c_double), parameter :: endval = 1145141919810.d0 ! end marker
  ! array sizes
  integer(c_int), parameter :: maxtab = 100    ! maximum lookup table length
  integer(c_int), parameter :: maxnlyr= 10     ! maximum number of fuel/htstr layers
  integer(c_int), parameter :: nnmlst = 11     ! total number of namelists
  integer(c_int), parameter :: nsolvr = 4      ! total number of solvers
  integer(c_int), parameter :: nfluids= 4      ! total number of fluids
  ! integers relative to file units
  integer(c_int), parameter :: inpfu  = 5      ! input file unit
  integer(c_int), parameter :: otpfu  = 10     ! output file unit baseline
  integer(c_int), parameter :: lbefu  = 10000  ! LBE property data file unit
  integer(c_int), parameter :: na1fu  = 10001  ! Na property data file unit (single phase)
  integer(c_int), parameter :: h2ofu  = 10005  ! H2O property data file unit
  integer(c_int), parameter :: h2osu  = 10006  ! H2O saturated data file unit
  ! indexes of namelist titles at array nmlstr(:)
  integer(c_int), parameter :: idxtime= 1      ! index of &XTIME
  integer(c_int), parameter :: idxsolv= 2      ! index of &XSOLV
  integer(c_int), parameter :: idxpipe= 3      ! index of &XPIPE
  integer(c_int), parameter :: idxjunc= 4      ! index of &XJUNC
  integer(c_int), parameter :: idxcore= 5      ! index of &XCORE
  integer(c_int), parameter :: idxfuel= 6      ! index of &XFUEL
  integer(c_int), parameter :: idxclad= 7      ! index of &XCLAD
  integer(c_int), parameter :: idxfrod= 8      ! index of &XFROD
  integer(c_int), parameter :: idxthbc= 9      ! index of &XTHBC
  integer(c_int), parameter :: idxhstr= 10     ! index of &XHSTR
  integer(c_int), parameter :: idxsgnl= 11     ! index of &XSGNL
  ! indexes of fluid name at fluid counter array numfld(:)
  integer(c_int), parameter :: idxlbe = 1      ! index of LBE fluid      
  integer(c_int), parameter :: idxna1 = 2      ! index of Na fluid (single phase)
  integer(c_int), parameter :: idxna2 = 3      ! index of Na fluid (two phase)
  integer(c_int), parameter :: idxh2o = 4      ! index of H2O fluid
  !
  !== Global Variables ==
  ! logical variables
  logical :: solve_htstr  = .FALSE.
  logical :: solve_fuelrod= .FALSE.
  logical :: solve_fluid  = .FALSE.
  logical :: solve_pkcore = .FALSE.
  logical :: set_dtprint  = .FALSE.
  logical :: set_corepow  = .FALSE.
  ! LBE thermal property data
  real(c_double), allocatable :: tpf_lbe(:,:)
  ! Na1 thermal property data
  real(c_double), allocatable :: tpf_na1(:,:)
  ! H2O thermal property data
  real(c_double), allocatable :: p_h2o(:)
  real(c_double), allocatable :: h_h2o(:)
  real(c_double), allocatable :: tpf_h2o(:,:,:)
  real(c_double), allocatable :: sat_h2o(:,:)
  ! fluid index data
  integer(c_int), allocatable :: findindex(:,:) ! findindex(j,k) means the j-th pipe, k-th node's 'i' in 1-d array
  integer(c_int) :: pindex   ! index of the first pressure equation in Bx=b
  integer(c_int) :: neq_mdot = 0 ! number of ode equations of mdots (jun-head)
  integer(c_int) :: neq_flen = 0 ! number of ode equations of freelevel lengths 
  integer(c_int) :: neq_enth = 0 ! number of ode equations of fluid enthalpies
  ! solid index data
  integer(c_int) :: neq_hstr = 0 ! number of ode equations of heat structure node temperature
  integer(c_int) :: neq_frod = 0 ! number of ode equations of fuelrod node temperature
  ! total ode number
  integer(c_int) :: neq_tot = 0
  ! selection for sundials linear solver
  integer(c_int) :: solopt  = 1
  ! signal-related data
  character(len=20), allocatable :: sgnids_arr(:) ! global array storing signal id strings
  real(c_double), allocatable    :: sgnval_arr(:) ! global array storing signal real values
  ! variables for getting CPU time
  real(c_double) :: time_rhs = 0.d0
  ! parallel calculation data
  integer(c_int) :: num_of_threads = 1 ! number of omp threads
  ! tolerances
  real(c_double) :: rtol
  real(c_double), allocatable :: avtol(:)
  
end module Glodata_mod