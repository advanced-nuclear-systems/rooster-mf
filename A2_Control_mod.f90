! ----------------------------------------------------------------------
! Programmer(s): Yutong CHEN @XJTU, MER Dr. Konstantin MIKITYUK @PSI
! ---------------------------------------------------------------------- 
!
! A module for defining Control module
! 
! ----------------------------------------------------------------------
module Control_mod
  !======= Inclusions ===========
  use Reactor_mod
  use omp_lib
  !======= Declarations =========
  implicit none
  !== SHARED variables ==
  character(len=20) :: id, mat, cate, pid, sigid
  real(c_double) :: ri, ro, tmp0, p0
  integer(c_int) :: nr, pnode, mltpl, nomp
  !== UNIQUE variables ==
  ! XTIME variables-----------------------------------------------------
  real(c_double) :: tstart= 0.0d0  , tend  = 1.d+02, dtout = 1.d+01
  real(c_double) :: dtini = 1.0d-08, dtmin = 1.d-08, dtmax = 1.d-1
  ! XSOLV variables-----------------------------------------------------
  character(len=20) :: solver(nsolvr) = 'end'
  real(c_double) :: rt    = 1.0d-4
  real(c_double) :: at(4) = 1.0d-6
  integer(c_long):: mxsteps= 1000000
  integer(c_int) :: maxnef = 7, maxcor = 3
  ! XPIPE variables-----------------------------------------------------
  real(c_double) :: dhyd, len, dir, areaz, p2d, sb2st
  integer(c_int) :: nz
  ! XJUNC variables-----------------------------------------------------
  character(len=20) :: pid_f, pid_t
  real(c_double) :: kfac = 0.d0
  ! XCORE variables-----------------------------------------------------
  real(c_double) :: power0, lambda, beta, life
  ! XFUEL variables-----------------------------------------------------
  
  ! XCLAD variables-----------------------------------------------------
  
  ! XFROD variables-----------------------------------------------------
  character(len=20) :: fid, cid
  real(c_double) :: hgap, kr, kz
  ! XTHBC variables-----------------------------------------------------
  real(c_double) :: bcval(2)
  ! XHSTR variables-----------------------------------------------------
  character(len=20) :: bcl, bcr
  ! XSGNL variables-----------------------------------------------------
  real(c_double) :: sgnval = endval
  real(c_double) :: sgntab(2, maxtab) = endval
  ! --------------------------------------------------------------------
  !== NAMELIST declarations ==
  ! XTIME namelist  1.--------------------------------------------------
  NAMELIST /XTIME/ tstart, tend  , dtout , dtini , dtmin , dtmax 
  ! XSOLV namelist  2.--------------------------------------------------
  NAMELIST /XSOLV/ solver, nomp  , rt    , at    , maxnef, maxcor, &
                  mxsteps, solopt
  ! XPIPE namelist  3.--------------------------------------------------
  NAMELIST /XPIPE/ id    , mat   , cate  , dhyd  , len   , dir   , &
                   areaz , p2d   , sb2st , nz    , tmp0  , p0    , sigid
  ! XJUNC namelist  4.--------------------------------------------------
  NAMELIST /XJUNC/ cate  , pid_f , pid_t , sigid , kfac  
  ! XCORE namelist  5.--------------------------------------------------
  NAMELIST /XCORE/ power0, sigid , lambda, beta  , life  
  ! XFUEL namelist  6.--------------------------------------------------
  NAMELIST /XFUEL/ id    , mat   , ri    , ro    , nr    , tmp0  
  ! XCLAD namelist  7.--------------------------------------------------
  NAMELIST /XCLAD/ id    , mat   , ri    , ro    , nr    , tmp0  
  ! XFROD namelist  8.--------------------------------------------------
  NAMELIST /XFROD/ id    , fid   , hgap  , cid   , pid   , pnode , &
                   mltpl , kr    , kz    
  ! XTHBC namelist  9.--------------------------------------------------
  NAMELIST /XTHBC/ id    , cate  , bcval , pid   , pnode 
  ! XHSTR namelist 10.--------------------------------------------------
  NAMELIST /XHSTR/ id    , mat   , bcl   , bcr   , ri    , ro    , &
                   nr    , tmp0  , mltpl
  ! XSGNL namelist 11.--------------------------------------------------
  NAMELIST /XSGNL/ id    , cate  , sgnval, sgntab
  !== NAMELIST information ==
  ! namelist titles
  character(len=6) :: nmlstr(nnmlst)
  ! array storing namelist counts
  integer(c_int)   :: nmlnum(nnmlst)
  !== Fluid information ==
  ! array storing fluid show-up counts
  integer(c_int)   :: numfld(nfluids)
contains
! ----------------------------------------------------------------------
! construct_input: a subroutine that read the input file and store all 
! input data into object type(Reactor) :: Reactor_obj
! ----------------------------------------------------------------------
  subroutine construct_input()
    !======= Declarations =========
    implicit none
    !== I/O variables ==
    !== local variables ==
    character(len=6) :: nml_str
    real(c_double) :: dz
    real(c_double),allocatable :: inp_pipe(:)
    integer(c_int) :: i, j, thread_limit
    ! the objects
    type(Pipe),    allocatable :: Pipe_arr(:)
    type(Junc),    allocatable :: Junc_arr(:)
    type(Fuel),    allocatable :: Fuel_arr(:)
    type(Clad),    allocatable :: Clad_arr(:)
    type(Fuelrod), allocatable :: Frod_arr(:)
    type(Thermbc), allocatable :: Thermbc_arr(:)
    type(Htstr),   allocatable :: Htstr_arr(:)
    type(Signal),  allocatable :: Signal_arr(:)
    type(Fuel)    :: fuel_obj
    type(Gasgap)  :: ggap_obj
    type(Clad)    :: clad_obj
    type(Thermbc) :: bcl_obj, bcr_obj
    !======= Internals ============
    ! assign the namelist titles
    nmlstr(idxtime)= '&XTIME'
    nmlstr(idxsolv)= '&XSOLV'
    nmlstr(idxpipe)= '&XPIPE'
    nmlstr(idxjunc)= '&XJUNC'
    nmlstr(idxcore)= '&XCORE'
    nmlstr(idxfuel)= '&XFUEL'
    nmlstr(idxclad)= '&XCLAD'
    nmlstr(idxfrod)= '&XFROD'
    nmlstr(idxthbc)= '&XTHBC'
    nmlstr(idxhstr)= '&XHSTR'
    nmlstr(idxsgnl)= '&XSGNL'
    ! open input file
    open(unit=inpfu, status='old', file='input', form='formatted')
    ! initialize the counters
    nmlnum = 0
    numfld = 0
    ! cycle to count the namelists
    do while (.TRUE.)
      read(inpfu, '(1X,A6)', end=100) nml_str
      do i = 1, nnmlst
        if(nml_str == nmlstr(i)) then
          nmlnum(i) = nmlnum(i) + 1
        end if
      end do
    end do
100 continue

    ! /XTIME/-----------------------------------------------------------
    if(nmlnum(idxtime)>=1) then
      rewind(inpfu)
      write(0,*) '===/XTIME/==='
      ! only read the first /XTIME/
      read(inpfu, nml=XTIME)

    else
      print*, 'no namelist /XTIME/ detected'
      stop 'program stop at Construct_input'
    end if

    ! /XSOLV/-----------------------------------------------------------
    if(nmlnum(idxsolv)>=1) then
      rewind(inpfu)
      write(0,*) '===/XSOLV/==='
      ! only read the first /XSOLV/
      read(inpfu, nml=XSOLV)
      ! cycle the solver array
      do i = 1, nsolvr
        select case (solver(i))
          case('htstr')
            solve_htstr  = .TRUE.
            print*, 'solves htstr'
          case('fuelrod')
            solve_fuelrod= .TRUE.
            print*, 'solves fuelrod'
          case('fluid')
            solve_fluid  = .TRUE.
            print*, 'solves fluid'
          case('pointcore')
            solve_pkcore = .TRUE.
            print*, 'solves pointcore'
          case default
            if(solver(i)=='end') then
              print*, 'all involved solvers are presented above...'
              exit
            else
              print*, 'no match for solver = ', solver(i)
              stop 'program stopped at Construct_input'
            end if
        end select
      end do
      ! set omp parameters
      thread_limit = omp_get_max_threads()
      write(0,'(1X,A27,I3)') 'maximum available threads =', thread_limit
      num_of_threads = min(max(nomp, 1), thread_limit)
      ! print SUNLinSol information
      select case (solopt)
        case(1)
          print*, 'Dense direct linear solver (internal) selected'
        case(2)
          print*, 'PCG iterative solver selected'
        case(3)
          print*, 'SPBCGS solver selected'
        case(4)
          print*, 'SPFGMR solver selected'
        case(5)
          print*, 'SPGMR solver selected'
        case(6)
          print*, 'SPTFQMR solver selected'
        case default
          print*, 'no match for solopt = ', solopt
          stop 'program stopped at Construct_input'
      end select
    else
      print*, 'no namelist /XSOLV/ detected'
      stop 'program stop at Construct_input'
    end if

    if(solve_fluid) then
    ! /XPIPE/-----------------------------------------------------------
      if(nmlnum(idxpipe)>=1) then
        allocate(Pipe_arr(nmlnum(idxpipe)))
        rewind(inpfu)
        write(0,'(1X,A15,I3)') '===/XPIPE/=== *', nmlnum(idxpipe)
        do i = 1, nmlnum(idxpipe)
          read(inpfu, nml=xpipe)
          ! read different categories of pipes seperately
          select case (cate)
            case('round')
              allocate(inp_pipe(6))
              inp_pipe(1) = dhyd
              inp_pipe(2) = len
              inp_pipe(3) = areaz
              inp_pipe(4) = dir
              inp_pipe(5) = tmp0
              inp_pipe(6) = p0
              print*, id, tmp0
            case('brods')
              allocate(inp_pipe(8))
              inp_pipe(1) = dhyd
              inp_pipe(2) = len
              inp_pipe(3) = areaz
              inp_pipe(4) = dir
              inp_pipe(5) = p2d
              inp_pipe(6) = sb2st
              inp_pipe(7) = tmp0
              inp_pipe(8) = p0
            case('wrods')
              allocate(inp_pipe(8))
              inp_pipe(1) = dhyd
              inp_pipe(2) = len
              inp_pipe(3) = areaz
              inp_pipe(4) = dir
              inp_pipe(5) = p2d
              inp_pipe(6) = sb2st
              inp_pipe(7) = tmp0
              inp_pipe(8) = p0
            case('freelevel')
              nz = 1
              allocate(inp_pipe(5))
              inp_pipe(1) = dhyd
              inp_pipe(2) = len
              inp_pipe(3) = areaz
              inp_pipe(4) = tmp0
              inp_pipe(5) = p0
            case('fill')
              nz = 1
              allocate(inp_pipe(2))
              inp_pipe(1) = tmp0
              inp_pipe(2) = p0
              Pipe_arr(i)%sigid = sigid
            case('break')
              nz = 1
              allocate(inp_pipe(2))
              inp_pipe(1) = tmp0
              inp_pipe(2) = p0
            case default
              print*, 'no match for pipe category = ', cate
              stop 'program stopped at Construct_inputs'
          end select
          call Pipe_arr(i)%init_pipe(id,mat,cate,nz,inp_pipe)
          deallocate(inp_pipe)
          ! count the fluid show-ups
          select case (mat)
            case ('LBE')
              numfld(idxlbe) = numfld(idxlbe)+ 1
            case ('Na1')
              numfld(idxna1) = numfld(idxna1)+ 1
            case ('H2O')
              numfld(idxh2o) = numfld(idxh2o)+ 1
            case default
              print*, 'no match for fluid = ', mat
              stop 'program stopped at Construct_inputs'
          end select
        end do
        ! assign the Pipe_arr to Reactor_obj
        call Reactor_glo%fluid%init_fld_p(Pipe_arr)
        ! get the fluid thermophysical data
        if(numfld(idxlbe)>0) call get_lbe_dat()
        if(numfld(idxna1)>0) call get_na1_dat()
        if(numfld(idxh2o)>0) call get_h2o_dat()
      else
        print*, 'missing /XPIPE/ input for fluid solver'
        stop 'program stopped at Construct_input'
      end if
      deallocate(Pipe_arr)
    ! /XJUNC/-----------------------------------------------------------
      if(nmlnum(idxjunc)>=1) then
        allocate(Junc_arr(nmlnum(idxjunc)))
        rewind(inpfu)
        write(0,'(1X,A15,I3)') '===/XJUNC/=== *', nmlnum(idxjunc)
        do i = 1, nmlnum(idxjunc)
          read(inpfu, nml=xjunc)
          ! read different categories of pipes seperately
          select case(cate)
            case('jun-norm') 
              ! normal junction whose mdot is dependent variable
              ! mdot is dependent
              call Junc_arr(i)%init_junc(cate, pid_f, pid_t)
            case('jun-head') 
              ! junciton with pump driving head
              ! mdot is independently solved by ode
              ! pump head controlled by signal sigid
              call Junc_arr(i)%init_junc(cate, pid_f, pid_t, sigid)
            case('jun-flow')
              ! junction with massflowrate signal
              ! mdot is independent, but controlled by signal sigid 
              call Junc_arr(i)%init_junc(cate, pid_f, pid_t, sigid)
            case default
              print*, 'no match for junction category = ', cate
              stop 'program stopped at Construct_inputs'
          end select
          ! local pressure loss coefficient
          Junc_arr(i)%kfac = kfac
        end do
        ! assign the Junc_arr to Reactor_obj
        call Reactor_glo%fluid%init_fld_j(Junc_arr)
      else
        print*, 'missing /XJUNC/ input for fluid solver'
        stop 'program stopped at Construct_input'
      end if
      ! junction mdot solved by ode
      neq_tot = neq_tot + neq_mdot
      ! freelevel length solved by ode
      neq_tot = neq_tot + neq_flen
      ! fluid enthalpy solved by ode
      neq_tot = neq_tot + neq_enth
      !
      call Reactor_glo%fluid%init_fluid()
    end if

    if(solve_fuelrod) then
      if(.not.solve_fluid) then
        print*, 'fuelrod simulation cannot proceed without fluid'
        stop 'program stopped at Construct_input'
      end if
    ! /XCORE/-----------------------------------------------------------
      if(nmlnum(idxcore)>=1) then
        rewind(inpfu)
        write(0,*) '===/XCORE/==='
        ! only read the first /XCORE/
        read(inpfu, nml=XCORE)
        Reactor_glo%core%power0 = power0
        if(solve_pkcore) then
          call Reactor_glo%core%init_core_pkcore(lambda, beta, life)
        end if
      else
        print*, 'missing /XCORE/ input for fuelrod solver'
        stop 'program stopped at Construct_input'
      end if

    ! /XFUEL/-----------------------------------------------------------
      if(nmlnum(idxfuel)>=1) then
        allocate(Fuel_arr(nmlnum(idxfuel)))
        rewind(inpfu)
        write(0,'(1X,A15,I3)') '===/XFUEL/=== *', nmlnum(idxfuel)
        do i = 1, nmlnum(idxfuel)
          read(inpfu, nml=XFUEL)
          call Fuel_arr(i)%init_fp(id, mat, ri, ro, nr, tmp0)
        end do
      else
        print*, 'missing /XFUEL/ input for fuelrod solver'
        stop 'program stopped at Construct_input'
      end if

    ! /XCLAD/-----------------------------------------------------------
      if(nmlnum(idxclad)>=1) then
        allocate(Clad_arr(nmlnum(idxclad)))
        rewind(inpfu)
        write(0,'(1X,A15,I3)') '===/XCLAD/=== *', nmlnum(idxclad)
        do i = 1, nmlnum(idxclad)
          read(inpfu, nml=XCLAD)
          call Clad_arr(i)%init_cd(id, mat, ri, ro, nr, tmp0)
        end do
      else
        print*, 'missing /XCLAD/ input for fuelrod solver'
        stop 'program stopped at Construct_input'
      end if

    ! /XFROD/-----------------------------------------------------------
      if(nmlnum(idxfrod)>=1) then
        allocate(Frod_arr(nmlnum(idxfrod)))
        rewind(inpfu)
        write(0,'(1X,A15,I3)') '===/XFROD/=== *', nmlnum(idxfrod)
        do i = 1, nmlnum(idxfrod)
          read(inpfu, nml=XFROD)
          ! loop for the matching fuel pellet object
          do j = 1, nmlnum(idxfuel)
            if(fid==Fuel_arr(j)%idstr) then
              fuel_obj = Fuel_arr(j)
              exit
            end if
          end do
          if(j>nmlnum(idxfuel)) then
            print*, 'no match for fuel = ', fid
            stop 'program stopped at Construct_input'
          end if
          ! initialize gasgap object
          call ggap_obj%init_gg(hgap)
          ! loop for the matching cladding object
          do j = 1, nmlnum(idxclad)
            if(cid==Clad_arr(j)%idstr) then
              clad_obj = Clad_arr(j)
              exit
            end if
          end do
          if(j>nmlnum(idxclad)) then
            print*, 'no match for clad = ', cid
            stop 'program stopped at Construct_input'
          end if
          ! loop for the matching pipe object
          do j = 1, nmlnum(idxpipe)
            if(pid==Reactor_glo%fluid%pipes(j)%idstr) then
              dz = Reactor_glo%fluid%pipes(j)%dz
              exit
            end if
          end do
          if(j>nmlnum(idxpipe)) then
            print*, 'no match for pipeid = ', pid
            stop 'program stopped at Construct_input'
          end if
          call Frod_arr(i)%init_frd(id, fuel_obj, ggap_obj, clad_obj, &
          kr, kz, pid, pnode, mltpl, dz)
        end do
      else
        print*, 'missing /XFROD/ input for fuelrod solver'
        stop 'program stopped at Construct_input'
      end if
      ! assign the Frod_arr to Reactor_obj (shallow copy)
      call Reactor_glo%Solid%init_sld_f(Frod_arr)
      neq_tot = neq_tot + neq_frod
      ! free the dynamic memory
      deallocate(Fuel_arr)
      deallocate(Clad_arr)
      deallocate(Frod_arr)
    end if

    if(solve_htstr) then
    ! /XTHBC/-----------------------------------------------------------
      if(nmlnum(idxthbc)>=1) then
        allocate(Thermbc_arr(nmlnum(idxthbc)))
        rewind(inpfu)
        write(0,'(1X,A15,I3)') '===/XTHBC/=== *', nmlnum(idxthbc)
        do i = 1, nmlnum(idxthbc)
          read(inpfu, nml=XTHBC)
          if(cate=='coupling') then
            ! here the bcval does not have practical meaning
            bcval = 0.d0
            call Thermbc_arr(i)%init_bc(id, cate, bcval, pid, pnode)
          else
            call Thermbc_arr(i)%init_bc(id, cate, bcval)
          end if
        end do
      else
        print*, 'missing /XTHBC/ input for htstr solver'
        stop 'program stopped at Construct_inputs'
      end if

    ! /XHSTR/-----------------------------------------------------------
      if(nmlnum(idxhstr)>=1) then
        allocate(Htstr_arr(nmlnum(idxhstr)))
        rewind(inpfu)
        write(0,'(1X,A15,I3)') '===/XHSTR/=== *', nmlnum(idxhstr)
        do i = 1, nmlnum(idxhstr)
          read(inpfu, nml=XHSTR)
          ! search for the left thermbc object
          do j = 1, nmlnum(idxthbc)
            if(Thermbc_arr(j)%bcname==bcl) then
              bcl_obj = Thermbc_arr(j)
              exit
            end if
          end do
          if(j>nmlnum(idxthbc)) then
            print*, 'no match for left thermbc name = ', bcl
            stop 'program stopped at Construct_inputs'
          end if
          ! search for the right thermbc object
          do j = 1, nmlnum(idxthbc)
            if(Thermbc_arr(j)%bcname==bcr) then
              bcr_obj = Thermbc_arr(j)
              exit
            end if
          end do
          if(j>nmlnum(idxthbc)) then
            print*, 'no match for right thermbc name = ', bcr
            stop 'program stopped at Construct_inputs'
          end if
          ! initialize the heat structure objects
          call Htstr_arr(i)%init_hs( id, mat, ri, ro, nr, tmp0, bcl_obj, bcr_obj, mltpl)
        end do
        ! assign the Htstr_arr to Reactor_obj (shallow copy)
        call Reactor_glo%Solid%init_sld_h(Htstr_arr)
        ! accumulate the htstr node number
        neq_tot = neq_tot + neq_hstr
        ! free the dynamic memory
        deallocate(Htstr_arr)
        deallocate(Thermbc_arr)
      else
        print*, 'missing /XHSTR/ input for htstr solver'
        stop 'program stopped at Construct_inputs'
      end if
    end if

    ! /XSGNL/-----------------------------------------------------------
    if(nmlnum(idxsgnl)>=1) then
      allocate(Signal_arr(nmlnum(idxsgnl)))
      allocate(Reactor_glo%signals(nmlnum(idxsgnl)))
      ! global arrarys storing signal data
      allocate(sgnids_arr(nmlnum(idxsgnl)))
      allocate(sgnval_arr(nmlnum(idxsgnl)))
      rewind(inpfu)
      write(0,'(1X,A15,I3)') '===/XSGNL/=== *', nmlnum(idxsgnl)
      do i = 1, nmlnum(idxsgnl)
        sgntab = endval
        read(inpfu, nml=XSGNL)
        sgnids_arr(i) = id
        ! echo the signals with fixed name and purpose
        select case (id)
          case('dtprint')
            print*, '*** system signal dtprint specified'
            set_dtprint = .TRUE.
          case('corepow')
            print*, '*** system signal corepow specified'
            set_corepow = .TRUE.
          case default
            print*, ADJUSTL(TRIM(id)), ' signal specified'
        end select
        !
        select case (cate)
          case('constant')
            call Signal_arr(i)%init_sgnl(id, cate)
            sgnval_arr(i) = sgnval
          case('lookup')
            ! overwrite the sgnval with the first column value
            sgnval = sgntab(2,1)
            ! j is the counter of actual input number of sgntab
            j = 1
            do while (.TRUE.)
              if(sgntab(1,j)/=endval .and. sgntab(2,j)/=endval) then
                j = j+1
              else
                if(sgntab(1,j)/=endval .or. sgntab(2,j)/=endval) then
                  print*, 'Warning in signal id = ', id, ', lookup table'
                  print*, 'the x-array and y-array sizes do not match'
                end if
                exit
              end if
            end do
            ! check the lookup table input
            if(j<3) then
              print*, 'lookup table needs at leat 2 columns for signal = ', id
              stop 'program stopped at Construct_inputs'
            else
              call Signal_arr(i)%init_sgnl(id, cate, sgntab(1:2,1:j-1))
            end if
          case default
            print*, 'no match for signal category = ', cate
            stop 'program stopped at Construct_inputs'
        end select
      end do
      Reactor_glo%signals = Signal_arr
      Reactor_glo%nsgnl   = size(Signal_arr)
      deallocate(Signal_arr)
    end if
    ! check input
    if((solve_fuelrod).and.(.not.solve_pkcore).and.(.not.set_corepow)) then
        print*, 'need pkcore solver or corepow signal for core calculation'
        stop 'program stopped at Construct_inputs'
    end if
    ! close input file
    close(inpfu)

  end subroutine construct_input
  
! ----------------------------------------------------------------------
! evaluate_signals: a subroutine that evaluate the siganls, update the 
! signal values during transient calculations
! ----------------------------------------------------------------------
  subroutine evaluate_signals(time)
    !======= Declarations =========
    implicit none
    !== I/O variables ==
    real(c_double) :: time
    !== local variables ==
    integer(c_int) :: i, j
    real(c_double) :: t_arr(1), y_arr(1)
    ! pointer for sliced array
    !======= Internals ============
    do i = 1, Reactor_glo%nsgnl
      ! update signal values array
      select case (Reactor_glo%signals(i)%cate)
        case ('constant')
          ! do not need further action
          continue
        case ('lookup')
          t_arr(1) = time
          call interp1d(Reactor_glo%signals(i)%sgntab(1,:), &
                        Reactor_glo%signals(i)%sgntab(2,:),t_arr,y_arr)
          sgnval_arr(i) = y_arr(1)
        case default
          print*, 'no match for signal type = ', Reactor_glo%signals(i)%cate
          stop 'program stopped at evaluate_signals'
      end select
    end do 
    ! assign the signal values relative to fluid
    if(solve_fluid) then
      do i = 1, Reactor_glo%fluid%npipe
        ! update fill pipe temperature
        if(Reactor_glo%fluid%pipes(i)%cate=='fill') then
          do j = 1, Reactor_glo%nsgnl
            if(Reactor_glo%fluid%pipes(i)%sigid==sgnids_arr(j)) then
              select case (Reactor_glo%fluid%pipes(i)%matid)
                case('LBE')
                  call ht_lbe(sgnval_arr(j:j), Reactor_glo%fluid%pipes(i)%enth)
                case('Na1')
                  call ht_na1(sgnval_arr(j:j), Reactor_glo%fluid%pipes(i)%enth)
                case('H2O')
                  call hpt_h2o(Reactor_glo%fluid%pipes(i)%pval,sgnval_arr(j:j), &
                               Reactor_glo%fluid%pipes(i)%enth)
                case default
                  print*, 'no match for fluid material name = ', Reactor_glo%fluid%pipes(i)%matid
                  stop 'program stopped at Fluid%init_fluid'
              end select
            end if
          end do
        end if
        if(Reactor_glo%fluid%pipes(i)%cate=='break') then
          ! update break pipe pressure
          do j = 1, Reactor_glo%nsgnl
            if(Reactor_glo%fluid%pipes(i)%sigid==sgnids_arr(j)) then
              Reactor_glo%fluid%pipes(i)%pval = sgnval_arr(j)
            end if
          end do
        end if
      end do
      ! update the mdot of jun-flow junctions
      do i = 1, Reactor_glo%fluid%njunc
        if(Reactor_glo%fluid%juncs(i)%cate=='jun-flow') then
          do j = 1, Reactor_glo%nsgnl
            if(Reactor_glo%fluid%juncs(i)%sigid==sgnids_arr(j)) then
              Reactor_glo%fluid%mdots(i) = sgnval_arr(j)
              exit
            end if
          end do
        end if
      end do
      ! update the phead of jun-head junctions
      do i = 1, Reactor_glo%fluid%njunc
        if(Reactor_glo%fluid%juncs(i)%cate=='jun-head') then
          do j = 1, Reactor_glo%nsgnl
            if(Reactor_glo%fluid%juncs(i)%sigid==sgnids_arr(j)) then
              Reactor_glo%fluid%juncs(i)%phead = sgnval_arr(j)
              exit
            end if
          end do
        end if
      end do
    end if
    ! assign the signal values relative to fuelrods
    if(solve_fuelrod.and.set_corepow) then
      do i = 1, Reactor_glo%nsgnl
        if(Reactor_glo%signals(i)%id=='corepow') then
          Reactor_glo%core%power = sgnval_arr(i)*Reactor_glo%core%power0
        end if
      end do
    end if 
  end subroutine evaluate_signals
  
! ----------------------------------------------------------------------
! open_output_files: a subroutine that open the output files, assigning 
! the file units and name strings
! ----------------------------------------------------------------------
  subroutine open_output_files()
    !======= Declarations =========
    implicit none
    !== I/O variables ==
    !== local variables ==
    ! output directory name
    character(len=20) :: rsltdir
    ! path for results
    character(len=50) :: path4results
    ! file id string
    character(len=250):: fileid
    ! column id string
    character(len=18) :: colid
    ! system command string
    character(len=100):: command
    ! loop counters
    integer(c_int) :: i, j, k
    ! status variable
    integer(c_int) :: status
    ! integers recording system time 
    integer(c_int) :: time_arr(8)
    integer(c_int) :: year, month, day, hour, minute, second
    
    !======= Internals ============
    ! the primary output dir 
    command = 'mkdir -p '//'output'
    call SYSTEM(command, status)
    ! the path for results
    call DATE_AND_TIME(VALUES = time_arr)
    year  = time_arr(1) ! 4
    month = time_arr(2) ! 2
    day   = time_arr(3) ! 2
    hour  = time_arr(5) ! 2
    minute= time_arr(6) ! 2
    second= time_arr(7) ! 2
    write(rsltdir, 100) year,'-',month,'-',day,'-',hour,'-',minute,'-',second
    path4results = 'output'//'/'//ADJUSTL(TRIM(rsltdir))
    command = 'mkdir '//path4results
    call SYSTEM(command, status)
    k = 1 ! k serves as the pointer for output file unit
    ! open fluid output files
    if(solve_fluid) then
      command = 'mkdir '//ADJUSTL(TRIM(path4results))//'/'//'fluid'
      call SYSTEM(command, status)
      if(status /= 0) then
        print*, 'Failed to creat output directory for fluid, check your permession.'
      else
        print*, 'Directory created successfully for fluid.'
      end if
      ! open output .dat files relative to fluid
      do i = 1, Reactor_glo%fluid%npipe
        ! fluid temperature
        write(fileid,110) 'pipe-temp-',ADJUSTL(TRIM(Reactor_glo%fluid%pipes(i)%idstr)),'.dat'
        open(unit=otpfu+k, file=ADJUSTL(TRIM(path4results))//'/'//'fluid'//'/'//ADJUSTL(TRIM(fileid)), &
        status='new', form='formatted')
        colid = 'time'
        write(otpfu+k,'(A18)',advance='no')colid
        do j = 1, Reactor_glo%fluid%pipes(i)%nz
          write(colid,120)'temp-',j
          write(otpfu+k,'(A18)',advance='no')colid
        end do
        write(otpfu+k,*)''
        k = k + 1
        ! fluid pressure
        write(fileid,110) 'pipe-p-',ADJUSTL(TRIM(Reactor_glo%fluid%pipes(i)%idstr)),'.dat'
        open(unit=otpfu+k, file=ADJUSTL(TRIM(path4results))//'/'//'fluid'//'/'//ADJUSTL(TRIM(fileid)), &
        status='new', form='formatted')
        colid = 'time'
        write(otpfu+k,'(A18)',advance='no')colid
        do j = 1, Reactor_glo%fluid%pipes(i)%nz
          write(colid,120)'p-',j
          write(otpfu+k,'(A18)',advance='no')colid
        end do
        write(otpfu+k,*)''
        k = k + 1
        ! fluid re
        write(fileid,110) 'pipe-re-',ADJUSTL(TRIM(Reactor_glo%fluid%pipes(i)%idstr)),'.dat'
        open(unit=otpfu+k, file=ADJUSTL(TRIM(path4results))//'/'//'fluid'//'/'//ADJUSTL(TRIM(fileid)), &
        status='new', form='formatted')
        colid = 'time'
        write(otpfu+k,'(A18)',advance='no')colid
        do j = 1, Reactor_glo%fluid%pipes(i)%nz
          write(colid,120)'re-',j
          write(otpfu+k,'(A18)',advance='no')colid
        end do
        write(otpfu+k,*)''
        k = k + 1
        ! fluid pr
        write(fileid,110) 'pipe-pr-',ADJUSTL(TRIM(Reactor_glo%fluid%pipes(i)%idstr)),'.dat'
        open(unit=otpfu+k, file=ADJUSTL(TRIM(path4results))//'/'//'fluid'//'/'//ADJUSTL(TRIM(fileid)), &
        status='new', form='formatted')
        colid = 'time'
        write(otpfu+k,'(A18)',advance='no')colid
        do j = 1, Reactor_glo%fluid%pipes(i)%nz
          write(colid,120)'pr-',j
          write(otpfu+k,'(A18)',advance='no')colid
        end do
        write(otpfu+k,*)''
        k = k + 1
        ! fluid pe
        write(fileid,110) 'pipe-pe-',ADJUSTL(TRIM(Reactor_glo%fluid%pipes(i)%idstr)),'.dat'
        open(unit=otpfu+k, file=ADJUSTL(TRIM(path4results))//'/'//'fluid'//'/'//ADJUSTL(TRIM(fileid)), &
        status='new', form='formatted')
        colid = 'time'
        write(otpfu+k,'(A18)',advance='no')colid
        do j = 1, Reactor_glo%fluid%pipes(i)%nz
          write(colid,120)'pe-',j
          write(otpfu+k,'(A18)',advance='no')colid
        end do
        write(otpfu+k,*)''
        k = k + 1
        ! fluid xe
        write(fileid,110) 'pipe-xe-',ADJUSTL(TRIM(Reactor_glo%fluid%pipes(i)%idstr)),'.dat'
        open(unit=otpfu+k, file=ADJUSTL(TRIM(path4results))//'/'//'fluid'//'/'//ADJUSTL(TRIM(fileid)), &
        status='new', form='formatted')
        colid = 'time'
        write(otpfu+k,'(A18)',advance='no')colid
        do j = 1, Reactor_glo%fluid%pipes(i)%nz
          write(colid,120)'xe-',j
          write(otpfu+k,'(A18)',advance='no')colid
        end do
        write(otpfu+k,*)''
        k = k + 1
        ! fluid htcmod
        write(fileid,110) 'pipe-hmod-',ADJUSTL(TRIM(Reactor_glo%fluid%pipes(i)%idstr)),'.dat'
        open(unit=otpfu+k, file=ADJUSTL(TRIM(path4results))//'/'//'fluid'//'/'//ADJUSTL(TRIM(fileid)), &
        status='new', form='formatted')
        colid = 'time'
        write(otpfu+k,'(A18)',advance='no')colid
        do j = 1, Reactor_glo%fluid%pipes(i)%nz
          write(colid,120)'hmod-',j
          write(otpfu+k,'(A18)',advance='no')colid
        end do
        write(otpfu+k,*)''
        k = k + 1
      end do
      ! fluid mdots
      write(fileid,110) 'mdots-','juncs','.dat'
      open(unit=otpfu+k, file=ADJUSTL(TRIM(path4results))//'/'//'fluid'//'/'//ADJUSTL(TRIM(fileid)), &
      status='new', form='formatted')
      colid = 'time'
      write(otpfu+k,'(A18)',advance='no')colid
      do i = 1, Reactor_glo%fluid%njunc
        write(colid,120)'mdot-',i
        write(otpfu+k,'(A18)',advance='no')colid
      end do
      write(otpfu+k,*)''
      k = k + 1
    end if
    ! open htstr output files
    if(solve_htstr) then
      command = 'mkdir '//ADJUSTL(TRIM(path4results))//'/'//'htstr'
      call SYSTEM(command, status)
      if(status /= 0) then
        print*, 'Failed to creat output directory for htstr, check your permession.'
      else
        print*, 'Directory created successfully for htstr.'
      end if
      ! open output .dat files relative to htstrs
      do i = 1, Reactor_glo%solid%nhtstr
        write(fileid,110) 'htstr-temp-',ADJUSTL(TRIM(Reactor_glo%solid%htstrs(i)%idstr)),'.dat'
        open(unit=otpfu+k, file=ADJUSTL(TRIM(path4results))//'/'//'htstr'//'/'//ADJUSTL(TRIM(fileid)), &
        status='new', form='formatted')
        colid = 'time'
        write(otpfu+k,'(A18)',advance='no')colid
        do j = 1, Reactor_glo%solid%htstrs(i)%nr
          write(colid,120)'temp-',j
          write(otpfu+k,'(A18)',advance='no')colid
        end do
        write(otpfu+k,*)''
        k = k + 1
      end do
    end if
    ! open fuelrod output files
    if(solve_fuelrod) then
      command = 'mkdir '//ADJUSTL(TRIM(path4results))//'/'//'fuelrod'
      call SYSTEM(command, status)
      if(status /= 0) then
        print*, 'Failed to creat output directory for fuelrod, check your permession.'
      else
        print*, 'Directory created successfully for fuelrod.'
      end if
      ! open output .dat files relative to fuelrods
      do i = 1, Reactor_glo%solid%nfrod
        write(fileid,110) 'frod-temp-',ADJUSTL(TRIM(Reactor_glo%solid%frods(i)%idstr)),'.dat'
        open(unit=otpfu+k, file=ADJUSTL(TRIM(path4results))//'/'//'fuelrod'//'/'//ADJUSTL(TRIM(fileid)), &
        status='new', form='formatted')
        colid = 'time'
        write(otpfu+k,'(A18)',advance='no')colid
        do j = 1, Reactor_glo%solid%frods(i)%nf
          write(colid,120)'fuel-',j
          write(otpfu+k,'(A18)',advance='no')colid
        end do
        do j = 1, Reactor_glo%solid%frods(i)%nc
          write(colid,120)'clad-',j
          write(otpfu+k,'(A18)',advance='no')colid
        end do
        write(otpfu+k,*)''
        k = k + 1
      end do
    end if

100 format(I4.4, A1, I2.2, A1, I2.2, A1, I2.2, A1, I2.2, A1, I2.2)
110 format(A, A, A)
120 format(A, I2.2)
  end subroutine open_output_files

  subroutine print_output_files(time)
    !======= Declarations =========
    implicit none
    !== I/O variables ==
    real(c_double) :: time
    !== local variables ==
    character(len=18) :: colstr
    integer(c_int) :: i, j, k

    !======= Internals ============
    k = 1 ! k serves as the pointer for output file unit
    ! fluid outputs
    if(solve_fluid) then
      ! write to output .dat files relative to fluid
      do i = 1, Reactor_glo%fluid%npipe
        ! fluid temperature
        write(colstr,'(ES18.8e2)')time
        write(otpfu+k,'(A18)',advance='no') ADJUSTL(colstr)
        do j = 1, Reactor_glo%fluid%pipes(i)%nz
          write(colstr,'(ES18.8e2)')Reactor_glo%fluid%pipes(i)%temp(j)
          write(otpfu+k,'(A18)',advance='no') ADJUSTL(colstr)
        end do
        write(otpfu+k,*)''
        k = k + 1
        ! fluid pressure
        write(colstr,'(ES18.8e2)')time
        write(otpfu+k,'(A18)',advance='no') ADJUSTL(colstr)
        do j = 1, Reactor_glo%fluid%pipes(i)%nz
          write(colstr,'(ES18.8e2)')Reactor_glo%fluid%pipes(i)%pval(j)
          write(otpfu+k,'(A18)',advance='no') ADJUSTL(colstr)
        end do
        write(otpfu+k,*)''
        k = k + 1
        ! fluid re
        write(colstr,'(ES18.8e2)')time
        write(otpfu+k,'(A18)',advance='no') ADJUSTL(colstr)
        do j = 1, Reactor_glo%fluid%pipes(i)%nz
          write(colstr,'(ES18.8e2)')Reactor_glo%fluid%pipes(i)%re(j)
          write(otpfu+k,'(A18)',advance='no') ADJUSTL(colstr)
        end do
        write(otpfu+k,*)''
        k = k + 1
        ! fluid pr
        write(colstr,'(ES18.8e2)')time
        write(otpfu+k,'(A18)',advance='no') ADJUSTL(colstr)
        do j = 1, Reactor_glo%fluid%pipes(i)%nz
          write(colstr,'(ES18.8e2)')Reactor_glo%fluid%pipes(i)%pr(j)
          write(otpfu+k,'(A18)',advance='no') ADJUSTL(colstr)
        end do
        write(otpfu+k,*)''
        k = k + 1
        ! fluid pe
        write(colstr,'(ES18.8e2)')time
        write(otpfu+k,'(A18)',advance='no') ADJUSTL(colstr)
        do j = 1, Reactor_glo%fluid%pipes(i)%nz
          write(colstr,'(ES18.8e2)')Reactor_glo%fluid%pipes(i)%pe(j)
          write(otpfu+k,'(A18)',advance='no') ADJUSTL(colstr)
        end do
        write(otpfu+k,*)''
        k = k + 1
        ! fluid xe
        write(colstr,'(ES18.8e2)')time
        write(otpfu+k,'(A18)',advance='no') ADJUSTL(colstr)
        do j = 1, Reactor_glo%fluid%pipes(i)%nz
          write(colstr,'(ES18.8e2)')Reactor_glo%fluid%pipes(i)%xe(j)
          write(otpfu+k,'(A18)',advance='no') ADJUSTL(colstr)
        end do
        write(otpfu+k,*)''
        k = k + 1
        ! fluid htcmod
        write(colstr,'(ES18.8e2)')time
        write(otpfu+k,'(A18)',advance='no') ADJUSTL(colstr)
        do j = 1, Reactor_glo%fluid%pipes(i)%nz
          write(colstr,'(I18)')Reactor_glo%fluid%pipes(i)%htcmod(j)
          write(otpfu+k,'(A18)',advance='no') ADJUSTL(colstr)
        end do
        write(otpfu+k,*)''
        k = k + 1
      end do
      ! fluid mdots
      write(colstr,'(ES18.8e2)')time
      write(otpfu+k,'(A18)',advance='no') ADJUSTL(colstr)
      do i = 1, Reactor_glo%fluid%njunc
        write(colstr,'(ES18.8e2)')Reactor_glo%fluid%mdots(i)
        write(otpfu+k,'(A18)',advance='no') ADJUSTL(colstr)
      end do
      write(otpfu+k,*)''
      k = k + 1
    end if
    ! htstr outputs
    if(solve_htstr) then
      ! write to output .dat files relative to htstrs
      do i = 1, Reactor_glo%solid%nhtstr
        write(colstr,'(ES18.8e2)')time
        write(otpfu+k,'(A18)',advance='no') ADJUSTL(colstr)
        do j = 1, Reactor_glo%solid%htstrs(i)%nr
          write(colstr,'(ES18.8e2)')Reactor_glo%solid%htstrs(i)%temp(j)
          write(otpfu+k,'(A18)',advance='no') ADJUSTL(colstr)
        end do
        write(otpfu+k,*)''
        k = k + 1
      end do
    end if
    ! fuelrod outputs
    if(solve_fuelrod) then
      ! write to output .dat files relative to fuelrods
      do i = 1, Reactor_glo%solid%nfrod
        write(colstr,'(ES18.8e2)')time
        write(otpfu+k,'(A18)',advance='no') ADJUSTL(colstr)
        do j = 1, Reactor_glo%solid%frods(i)%nf
          write(colstr,'(ES18.8e2)')Reactor_glo%solid%frods(i)%fuel%temp(j)
          write(otpfu+k,'(A18)',advance='no') ADJUSTL(colstr)
        end do
        do j = 1, Reactor_glo%solid%frods(i)%nc
          write(colstr,'(ES18.8e2)')Reactor_glo%solid%frods(i)%clad%temp(j)
          write(otpfu+k,'(A18)',advance='no') ADJUSTL(colstr)
        end do
        write(otpfu+k,*)''
        k = k + 1
      end do
    end if
  end subroutine print_output_files

end module Control_mod
