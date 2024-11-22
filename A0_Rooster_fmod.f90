! ----------------------------------------------------------------------
! Programmer(s): Yutong CHEN @XJTU, MER Dr. Konstantin MIKITYUK @PSI
! ----------------------------------------------------------------------
!
! A program for realizing Rooster framework using modern fortran
! 
! ----------------------------------------------------------------------
program Rooster
  !======= Inclusions ===========
  use Comprhs_mod
  use Control_mod
  use Coupling_mod

  use fcvode_mod                    ! Fortran interface to CVODE
  use fnvector_serial_mod           ! Fortran interface to serial N_Vector
  use fsunmatrix_dense_mod          ! Fortran interface to dense SUNMatrix
  
  use fsunlinsol_dense_mod          ! Fortran interface to dense SUNLinearSolver    1.
  use fsunlinsol_pcg_mod            ! Fortran interface to pcg SUNLinearSolver      2.
  use fsunlinsol_spbcgs_mod         ! Fortran interface to spbcgs SUNLinearSolver   3.
  use fsunlinsol_spfgmr_mod         ! Fortran interface to spfgmr SUNLinearSolver   4.
  use fsunlinsol_spgmr_mod          ! Fortran interface to spgmr SUNLinearSolver    5.
  use fsunlinsol_sptfqmr_mod        ! Fortran interface to sptfqmr SUNLinearSolver  6.

  !======= Declarations =========
  implicit none
  ! solution and tolerance vectors, size of which would be determined later
  real(c_double), allocatable :: yval(:)
  ! time variables
  real(c_double) :: tout, tret(1)
  ! CPU time variables
  real(c_double) :: tic, tac
  ! integers
  integer(c_int) :: retval
  integer(c_int) :: i
  integer(c_long):: neq
  ! sundials types
  type(N_Vector),           pointer :: sunvec_y      ! sundials solution vector
  type(N_Vector),           pointer :: sunvec_av     ! sundials tolerance vector
  type(SUNMatrix),          pointer :: sunmat_A      ! sundials matrix
  type(SUNLinearSolver),    pointer :: sunlinsol_LS  ! sundials linear solver
  type(c_ptr)                       :: cvode_mem     ! CVode memory
  type(c_ptr)                       :: sunctx        ! SUNDIALS simulation context
  
  !======= Internals ============
  call cpu_time(tic)

  retval = FSUNContext_Create(SUN_COMM_NULL, sunctx)
  ! initialize solution vectors and tolerances
  call construct_input()
  ! open output files
  call open_output_files()
  ! convert type of neq from c_int to c_double
  neq = int(neq_tot, kind=c_long)
  allocate(yval(neq))
  allocate(avtol(neq))
  ! set rtol
  rtol = rt
  ! set avtol
  ! mass flow rate [kg/s]
  i = 1
  if(neq_mdot>0) then
    avtol(i:i+neq_mdot-1) = at(1)
  end if
  ! ---
  ! freelevel length [m]
  i = i + neq_mdot
  if(neq_flen>0) then
    avtol(i:i+neq_flen-1) = at(2)
  end if
  ! ---
  ! fluid enthalpy [J/kg]
  i = i+neq_flen
  if(neq_enth>0) then
    avtol(i:i+neq_enth-1) = at(3)
  end if
  ! ---
  ! solid temperature [K]
  i = i+neq_enth
  if(neq_frod>0) then
    avtol(i:i+neq_frod-1) = at(4)
  end if
  i = i+neq_frod
  if(neq_hstr>0) then
    avtol(i:i+neq_hstr-1) = at(4)
  end if
  ! ---
  call Reactor_glo%write_to_y(yval)

  ! create serial vectors
  sunvec_y => FN_VMake_Serial(neq, yval, sunctx)
  if (.not. associated(sunvec_y)) then
    print *, 'ERROR: sunvec = NULL'
    stop 1
  end if
  sunvec_av => FN_VMake_Serial(neq, avtol, sunctx)
  if (.not. associated(sunvec_av)) then
    print *, 'ERROR: sunvec = NULL'
    stop 1
  end if
  
  ! Call FCVodeCreate and FCVodeInit to create and initialize CVode memory
  cvode_mem = FCVodeCreate(CV_BDF, sunctx)
  if (.not. c_associated(cvode_mem)) print *, 'ERROR: cvode_mem = NULL'

  retval = FCVodeInit(cvode_mem, c_funloc(fcnrob), tstart, sunvec_y)
  if (retval /= 0) then
     print *, 'Error in FCVodeInit, retval = ', retval, '; halting'
     stop 1
  end if

  retval = FCVodeSetMaxNumSteps(cvode_mem, mxsteps)
  if (retval /= 0) then
     print *, 'Error in FCVodeSetMaxNumSteps, retval = ', retval, '; halting'
     stop 1
  end if

  retval = FCVodeSetMaxErrTestFails(cvode_mem, maxnef)
  if (retval /= 0) then
    print *, 'Error in FCVodeSetMaxErrTestFails, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FCVodeSetMaxNonlinIters(cvode_mem, maxcor)
  if (retval /= 0) then
    print *, 'Error in FCVodeSetMaxNonlinIters, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FCVodeSetInitStep(cvode_mem, dtini)
  if (retval /= 0) then
    print *, 'Error in FCVodeSetInitStep, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FCVodeSetMinStep(cvode_mem, dtmin)
  if (retval /= 0) then
    print *, 'Error in FCVodeSetMinStep, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FCVodeSetMaxStep(cvode_mem, dtmax)
  if (retval /= 0) then
    print *, 'Error in FCVodeSetMaxStep, retval = ', retval, '; halting'
    stop 1
  end if

  ! Call FCVodeSVtolerances to set tolerances
  retval = FCVodeSVtolerances(cvode_mem, rtol, sunvec_av)
  if (retval /= 0) then
     print *, 'Error in FCVodeSVtolerances, retval = ', retval, '; halting'
     stop 1
  end if
  select case (solopt)
    case(1) ! Dense direct linear solver (internal) (Dense)
      ! Create dense SUNMatrix for use in linear solves
      sunmat_A => FSUNDenseMatrix(neq, neq, sunctx)
      if (.not. associated(sunmat_A)) then
         print *, 'ERROR: sunmat = NULL'
         stop 1
      end if
      ! Create dense SUNLinearSolver object
      sunlinsol_LS => FSUNLinSol_Dense(sunvec_y, sunmat_A, sunctx)
      if (.not. associated(sunlinsol_LS)) then
         print *, 'ERROR: sunlinsol = NULL, Dense'
         stop 1
      end if
    case(2) ! Preconditioned Conjugate Gradient iterative solver (PCG)
      ! SUNMatrix
      sunmat_A => null()
      ! Tell CVODE to use a PCG linear solver.
      sunlinsol_LS => FSUNLinSol_PCG(sunvec_y,SUN_PREC_NONE, 5, sunctx)
      if (.not. associated(sunlinsol_LS)) then
         print *, 'ERROR: sunls = NULL, PCG'
         stop 1
      end if
      !
    case(3) ! Scaled-Preconditioned BiCGtab iterative Solver (SPBCGS)
      ! SUNMatrix
      sunmat_A => null()
      ! Tell CVODE to use a SPBCGS linear solver
      sunlinsol_LS => FSUNLinSol_SPBCGS(sunvec_y,SUN_PREC_NONE, 5, sunctx)
      if (.not. associated(sunlinsol_LS)) then
         print *, 'ERROR: sunls = NULL, SPBCGS'
         stop 1
      end if
      !
    case(4) ! Scaled-Preconditioned FGMRES iterative solver (SPFGMR)
      ! SUNMatrix
      sunmat_A => null()
      ! Tell CVODE to use a SPFGMR linear solver.
      sunlinsol_LS => FSUNLinSol_SPFGMR(sunvec_y,SUN_PREC_NONE, 5, sunctx)
      if (.not. associated(sunlinsol_LS)) then
         print *, 'ERROR: sunls = NULL, SPFGMR'
         stop 1
      end if
      ! Set the type of Gram-Schmidt orthogonalization to use
      retval = FSUNLinSol_SPFGMRSetGSType(sunlinsol_LS,SUN_MODIFIED_GS)
      if (retval /= 0) then
         print *, 'Error in FSUNLinSol_SPGMRSetGSType'
         stop 1
      end if
    case(5) ! Scaled-Preconditioned GMRES iterative solver
      ! SUNMatrix
      sunmat_A => null()
      ! Tell CVODE to use a SPGMR linear solver.
      sunlinsol_LS => FSUNLinSol_SPGMR(sunvec_y,SUN_PREC_NONE, 5, sunctx)
      if (.not. associated(sunlinsol_LS)) then
         print *, 'ERROR: sunls = NULL, SPGMR'
         stop 1
      end if
      ! Set the type of Gram-Schmidt orthogonalization to use
      retval = FSUNLinSol_SPGMRSetGSType(sunlinsol_LS,SUN_MODIFIED_GS)
      if (retval /= 0) then
         print *, 'Error in FSUNLinSol_SPGMRSetGSType'
         stop 1
      end if
    case(6) ! Scaled-Preconditioned TFQMR iterative solver
      ! SUNMatrix
      sunmat_A => null()
      ! Tell CVODE to use a SPTFQMR linear solver.
      sunlinsol_LS => FSUNLinSol_SPTFQMR(sunvec_y,SUN_PREC_NONE, 5, sunctx)
      if (.not. associated(sunlinsol_LS)) then
         print *, 'ERROR: sunls = NULL, SPGMR'
         stop 1
      end if
    case default
      print*, 'no match for solopt = ', solopt
      stop 'program stopped at program Rooster'
  end select

  ! Attach the matrix and linear solver
  retval = FCVodeSetLinearSolver(cvode_mem, sunlinsol_LS, sunmat_A)
  if (retval /= 0) then
     print *, 'Error in FCVodeSetLinearSolver, retval = ', retval, '; halting'
     stop 1
  end if

  ! Output the initial values (read from input files)
  tout = tstart
  call print_output_files(tout)
  ! Push on the time step
  tout = tout + dtout

  ! Main integration loop
  do while(tout<=tend)
    retval = FCVode(cvode_mem, tout, sunvec_y, tret(1), CV_NORMAL)
    if (retval < 0) then
       print *, 'Error in FCVode, retval = ', retval, '; halting'
       stop 1
    end if
    ! Output the data
    call print_output_files(tout)
    print*, tout
    if (retval .eq. CV_SUCCESS) then
      tout = tout + dtout
    end if
 end do

  ! free memory
  call FCVodeFree(cvode_mem)
  retval = FSUNLinSolFree(sunlinsol_LS)
  call FSUNMatDestroy(sunmat_A)
  call FN_VDestroy(sunvec_y)
  call FN_VDestroy(sunvec_av)
  retval = FSUNContext_Free(sunctx)
  !
  call cpu_time(tac)
  print*, 'total CPU time = ', tac - tic
  print*, 'CPU time of composing RHS', time_rhs
end program Rooster