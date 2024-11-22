! ----------------------------------------------------------------------
! Programmer(s): Yutong CHEN @XJTU, MER Dr. Konstantin MIKITYUK @PSI
! ----------------------------------------------------------------------
!
! A module for defining ...
! 
! ----------------------------------------------------------------------
module Comprhs_mod
  !======= Inclusions =========
  use Control_mod
  use Coupling_mod
  use fsundials_core_mod
  !======= Declarations =========
  implicit none
  

contains
  ! ----------------------------------------------------------------
  ! fcnrob: The CVODE RHS operator function
  !
  ! Return values:
  !    0 = success,
  !    1 = recoverable error,
  !   -1 = non-recoverable error
  ! ----------------------------------------------------------------
  integer(c_int) function fcnrob(t, sunvec_y, sunvec_f, user_data) &
       result(ierr) bind(C,name='fcnrob')

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding

    !======= Declarations =========
    implicit none
    ! calling variables
    real(c_double), value :: t         ! current time
    type(N_Vector)        :: sunvec_y  ! solution N_Vector
    type(N_Vector)        :: sunvec_f  ! function N_Vector
    type(c_ptr),    value :: user_data ! user-defined data
    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer :: yval(:)
    real(c_double), pointer :: fval(:)
    ! local arrays/variables
    real(c_double), allocatable :: rhs(:)
    real(c_double), allocatable :: rhs_fld1(:)
    real(c_double), allocatable :: rhs_fld2(:)
    real(c_double), allocatable :: rhs_sld(:)
    real(c_double) :: t1, t2 ! CPU time
    integer, pointer :: f_user_data
    ! pointer
    !======= Internals ============
    call cpu_time(t1)
    ! upgate the signal values
    call evaluate_signals(t)
    
    call C_F_POINTER(user_data, f_user_data)
    
    ! get data arrays from SUNDIALS vectors
    yval => FN_VGetArrayPointer(sunvec_y)
    fval => FN_VGetArrayPointer(sunvec_f)

    ! update the Reactor_glo variables
    call Reactor_glo%read_from_y(yval)

    ! allocate rhs array
    allocate(rhs(neq_tot))
    rhs = 0.d0
    allocate(rhs_fld1(neq_mdot + neq_flen))
    allocate(rhs_fld2(neq_enth))
    allocate(rhs_sld(neq_frod + neq_hstr))
    ! fill residual vector
    if(solve_fluid) then
      call Reactor_glo%fluid%comp_rhs_fld1(rhs_fld1)
      call Reactor_glo%fluid%comp_rhs_fld2(rhs_fld2)
      if(solve_fuelrod .or. solve_htstr) then
        call Calc_coupling()
      end if
    end if
    if(solve_htstr.or.solve_fuelrod) then
      call Reactor_glo%solid%comp_rhs_sld(rhs_sld)
    end if
    rhs = (/rhs_fld1, rhs_fld2, rhs_sld/)
    fval = rhs
    deallocate(rhs)
    deallocate(rhs_fld1)
    deallocate(rhs_fld2)
    deallocate(rhs_sld)
    ! return success
    ierr = 0
    call cpu_time(t2)
    ! accumulate total CPU time of composing rhs
    time_rhs = time_rhs + t2 - t1
    return

  end function fcnrob
end module Comprhs_mod
