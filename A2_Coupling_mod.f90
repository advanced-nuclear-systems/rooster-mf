! ----------------------------------------------------------------------
! Programmer(s): Yutong CHEN @XJTU, MER Dr. Konstantin MIKITYUK @PSI
! ----------------------------------------------------------------------
!
! A module for defining ...
! 
! ----------------------------------------------------------------------
module Coupling_mod
  !======= Inclusions =========
  use Reactor_mod
  !======= Declarations =========
  implicit none
  
contains
  ! Calculate the coupling variables
  subroutine Calc_coupling()
    !======= Declarations =========
    implicit none
    !======= Internals ============
    if(solve_fluid) then
      if(solve_htstr.or.solve_fuelrod) then
        call qflux_w_to_f()
        if(solve_fuelrod) then
          call qvol_fuel()
        end if
      end if
    end if
  end subroutine Calc_coupling

  ! Get the wall_to_fluid convective heat flux for coupling calculation
  subroutine qflux_w_to_f()
    !======= Declarations =========
    implicit none
    !== local variables ==
    real(c_double) :: area_ht, qflux, ri, ro
    real(c_double) :: areaz, dz
    real(c_double) :: inp(10)
    integer(c_int) :: i, j, k
    integer(c_int) :: mltpl, nr
    type(Thermbc)  :: bcl, bcr
    type(Fuelrod)  :: frod
    !======= Internals ============
    ! initialize the pipe cell heating power
    do i = 1, Reactor_glo%fluid%npipe
      Reactor_glo%fluid%pipes(i)%qpow = 0.d0
    end do

    ! heat structure to fluid heat transfer
    do i = 1, Reactor_glo%solid%nhtstr
      bcl   = Reactor_glo%solid%htstrs(i)%bcleft
      bcr   = Reactor_glo%solid%htstrs(i)%bcright
      mltpl = Reactor_glo%solid%htstrs(i)%mltpl
      nr    = Reactor_glo%solid%htstrs(i)%nr
      ri    = Reactor_glo%solid%htstrs(i)%ri
      ro    = Reactor_glo%solid%htstrs(i)%ro
      ! treat left boundary
      if(bcl%bctype=='coupling') then
        do j = 1, Reactor_glo%fluid%npipe
          if(Reactor_glo%fluid%pipes(j)%idstr==bcl%pipeid) then
            exit
          end if
        end do
        areaz = Reactor_glo%fluid%pipes(j)%areaz
        dz    = Reactor_glo%fluid%pipes(j)%dz
        k = bcl%pipenode
        inp(1) = Reactor_glo%fluid%pipes(j)%re(k)
        inp(2) = Reactor_glo%fluid%pipes(j)%pr(k)
        inp(3) = Reactor_glo%fluid%pipes(j)%pe(k)
        inp(4) = Reactor_glo%fluid%pipes(j)%temp(k)
        inp(5) = Reactor_glo%solid%htstrs(i)%temp(1) ! left-most node temperature
        inp(6) = Reactor_glo%fluid%pipes(j)%dhyd
        inp(7) = Reactor_glo%fluid%pipes(j)%kpa(k)
        inp(8) = Reactor_glo%fluid%pipes(j)%pval(k)
        ! additional inputs for possible boiling issues
        if(Reactor_glo%fluid%pipes(j)%matid=='H2O') then
          inp(9) = Reactor_glo%fluid%pipes(j)%xe(k)
          inp(10)= Reactor_glo%fluid%pipes(j)%mdot/Reactor_glo%fluid%pipes(j)%areaz ! G
        end if
        call qw_cal(Reactor_glo%fluid%pipes(j)%matid, inp, qflux, &
                    Reactor_glo%fluid%pipes(j)%htcmod(k))
        Reactor_glo%solid%htstrs(i)%qflux_i =-qflux
        area_ht = 2.d0*pi*ri*dz*mltpl
        Reactor_glo%fluid%pipes(j)%qpow(k) = Reactor_glo%fluid%pipes(j)%qpow(k)+area_ht*qflux
      end if
      ! treat right boundary
      if(bcr%bctype=='coupling') then
        do j = 1, Reactor_glo%fluid%npipe
          if(Reactor_glo%fluid%pipes(j)%idstr==bcr%pipeid) then
            exit
          end if
        end do
        areaz = Reactor_glo%fluid%pipes(j)%areaz
        dz    = Reactor_glo%fluid%pipes(j)%dz
        k = bcr%pipenode
        inp(1) = Reactor_glo%fluid%pipes(j)%re(k)
        inp(2) = Reactor_glo%fluid%pipes(j)%pr(k)
        inp(3) = Reactor_glo%fluid%pipes(j)%pe(k)
        inp(4) = Reactor_glo%fluid%pipes(j)%temp(k)
        inp(5) = Reactor_glo%solid%htstrs(i)%temp(nr) ! right-most node temperature
        inp(6) = Reactor_glo%fluid%pipes(j)%dhyd
        inp(7) = Reactor_glo%fluid%pipes(j)%kpa(k)
        inp(8) = Reactor_glo%fluid%pipes(j)%pval(k)
        ! additional inputs for possible boiling issues
        if(Reactor_glo%fluid%pipes(j)%matid=='H2O') then
          inp(9) = Reactor_glo%fluid%pipes(j)%xe(k)
          inp(10)= Reactor_glo%fluid%pipes(j)%mdot/areaz ! G
        end if
        call qw_cal(Reactor_glo%fluid%pipes(j)%matid, inp, qflux, &
                    Reactor_glo%fluid%pipes(j)%htcmod(k))
        Reactor_glo%solid%htstrs(i)%qflux_o = qflux
        area_ht = 2.d0*pi*ro*dz*mltpl
        Reactor_glo%fluid%pipes(j)%qpow(k) = Reactor_glo%fluid%pipes(j)%qpow(k)+area_ht*qflux
      end if
    end do
    ! fuelrod to fluid heat transfer
    do i = 1, Reactor_glo%solid%nfrod
      frod  = Reactor_glo%solid%frods(i)
      mltpl = frod%mltpl
      nr    = frod%clad%nr
      ro    = frod%clad%ro
      do j = 1, Reactor_glo%fluid%npipe
        if(Reactor_glo%fluid%pipes(j)%idstr==frod%pipeid) then
          exit
        end if
      end do
      areaz = Reactor_glo%fluid%pipes(j)%areaz
      dz    = Reactor_glo%fluid%pipes(j)%dz
      k = frod%pipenode
      inp(1) = Reactor_glo%fluid%pipes(j)%re(k)
      inp(2) = Reactor_glo%fluid%pipes(j)%pr(k)
      inp(3) = Reactor_glo%fluid%pipes(j)%pe(k)
      inp(4) = Reactor_glo%fluid%pipes(j)%temp(k)
      inp(5) = frod%clad%temp(nr)
      inp(6) = Reactor_glo%fluid%pipes(j)%dhyd
      inp(7) = Reactor_glo%fluid%pipes(j)%kpa(k)
      inp(8) = Reactor_glo%fluid%pipes(j)%pval(k)
      ! additional inputs for possible boiling issues
      if(Reactor_glo%fluid%pipes(j)%matid=='H2O') then
        inp(9) = Reactor_glo%fluid%pipes(j)%xe(k)
        inp(10)= Reactor_glo%fluid%pipes(j)%mdot/areaz ! G
      end if
      call qw_cal(Reactor_glo%fluid%pipes(j)%matid, inp, qflux, &
      Reactor_glo%fluid%pipes(j)%htcmod(k))
      Reactor_glo%solid%frods(i)%qflux_co = qflux
      area_ht = 2.d0*pi*ro*dz*mltpl
      Reactor_glo%fluid%pipes(j)%qpow(k) = Reactor_glo%fluid%pipes(j)%qpow(k)+area_ht*qflux
    end do
  end subroutine qflux_w_to_f

  subroutine qvol_fuel()
    !======= Declarations =========
    implicit none
    !== local variables ==
    real(c_double) :: qvavg, kr, kz
    integer(c_int) :: i
    !======= Internals ============
    qvavg = Reactor_glo%core%power/Reactor_glo%solid%fvol_tot
    do i = 1, Reactor_glo%solid%nfrod
      kr = Reactor_glo%solid%frods(i)%kr
      kz = Reactor_glo%solid%frods(i)%kz
      Reactor_glo%solid%frods(i)%qv = qvavg*kr*kz
    end do
  end subroutine qvol_fuel
end module Coupling_mod