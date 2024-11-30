! ----------------------------------------------------------------------
! Programmer(s): Yutong CHEN @XJTU, MER Dr. Konstantin MIKITYUK @PSI
! ----------------------------------------------------------------------
!
! A module for defining auxiliary correlations
! 
! ----------------------------------------------------------------------
module Correlate_mod
  !======= Inclusions =========
  use Interp_mod
  use Getdata_mod
  use Matrix_mod
  !======= Declarations =========
  implicit none
  
contains
! ----------------------------------------------------------------------
! matpro_s: a subroutine for solid material properties calculation
! input  - solid material id and temperature
! output - solid material properties
! all data need to be arranged into array format
! ----------------------------------------------------------------------
  subroutine matpro_s(matid, temp, rho, cp, kpa)
    !======= Declarations =========
    implicit none
    !== input variables ==
    ! material id
    character(len=20) :: matid
    ! temperature
    real(c_double), intent(in) :: temp(:)
    !== output variables ==
    ! density, heat_capacity, heat_conductivity
    real(c_double), intent(out) :: rho(:), cp(:), kpa(:)
    !== local variables ==
    real(c_double), allocatable :: tc(:), tf(:)
    integer(c_int) :: i, nlen
    ! coefficients for correlations
    real(c_double) :: a0, a1, a2, a3, a4
    real(c_double) :: b1, b2, b3, b4, b5, b6
    real(c_double) :: c1, c2, c3, c4, c5, c6, c7
    
    !======= Internals ============
    ! check if input is correct
    nlen = size(temp)
    if(size(rho)/=nlen .or. size(cp)/=nlen .or. size(kpa)/=nlen) then
      stop 'Array size does not match in matpro_s'
    end if
    ! alllcate temporary arrays
    allocate(tc(nlen))
    allocate(tf(nlen))
    tc = temp - 273.15d0   ! Celsius temperature array
    tf = tc*1.8d0 + 32.d0  ! Fahrenheit temperature array
    select case (matid)
      case('ss316') ! Stainless Steel, Type 316, developed based on TRACE manual
        ! density
        rho = 7954.d0
        ! heat capacity
        do i = 1, nlen
          ! TRACE model limitations
          tf(i) = MIN(tf(i), 2400.d0)
        end do
        a0 = 426.17d0
        a1 = 4.3816d-1
        a2 =-6.3759d-4
        a3 = 4.4803d-7
        a4 =-1.0729d-10
        cp  = a0 + a1*tf + a2*tf**(2.d0) + a3*tf**(3.d0) + a4*tf**(4.d0)
        ! heat conductivity
        kpa = 9.248d0 + 1.571d-2*temp
      case('uo2') ! UO2 Fuel, developed based on TRACE manual
        ! density
        rho = 10980.d0
        ! heat capacity
        b1 = 19.1450d0
        b2 = 7.84730d-4
        b3 = 5.64370d6
        b4 = 535.285d0
        b5 = 37694.6d0
        b6 = 1.98700d0
        cp  = b1*b4**(2.d0)*EXP(b4/temp)/(temp*(EXP(b4/temp)-1.d0))**(2.d0)
        cp  = cp + 2.d0*b2*temp + b3*b5*EXP(-b5/b6/temp)/b6/temp**(2.d0)
        cp  = cp*15.496d0
        ! heat conductivity - Legacy Model of TRACE manual
        c1 = 40.40d0
        c2 = 464.0d0
        c3 = 1.216d-4
        c4 = 1.867d-3
        c5 = 1.910d-2
        c6 = 2.580d0
        c7 =-5.800d-4
        do i = 1, nlen
          kpa(i) = c1/(c2 + MIN(tc(i), 1650.d0)) + c3*EXP(c4*tc(i))
        end do
        kpa = 100.d0*kpa
      case('bn') ! Boron Nitride insulators, developed based on TRACE manual
        ! density
        rho = 2002.0d0
        ! heat capacity
        a0 = 760.59d0
        a1 = 1.7955d0
        a2 =-8.6704d-4
        a3 = 1.5896d-7
        cp  = a0 + a1*tf + a2*tf**(2.d0) + a3*tf**(3.d0)
        ! heat conductivity
        kpa = 25.27d0 - 1.365d-3*tf
      case('nicr') ! Constantan/Nichrome heater coil, developed based on TRACE manual
        ! density
        rho = 8393.4d0
        ! heat capacity
        cp  = 110.d0/tf**(0.2075d0)
        ! heat conductivity
        kpa = 29.18d0 + 2.683d-3*(tf - 100.d0)
      case('cu') ! Copper, developed based on NACIE-UP benchmark, NA-I-R-542, Feb. 2023
        ! density
        rho = 8933.d0
        ! heat capacity
        cp  = 385.0d0
        ! heat conductivity
        kpa = 401.0d0
      case('powder') ! SS powser, developed based on NACIE-UP benchmark, NA-I-R-542, Feb. 2023
        ! density
        rho = 7954.d0
        ! heat capacity
        cp  = (6.181d0 + 1.788d-3*temp)*10.165d0*4.184d0
        ! heat conductivity
        kpa = max(0.3d0 + 5.d-3*(tc - 200.d0), 0.3d0)
      case('mgo') ! Mgo material, constants from google
        ! density
        rho = 3580.d0
        ! heat capacity
        cp  = 880.0d0
        ! heat conductivity
        kpa = 3.d0 ! 8-30 W/(m*K)
      case default
        print*, 'no match for material name = ', matid
        stop 'program stopped at matpro_s'
    end select  
    ! using 'deallocate' is crucial for dynamic memory management, since fortran dose not have 
    ! automatic garbage collection mechanism like high-level languages such as Java or Python
    deallocate(tc)
    deallocate(tf)
  end subroutine matpro_s
! ----------------------------------------------------------------------
! fh_lbe: a subroutine for calculating thermal properties of LBE fluid
! input  - LBE enthalpy
! output - LBE thermal-physical properties
! all data need to be arranged into array format
! ----------------------------------------------------------------------
  subroutine fh_lbe(enth, temp, rho, miu, niu, cp, kpa)
    !======= Declarations =========
    implicit none
    !== input variables ==
    ! LBE enthalpy
    real(c_double), intent(in) :: enth(:)
    !== output variables ==
    ! LBE temperature, density, 
    real(c_double), intent(out) :: temp(:), rho(:)
    ! dynamic/kinetic viscosity, heat capacity, heat conductivity
    real(c_double), intent(out) :: miu(:), niu(:), cp(:), kpa(:)
    
    !======= Internals ============
    ! interpolate to get thermal-physical data          !col1 - enthalpy
    call interp1d(tpf_lbe(:,1),tpf_lbe(:,2),enth,temp)  !col2 - temperature
    call interp1d(tpf_lbe(:,1),tpf_lbe(:,3),enth,rho)   !col3 - density
    call interp1d(tpf_lbe(:,1),tpf_lbe(:,4),enth,miu)   !col4 - dynamic viscosity
    call interp1d(tpf_lbe(:,1),tpf_lbe(:,5),enth,niu)   !col5 - kinetic viscosity
    call interp1d(tpf_lbe(:,1),tpf_lbe(:,6),enth,cp)    !col6 - heat capacity
    call interp1d(tpf_lbe(:,1),tpf_lbe(:,7),enth,kpa)   !col7 - heat conductivity
  end subroutine fh_lbe
! ----------------------------------------------------------------------
! ht_lbe: a subroutine for calculating enthalpy of LBE fluid
! input  - LBE temperature
! output - LBE enthalpy
! all data need to be arranged into array format
! ----------------------------------------------------------------------
  subroutine ht_lbe(temp, enth)
    !======= Declarations =========
    implicit none
    !== input variables ==
    ! fluid temperature
    real(c_double), intent(in) :: temp(:)
    !== output variables ==
    ! fluid enthalpy
    real(c_double), intent(out) :: enth(:)
    
    !======= Internals ============
    ! interpolate to get thermal-physical data
    call interp1d(tpf_lbe(:,2),tpf_lbe(:,1),temp,enth)  ! temperature
    
  end subroutine ht_lbe
! ----------------------------------------------------------------------
! fh_na1: a subroutine for calculating thermal properties of Na1 fluid
! input  - Na1 enthalpy
! output - Na1 thermal-physical properties
! all data need to be arranged into array format
! ----------------------------------------------------------------------
  subroutine fh_na1(enth, temp, rho, miu, niu, cp, kpa)
    !======= Declarations =========
    implicit none
    !== input variables ==
    ! Na1 enthalpy
    real(c_double), intent(in) :: enth(:)
    !== output variables ==
    ! Na1 temperature, density, 
    real(c_double), intent(out) :: temp(:), rho(:)
    ! dynamic/kinetic viscosity, heat capacity, heat conductivity
    real(c_double), intent(out) :: miu(:), niu(:), cp(:), kpa(:)
    
    !======= Internals ============
    ! interpolate to get thermal-physical data          !col1 - enthalpy
    call interp1d(tpf_na1(:,1),tpf_na1(:,2),enth,temp)  !col2 - temperature
    call interp1d(tpf_na1(:,1),tpf_na1(:,3),enth,rho)   !col3 - density
    call interp1d(tpf_na1(:,1),tpf_na1(:,4),enth,miu)   !col4 - dynamic viscosity
    call interp1d(tpf_na1(:,1),tpf_na1(:,5),enth,niu)   !col5 - kinetic viscosity
    call interp1d(tpf_na1(:,1),tpf_na1(:,6),enth,cp)    !col6 - heat capacity
    call interp1d(tpf_na1(:,1),tpf_na1(:,7),enth,kpa)   !col7 - heat conductivity
  end subroutine fh_na1
! ----------------------------------------------------------------------
! ht_na1: a subroutine for calculating enthalpy of Na1 fluid
! input  - Na1 temperature
! output - Na1 enthalpy
! all data need to be arranged into array format
! ----------------------------------------------------------------------
  subroutine ht_na1(temp, enth)
    !======= Declarations =========
    implicit none
    !== input variables ==
    ! fluid temperature
    real(c_double), intent(in) :: temp(:)
    !== output variables ==
    ! fluid enthalpy
    real(c_double), intent(out) :: enth(:)
    
    !======= Internals ============
    ! interpolate to get thermal-physical data
    call interp1d(tpf_na1(:,2),tpf_na1(:,1),temp,enth)  ! temperature
    
  end subroutine ht_na1
! ----------------------------------------------------------------------
! fph_h2o: a subroutine for calculating thermal properties of H2O fluid
! input  - H2O pressure and enthalpy
! output - H2O thermal-physical properties
! all data need to be arranged into array format
! ----------------------------------------------------------------------
  subroutine fph_h2o(pval, enth, rho, niu, miu, kpa, cp, temp, xe)
    !======= Declarations =========
    implicit none
    !== input variables ==
    ! water pressure, enthalpy
    real(c_double), intent(in) :: pval(:), enth(:)
    !== output variables ==
    ! water density, kinetic/dynamic viscosity
    real(c_double), intent(out):: rho(:), niu(:), miu(:)
    ! heat conductivity, heat capacity, temperature, gas quality
    real(c_double), intent(out):: kpa(:), cp(:), temp(:), xe(:)
    
    !======= Internals ============
    ! interpolate to get thermal-physical data          
    call interp2d(p_h2o,h_h2o,tpf_h2o(1,:,:),pval,enth,rho)    !col1 - density
    call interp2d(p_h2o,h_h2o,tpf_h2o(2,:,:),pval,enth,niu)    !col2 - niu
    call interp2d(p_h2o,h_h2o,tpf_h2o(3,:,:),pval,enth,miu)    !col3 - miu
    call interp2d(p_h2o,h_h2o,tpf_h2o(4,:,:),pval,enth,kpa)    !col4 - kpa
    call interp2d(p_h2o,h_h2o,tpf_h2o(5,:,:),pval,enth,cp)     !col5 - cp
    call interp2d(p_h2o,h_h2o,tpf_h2o(6,:,:),pval,enth,temp)   !col6 - temp
    call interp2d(p_h2o,h_h2o,tpf_h2o(7,:,:),pval,enth,xe)     !col7 - xe
    
  end subroutine
! ----------------------------------------------------------------------
! hpt_h2o: a subroutine for calculating enthalpy of H2O fluid
! input  - H2O pressure and temperature
! output - H2O enthalpy
! all data need to be arranged into array format
! ----------------------------------------------------------------------
  subroutine hpt_h2o(pval, temp, enth)
    !======= Declarations =========
    implicit none
    !== input variables ==
    ! water pressure, temperature
    real(c_double), intent(in) :: pval(:), temp(:)
    !== output variables ==
    ! water enthalpy
    real(c_double), intent(out) :: enth(:)
    !== local variables ==
    real(c_double), allocatable :: tsat(:)
    real(c_double) :: ptmp(2), etmp(2), ttmp(2), ctmp(1)
    real(c_double) :: res1, res2
    integer(c_int) :: i, j, nlen

    !======= Internals ============
    ! interpolate to get thermal-physical data
    nlen = size(pval)
    allocate(tsat(nlen))
    call interp1d(p_h2o,sat_h2o(9,:),pval,tsat)
    ! loop
    do i = 1, nlen
      ! initial guess of enthalpy
      if(temp(i)<=tsat(i)) then
        ! subcooled water
        call interp1d(sat_h2o(9,:),sat_h2o(10,:),temp(i:i),enth(i:i))
      else
        ! superheated steam
        call interp1d(p_h2o,sat_h2o(11,:),pval(i:i),enth(i:i))
        call interp1d(p_h2o,sat_h2o( 6,:),pval(i:i),ctmp)
        ! initial guess of enthalpy
        enth(i) = enth(i) + 0.1d0*ctmp(1)*(temp(i)-tsat(i))
      end if
      ! obtain true enth value with an iterative approach
      ptmp = pval(i)
      do j = 1, 10
        ! Newton algorithm
        etmp(1) = enth(i)
        etmp(2) = enth(i) + 1.d0
        call interp2d(p_h2o,h_h2o,tpf_h2o(6,:,:),ptmp,etmp,ttmp)
        res1 = temp(i) - ttmp(1)
        res2 = temp(i) - ttmp(2)
        enth(i) = enth(i)-res1/((res2-res1)/1.d0)
      end do
    end do
    deallocate(tsat)
  end subroutine hpt_h2o
! ----------------------------------------------------------------------
! fp_h2o: a subroutine for calculating saturated properties of h2o fluid
! input  - H2O pressure
! output - H2O saturated properties
! all data need to be arranged into array format
! ----------------------------------------------------------------------
  subroutine fp_h2o(pval, satpp)
    !======= Declarations =========
    implicit none
    !== input variables ==
    ! water pressure
    real(c_double), intent(in) :: pval
    !== output variables ==
    ! saturated properties
    real(c_double), intent(out) :: satpp(13)
    !== local variables ==
    real(c_double) :: parr(1)
    !======= Internals ============
    parr(1) = pval
    call interp1d(p_h2o,sat_h2o( 1,:),parr,satpp( 1: 1)) !  1-rhol
    call interp1d(p_h2o,sat_h2o( 2,:),parr,satpp( 2: 2)) !  2-rhog
    call interp1d(p_h2o,sat_h2o( 3,:),parr,satpp( 3: 3)) !  3-miul
    call interp1d(p_h2o,sat_h2o( 4,:),parr,satpp( 4: 4)) !  4-miug
    call interp1d(p_h2o,sat_h2o( 5,:),parr,satpp( 5: 5)) !  5-cpl
    call interp1d(p_h2o,sat_h2o( 6,:),parr,satpp( 6: 6)) !  6-cpg
    call interp1d(p_h2o,sat_h2o( 7,:),parr,satpp( 7: 7)) !  7-kl
    call interp1d(p_h2o,sat_h2o( 8,:),parr,satpp( 8: 8)) !  8-kg
    call interp1d(p_h2o,sat_h2o( 9,:),parr,satpp( 9: 9)) !  9-ts
    call interp1d(p_h2o,sat_h2o(10,:),parr,satpp(10:10)) ! 10-hls
    call interp1d(p_h2o,sat_h2o(11,:),parr,satpp(11:11)) ! 11-hgs
    call interp1d(p_h2o,sat_h2o(12,:),parr,satpp(12:12)) ! 12-hgl
    call interp1d(p_h2o,sat_h2o(13,:),parr,satpp(13:13)) ! 13-sgm
  end subroutine 
! ----------------------------------------------------------------------
! qw_cal: a subroutine for calculating wall-to-fluid heat flux with 
! classical Dittus-Boelter correlation
! input  - 
! output - 
! all data need to be arranged into array format
! ----------------------------------------------------------------------
  subroutine qw_cal(matfld, inp, qflux, htcmod)
    !======= Declarations =========
    implicit none
    !== input variables ==
    character(len=20), intent(in):: matfld  ! material id
    real(c_double), intent(in)   :: inp(:)  ! input array
    !== output variables ==
    real(c_double), intent(out)  :: qflux   ! heat flux
    integer(c_int), intent(out)  :: htcmod  ! flag for heat transfer mode
    !== local variables ==
    real(c_double) :: re, pr, pe, tf, tw, dh, kf, pp, xe, gg
    real(c_double) :: nu, htc
    real(c_double) :: satval(13)
    real(c_double) :: rhog, rhol, miul, miug, cpl, cpg, kl, kg, ts
    real(c_double) :: hls, hgs, hgl, sgm
    real(c_double) :: alfa, reg, prg, xx, xcr, yy, phy
    real(c_double) :: tchf, tmfb
    real(c_double) :: qsp, qnb, qchf, qmfb, qdry
    real(c_double) :: ntmp
    !======= Internals ============
    re = inp(1) ! Reynolds number    [-]
    pr = inp(2) ! Prandtl number     [-]
    pe = inp(3) ! Prandtl number     [-]
    tf = inp(4) ! Fluid temperature  [K]
    tw = inp(5) ! Wall temperature   [K]
    dh = inp(6) ! Hydraulic diameter [m]
    kf = inp(7) ! Fluid heat conductivity  [W/(m⋅K)]
    pp = inp(8) ! Cell pressure      [Pa]
    
    select case (matfld)
      case('LBE','Na1')
        ! liquid metal, round tube
        htcmod = 0  ! Default - Single phase liquid convection (SPL)
        nu = 4.8d0 + 0.025d0*(pe)**(0.8d0) ! Lyon-Martinelli correlation (TRACE default)
        htc = kf*nu/dh      ! heat transfer coefficient
        qflux = htc*(tw-tf) ! wall-to-fluid heat flux
      case('H2O')
        xe = inp(9) ! Gas quality [-]
        if(xe>=1.d0) then
          htcmod = 7 ! SPV regime
          nu = 0.023d0*(re)**(0.8d0)*(pr)**(0.33d0)*(tf/tw)**(0.5d0)
          nu = Max(4.33d0, nu) ! heat conductivity effect
          htc = kf*nu/dh     ! heat transfer coefficient
          qflux = htc*(tw-tf) ! wall-to-fluid heat flux
        else
          ! xe<1.d0
          call fp_h2o(pp,satval)
          rhol = satval( 1) 
          rhog = satval( 2) 
          miul = satval( 3) 
          miug = satval( 4) 
          cpl  = satval( 5) 
          cpg  = satval( 6) 
          kl   = satval( 7) 
          kg   = satval( 8) 
          ts   = satval( 9) 
          hls  = satval(10) 
          hgs  = satval(11) 
          hgl  = satval(12) 
          sgm  = satval(13)
          if(xe<0.d0) then 
            ! Subcooled flow
            if(tw<=ts) then
              htcmod = 0 ! SPL regime
              nu = Max(0.023d0*(re)**(0.8d0)*(pr)**(0.33d0),4.33d0)
              htc = kf*nu/dh
              qflux = htc*(tw-tf) ! wall-to-fluid heat flux
            else
              qnb = 2.d3*exp(pp/4.34d6)*(tw-ts)**(2.0)
              nu = Max(0.023d0*(re)**(0.8d0)*(pr)**(0.33d0),4.33d0)
              htc = kf*nu/dh     ! heat transfer coefficient
              qsp = htc*(tw-tf) ! wall-to-fluid heat flux
              if(qsp>qnb) then
                htcmod = 0 ! SPL regime
                qflux = qsp
              else
                ! tw >= tonb
                qchf = 0.14*hgl*(sgm*9.8*rhog**(2.0)*(rhol-rhog))**(0.25)
                tchf = ts + (qchf*exp(-pp/4.34d6)/2.d3)**(0.50)
                if(tw<tchf) then
                  htcmod = 1 ! SCNB regime
                  qflux = qnb
                else
                  ! tw >= tchf
                  if(pp<9.d6) then
                    tmfb = 557.90 + 44.1*(pp/1.d6) - 3.72*(pp/1.d6)**(2.0)
                  else
                    tmfb = 647.09 + 0.71*(pp/1.d6)
                  end if
                  if(tw<tmfb) then
                    htcmod = 2 ! SCTRAN regime
                    qmfb = 0.14*hgl*(kg**(2.0)*cpg*rhog*9.8*(rhol-rhog)/miug)**(1.0/3.0)*(tmfb-ts)
                    ntmp = log((tmfb-ts)/(tw-ts))/log((tmfb-ts)/(tchf-ts))
                    qflux = qmfb*(qchf/qmfb)**(ntmp)
                  else
                    htcmod = 3 ! SCIAFB regime
                    qflux = 0.14*hgl*(kg**(2.0)*cpg*rhog*9.8*(rhol-rhog)/miug)**(1.0/3.0)*(tw-ts)
                  end if
                end if 
              end if
            end if
          else
            ! here 0.0<xe<1.d0
            gg = inp(10)
            alfa = rhol*xe/(rhol*xe+rhog*(1.0-xe)) ! void fraction
            reg= gg*dh/miug
            prg= miug*cpg/kg
            xx = pp/1.d6/9.8d0
            xcr= (0.39+1.57*xx-2.04*xx**(2.0)+0.68*xx**(3.0))*(gg/1.d3)**(-0.5)*(dh/8.d-3)**(0.15)
            xcr= Min(xcr,0.99d0)
            if(xe<xcr) then
              ! here 0.0<xe<xcr, no dry out
              qchf = 0.14*hgl*(sgm*9.8*rhog**(2.0)*(rhol-rhog))**(0.25)*(1.d0-alfa)
              tchf=ts+(qchf*exp(-pp/4.34d6)/2.d3)**(0.50)
              if(tw<tchf) then
                htcmod = 4 ! NB regime
                qflux = 2.d3*exp(pp/4.34d6)*(tw-ts)**(2.0)
              else
                if(pp<9.d6) then
                  tmfb = 557.90 + 44.1*(pp/1.d6) - 3.72*(pp/1.d6)**(2.0)
                else
                  tmfb = 647.09 + 0.71*(pp/1.d6)
                end if
                if(tw<tmfb) then
                  htcmod = 5 ! TRAN regime
                  qmfb = 0.14*hgl*(kg**(2.0)*cpg*rhog*9.8*(rhol-rhog)/miug)**(1.0/3.0)*(tmfb-ts)
                  ntmp = log((tmfb-ts)/(tw-ts))/log((tmfb-ts)/(tchf-ts))
                  qflux= qmfb*(qchf/qmfb)**(ntmp)
                else
                  htcmod = 6 ! IAFB regime
                  qflux= 0.14*hgl*(kg**(2.0)*cpg*rhog*9.8*(rhol-rhog)/miug)**(1.0/3.0)*(tmfb-ts)
                end if
                ! heat flux at critical dry out point
                yy = 1.0-0.1*((rhol/rhog-1.0)*(1.0-xcr))**(0.40)
                alfa = rhol*xcr/(rhol*xcr+rhog*(1.0-xcr))
                nu = 0.023d0*(reg*(xcr/alfa))**(0.8d0)*prg**(0.33)*yy*(tf/tw)**(0.50)
                htc= nu*kg/dh
                qdry = htc*(tw-tf)
                phy= exp(1.d0*xe/(xe-xcr))
                ! interpolation for robustness of calculation
                qflux = qflux*phy + (1.0-phy)*qdry
              end if
            else
              ! here xcr<=xe, dry out
              htcmod = 7
              yy = 1.0-0.1*((rhol/rhog-1.0)*(1.0-xe))**(0.40)
              nu = 0.023d0*(reg*(xe/alfa))**(0.8d0)*prg**(0.33)*yy*(tf/tw)**(0.50)
              htc= kg*nu/dh
              qflux = htc*(tw-tf)
            end if
          end if
        end if
        
        
        
      case default
        print*, 'no match for fluid = ', matfld
        stop 'program stopped at qw_cal'
    end select
    
    
  end subroutine qw_cal
  
  subroutine fric_cal(cate, re0, fric, inp)
    !======= Declarations =========
    implicit none
    !== input variables ==
    character(len=20), intent(in) :: cate     ! pipe category
    real(c_double),    intent(in) :: re0(:) ! input re array (not editable)
    !== output variables ==
    real(c_double), intent(out)   :: fric(:)  ! friction factor
    !== optional variables ==
    real(c_double), intent(in), optional :: inp(:) ! input array
    !== local variables ==
    real(c_double), allocatable :: re(:)
    real(c_double), allocatable :: a(:), b(:)
    integer(c_int) :: nnodes, i
    !======= Internals ============
    nnodes = size(re0)
    allocate(re(nnodes))
    do i = 1, nnodes
      re(i) = max(re0(i), 1.d0) ! to avoid numerical error
    end do
    select case(cate)
      case('round')
        allocate(a(nnodes))
        allocate(b(nnodes))
        ! ChurchillS
        a = (2.457d0*LOG(1.d0/((7.d0/re)**(0.9d0))))**(16.d0)
        b = (37530.d0/re)**(16.d0)
        fric = 8.d0*((8.d0/re)**(12.d0)+1.d0/(a+b)**1.5d0)**(1.d0/12.d0)
        deallocate(a)
        deallocate(b)
      case('brods')
        if(.not.present(inp)) print*, 'lacking input'
        print*, 'brods to be developed'
        stop 'program stopped at fric_cal'
      case('wrods')
        if(.not.present(inp)) print*, 'lacking input'
        print*, 'wrods to be developed'
        stop 'program stopped at fric_cal'
      case('freelevel','fill','break')
        print*, 'fric_cal not needed for cate = ', cate
        stop 'program stopped at fric_cal'
      case default
        print*, 'no match for cate = ', cate
        stop 'program stopped at fric_cal'
    end select
  end subroutine fric_cal

end module Correlate_mod
