! ----------------------------------------------------------------------
! Programmer(s): Yutong CHEN @XJTU, MER Dr. Konstantin MIKITYUK @PSI
! ----------------------------------------------------------------------
!
! A module for defining Subroutines reading fluid data from .dat files
! 
! ----------------------------------------------------------------------
module Getdata_mod
  !======= Inclusions ===========
  use Glodata_mod
  
contains
! ----------------------------------------------------------------------
! get_lbe_dat: a subroutine reading thermal property datapoints for LBE
! fluid, filling global tpf_lbe(:,:) array containing list of enthalpy, 
! temperature, density, dynamic_viscosity, kinetic_viscosity, 
! heat_capacity and heat_conductivity
! ----------------------------------------------------------------------
  subroutine get_lbe_dat()
    !======= Declarations =========
    implicit none
    !== input variables ==
    !== output variables ==
    !== local variables ==
    real(c_double) :: c_real
    integer(c_int) :: i, filelen

    !======= Internals ============
    print*, 'Reading LBE property data from tpf_lbe.dat file...'
    ! open lbe fluid data file
    open(unit=lbefu, status='old', file='tpf_lbe.dat', form='formatted')
    filelen = 0
    ! loop over the .dat file to get array length
    do while(.TRUE.)
      read(lbefu, '(F12.5)', end=100) c_real
      filelen = filelen + 1
    end do
    100 continue
    ! read data into global tpf_lbe(:,:) array
    ! column 1-enthalpy, 2-temperature, 3-density, 4-dynamic_viscosity,
    ! 5-kinetic_viscosity, 6-heat_capacity, 7-heat_conductivity
    allocate(tpf_lbe(filelen,7))
    rewind(lbefu)
    do i = 1, filelen
      read(lbefu, '(7F12.5)') tpf_lbe(i,1),tpf_lbe(i,2),tpf_lbe(i,3), &
      tpf_lbe(i,4),tpf_lbe(i,5),tpf_lbe(i,6),tpf_lbe(i,7)
    end do
    close(lbefu)
    print*, 'Done.'
  end subroutine get_lbe_dat
! ----------------------------------------------------------------------
! get_na1_dat: a subroutine reading thermal property datapoints for Na1
! fluid, filling global tpf_na1(:,:) array containing list of enthalpy, 
! temperature, density, dynamic_viscosity, kinetic_viscosity, 
! heat_capacity and heat_conductivity
! ----------------------------------------------------------------------
  subroutine get_na1_dat()
    !======= Declarations =========
    implicit none
    !== input variables ==
    !== output variables ==
    !== local variables ==
    real(c_double) :: c_real
    integer(c_int) :: i, filelen

    !======= Internals ============
    print*, 'Reading Na1 property data from tpf_na1.dat file...'
    ! open na1 fluid data file
    open(unit=na1fu, status='old', file='tpf_na1.dat', form='formatted')
    filelen = 0
    ! loop over the .dat file to get array length
    do while(.TRUE.)
      read(na1fu, '(F12.5)', end=133) c_real
      filelen = filelen + 1
    end do
    133 continue
    ! read data into global tpf_na1(:,:) array
    ! column 1-enthalpy, 2-temperature, 3-density, 4-dynamic_viscosity,
    ! 5-kinetic_viscosity, 6-heat_capacity, 7-heat_conductivity
    allocate(tpf_na1(filelen,7))
    rewind(na1fu)
    do i = 1, filelen
      read(na1fu, '(7F12.5)') tpf_na1(i,1),tpf_na1(i,2),tpf_na1(i,3), &
      tpf_na1(i,4),tpf_na1(i,5),tpf_na1(i,6),tpf_na1(i,7)
    end do
    close(na1fu)
    print*, 'Done.'
  end subroutine get_na1_dat
! ----------------------------------------------------------------------
! get_h2o_dat: a subroutine reading thermal property datapoints for h2o
! fluid, filling global tpf_h2o(:,:) array containing list of enthalpy, 
!  
!  
!  
! ----------------------------------------------------------------------
  subroutine get_h2o_dat()
    !======= Declarations =========
    implicit none
    !== input variables ==
    !== output variables ==
    !== local variables ==
    real(c_double) :: c_real
    real(c_double) :: p1, p2
    integer(c_int) :: i, j, filelen
    integer(c_int) :: lenp, lenh     ! p, h array length
    !======= Internals ============
    print*, 'Reading H2O property data from tpf_h2o.dat & sat_h2o.dat file...'
    ! open h2o fluid data file
    open(unit=h2ofu, status='old', file='tpf_h2o.dat', form='formatted')
    open(unit=h2osu, status='old', file='sat_h2o.dat', form='formatted')
    filelen = 0
    ! loop over the .dat file to get array length
    do while(.TRUE.)
      read(h2ofu, '(F12.5)', end=200) c_real
      filelen = filelen + 1
    end do
    200 continue
    ! get the length of h array
    lenh = 0 
    rewind(h2ofu)
    do while(.True.)
      read(h2ofu,'(F12.5)') p1
      if(lenh == 0) then
        p2 = p1
        lenh = lenh + 1
      else
        if(p1 == p2) then
          lenh = lenh + 1
        else
          exit
        end if
      end if
    end do
    ! the length of p array
    lenp = filelen/lenh
    ! allocate data arrays
    allocate(p_h2o(lenp))
    allocate(h_h2o(lenh))
    allocate(tpf_h2o( 7,lenp,lenh))
    allocate(sat_h2o(13,lenp))
    ! read data into global tpf_h2o(7,lenp,lenh) and sat_h2o(11,lenp) array
    ! tpf_h2o(7,lenp,lenh):
    ! column 1-density, 2-kinetic_viscosity, 3-dynamic_viscosity, 4-heat_conductivity,
    ! 5-heat_capacity, 6-temperature, 7-gas_mass_quality
    ! sat_h2o(11,lenp):
    ! column 1-rhol, 2-rhog, 3-miul, 4-miug, 5-cpl, 6-cpg, 7-kpal, 8-kpag, 9-tsat, 
    ! 10-hls, 11-hgs, 12-hgl, 13-sigma
    rewind(h2ofu)
    rewind(h2osu)
    do i = 1, lenp
      read(h2osu,'(13F12.5)') sat_h2o(1,i), sat_h2o(2,i), sat_h2o(3,i),     &
      sat_h2o(4,i), sat_h2o(5,i), sat_h2o(6,i), sat_h2o(7,i), sat_h2o(8,i), &
      sat_h2o(9,i), sat_h2o(10,i),sat_h2o(11,i),sat_h2o(12,i),sat_h2o(13,i)
      do j = 1, lenh
        read(h2ofu,'(9F12.5)')p_h2o(i),h_h2o(j),tpf_h2o(1,i,j),&
        tpf_h2o(2,i,j), tpf_h2o(3,i,j), tpf_h2o(4,i,j), tpf_h2o(5,i,j),  &
        tpf_h2o(6,i,j), tpf_h2o(7,i,j)
      end do
    end do
    close(h2ofu)
    close(h2osu)
    print*, 'Done.'
  end subroutine get_h2o_dat
end module Getdata_mod
