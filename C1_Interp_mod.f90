! ----------------------------------------------------------------------
! Programmer(s): Yutong CHEN @XJTU, MER Dr. Konstantin MIKITYUK @PSI
! ----------------------------------------------------------------------
!
! A module for defining Subroutines for interpolation
! 
! ----------------------------------------------------------------------
module Interp_mod
  !======= Inclusions ===========
  use Glodata_mod
  
contains
! ----------------------------------------------------------------------
! interp1d: a subroutine for 1-D linear interpolation
! output yval array corresponding to input xval array based on 
! data arrays xdat and ydat 
! ----------------------------------------------------------------------
  subroutine interp1d(xdat,ydat,xval,yval)
    !======= Declarations =========
    implicit none
    !== input variables ==
    real(c_double),  intent(in) :: xdat(:), ydat(:)   ! Data points
    real(c_double),  intent(in) :: xval(:)            ! Input array
    !== output variables ==
    real(c_double),  intent(out):: yval(:)            ! Output array
    !== local variables ==
    integer(c_int) :: i, n, nlen
    integer(c_int) :: l_indx, r_indx                  ! Left/Right points
    real(c_double) :: dx, xdmin, xdmax, xtemp
    !======= Internals ============
    ! check if input is correct
    nlen = size(xval)
    if(size(yval)/=nlen) then
      stop 'Array size does not match in interp1d'
    end if
    ! main loop
    do i = 1, nlen
      xdmin = xdat(1)+small
      xdmax = xdat(size(xdat))-small
      xtemp = MIN(MAX(xval(i),xdmin), xdmax)
      n = MINLOC(abs(xdat(:)-xtemp), dim=1)
      if(xdat(n)<xtemp)then
        l_indx = n
        r_indx = n+1
      else
        l_indx = n-1
        r_indx = n
      end if
        dx = xdat(r_indx)-xdat(l_indx)
        yval(i) = (xtemp-xdat(l_indx))/dx*ydat(r_indx) + (xdat(r_indx)-xtemp)/dx*ydat(l_indx)
    end do
  end subroutine interp1d

! ----------------------------------------------------------------------
! interp2d: a subroutine for 2-D linear interpolation
! output fval array corresponding to input xval and yval array 
! based on data arrays xdat, ydat and fdat
! ----------------------------------------------------------------------
  subroutine interp2d(xdat, ydat, fdat, xval, yval, fval)
    !======= Declarations =========
    implicit none
    !== input variables ==
    real(c_double), intent(in) :: xdat(:), ydat(:), fdat(:,:)   ! Data arrays
    real(c_double), intent(in) :: xval(:), yval(:)              ! Input
    !== output variables ==
    real(c_double), intent(out):: fval(:)                     ! Output
    !== local variables ==
    real(c_double) :: dx, dy, xdmin, xdmax, xtemp, ydmin, ydmax, ytemp
    real(c_double) :: f_lb, f_rb, f_lt, f_rt
    real(c_double) :: w_lb, w_rb, w_lt, w_rt
    integer(c_int) :: i, n, nlen
    integer(c_int) :: l_indx, r_indx, b_indx, t_indx   ! Left/Right/Bottom/Top points
    !======= Internals ============
    ! check if input is correct
    nlen = size(xval)
    if(size(yval)/=nlen .OR. size(fval)/=nlen) then
      stop 'Array size does not match in interp2d'
    end if
    ! main loop
    do i = 1, nlen
      ! find interval in x direction
      xdmin = xdat(1)+small
      xdmax = xdat(size(xdat))-small
      xtemp = MIN(MAX(xval(i),xdmin), xdmax)
      n = MINLOC(abs(xdat(:)-xtemp), dim=1)
      if(xdat(n)<xtemp)then
        l_indx = n
        r_indx = n+1
      else
        l_indx = n-1
        r_indx = n
      end if
      dx = xdat(r_indx)-xdat(l_indx)
      ! find interval in y direction
      ydmin = ydat(1)+small
      ydmax = ydat(size(ydat))-small
      ytemp = MIN(MAX(yval(i),ydmin), ydmax)
      n = MINLOC(abs(ydat(:)-ytemp), dim=1)
      if(ydat(n)<ytemp)then
        b_indx = n
        t_indx = n+1
      else
        b_indx = n-1
        t_indx = n
      end if
      dy = ydat(t_indx)-ydat(b_indx)
      ! neighbouring data points
      f_lb = fdat(l_indx,b_indx)
      f_rb = fdat(r_indx,b_indx)
      f_lt = fdat(l_indx,t_indx)
      f_rt = fdat(r_indx,t_indx)
      ! weight factors
      w_lb = (xdat(r_indx)-xtemp)*(ydat(t_indx)-ytemp)/dx/dy
      w_rb = (xtemp-xdat(l_indx))*(ydat(t_indx)-ytemp)/dx/dy
      w_lt = (xdat(r_indx)-xtemp)*(ytemp-ydat(b_indx))/dx/dy
      w_rt = (xtemp-xdat(l_indx))*(ytemp-ydat(b_indx))/dx/dy
      ! interpolation
      fval(i) = f_lb*w_lb + f_rb*w_rb + f_lt*w_lt + f_rt*w_rt
    end do
  end subroutine interp2d

end module Interp_mod