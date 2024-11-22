! ----------------------------------------------------------------------
! Programmer(s): Yutong CHEN @XJTU, MER Dr. Konstantin MIKITYUK @PSI
! ----------------------------------------------------------------------
!
! A module for defining Subroutines for matrix operations
! 
! ----------------------------------------------------------------------
module Matrix_mod
  !======= Inclusions ===========
  use Glodata_mod
contains
! ----------------------------------------------------------------------
! bsdet: a subroutine for determinant calculation
! input square matrix a(n,n)
! output corresponding determinant value det_a
! ----------------------------------------------------------------------
  subroutine bsdet(a, det_a)
    !======= Declarations =========
    implicit none
    !== input variables ==
    real(c_double),  intent(in)  :: a(:,:)   ! Input square matrix
    !== output variables ==
    real(c_double),  intent(out) :: det_a    ! Output determinant value
    !== local variables ==
    real(c_double), allocatable :: aa(:,:)   ! temporary copy of matrix a 
    real(c_double) :: f, d, q
    integer(c_int) :: i, j, k, n
    integer(c_int) :: is, js
    !======= Internals ============
    ! initialization
    f = 1.d0
    det_a = 1.d0
    is = 0
    js = 0
    ! check input
    n = size(a(1,:))
    if(size(a(1,:))==size(a(:,1))) then
      n = size(a(1,:))
      allocate(aa(n,n))
      aa = a
    else
      print*, 'input matrix a is not a square matrix, please check'
      stop 'program stop at bsdet'
    end if
    ! main loop
    do k = 1, n-1
      q = 0.d0
      do i = k, n
        do j = k, n
          if(ABS(aa(i,j))>q) then
            q  = ABS(aa(i,j))
            is = i
            js = j
          end if
        end do
      end do
      if(q+1.d0==1.d0) then
        det_a = 0.d0
        return
      end if
      if(is/=k) then
        f =-f
        do j = k, n
          d = aa(k,j)
          aa(k,j) = aa(is,j)
          aa(is,j)= d
        end do
      end if
      if(js/=k) then
        f =-f
        do i = k, n
          d = aa(i,js)
          aa(i,js)= aa(i,k)
          aa(i,k) = d
        end do
      end if
      det_a = det_a*aa(k,k)
      do i = k+1, n
        d = aa(i,k)/aa(k,k)
        do j = k+1, n
          aa(i,j) = aa(i,j)-d*aa(k,j)
        end do
      end do
    end do
    det_a = f*det_a*aa(n,n)
    return
    deallocate(aa)
  end subroutine bsdet
! ----------------------------------------------------------------------
! brinv: a subroutine for inverse matrix computation
! input square matrix a(n,n)
! output corresponding inverse matrix inv_a(n,n)
! ----------------------------------------------------------------------
  subroutine brinv(a, inv_a)
    !======= Declarations =========
    implicit none
    !== input variables ==
    real(c_double),  intent(in)  :: a(:,:)     ! Input square matrix
    !== output variables ==
    real(c_double),  intent(out) :: inv_a(:,:) ! Output inverse matrix
    !== local variables ==
    real(c_double) :: t, d, det
    integer(c_int) :: i, j, k, n
    integer(c_int), allocatable :: is(:), js(:)
    !======= Internals ============
    call bsdet(a, det)
    if(det+1.d0==1.d0) then
      print*, 'input matrix a is singular, please check'
      stop 'program stop at brinv'
    end if
    n = size(a(1,:))
    if(size(inv_a(1,:))/=n .or. size(inv_a(:,1))/=n) then
      print*, 'matrix inv_a size does not match, please check'
      stop 'program stop at brinv'
    end if
    allocate(is(n))
    allocate(js(n))
    inv_a = a
    ! main loop
    do k = 1, n
      d = 0.d0
      do i = k, n
        do j = k, n
          if(ABS(inv_a(i,j))>d) then
            d = ABS(inv_a(i,j))
            is(k) = i
            js(k) = j
          end if
        end do
      end do
      do j = 1, n
        t = inv_a(k,j)
        inv_a(k,j) = inv_a(is(k),j)
        inv_a(is(k),j) = t
      end do
      do i = 1, n
        t = inv_a(i,k)
        inv_a(i,k) = inv_a(i,js(k))
        inv_a(i,js(k)) = t
      end do
      inv_a(k,k) = 1.d0/inv_a(k,k)
      do j = 1, n
        if(j/=k) then
          inv_a(k,j) = inv_a(k,j)*inv_a(k,k)
        end if
      end do
      do i = 1, n
        if(i/=k) then
          do j = 1, n
            if(j/=k) then
              inv_a(i,j) = inv_a(i,j)-inv_a(i,k)*inv_a(k,j)
            end if
          end do
        end if
      end do
      do i = 1, n
        if(i/=k) then
          inv_a(i,k) = -inv_a(i,k)*inv_a(k,k)
        end if
      end do
    end do
    do k = n, 1, -1
      do j = 1, n
        t = inv_a(k,j)
        inv_a(k,j)=inv_a(js(k),j)
        inv_a(js(k),j)=t
      end do
      do i = 1, n
        t = inv_a(i,k)
        inv_a(i,k) = inv_a(i,is(k))
        inv_a(i,is(k)) = t
      end do
    end do
    return
  end subroutine brinv
end module Matrix_mod