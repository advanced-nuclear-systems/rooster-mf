! ----------------------------------------------------------------------
! Programmer(s): Yutong CHEN @XJTU, MER Dr. Konstantin MIKITYUK @PSI
! ----------------------------------------------------------------------
!
! A module for defining type fluid
! 
! ----------------------------------------------------------------------
module Fluid_mod
  !======= Inclusions =========
  use Pipe_mod
  use Junction_mod
  
  !======= Declarations =========
  implicit none
  
  type :: Fluid
    !======= Member Variables =========
    ! matrix for solving junction mass flow rate
    real(c_double), allocatable :: inv_A(:,:)
    ! array 'mdots' stores mdot of all junctions
    real(c_double), allocatable :: mdots(:)
    ! matrix for solving cell pressure
    real(c_double), allocatable :: inv_B(:,:)
    ! number of pipes, external junctions
    integer(c_int) :: npipe, njunc
    ! number of pipes whose h solved by ode (pipes except break,fill)
    integer(c_int) :: npipe_h = 0
    ! number of internal junctions (pipes except break,fill,freelevel)
    integer(c_int) :: njint = 0
    ! total number of pipe nodes (number of scalar nodes)
    integer(c_int) :: pntot = 0
    ! array storing index of pipes whose enthalpy calculated by ode
    integer(c_int), allocatable :: pidx_h(:)
    ! array storing pipe start/end index in enthalpy rhs array
    integer(c_int), allocatable :: indx1(:), indx2(:)
    ! array storing pipe objects
    type(Pipe),     allocatable :: pipes(:)
    ! array storing junction objects
    type(Junc), allocatable :: juncs(:)

  contains
    !======= Member Procedures ========
    procedure :: init_fld_p      ! constructor for setting pipes array
    procedure :: init_fld_j      ! constructor for setting junctions array
    procedure :: init_fluid      ! check inputs and initialize fluid field
    procedure :: comp_rhs_fld1   ! update mdot and pressure field, calculate dmdot/dt
    procedure :: comp_rhs_fld2   ! calculate fluid dh/dt
    !compose rhs array
  end type Fluid

contains
  
  subroutine init_fld_p(self, pipe_array)
    !======= Declarations =========
    implicit none
    !== I/O variables ==
    class(Fluid), intent(inout) :: self
    class(Pipe), intent(in)     :: pipe_array(:)
    !== local variables ==
    integer(c_int) :: i, j
    !======= Internals ============
    self%npipe = size(pipe_array)
    ! allocate and assign the arrays storing pipe objects
    allocate(self%pipes(self%npipe))
    self%pipes = pipe_array

    ! loop through the pipes
    do i = 1, self%npipe
      ! calculate the total number of pipe nodes solved by ode
      self%pntot = self%pntot + self%pipes(i)%nz
      ! get neq_flen
      if(self%pipes(i)%cate=='freelevel') then
        neq_flen = neq_flen+1
      end if
      ! get fluid%njint
      if(self%pipes(i)%cate/='freelevel'.and.self%pipes(i)%cate/='break'.and. &
         self%pipes(i)%cate/='fill') then
        self%njint = self%njint+1
      end if
      ! get neq_enth and npipe_h
      if(self%pipes(i)%cate/='break'.and.self%pipes(i)%cate/='fill') then
        neq_enth = neq_enth + self%pipes(i)%nz
        self%npipe_h  = self%npipe_h + 1
      end if
      ! get fluid%indx1 and fluid%indx2
    end do
    ! allocate fluid%pidx_h, fluid%indx1 and fluid%indx2
    allocate(self%pidx_h(self%npipe_h))
    allocate(self%indx1(self%npipe_h))
    allocate(self%indx2(self%npipe_h))
    j = 1
    do i = 1, self%npipe
      if(self%pipes(i)%cate/='break'.and.self%pipes(i)%cate/='fill') then
        self%pidx_h(j) = i
        if(j==1) then
          self%indx1(j) = 1
        else
          self%indx1(j) = self%indx2(j-1) + 1
        end if
        self%indx2(j) = self%indx1(j) + self%pipes(i)%nz - 1
        j = j + 1
      end if
    end do
  end subroutine init_fld_p

  subroutine init_fld_j(self, junc_array)
    !======= Declarations =========
    implicit none
    !== I/O variables ==
    class(Fluid), intent(inout) :: self
    class(Junc), intent(in)     :: junc_array(:)
    !== local variables ==
    integer(c_int) :: i
    !======= Internals ============
    self%njunc = size(junc_array)
    ! allocate and assign the arrays storing junction objects
    allocate(self%juncs(self%njunc))
    ! allocate the array storing mdots of all junctions
    allocate(self%mdots(self%njunc))
    self%mdots = 0.d0
    self%juncs = junc_array
    ! calculate the total number of junctions solved by ode
    do i = 1, self%njunc
      if(self%juncs(i)%cate=='jun-head') neq_mdot = neq_mdot+1
    end do
  end subroutine init_fld_j

  subroutine init_fluid(self)
    !======= Declarations =========
    implicit none
    !== I/O variables ==
    class(Fluid), intent(inout) :: self
    !== local variables ==
    character(len=20) :: mid_f, mid_t
    real(c_double), allocatable :: A(:,:), B(:,:)
    real(c_double) :: l1, a1, l2, a2
    real(c_double) :: l_over_a
    integer(c_int) :: i, j, k, l, m, n, nzmax
    !======= Internals ============
    ! loop through the pipes
    nzmax = 0
    do i = 1, self%npipe
      ! the 'fill' pipe can only serve as 'from' of junctions
      if(self%pipes(i)%cate=='fill') then
        do j = 1, self%njunc
          if(self%juncs(j)%pid_t==self%pipes(i)%idstr) then
            print*, 'error in fill pipe ', self%pipes(i)%idstr
            print*, 'the fill pipe cannot serve as "to" of juncs'
            stop 'program stopped at Fluid%init_fluid()'
          end if
          if(self%juncs(j)%pid_f==self%pipes(i)%idstr) then
            if(self%juncs(j)%cate/='jun-flow') then
              print*, 'error in fill pipe ', self%pipes(i)%idstr
              print*, 'the juncs linking fill must be jun-flow type'
              stop 'program stopped at Fluid%init_fluid()'
            end if
          end if
        end do
      ! the 'break' pipe can only serve as 'to' of junctions
      else if(self%pipes(i)%cate=='break') then
        do j = 1, self%njunc
          if(self%juncs(j)%pid_f==self%pipes(i)%idstr) then
            print*, 'error in break pipe ', self%pipes(i)%idstr
            print*, 'the fill pipe cannot serve as "from" of juncs'
            stop 'program stopped at Fluid%init_fluid()'
          end if
        end do
      end if
      ! get the maximum pipe%nz
      if(self%pipes(i)%nz>nzmax) nzmax = self%pipes(i)%nz
    end do
    ! loop through the junctions
    do i = 1, self%njunc
      ! the "from" and "to" cannot be the same pipe
      if(self%juncs(i)%pid_f==self%juncs(i)%pid_t) then
        print*, 'error in junction ', i
        print*, '"from" and "to" cannot be the same pipe'
        stop 'program stopped at Fluid%init_fluid()'
      end if
      ! the fluid material of "from" and "to" cannot be different
      do j = 1, self%npipe
        if(self%pipes(j)%idstr==self%juncs(i)%pid_f) then
          mid_f = self%pipes(j)%matid
        else if(self%pipes(j)%idstr==self%juncs(i)%pid_t) then
          mid_t = self%pipes(j)%matid
        end if
      end do
      if(mid_f/=mid_t) then
        print*, 'error in junction ', i
        print*, 'the fluid material of "from" and "to" cannot be different'
        stop 'program stopped at Fluid%init_fluid()'
      end if
    end do
    ! calculate inv_A matrix, A.*[mdot] = [vec1]
    allocate(A(self%njunc,self%njunc))
    A = 0.d0
    k = 1 ! matrix line counter (equation number)
    ! juncs whose mdot obtained from ode/signal 
    do i = 1, self%njunc
      if(self%juncs(i)%cate/='jun-norm') then
        A(k,i) = 1.d0
        k = k+1
      end if
    end do
    ! normal dependent junctions
    do i = 1, self%npipe
      if(self%pipes(i)%cate/='freelevel'.and. &
         self%pipes(i)%cate/='break'.and.self%pipes(i)%cate/='fill') then
        do j = 1, self%njunc
          if(self%juncs(j)%pid_t==self%pipes(i)%idstr) then
            A(k,j) = 1.d0
          else if (self%juncs(j)%pid_f==self%pipes(i)%idstr) then
            A(k,j) =-1.d0
          end if
        end do
        k = k+1
      end if
    end do
    allocate(self%inv_A(self%njunc,self%njunc))
    ! calculate inv_A
    call brinv(A,self%inv_A)
    deallocate(A)
    ! ---
    ! construct the findindex array
    j = 1 ! the pipe index
    k = 1 ! the pipe node index
    allocate(findindex(self%npipe,nzmax))
    findindex = 0
    do i = 1, self%pntot
      findindex(j,k) = i
      k = k+1
      if(k>self%pipes(j)%nz) then
        k = 1
        j = j+1
      end if
    end do
    ! ---
    ! constructe the B matrix, B.*x = b
    allocate(B(self%njunc+self%njint+self%pntot, self%njunc+self%njint+self%pntot))
    allocate(self%inv_B(self%njunc+self%njint+self%pntot, self%njunc+self%njint+self%pntot))
    B = 0.d0
    ! equations for getting external junctions dmdot/dt
    l = 1 ! matrix line counter (equation number)
    do i = 1, self%njunc
      if(self%juncs(i)%cate=='jun-flow') then
        B(l,i) = 1.d0
        l = l+1
      end if
    end do
    do i = 1, self%npipe
      if(self%pipes(i)%cate/='freelevel'.and. &
         self%pipes(i)%cate/='break'.and.self%pipes(i)%cate/='fill') then
        do j = 1, self%njunc
          if(self%juncs(j)%pid_t==self%pipes(i)%idstr) then
            B(l,j) = 1.d0
          else if (self%juncs(j)%pid_f==self%pipes(i)%idstr) then
            B(l,j) =-1.d0
          end if
        end do
        l = l+1
      end if
    end do
    ! ---
    ! equations for getting internal junctions dmdot/dt (incompressible assumption)
    k = 1
    do i = 1, self%npipe
      if(self%pipes(i)%cate/='freelevel'.and. &
         self%pipes(i)%cate/='break'.and.self%pipes(i)%cate/='fill') then ! pipe containing internal junction
        do j = 1, self%njunc
          if(self%juncs(j)%pid_t==self%pipes(i)%idstr) then
            B(l,j) =-1.d0
          end if
        end do
        B(l, self%njunc+k) = 1.d0
        k = k+1
        l = l+1
      end if
    end do
    pindex = l
    ! construct pressure equations
    do i = 1, self%njunc ! loop through the external junctions
      do j = 1, self%npipe
        if(self%pipes(j)%idstr==self%juncs(i)%pid_f) exit
      end do
      k = self%pipes(j)%nz
      if((self%juncs(i)%cate=='jun-flow').and.(self%pipes(j)%cate/='fill')) then
        continue
      else
        ! here (j, k) denotes the 'from' node index of i-th junction
        do m = 1, self%npipe
          if(self%pipes(m)%idstr==self%juncs(i)%pid_t) exit
        end do
        ! here (m, 1) denotes the 'to' node index of i-th junction
        B(l, self%njunc+self%njint+findindex(j,k)) = 1.d0
        B(l, self%njunc+self%njint+findindex(m,1)) =-1.d0
        l_over_a = 0.d0
        if(self%pipes(j)%cate/='freelevel'.and. &
           self%pipes(j)%cate/='break'.and.self%pipes(j)%cate/='fill') then
          l1 = self%pipes(j)%dz
          a1 = self%pipes(j)%areaz
          l_over_a = l_over_a + l1/a1/2.d0
        end if
        if(self%pipes(m)%cate/='freelevel'.and. &
           self%pipes(m)%cate/='break'.and.self%pipes(m)%cate/='fill') then
          l2 = self%pipes(m)%dz
          a2 = self%pipes(m)%areaz
          l_over_a = l_over_a + l2/a2/2.d0
        end if
        B(l, i) = -1.d0*l_over_a
        l = l+1
      end if
    end do
    do j = 1, self%npipe ! loop for pbnds and internal junctions
      if(self%pipes(j)%cate=='freelevel'.or.self%pipes(j)%cate=='break') then ! pressure boundary
        B(l,self%njunc+self%njint+findindex(j,1)) = 1.d0
        l = l+1
      else if(self%pipes(j)%nz>1) then
        do k = 1, self%pipes(j)%nz-1
          B(l,self%njunc+self%njint+findindex(j,k))  = 1.d0
          B(l,self%njunc+self%njint+findindex(j,k+1))=-1.d0
          ! loop to find the internal junction index of the j-th pipe
          n = j
          do m = 1, j 
            if(self%pipes(m)%cate=='freelevel'.or. &
            self%pipes(m)%cate=='break'.or.self%pipes(m)%cate=='fill') then
              n = n-1
            end if
          end do 
          ! now n is the internal junction index
          l1 = self%pipes(j)%dz
          a1 = self%pipes(j)%areaz
          l_over_a = l1/a1
          B(l,self%njunc+n) = -1.d0*l_over_a
          l = l+1
        end do
      end if
    end do
    ! calculate inv_B
    call brinv(B,self%inv_B)
    ! fluid initialization
    do i = 1, self%npipe
      ! initialize the pipe fluid enthalpies
      select case (self%pipes(i)%matid)
        case('LBE')
          call ht_lbe(self%pipes(i)%temp, self%pipes(i)%enth)
        case('Na1')
          call ht_na1(self%pipes(i)%temp, self%pipes(i)%enth)
        case('H2O')
          call hpt_h2o(self%pipes(i)%pval, self%pipes(i)%temp, self%pipes(i)%enth)
        case default
          print*, 'no match for fluid material name = ', self%pipes(i)%matid
          stop 'program stopped at Fluid%init_fluid'
      end select
    ! initialize the pipe fluid properties
      call self%pipes(i)%calc_prop()
    end do
  end subroutine init_fluid

  subroutine comp_rhs_fld1(self, rhs)
    !======= Declarations =========
    implicit none
    !== I/O variables ==
    class(Fluid), intent(inout) :: self
    ! rhs array
    real(c_double), intent(out) :: rhs(:)
    !== local variables ==
    real(c_double), allocatable :: b(:), x(:) ! x = [dmdot/dt_junc, dmdot/dt_jint, p]
    real(c_double), allocatable :: dpf(:), dpg(:), dpk(:)
    real(c_double), allocatable :: fric(:)
    real(c_double) :: rr, gg, aa
    integer(c_int) :: i, j, k, l, m
    !======= Internals ============
    ! update the fluid properties
    do i = 1, self%npipe
      call self%pipes(i)%calc_prop()
    end do 
    ! calculate the junction mass flowrates
    allocate(b(self%njunc))
    b = 0.d0
    k = 1 ! matrix line number (equation number)
    do i =1, self%njunc
      if(self%juncs(i)%cate/='jun-norm') then
        select case (self%juncs(i)%cate)
        case('jun-head','jun-flow')
          b(k) = self%mdots(i)
        case default
          continue
        end select
        k = k+1
      end if
    end do
    self%mdots = MATMUL(self%inv_A, b)
    deallocate(b)
    
    ! calculate mdot of internal junctions (pipe internal)
    do j = 1, self%npipe
      self%pipes(j)%mdot = 0.d0
      do i = 1, self%njunc
        if(self%juncs(i)%pid_t==self%pipes(j)%idstr) then
          self%pipes(j)%mdot = self%pipes(j)%mdot + self%mdots(i)
        end if
      end do
      if(self%pipes(j)%cate/='break'.and.self%pipes(j)%cate/='fill') then
        !                     |              mdot|*              dhyd/              miu/              areaz
        self%pipes(j)%re = abs(self%pipes(j)%mdot)*self%pipes(j)%dhyd/self%pipes(j)%miu/self%pipes(j)%areaz
        self%pipes(j)%pe = self%pipes(j)%re*self%pipes(j)%pr 
      else
        self%pipes(j)%re = 0.d0
        self%pipes(j)%pe = 0.d0
      end if
    end do
    ! vector storing dm/dt and node pressure 
    allocate(x(self%njunc+self%njint+self%pntot))
    ! vector storing rhs of pressure equations
    allocate(b(self%njunc+self%njint+self%pntot))
    b = 0.d0
    ! calculate the frictional pressure drop
    allocate(dpf(self%pntot))
    i = 1
    do j = 1, self%npipe
      allocate(fric(self%pipes(j)%nz))
      select case (self%pipes(j)%cate)
        case('round')
          call fric_cal(self%pipes(j)%cate,self%pipes(j)%re, fric)
          dpf(i:i+self%pipes(j)%nz-1) = fric*(self%pipes(j)%dz/self%pipes(j)%dhyd)* &
          0.5d0*self%pipes(j)%mdot*abs(self%pipes(j)%mdot)/self%pipes(j)%rho/(self%pipes(j)%areaz)**(2.d0)
          i = i + self%pipes(j)%nz
        case('brods')
          stop 'TBD'
        case('wrods')
          stop 'TBD'
        case('freelevel','fill','break')
          dpf(i) = 0.d0
          i = i+1
        case default
          print*, 'no match for cate = ', self%pipes(j)%cate
          stop 'program stopped at fluid%comp_rhs_fld'
      end select
      deallocate(fric)
    end do
    ! calculate the gravational pressure drop array of all fluid nodes
    allocate(dpg(self%pntot))
    i = 1
    do j = 1, self%npipe
      do k = 1, self%pipes(j)%nz
        if(self%pipes(j)%cate/='freelevel'.and.self%pipes(j)%cate/='break'.and.& 
           self%pipes(j)%cate/='fill') then
          !                         rho*   g*              dz*              dir
          dpg(i) = self%pipes(j)%rho(k)*grav*self%pipes(j)%dz*self%pipes(j)%dir
        else
          ! freelevel/break/fill pipe
          dpg(i) = 0.d0
        end if
        i = i + 1
      end do
    end do
    ! calculate the local pressure drop array of external junctions
    allocate(dpk(self%njunc))
    do i = 1, self%njunc
      ! here (j, k) denotes the upwind node index of i-th junction
      if(self%mdots(i)>0.d0) then
        do j = 1, self%npipe
          if(self%pipes(j)%idstr==self%juncs(i)%pid_f) exit
        end do
        k = self%pipes(j)%nz
      else
        do j = 1, self%npipe
          if(self%pipes(j)%idstr==self%juncs(i)%pid_t) exit
        end do
        k = 1
      end if
      rr = self%pipes(j)%rho(k)
      if(self%pipes(j)%cate/='fill'.and.self%pipes(j)%cate/='break') then
        aa = self%pipes(j)%areaz
      else
        if(self%mdots(i)>0.d0) then
          do m = 1, self%npipe
            if(self%pipes(m)%idstr==self%juncs(i)%pid_t) exit
          end do
        else
          do m = 1, self%npipe
            if(self%pipes(m)%idstr==self%juncs(i)%pid_f) exit
          end do
        end if
        aa = self%pipes(m)%areaz
      end if
      gg = self%mdots(i)/aa
      dpk(i) = self%juncs(i)%kfac*gg*abs(gg)/rr
    end do
    ! construct the b vector
    l = pindex
    do i = 1, self%njunc! loop through the external junctions
      do j = 1, self%npipe
        if(self%pipes(j)%idstr==self%juncs(i)%pid_f) exit
      end do
      k = self%pipes(j)%nz
      ! here (j, k) denotes the 'from' node index of i-th junction
      if((self%juncs(i)%cate=='jun-flow').and.(self%pipes(j)%cate/='fill')) then
        continue
      else
        do m = 1, self%npipe
          if(self%pipes(m)%idstr==self%juncs(i)%pid_t) exit
        end do
        ! here (m, 1) denotes the 'to' node index of i-th junction
        b(l) = b(l) + (dpg(findindex(j,k))+ dpg(findindex(m,1)))/2.d0 ! gravational pressure drop
        b(l) = b(l) + (dpf(findindex(j,k))+ dpf(findindex(m,1)))/2.d0 ! frictionla pressure drop
        b(l) = b(l) + dpk(i)
        if(self%juncs(i)%cate=='jun-head') then
          b(l) = b(l) - self%juncs(i)%phead
        end if
        l = l+1
      end if
    end do
    do j = 1, self%npipe ! loop through the internal junctions
      if(self%pipes(j)%cate=='freelevel'.or.self%pipes(j)%cate=='break') then ! pressure boundary
        b(l) = self%pipes(j)%pval(1)
        l = l+1
      else if(self%pipes(j)%nz>1) then
        do k = 1, self%pipes(j)%nz-1
          b(l) = b(l) + (dpg(findindex(j,k))+ dpg(findindex(j,k+1)))/2.d0
          b(l) = b(l) + (dpf(findindex(j,k))+ dpf(findindex(j,k+1)))/2.d0
          l = l+1
        end do
      end if
    end do
    x = MATMUL(self%inv_B,b)
    ! read pressure from p vector
    i = self%njunc + self%njint + 1
    do j = 1, self%npipe
      do k = 1, self%pipes(j)%nz
        self%pipes(j)%pval(k) = x(i)
        i = i+1
      end do
    end do
    j = 1
    ! time derivatives of mdot
    if(neq_mdot>0) then
      do i = 1, self%njunc
        if(self%juncs(i)%cate=='jun-head') then
          rhs(j) = x(j)
          j = j+1
        end if
      end do
    end if
    ! time derivatives of freelevel volume length
    if(neq_flen>0) then
      do i = 1, self%npipe
        if(self%pipes(i)%cate=='freelevel') then
          rhs(j) = 0.d0
          do k = 1, self%njunc
            if(self%juncs(k)%pid_t==self%pipes(i)%idstr) then
              rhs(j) = rhs(j) + self%mdots(k)
            else if(self%juncs(k)%pid_f==self%pipes(i)%idstr) then
              rhs(j) = rhs(j) - self%mdots(k)
            end if
          end do
          rhs(j) = rhs(j)/self%pipes(i)%rho(1)/self%pipes(i)%areaz
          j = j + 1
        end if
      end do
    end if
    ! deallocate the temporary arrays
    deallocate(b)
    deallocate(x)
    deallocate(dpf)
    deallocate(dpg)
    ! time derivatives of fluid enthalpies
  end subroutine comp_rhs_fld1

  subroutine comp_rhs_fld2(self, rhs)
    !======= Declarations =========
    implicit none
    !== I/O variables ==
    class(Fluid), intent(inout) :: self
    ! rhs array
    real(c_double), intent(out) :: rhs(:)
    !== local variables ==
    real(c_double) :: mh
    integer(c_int) :: i, j, k, n
    !======= Internals ============
    ! calculate inter-pipe enthalpy convection
    do i = 1, self%npipe ! loop to initialize
      self%pipes(i)%mh_i = 0.d0
      self%pipes(i)%mh_o = 0.d0
    end do
    do i = 1, self%njunc ! loop to accumulate mh_i and mh_o
      do j = 1, self%npipe
        if(self%pipes(j)%idstr==self%juncs(i)%pid_f) then
          exit
        end if
      end do
      do k = 1, self%npipe
        if(self%pipes(k)%idstr==self%juncs(i)%pid_t) then
          exit
        end if
      end do
      ! j-th, k-th pipe are the 'from', 'to' of i-th junction respectively
      if(self%mdots(i)>0.d0) then
        n  = self%pipes(j)%nz
        mh = self%mdots(i)*self%pipes(j)%enth(n)
      else ! self%mdots(i)<=0.d0
        mh = self%mdots(i)*self%pipes(k)%enth(1)
      end if
      self%pipes(j)%mh_o = self%pipes(j)%mh_o + mh
      self%pipes(k)%mh_i = self%pipes(k)%mh_i + mh
    end do
    ! ---
    ! calculate the rhs array
    ! the pointer manipulation maybe inappropriate for omp parallel use
    do i = 1, self%npipe_h
      j = self%pidx_h(i)
      call self%pipes(j)%calc_rhs_pipe(rhs(self%indx1(i):self%indx2(i)))
    end do
  end subroutine comp_rhs_fld2
end module Fluid_mod