! rotate_alm testbed
!
! Modifications: Dag Sverre Seljebotn
!
! Original HEALPix copyright:
!
!  Copyright (C) 1997-2008 Krzysztof M. Gorski, Eric Hivon, 
!                          Benjamin D. Wandelt, Anthony J. Banday, 
!                          Matthias Bartelmann, Hans K. Eriksen, 
!                          Frode K. Hansen, Martin Reinecke

  !========================================================
  subroutine rotate_alm_single_KLOAD(lmax, alm, psi, theta, phi)
    !=========================================================
    !Input: Complex array alm(p,l,m) with (l,m) in [0,lmax]^2, and p in [1,nd]
    !Euler rotation angles psi, theta, phi in radians
    !Output: Rotated array alm(p, l,m)
    !
    ! Euler angle convention  is right handed, active rotation
    ! psi is the first rotation about the z-axis (vertical), in [-2pi,2pi]
    ! then theta about the ORIGINAL (unrotated) y-axis, in [-2pi,2pi]
    ! then phi  about the ORIGINAL (unrotated) z-axis (vertical), in [-2pi,2pi]
    !
    ! Equivalently
    ! phi is the first rotation about the z-axis (vertical)
    ! then theta  about the NEW   y-axis (line of nodes)
    ! then psi    about the FINAL z-axis (figure axis)
    ! ---
    ! the recursion on the Wigner d matrix is inspired from the very stable
    ! double sided one described in Risbo (1996, J. of Geodesy, 70, 383)
    ! based on equation (4.4.1) in Edmonds (1957).
    ! the Risbo's scatter scheme has been repladed by a gather scheme for
    ! better computing efficiency
    ! the size of the matrix is divided by 2 using Edmonds Eq.(4.2.5) 
    ! to speed up calculations
    ! the loop on j has been unrolled for further speed-up
    ! EH, March--April 2005
    !=========================================================
    integer(I4B),   intent(in) :: lmax
    complex(KALMC), intent(inout), dimension(1:,0:,0:) :: alm
    real(DP),       intent(in) :: psi, theta, phi
    ! local variables
    complex(DPC), dimension(0:lmax) :: exppsi, expphi
    complex(DPC), dimension(:,:), allocatable :: alm1, alm2
    real(DP),     dimension(:,:), allocatable :: d, dd
    real(DP),     dimension(:),   allocatable :: sqt, rsqt
    real(DP),     dimension(:),   allocatable :: tsign
    integer(I4B) :: status
    integer(I4B) :: mm, ll, na1, na2, nd
    integer(I4B) :: i, j, k, kd, hj
    real(DP)     :: p, q, pj, qj, fj, temp
    character(len=*), parameter :: code = 'ROTATE_ALM'
    
    interface
       subroutine omp_set_num_threads(set)
         integer, intent(in) :: set 
       end subroutine omp_set_num_threads
    end interface
    integer :: omp_get_num_threads, omp_get_max_threads, omp_get_thread_num, omp_get_level

    !==========================================================
    
!    call omp_set_num_threads(2)

    if (abs(psi) > 2.d0*PI .or. abs(phi) > 2.d0*PI .or. abs(theta) > 2.d0*PI) then
       write(*,'(a,3(g12.4))') code,psi,theta,phi
       call fatal_error(code//': angles should be in Radians')
    endif

    nd = size(alm,1)
    na1 = size(alm,2)
    na2 = size(alm,3)
    if (na1 < (lmax+1) .or. na2 < (lmax+1)) then
       call fatal_error(code//': unconsistent alm array size and lmax')
    endif

    allocate(d (-1:2*lmax,   -1:lmax),   stat = status)
    call assert_alloc(status,code,'d')
    allocate(dd(-1:2*lmax, -1:lmax), stat = status)
    call assert_alloc(status,code,'dd')
    call assert_alloc(status,code,'alm1 & alm2')
    allocate(alm1(1:nd,0:lmax), alm2(1:nd,0:lmax), stat = status)
    
    d  = 0.0_dp ! very important for gather scheme
    dd = 0.0_dp

!$OMP parallel default(none) &
!$OMP   shared(d, dd, alm, alm1, alm2, hj) &
!$OMP   firstprivate(nd, lmax, theta, psi, phi) &
!$OMP   private(rsqt, tsign, sqt, exppsi, expphi, p, q, k, j, fj, qj, pj, status)

    allocate(sqt(0:2*lmax), rsqt(0:2*lmax), stat = status)
    call assert_alloc(status,code,'sqt & rsqt')
    allocate(tsign(0:lmax+1), stat = status)
    call assert_alloc(status,code,'tsign')

    do i=0, lmax,2
       tsign(i)   =  1.0_dp
       tsign(i+1) = -1.0_dp
    enddo
    !     initialization of square-root  table
    do i=0,2*lmax
       sqt(i) = SQRT(DBLE(i))
    enddo

    ! initialisation of exponential table
    exppsi(0)=cmplx(1, 0, kind=DPC)
    expphi(0)=cmplx(1, 0, kind=DPC)

    do i=1,lmax
       exppsi(i)= cmplx(cos(psi*i), -sin(psi*i), kind=DPC)
       expphi(i)= cmplx(cos(phi*i), -sin(phi*i), kind=DPC)
    enddo

    ! Note: theta has the correct sign.
    p = sin(theta/2.d0)
    q = cos(theta/2.d0)


    write(*,*) omp_get_level(), omp_get_thread_num(), omp_get_num_threads()
    do ll=0,lmax

       ! ------ build d-matrix of order l ------
       if (ll == 0) then
!$OMP single
          d(0,0) = 1.d0
!$OMP end single
          goto 2000
       endif
       if (ll == 1) then
          !     initialize d-matrix degree 1/2
!$OMP single
          dd(0,0)  =  q
          dd(1,0)  = -p
          dd(0,1)  =  p
          dd(1,1)  =  q
!$OMP end single
          goto 1000
       endif

       !  l - 1 --> l - 1/2

       j = 2*ll - 1
       rsqt(0:j) = sqt(j:0:-1)
       fj = DBLE(j)
       qj = q / fj
       pj = p / fj

!$OMP do schedule(dynamic, 32)
       do k = 0, j/2 ! keep only m' <= -1/2

          dd(0:j,k) = rsqt(0:j) * ( d(0:j,k)      * (sqt(j-k)  * qj)   &
               &                  + d(0:j,k-1)    * (sqt(k)    * pj) ) &
               &    +  sqt(0:j) * ( d(-1:j-1,k-1) * (sqt(k)    * qj)   &
               &                  - d(-1:j-1,k)   * (sqt(j-k)  * pj) )
       enddo ! loop on k
!$OMP end do


!$OMP single
       hj = ll-1
       ! l=half-integer, reconstruct m'= 1/2 by symmetry
       if (mod(ll,2) == 0) then
          do k = 0, j-1, 2
             dd(k,   ll) =   dd(j-k,   hj)
             dd(k+1, ll) = - dd(j-k-1, hj)
          enddo
       else
          do k = 0, j-1, 2
             dd(k,   ll) = - dd(j-k,   hj)
             dd(k+1, ll) =   dd(j-k-1, hj)
          enddo
       endif
!$OMP end single

1000   continue

       !  l - 1/2 --> l
       j = 2*ll
       rsqt(0:j) = sqt(j:0:-1)
       fj = DBLE(j)
       qj = q / fj
       pj = p / fj

!$OMP do schedule(dynamic,32)
       do k = 0, j/2 ! keep only m' <= 0
          d (0:j,k) = rsqt(0:j) * ( dd(0:j,k)      * (sqt(j-k)  * qj)   &
               &                  + dd(0:j,k-1)    * (sqt(k)    * pj) ) &
               &    +  sqt(0:j) * ( dd(-1:j-1,k-1) * (sqt(k)    * qj)   &
               &                  - dd(-1:j-1,k)   * (sqt(j-k)  * pj) )
       enddo ! loop on k
!$OMP end do

2000   continue
       ! ------- apply rotation matrix -------
!$OMP single
       do kd = 1, nd
          alm1(kd,0:ll)  = alm(kd,ll,0:ll) * exppsi(0:ll)
       enddo

       ! m = 0
       do kd = 1, nd
          alm2(kd,0:ll) = alm1(kd,0) * d(ll:2*ll,ll)
       enddo
!$OMP end single

!$OMP do schedule(dynamic,32)
       do mm = 0, ll
          do kd = 1, nd
             alm2(kd, mm) = alm2(kd,mm) + sum(alm1(kd,1:ll) *                d(ll-1:0:-1,ll-mm)) &
                  &                +conjg(sum(alm1(kd,1:ll) * (tsign(1:ll) * d(ll+1:2*ll,ll-mm))))
          enddo
       enddo
!$OMP end do

       ! new alm for ll
!$OMP single
       do kd = 1,nd
          alm(kd,ll,0:ll) = alm2(kd,0:ll)*expphi(0:ll)
       enddo
!$OMP end single

    enddo ! loop on ll

    deallocate(sqt, rsqt)
    deallocate(tsign)
!$OMP end parallel

    deallocate(d)
    deallocate(dd)
    deallocate(alm1, alm2)

  end subroutine rotate_alm_single_KLOAD

