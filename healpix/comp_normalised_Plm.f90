  !  THIS routine returns the value of normalised P_lm(theta) such that
  !  2.*PI*Integral_-1^+1 dx P_lm(x)^2 = 1., where x = cos(theta)
  !  
  !  modified P_lm generating routine from HEALPix, by K.M.G. 4, Sept. 2000
  
  !=======================================================================
  subroutine comp_normalised_Plm(nlmax, m, theta, plm)
    use healpix_types
    !=======================================================================
    !nlmax (integer) = maximum l
    !theta in radians (double precision)
    !lambda (double precision) = modulus of the complex spherical harmonics
    !  contains lambda(l,m) with l,m in [0,nlmax]
    !  contains lambda(l-m) with l-m in [0,nlmax-m]
    !=======================================================================
    IMPLICIT none
    !

    
    INTEGER(I4B), INTENT(IN)  :: m
    INTEGER(I4B), INTENT(IN)  :: nlmax
    REAL(DP),     INTENT(IN)  :: theta
    
    
    REAL(DP) :: lambda(0:nlmax)
    REAL(DP),     DIMENSION(0:nlmax), INTENT(OUT) :: plm
    
    INTEGER(I4B) :: nmmax
    INTEGER(I4B) l, ith, indl, mm               !, m ...  alm related
    
!    REAL(DP) sq4pi_inv
    REAL(DP) cth, sth
    REAL(DP) a_rec, lam_mm, lam_lm, lam_0, lam_1, lam_2, par_lm
    REAL(DP) f2m, fm2, fl2
    
    Character(LEN=7), PARAMETER :: code = 'ALM2MAP'
    INTEGER(I4B) :: status
    
    REAL(DP), PARAMETER :: bignorm = 1.d-20*max_dp
    !=======================================================================
    
    
    !      write(*,*)'   PI   =    ',PI
    
    nmmax = nlmax
    
    LAMBDA = 0.0d0
    
    !     --------------------------------------------
!    sq4pi_inv = 1.D0 / SQRT(4.D0*PI)
    
    cth = COS(theta)
    sth = SIN(theta)
    
    plm=0.d0
    
    !      write(*,*)cth,sth
    
    !        lambda_mm tends to go down when m increases (risk of underflow)
    !        lambda_lm tends to go up   when l increases (risk of overflow)
    
    lam_mm = sq4pi_inv * bignorm ! lambda_0_0 * big number --> to avoid underflow
    
    !      do m = 0, nmmax
    
    fm2 = DFLOAT(m) **2
    
    !           ---------- l = m ----------
    par_lm = 1.d0  ! = (-1)^(l+m)
    if (m .ge. 1) then ! lambda_0_0 for m>0
       do mm = 1, m
          f2m = 2.d0 * mm
          lam_mm = - lam_mm * sth * DSQRT( (f2m+1.d0)/ f2m )
       enddo
    endif
    
    lam_lm = lam_mm / bignorm ! actual lambda_mm
    
    LAMBDA(M) = LAM_LM
    
    
    !           ---------- l > m ----------
    lam_0 = 0.d0
    lam_1 = 1.d0 / bignorm    ! small number --> to avoid overflow
    fl2 = DFLOAT(m+1) **2
    a_rec = DSQRT( (4.d0 * fl2 - 1.d0) / (fl2-fm2) )
    lam_2 = cth * lam_1 * a_rec
    do l = m+1, nmmax
       par_lm = - par_lm  ! = (-1)^(l+m)
       
       lam_lm = lam_2 * lam_mm ! actual lambda_lm (small and big numbers cancel out)
       
       !            lambda(l,m) = lam_lm
       !            lambda(l-m) = lam_lm
       
       LAMBDA(L) = LAM_LM
       
       lam_0 = lam_1 / a_rec
       lam_1 = lam_2
       fl2 = DFLOAT(l+1) **2
       a_rec = DSQRT( (4.d0 * fl2 - 1.d0) / (fl2-fm2) )
       lam_2 = (cth * lam_1 - lam_0) * a_rec
    enddo
    
    !      enddo
    
    
    plm = lambda
    
    return
  end subroutine comp_normalised_Plm
