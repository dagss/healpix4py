
  ! Eventually, this should use fwrap or at least ISO_C_BINDING,
  ! but until I manage to get gfortran 4.3 built locally, this is
  ! what we're stuck with

  subroutine cywrap_nside2npix(nside, out)
    use healpix_types
    use pix_tools
    implicit none

    integer :: nside, out

    out = nside2npix(nside)
  end subroutine cywrap_nside2npix

  subroutine cywrap_npix2nside(npix, out)
    use healpix_types
    use pix_tools
    implicit none

    integer :: npix, out

    out = npix2nside(npix)
  end subroutine cywrap_npix2nside

  subroutine cywrap_alm2map_sc_d(nsmax, nlmax, nmmax, alm, map)
    use healpix_types
    use alm_tools
    implicit none

    integer(I4B), intent(IN)                   :: nsmax, nlmax, nmmax
    complex(DPC), intent(IN),  dimension(1:1,0:nlmax,0:nmmax) :: alm
    real(DP),   intent(OUT), dimension(0:(12_i8b*nsmax)*nsmax-1) :: map

    call alm2map(nsmax, nlmax, nmmax, alm, map)
  end subroutine cywrap_alm2map_sc_d

  subroutine cywrap_map2alm_sc_d(nsmax, nlmax, nmmax, map, alm, zbounds, w8ring)
    use healpix_types
    use alm_tools
    implicit none

    integer(I4B), intent(IN)                    :: nsmax, nlmax, nmmax
    real(DP),   intent(IN),  dimension(0:(12_i8b*nsmax)*nsmax-1) :: map
    complex(DPC), intent(OUT), dimension(1:1,0:nlmax,0:nmmax) :: alm
    real(DP),     intent(IN),  dimension(1:2) :: zbounds
    real(DP),     intent(IN),  dimension(1:2*nsmax,1) :: w8ring

    call map2alm(nsmax, nlmax, nmmax, map, alm, zbounds, w8ring)
  end subroutine cywrap_map2alm_sc_d

  subroutine cywrap_output_map_d(map, outfile, shapes)
    use healpix_types
    use fitstools, only: output_map

    integer, dimension(1:3) :: shapes
    real(dp), intent(inout), dimension(101:shapes(1)+100, 1:shapes(2)) :: map
    character(len=shapes(3)), intent(in) :: outfile

    character(len=80), dimension(1:100)           :: header

    header(:) = ''
    call output_map(map, header, outfile)
  end subroutine cywrap_output_map_d

  subroutine cywrap_convert_ring2nest_d(nside, map, nd)
    use healpix_types
    use pix_tools, only: convert_ring2nest
    integer(i4b), intent(in) :: nside, nd
    real(dp), intent(inout), dimension(0:12*nside**2-1, 1:nd) :: map
    call convert_ring2nest(nside, map)
  end subroutine cywrap_convert_ring2nest_d

  subroutine cywrap_convert_nest2ring_d(nside, map, nd)
    use healpix_types
    use pix_tools, only: convert_nest2ring
    integer(i4b), intent(in) :: nside, nd
    real(dp), intent(inout), dimension(0:12*nside**2-1, 1:nd) :: map
    call convert_nest2ring(nside, map)
  end subroutine cywrap_convert_nest2ring_d

  subroutine cywrap_pix2vec_nest(nside, ipix, vector)
    use healpix_types
    use fitstools
    use pix_tools
    INTEGER(KIND=I4B), INTENT(IN) :: nside, ipix
    REAL(KIND=DP), INTENT(OUT), dimension(3) :: vector
    call pix2vec_nest(nside, ipix, vector)
  end subroutine cywrap_pix2vec_nest

  subroutine cywrap_pix2vec_ring(nside, ipix, vector)
    use healpix_types
    use fitstools
    use pix_tools
    INTEGER(KIND=I4B), INTENT(IN) :: nside, ipix
    REAL(KIND=DP), INTENT(OUT), dimension(3) :: vector
    call pix2vec_ring(nside, ipix, vector)
  end subroutine cywrap_pix2vec_ring

  subroutine cywrap_vec2ang(vector, theta, phi)
    !need to specify end-point of vector..., not assumed-shape
    use healpix_types
    use pix_tools
    REAL(KIND=DP), INTENT(IN), dimension(1:3) :: vector
    REAL(KIND=DP), INTENT(OUT) :: theta, phi
    call vec2ang(vector, theta, phi)
  end subroutine cywrap_vec2ang

  subroutine cywrap_vec2pix_ring(nside, vector, ipix)
    use healpix_types
    use pix_tools
    INTEGER(KIND=I4B), INTENT(IN) :: nside
    INTEGER(KIND=I4B), INTENT(OUT) :: ipix
    REAL(KIND=DP), INTENT(IN), dimension(1:3) :: vector
    call vec2pix_ring(nside, vector, ipix)
  end subroutine cywrap_vec2pix_ring

  subroutine cywrap_sub_udgrade_nest_d(map_in, nside_in, map_out, nside_out)
    use healpix_types
    use udgrade_nr
    INTEGER(I4B), INTENT(IN) :: nside_in, nside_out
    REAL(DP),     INTENT(IN),  DIMENSION(0:12*nside_in**2-1), target :: map_in
    REAL(DP),     INTENT(OUT), DIMENSION(0:12*nside_out**2-1) :: map_out
!    REAL(DP),     INTENT(IN), OPTIONAL :: fmissval
!    LOGICAL(LGT), INTENT(IN), OPTIONAL :: pessimistic
    call sub_udgrade_nest(map_in, nside_in, map_out, nside_out) !, fmissval, pessimistic)
  end subroutine cywrap_sub_udgrade_nest_d

  subroutine cywrap_read_unformatted_2d_complex_d(data, n, m, filename, filename_len)
    use healpix_types
    complex(dpc), intent(out), dimension(1:n,1:m)   :: data
    integer(i4b), intent(in) :: filename_len
    character(len=filename_len),  intent(in)    :: filename
    open(10, FILE=trim(filename), FORM='UNFORMATTED')
    read(10) data
    close(10)
  end subroutine cywrap_read_unformatted_2d_complex_d

  subroutine cywrap_read_unformatted_1d_real_d(data, n, filename, filename_len)
    use healpix_types
    real(dp), intent(out), dimension(1:n)   :: data
    integer(i4b), intent(in) :: filename_len
    character(len=filename_len),  intent(in)    :: filename
    open(10, FILE=trim(filename), FORM='UNFORMATTED')
    read(10) data
    close(10)
  end subroutine cywrap_read_unformatted_1d_real_d


  subroutine cywrap_rotate_alm_d(lmax, alm, psi, theta, phi, alm_s0, alm_s1, alm_s2)
    use healpix_types
    use alm_tools
    implicit none
    integer(I4B),   intent(in) :: lmax, alm_s0, alm_s1, alm_s2
    complex(DPC), intent(inout), dimension(1:alm_s0,0:alm_s1-1,0:alm_s2-1) :: alm
    real(DP),       intent(in) :: psi, theta, phi
    call rotate_alm(lmax, alm, psi, theta, phi)
  end subroutine cywrap_rotate_alm_d

  subroutine cywrap_remove_dipole_double(nside, map, ordering, degree, multipoles, mask)
    use healpix_types
    use pix_tools

    integer(kind=i4b),                  intent(in)    :: nside
    integer(kind=i4b),                  intent(in)    :: ordering, degree
    real   (kind=DP),   dimension(0:degree*degree-1),  intent(out)   :: multipoles
    real   (kind=DP),   dimension(0:12*nside*nside-1),  intent(inout) :: map
    real   (kind=DP),   dimension(0:12*nside*nside-1),  intent(in) :: mask
    call remove_dipole(nside, map, ordering, degree, multipoles, (/ -1.0_DP, 1.0_DP /), &
 & mask=mask)
  end subroutine
