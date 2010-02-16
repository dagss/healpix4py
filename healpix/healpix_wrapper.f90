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

  subroutine cywrap_alms2fits(filename, nalms, alms, ncl, next, shapes)
    use healpix_types
    use fitstools
    integer(i4b), intent(in), dimension(1:1) :: shapes
    character(len=shapes(1)),  intent(in)    :: filename
    integer(I4B),      intent(in)            :: nalms, next, ncl
    real(DP),          intent(in), dimension(1:nalms,1:ncl+1,1:next), target :: alms

    character(len=80), dimension(1:1,1:next) :: header

    header(:,:) = ''

    call alms2fits(filename, nalms, alms, ncl, header, 1, next)
  end subroutine cywrap_alms2fits

  subroutine cywrap_fits2alms(filename, nalms, alms, ncl, next, shapes)
    use healpix_types
    use fitstools
    integer(i4b), intent(in), dimension(1:1) :: shapes
    character(len=shapes(1)),  intent(in)    :: filename
    integer(I4B),      intent(in)            :: nalms, next, ncl
    real(DP),          intent(inout), dimension(1:nalms,1:ncl+1,1:next), target :: alms

    character(len=80), dimension(1:10,1:next) :: header

    call fits2alms(filename, nalms, alms, ncl, header, 10, next)
  end subroutine cywrap_fits2alms

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
