module rotate_alm

  use healpix_types
  use misc_utils, only: assert_alloc, assert_present, assert, fatal_error, string, strupcase !, wall_clock_time

!  use omp

  implicit none

  public :: rotate_alm_single

  interface rotate_alm_single
     module procedure rotate_alm_single_d, rotate_alm_single_s
  end interface

contains

  include 'rotate_alm_dd_inc.f90'
  include 'rotate_alm_ss_inc.f90'

end module rotate_alm
