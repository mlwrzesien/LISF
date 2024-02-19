!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !MODULE: AMSREcnnSnow_obsMod
! \label(AMSREcnnSnow_obsMod)
!
! !INTERFACE:
module AMSREcnnSnow_obsMod
!
! !USES:
  use ESMF

  implicit none
  PRIVATE
!
! !INPUT PARAMETERS:
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
!  This module handles the observation plugin for the
!  University of Arizona (UA) SWE/Snow Depth data v01.
!  The UA SNOW data is provided in the NAD 1983 grid
!  with ~4-km resolution. The domain extents are from
!  approximately (24N, -125W) to (50N, -66.5W).
!  The data entries are 16-bit signed integers.
!
!  Temporal coverage is from 1 Oct 1981 - 30 Sep 2017.
!  The data is organized by Water Years.
!
! !FILES USED:
!
! !REVISION HISTORY:
!  28 May 2019: Rhae Sung Kim, Initial Specification!
!  19 Jun 2019: David Mocko, Set valid time of data to 12Z
!
!EOP

  PUBLIC :: AMSREcnnSnow_obsinit
  PUBLIC :: amsrecnnsnowobs

  type, public :: amsrecnnsnowobsdec
     character*100        :: odir
     integer              :: nc, nr

     real,    allocatable     :: rlat(:)
     real,    allocatable     :: rlon(:)
     integer, allocatable     :: n11(:)
     integer, allocatable     :: n12(:)
     integer, allocatable     :: n21(:)
     integer, allocatable     :: n22(:)
     real,    allocatable     :: w11(:)
     real,    allocatable     :: w12(:)
     real,    allocatable     :: w21(:)
     real,    allocatable     :: w22(:)

  end type amsrecnnsnowobsdec

  type(amsrecnnsnowobsdec), allocatable :: amsrecnnsnowobs(:)

contains

!BOP
!
! !ROUTINE: AMSREcnnSnow_obsinit
! \label{AMSREcnnSnow_obsinit}
!
! !INTERFACE:
  subroutine AMSREcnnSnow_obsinit(i)
!
! !USES:
    use LVT_coreMod
    use LVT_histDataMod
    use LVT_timeMgrMod
    use LVT_logMod

    implicit none
!
! !INPUT PARAMETERS:
    integer,   intent(IN) :: i
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
! This subroutine initializes and sets up the data structures
! required for reading UA SNOW data.
!
! !FILES USED:
!
!EOP

    real               :: gridDesci(50)
    integer            :: status

    if(.not.allocated(amsrecnnsnowobs)) then
       allocate(amsrecnnsnowobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_config, amsrecnnsnowobs(i)%odir, &
         label='AMSRE CNN snow observation directory:',rc=status)
    call LVT_verify(status, 'AMSRE CNN snow observation directory: not defined')

    gridDesci = 0
    call LVT_update_timestep(LVT_rc, 3600)

    allocate(amsrecnnsnowobs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
    allocate(amsrecnnsnowobs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))

    allocate(amsrecnnsnowobs(i)%w11(LVT_rc%lnc*LVT_rc%lnr))
    allocate(amsrecnnsnowobs(i)%w12(LVT_rc%lnc*LVT_rc%lnr))
    allocate(amsrecnnsnowobs(i)%w21(LVT_rc%lnc*LVT_rc%lnr))
    allocate(amsrecnnsnowobs(i)%w22(LVT_rc%lnc*LVT_rc%lnr))

    allocate(amsrecnnsnowobs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
    allocate(amsrecnnsnowobs(i)%n12(LVT_rc%lnc*LVT_rc%lnr))
    allocate(amsrecnnsnowobs(i)%n21(LVT_rc%lnc*LVT_rc%lnr))
    allocate(amsrecnnsnowobs(i)%n22(LVT_rc%lnc*LVT_rc%lnr))

    amsrecnnsnowobs(i)%nc = 3599
    amsrecnnsnowobs(i)%nr = 700

    gridDesci(1) = 0
    gridDesci(2) = 3599
    gridDesci(3) = 700
    gridDesci(4) = 0.025
    gridDesci(5) = -179.975
    gridDesci(6) = 128
    gridDesci(7) = 69.925
    gridDesci(8) = -179.825
    gridDesci(9) = 0.1
    gridDesci(10) = 0.1
    gridDesci(20) = 64

    call bilinear_interp_input(gridDesci,LVT_rc%gridDesc,&
         LVT_rc%lnc*LVT_rc%lnr, &
         amsrecnnsnowobs(i)%rlat,  amsrecnnsnowobs(i)%rlon, &
         amsrecnnsnowobs(i)%n11, amsrecnnsnowobs(i)%n12,   &
         amsrecnnsnowobs(i)%n21, amsrecnnsnowobs(i)%n22,   &
         amsrecnnsnowobs(i)%w11, amsrecnnsnowobs(i)%w12,   &
         amsrecnnsnowobs(i)%w21, amsrecnnsnowobs(i)%w22)


  end subroutine AMSREcnnSnow_obsinit


end module AMSREcnnSnow_obsMod
