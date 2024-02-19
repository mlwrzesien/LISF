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
! !MODULE: AMSRcnnSnow_obsMod
! \label(AMSRcnnSnow_obsMod)
!
! !INTERFACE:
module AMSRcnnSnow_obsMod
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

  PUBLIC :: AMSRcnnSnow_obsinit
  PUBLIC :: amsrcnnsnowobs

  type, public :: amsrcnnsnowobsdec
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

  end type amsrcnnsnowobsdec

  type(amsrcnnsnowobsdec), allocatable :: amsrcnnsnowobs(:)

contains

!BOP
!
! !ROUTINE: AMSRcnnSnow_obsinit
! \label{AMSRcnnSnow_obsinit}
!
! !INTERFACE:
  subroutine AMSRcnnSnow_obsinit(i)
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

    if(.not.allocated(amsrcnnsnowobs)) then
       allocate(amsrcnnsnowobs(LVT_rc%nDataStreams))
    endif

    call ESMF_ConfigGetAttribute(LVT_config, amsrcnnsnowobs(i)%odir, &
         label='AMSR CNN snow observation directory:',rc=status)
    call LVT_verify(status, 'AMSR CNN snow observation directory: not defined')

    gridDesci = 0
    call LVT_update_timestep(LVT_rc, 86400)

    allocate(amsrcnnsnowobs(i)%rlat(LVT_rc%lnc*LVT_rc%lnr))
    allocate(amsrcnnsnowobs(i)%rlon(LVT_rc%lnc*LVT_rc%lnr))

    allocate(amsrcnnsnowobs(i)%w11(LVT_rc%lnc*LVT_rc%lnr))
    allocate(amsrcnnsnowobs(i)%w12(LVT_rc%lnc*LVT_rc%lnr))
    allocate(amsrcnnsnowobs(i)%w21(LVT_rc%lnc*LVT_rc%lnr))
    allocate(amsrcnnsnowobs(i)%w22(LVT_rc%lnc*LVT_rc%lnr))

    allocate(amsrcnnsnowobs(i)%n11(LVT_rc%lnc*LVT_rc%lnr))
    allocate(amsrcnnsnowobs(i)%n12(LVT_rc%lnc*LVT_rc%lnr))
    allocate(amsrcnnsnowobs(i)%n21(LVT_rc%lnc*LVT_rc%lnr))
    allocate(amsrcnnsnowobs(i)%n22(LVT_rc%lnc*LVT_rc%lnr))

    amsrcnnsnowobs(i)%nc = 3599
    amsrcnnsnowobs(i)%nr = 700

    gridDesci(1) = 0
    gridDesci(2) = 3599
    gridDesci(3) = 700
    gridDesci(4) = 0.025
    gridDesci(5) = -179.975
    gridDesci(6) = 128
    gridDesci(7) = 69.925
    gridDesci(8) = 179.825
    gridDesci(9) = 0.1
    gridDesci(10) = 0.1
    gridDesci(20) = 64

    call bilinear_interp_input(gridDesci,LVT_rc%gridDesc,&
         LVT_rc%lnc*LVT_rc%lnr, &
         amsrcnnsnowobs(i)%rlat,  amsrcnnsnowobs(i)%rlon, &
         amsrcnnsnowobs(i)%n11, amsrcnnsnowobs(i)%n12,   &
         amsrcnnsnowobs(i)%n21, amsrcnnsnowobs(i)%n22,   &
         amsrcnnsnowobs(i)%w11, amsrcnnsnowobs(i)%w12,   &
         amsrcnnsnowobs(i)%w21, amsrcnnsnowobs(i)%w22)


  end subroutine AMSRcnnSnow_obsinit


end module AMSRcnnSnow_obsMod
