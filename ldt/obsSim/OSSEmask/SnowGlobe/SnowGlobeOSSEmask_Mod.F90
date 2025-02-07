!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

module SnowGlobeOSSEmask_Mod
!BOP
! 
! !MODULE: SnowGlobeOSSEmask_Mod
! 
! !DESCRIPTION: 
!  This module handles the use of external orbital masks from
!  SnowGlobe instrument to be applied to the observations created
!  for OSSEs
!
! !REVISION HISTORY: 
!  02 Jul 2021    Sujay Kumar  Initial Specification
!
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: SnowGlobeOSSEmask_init
!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  PUBLIC :: SnowGlobeOSSEmaskData
!
!EOP
  
  type, public :: SnowGlobeOSSEmaskDatadec
     integer       :: nvars
     integer       :: nest
     integer       :: nc,nr

     character(len=LDT_CONST_PATH_LEN) :: odir

!--------------------------------------------------------
!  interpolation/upscaling weights
!--------------------------------------------------------
     integer, allocatable   :: n11(:)
     integer, allocatable   :: n12(:)
     integer, allocatable   :: n21(:)
     integer, allocatable   :: n22(:)
     real,  allocatable     :: w11(:)
     real,  allocatable     :: w12(:)
     real,  allocatable     :: w21(:)
     real,  allocatable     :: w22(:)

  end type SnowGlobeOSSEmaskDatadec

  type(SnowGlobeOSSEmaskDatadec)  :: SnowGlobeOSSEmaskData

contains

!BOP
! !ROUTINE: SnowGlobeOSSEmask_init
! \label{SnowGlobeOSSEmask_init}
! 
! !INTERFACE: 
  subroutine SnowGlobeOSSEmask_init()
! !USES: 
    use ESMF
    use LDT_coreMod
    use LDT_DAobsDataMod
    use LDT_logMod

    implicit none
! 
! !DESCRIPTION: 
! This routine initializes the structures required for the handling of 
! SnowGlobe observation masks
! 
!EOP

    integer                 :: n
    integer                 :: rc
    real                    :: datares
    real                    :: run_dd(6)
    real                    :: gridDesci(20)
    real                    :: cornerlat1, cornerlat2
    real                    :: cornerlon1, cornerlon2

    n = 1


    call ESMF_ConfigGetAttribute(LDT_config,SnowGlobeOSSEmaskData%odir, &
         label="SnowGlobe OSSE mask directory:",rc=rc)
    call LDT_verify(rc,'SnowGlobe OSSE mask directory: not defined')

    
    SnowGlobeOSSEmaskData%nc = 9900
    SnowGlobeOSSEmaskData%nr = 9240

    gridDesci = 0
    gridDesci(1) = 3 
    gridDesci(2) = SnowGlobeOSSEmaskData%nc
    gridDesci(3) = SnowGlobeOSSEmaskData%nr
    gridDesci(4) = 31.005
    gridDesci(5) = -120.995
    gridDesci(6) = 8
    gridDesci(7) = 50.0
    gridDesci(8) = 0.25
    gridDesci(9) = 0.25
    gridDesci(10) = 30.0
    gridDesci(11) = -97.9
    gridDesci(20) = 0


!-------------------------------------------------------------------
!  if the LIS output (obs) is at a coarser resolution than the 
!  LDT grid, then setup the weights for interpolation. Else 
!  setup the weights for upscaling. 
!-------------------------------------------------------------------

#if 0 
    allocate(SnowGlobeOSSEmaskData%n11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
    allocate(SnowGlobeOSSEmaskData%n12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
    allocate(SnowGlobeOSSEmaskData%n21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
    allocate(SnowGlobeOSSEmaskData%n22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
    
    allocate(SnowGlobeOSSEmaskData%w11(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
    allocate(SnowGlobeOSSEmaskData%w12(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
    allocate(SnowGlobeOSSEmaskData%w21(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
    allocate(SnowGlobeOSSEmaskData%w22(LDT_rc%lnc(n)*LDT_rc%lnr(n)))
    
    call bilinear_interp_input(n, gridDesci, &
         SnowGlobeOSSEmaskData%n11, &
         SnowGlobeOSSEmaskData%n12, SnowGlobeOSSEmaskData%n21, &
         SnowGlobeOSSEmaskData%n22, SnowGlobeOSSEmaskData%w11, &
         SnowGlobeOSSEmaskData%w12, SnowGlobeOSSEmaskData%w21, &
         SnowGlobeOSSEmaskData%w22)
#endif    
  end subroutine SnowGlobeOSSEmask_init
  
end module SnowGlobeOSSEmask_Mod
