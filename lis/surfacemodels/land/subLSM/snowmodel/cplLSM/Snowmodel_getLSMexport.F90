!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: snowmodel_getLSMexport
! \label{snowmodel_getLSMexport}
!
! !REVISION HISTORY:
! 19 Sep 2020: Sujay Kumar; Initial Specification
!  9 Dec 2020: Mahdi Navari; edited to take into account the Crocus slope correction 
!  2 Aug 2021: Kristi Arsenault; Edited to support SnowModel implementation
!
! !INTERFACE:
subroutine snowmodel_getLSMexport(n, SubLSM2LSM_State)
! !USES:

  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use snowmodel_lsmMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: n
  type(ESMF_State)       :: SubLSM2LSM_State
! 
! !DESCRIPTION:
! 
! 
!EOP

  type(ESMF_Field)   :: snwdField, sweField
  real, pointer      :: swe(:), snwd(:)
  integer            :: t
  integer            :: status
  real               :: tmp_SLOPE

 
  call ESMF_StateGet(SubLSM2LSM_State,"Total SWE",sweField,rc=status)
  call LIS_verify(status)
  call ESMF_StateGet(SubLSM2LSM_State,"Total snowdepth",snwdField,rc=status)
  call LIS_verify(status)

  call ESMF_FieldGet(sweField,localDE=0,farrayPtr=swe,rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(snwdField,localDE=0,farrayPtr=snwd,rc=status)
  call LIS_verify(status)

  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
!     tmp_SLOPE = snowmodel_struc(n)%crocus81(t)%SLOPE
!     swe(t) = snowmodel_struc(n)%crocus81(t)%SWE_1D / COS(tmp_SLOPE)
!     snwd(t) = snowmodel_struc(n)%crocus81(t)%SD_1D / COS(tmp_SLOPE)
!  SnowModel
     swe(t)  = snowmodel_struc(n)%sm(t)%swe_depth 
     snwd(t) = snowmodel_struc(n)%sm(t)%snow_depth
  enddo

end subroutine snowmodel_getLSMexport


