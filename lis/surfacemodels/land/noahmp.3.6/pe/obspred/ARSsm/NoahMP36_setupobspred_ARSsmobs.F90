!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: NoahMP36_setupobspred_ARSsmobs
!  \label{NoahMP36_setupobspred_ARSsmobs}
!
! !REVISION HISTORY:
! 2 Feb 2018: Soni Yatheendradas; Initial Specification
!
! !INTERFACE:
subroutine NoahMP36_setupobspred_ARSsmobs(OBSPred)
! !USES:
  use ESMF
  use LIS_coreMod,      only : LIS_rc, LIS_vecPatch
  use LIS_logMod,       only : LIS_verify

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)       :: OBSPred
!
! !DESCRIPTION:
!  
!  This routine creates an entry in the Obs pred object used for 
!  parameter estimation
! 
!EOP
  integer                :: n
  type(ESMF_ArraySpec)   :: realarrspec
  type(ESMF_Field)       :: smcField
!  type(ESMF_Field)       :: smstdField
  integer                :: status

  n = 1
  call ESMF_ArraySpecSet(realarrspec, rank=1,typekind=ESMF_TYPEKIND_R4,&
       rc=status)
  call LIS_verify(status)

  smcField = ESMF_FieldCreate(arrayspec=realarrspec, grid=LIS_vecPatch(n,LIS_rc%lsm_index), &
       name="ARS_sm", rc=status)
  call LIS_verify(status)
  
  call ESMF_StateAdd(OBSPred,(/smcField/),rc=status)
  call LIS_verify(status)

!  smstdField = ESMF_FieldCreate(arrayspec=realarrspec, grid=LIS_vecPatch(n,LIS_rc%lsm_index), &
!       name="ARSsm standard deviation of soil moisture", rc=status)
!  call LIS_verify(status)
!  
!  call ESMF_StateAdd(OBSPred,(/smstdField/),rc=status)
!  call LIS_verify(status)
end subroutine NoahMP36_setupobspred_ARSsmobs

