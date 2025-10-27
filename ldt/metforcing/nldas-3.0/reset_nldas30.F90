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
! !MODULE: reset_nldas30
! \label{reset_nldas30}
! 
! !REVISION HISTORY: 
! 18 Mar 2015: James Geiger, initial code (based on merra-land)
! 12 Nov 2015: KR Arsenault, added to LDT
! 23 Oct 2025: M Wrzesien, NLDAS-3 code based on MERRA2
! 
! !INTERFACE:
subroutine reset_nldas30
! !USES:
  use LDT_coreMod,  only : LDT_rc
  use LDT_timeMgrMod, only : LDT_date2time
  use nldas30_forcingMod
!
! !DESCRIPTION:
!  Routine to cleanup allocated structures for nldas30 forcing. 
!
!EOP  
  implicit none
  integer :: n 

  do n=1,LDT_rc%nnest
     nldas30_struc(n)%startFlag = .true. 
     nldas30_struc(n)%dayFlag = .true. 
     nldas30_struc(n)%nldas30time1 = 3000.0  ! original value
!     nldas30_struc(n)%nldas30time1 = 0.0  ! alternate value (KRA)
     nldas30_struc(n)%nldas30time2 = 0.0
     nldas30_struc(n)%ringtime = 0.0
     nldas30_struc(n)%reset_flag = .true.
  enddo
end subroutine reset_nldas30
