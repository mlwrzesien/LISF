!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: snowmodel_reset
! \label{snowmodel_reset}
!
! !REVISION HISTORY:
!  14 Apr 2020: Kristi Arsenault; Add G. Liston's SnowModel
! 
! !INTERFACE:
subroutine snowmodel_reset()
! !USES:
  use LIS_coreMod,    only : LIS_rc
  use snowmodel_lsmMod
  use LIS_surfaceModelDataMod, only : LIS_sfmodel_struc
  use LIS_timeMgrMod,   only : LIS_clock, LIS_calendar, &
       LIS_update_timestep
  use LIS_logMod,       only : LIS_verify, LIS_logunit
!
! !DESCRIPTION: 
! 
!  This routine is the reset point for parameters and variables
!  required for SnowModel. 
!  
! The routines invoked are: 
! \begin{description}
! \item[snowmodel\_resetvegparms](\ref{snowmodel_resetvegparms}) \newline
!   initializes the vegetation-related parameters in SnowModel
! \end{description}
!EOP
  implicit none
  integer                 :: i,n
  integer                 :: status


end subroutine snowmodel_reset
