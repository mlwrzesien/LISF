!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: snowmodel_setup
! \label{snowmodel_setup}
!
! !REVISION HISTORY:
!  14 Apr 2020: Kristi Arsenault; Add G. Liston's SnowModel 
! 
! !INTERFACE:
subroutine snowmodel_setup()
! !USES:
  use LIS_coreMod,    only : LIS_rc
  use snowmodel_lsmMod
!
! !DESCRIPTION: 
! 
!  This routine is the entry point to set up the parameters
!  required for Liston's SnowModel. These include vegetation
!  topography, and initialization of state variables in SnowModel.
!  
! The routines invoked are: 
! \begin{description}
! \item[snowmodel\_setvegparms](\ref{snowmodel_setvegparms}) \newline
!   initializes the vegetation-related parameters in SnowModel
! \end{description}
!EOP

  implicit none

!  call snowmodel_setvegparms(LIS_rc%lsm_index)

end subroutine snowmodel_setup

