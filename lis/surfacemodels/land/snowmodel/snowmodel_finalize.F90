!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: snowmodel_finalize
! \label{snowmodel_finalize}
!
! !REVISION HISTORY:
!  14 Apr 2020: Kristi Arsenault; Add G. Liston's SnowModel 
!
! !INTERFACE:
subroutine snowmodel_finalize()
! !USES:
  use LIS_coreMod,  only : LIS_rc
  use snowmodel_lsmMod
!
! !DESCRIPTION:
!  
!  This routine cleans up the allocated memory structures in SnowModel 
!  
!EOP
  implicit none

  integer :: t,n

  ! Print a banner when the model run is finished.
  print *
  print *,&
       & 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
  print *,&
       & 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
  print *,&
       & '                 The SnowModel Run Has Finished                '
  print *,&
       & 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
  print *,&
       & 'ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc'
  print *


end subroutine snowmodel_finalize


