!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: snowmodel_readrst
! \label{snowmodel_readrst}
!
! !REVISION HISTORY:
!  14 Apr 2020: Kristi Arsenault; Add G. Liston's SnowModel 
!
! !INTERFACE:
subroutine snowmodel_readrst()
! !USES:
  use LIS_coreMod,    only : LIS_rc, LIS_masterproc
  use LIS_historyMod, only : LIS_readvar_restart
  use LIS_logMod,     only : LIS_logunit, LIS_endrun, &
       LIS_getNextUnitNumber, LIS_releaseUnitNumber,&
       LIS_verify
  use snowmodel_lsmMod
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
!
! !DESCRIPTION:
!  This program reads restart files for SnowModel.  This
!  includes all relevant water/energy storages and tile information. 
!  The following is the list of variables specified in the SnowModel 
!  restart file: 
!
!  \begin{verbatim}
!   nc,nr,ntiles    - grid and tile space dimensions 
!   snowh        - SnowModel snow depth
!   sneqv        - SnowModel snow water equivalent
!  \end{verbatim}
!
!  The routines invoked are: 
! \begin{description}
! \item[LIS\_readvar\_restart](\ref{LIS_readvar_restart}) \newline
!  reads a variable from the restart file
! \item[snowmodel\_coldstart](\ref{snowmodel_coldstart}) \newline
!   initializes the SnowModel state variables
! \end{description}
!EOP

  implicit none

  character*20      :: wformat

  wformat="netcdf"

  write(LIS_logunit,*) '[INFO] Call to the SnowModel Read restart routine ...'

end subroutine snowmodel_readrst

