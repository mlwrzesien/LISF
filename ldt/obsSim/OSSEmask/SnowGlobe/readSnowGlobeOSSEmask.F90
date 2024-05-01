!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

#include "LDT_misc.h"
!BOP
! 
! !ROUTINE: readSnowGlobeOSSEmask
! \label{readSnowGlobeOSSEmask}
!
! !INTERFACE: 
subroutine readSnowGlobeOSSEmask(n)
! !USES:
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
#if (defined USE_GRIBAPI)
  use grib_api
#endif
  use LDT_coreMod
  use LDT_obsSimMod
  use LDT_historyMod
  use LDT_logMod
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use SnowGlobeOSSEmask_Mod
!
! !DESCRIPTION: 
!  This routine reads the SnowGlobe orbital masks and interpolates to the
!  LDT domain
!EOP
  implicit none

  integer,   intent(in) :: n
  integer,   parameter  :: ndays = 16
  character(len=LDT_CONST_PATH_LEN)         :: fname
  integer               :: c,r,k
  integer               :: maskid
  logical               :: file_exists
  integer               :: ftn
  integer               :: iret
  integer               :: tindex
  real                  :: mask(SnowGlobeOSSEmaskData%nc,&
       SnowGlobeOSSEmaskData%nr)
!  real                  :: mask1d(SnowGlobeOSSEmaskData%nc*SnowGlobeOSSEmaskData%nr)
!  real                  :: mask_2d(LDT_rc%lnc(n),LDT_rc%lnr(n))

  !index 1 in the data is assumed to be jan 1. The data repeats
  !every 25 days.
  !ignoring leap years for now     
  
  if(mod(LDT_rc%doy,25).ne.0) then 
     tindex = mod(LDT_rc%doy,25)
  else
     tindex = 25
  endif
  
  ! create SnowGlobe filename
  call create_SnowGlobe_ossemask_filename(&
       SnowGlobeOSSEmaskData%odir,&
       tindex,&
       fname)

! read, subset and interpolate the data
  
  inquire(file=trim(fname),exist=file_exists)

  if(file_exists) then 
     write(LDT_logunit,*) '[INFO] reading SnowGlobe mask ',trim(fname)

#if(defined USE_NETCDF3 || defined USE_NETCDF4) 
        
     iret = nf90_open(path=trim(fname),mode=nf90_nowrite, ncid=ftn)
     call LDT_verify(iret, 'Error opening file '//trim(fname))
     
     call LDT_verify(nf90_inq_varid(ftn, "track", maskid),&
          'nf90_inq_varid failed in readSnowGlobeOSSEmask')

     call LDT_verify(nf90_get_var(ftn,maskid,mask),&
          'nf90_get_Var failed in readSnowGlobeOSSEmask')

     
!     mask1d(:) = 0.0

!     do r=1,SnowGlobeOSSEmaskData%nr
!        do c=1, SnowGlobeOSSEmaskData%nc
!           mask1d(c+(r-1)*SnowGlobeOSSEmaskData%nc) = &
!                mask(r,c)
!        enddo
!     enddo

!     call convertSnowGlobeOSSEmaskToLDTgrid(n,mask1d(:),mask_2d(:,:))

#endif
  else
     write(LDT_logunit,*) '[WARN] file '//trim(fname)
     write(LDT_logunit,*) '[WARN] not found ...'
     mask = LDT_rc%udef
  endif

  !  LDT_obsSim_struc%datamask = mask_2d
    LDT_obsSim_struc%datamask = mask

end subroutine readSnowGlobeOSSEmask

!BOP
!
! !ROUTINE: create_SnowGlobe_ossemask_filename
! \label{create_SnowGlobe_ossemask_filename}
!
! !INTERFACE:
subroutine create_SnowGlobe_ossemask_filename(odir,tindex, fname)
  !USES:
   use LDT_coreMod,  only : LDT_rc
   use LDT_logMod

   implicit none 
! !ARGUMENTS:

   character(len=*)             :: odir
   real                         :: tindex
   character(len=*)             :: fname
! 
! !DESCRIPTION:  
!  Create the file name for the SnowGlobe mask files
!
!
!EOP

   character*10    :: temp

   write(unit=temp,fmt='(i2.2)') tindex
   fname = trim(odir)//'/SnowGlobe_mask_'//trim(temp)//'.nc'
   
 end subroutine create_SnowGlobe_ossemask_filename

