!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LVT_misc.h"
!BOP
!
! !ROUTINE: readAMSRcnnSnowObs
! \label{readAMSRcnnSnowObs}
!
! !INTERFACE:
subroutine readAMSRcnnSnowObs(source)
!
! !USES:
   use ESMF
   use LVT_histDataMod
   use LVT_coreMod,    only : LVT_rc
   use LVT_timeMgrMod, only : LVT_calendar
   use LVT_logMod,     only : LVT_logunit, LVT_verify
   use AMSRcnnSnow_obsMod,  only : amsrcnnsnowobs
   use map_utils

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
   use netcdf
#endif

   implicit none
!
! !INPUT PARAMETERS:
   integer,  intent(in) :: source
!
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
! This subroutine reads and processes the gridded UA SNOW
! data.  The data for a given year is read into memory at
! the start of a year and indexed into during each day.
!
! !REVISION HISTORY:
! 28 May 2019: Rhae Sung Kim, Initial Specification
! 19 Jun 2019: David Mocko, Set valid time of data to 12Z
! 25 Jun 2019: David Mocko, Only read one file per call
!
!EOP

   integer             :: ftn
   character*100       :: fname
   logical             :: file_exists
   real, allocatable   :: snwd1(:,:)
   real                :: snwd_in(amsrcnnsnowobs(source)%nc*amsrcnnsnowobs(source)%nr)
   logical*1           :: lb(amsrcnnsnowobs(source)%nc*amsrcnnsnowobs(source)%nr)
   logical*1           :: lo(LVT_rc%lnc*LVT_rc%lnr)
   real                :: snwd_out(LVT_rc%lnc*LVT_rc%lnr)
   real                :: snwd_final(LVT_rc%lnc, LVT_rc%lnr)
   integer             :: t,c,r,k,nt
   integer             :: nid,varid
   integer             :: iret
   type(ESMF_Time)     :: uatime1
   integer             :: status
   real                :: timenow
   logical             :: alarmCheck

   snwd_out = LVT_rc%udef
   snwd_final = LVT_rc%udef

   timenow = float(LVT_rc%dhr(source)*3600 + &
        LVT_rc%dmn(source)*60 + LVT_rc%dss(source))
   
! The UA dataset is a daily product, valid approximately
! at 12Z, which is the observation time of many of the
! surface products used to generate the data (personal
! communication with the UA SNOW product generators).
   alarmcheck = (mod(timenow, 86400.0).eq.0)

   if (alarmcheck) then
      LVT_rc%resetFlag(source) = .false.

      call create_AMSRcnn_filename(amsrcnnsnowobs(source)%odir, &
           LVT_rc%dyr(source), LVT_rc%dmo(source), LVT_rc%dda(source),&
           fname)

      inquire(file=trim(fname),exist=file_exists)

      if (file_exists) then
         write(LVT_logunit,*) '[INFO] Reading AMSR CNN data ',&
                                trim(fname)
         
         allocate(snwd1(amsrcnnsnowobs(source)%nc,   &
              amsrcnnsnowobs(source)%nr))

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
         
         iret = nf90_open(path=trim(fname),mode=NF90_NOWRITE,ncid=nid)
         call LVT_verify(iret, 'Error opening file'//trim(fname))

         iret = nf90_inq_varid(nid, 'SD',varid)
         call LVT_verify(iret, 'Error nf90_inq_varid: SD')

         iret = nf90_get_var(nid,varid, SNWD1)
         call LVT_verify(iret, 'Error nf90_get_var: SD')

         iret = nf90_close(nid)
         call LVT_verify(iret, 'Error nf90_close')
#endif
         do c=1,amsrcnnsnowobs(source)%nc
            do r=1,amsrcnnsnowobs(source)%nr
               if(snwd1(c,r).lt.0.or.snwd1(c,r).gt.1000) then
                  snwd_in(c+(amsrcnnsnowobs(source)%nr-r)*&
                       amsrcnnsnowobs(source)%nc) = -9999.0
               else
                  snwd_in(c+(amsrcnnsnowobs(source)%nr-r)*&
                       amsrcnnsnowobs(source)%nc) = snwd1(c,r)
               endif
            enddo
         enddo
         
!         print*, amsrcnnsnowobs(source)%nc,amsrcnnsnowobs(source)%nr
!         open(100,file='test.bin',form='unformatted')
!         write(100) snwd_in
!         close(100)
!         stop
         deallocate(snwd1)        

      else
         write(LVT_logunit,*) '[WARN] AMSR CNN file not found: ',&
              trim(fname)
         snwd_in = -9999.0
      endif

      lb = .false.
      
      do t=1,amsrcnnsnowobs(source)%nc*amsrcnnsnowobs(source)%nr
         if(snwd_in(t).ne.-9999.0) then 
            lb(t) = .true.
         endif
      enddo


      call bilinear_interp(LVT_rc%gridDesc,lb,snwd_in, &
           lo,snwd_out, &
           amsrcnnsnowobs(source)%nc*amsrcnnsnowobs(source)%nr, &
           LVT_rc%lnc*LVT_rc%lnr,  &
           amsrcnnsnowobs(source)%rlat, &
           amsrcnnsnowobs(source)%rlon, &
           amsrcnnsnowobs(source)%w11,  &
           amsrcnnsnowobs(source)%w12,  &
           amsrcnnsnowobs(source)%w21,  &
           amsrcnnsnowobs(source)%w22,  &
           amsrcnnsnowobs(source)%n11,  &
           amsrcnnsnowobs(source)%n12,  &
           amsrcnnsnowobs(source)%n21,  &
           amsrcnnsnowobs(source)%n22,  &
           LVT_rc%udef, iret)

      do r=1,LVT_rc%lnr
         do c=1, LVT_rc%lnc
            snwd_final(c,r) = snwd_out(c+(r-1)*LVT_rc%lnc)
         enddo
      enddo

      do r=1,LVT_rc%lnr
         do c=1,LVT_rc%lnc

            if(snwd_final(c,r).ge.0) then
               snwd_final(c,r) = snwd_final(c,r)/100.0 ! Convert cm to m
            else
               snwd_final(c,r) = LVT_rc%udef
            endif
         enddo
      enddo
!      print*, 'here'
!      open(100,file='test.bin',form='unformatted')
!      write(100) snwd_final
!      close(100)
!      stop
   endif

   call LVT_logSingleDataStreamVar(LVT_MOC_SNOWDEPTH,source,snwd_final,vlevel=1,units="m")

end subroutine readAMSRcnnSnowObs

!BOP
!
! !ROUTINE: create_AMSRcnn_filename
! \label(create_AMSRcnn_filename)
!
! !INTERFACE:
subroutine create_AMSRcnn_filename(odir,yr,mo,da,fname)
!
! !USES:
   use LVT_String_Utility
   implicit none
!
! !INPUT PARAMETERS:
   character(len=*), intent(in)  :: odir
   integer,          intent(in)  :: yr,mo,da
!
! !OUTPUT PARAMETERS:
   character(len=*), intent(out) :: fname
!
   character*4                   :: fyr
   character*2                   :: fmo
   character*2                   :: fda

   write(fyr, '(i4.4)' ) yr
   write(fmo, '(i2.2)' ) mo
   write(fda, '(i2.2)' ) da

   fname = trim(odir)//'/'//trim(fyr)//&
        '/AMSR_CNN_NH_'//trim(fyr)//'-'//&
        trim(fmo)//'-'//trim(fda)//'.nc'

end subroutine create_AMSRcnn_filename
