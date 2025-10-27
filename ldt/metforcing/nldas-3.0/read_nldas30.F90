!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LDTF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: read_nldas30
! \label{read_nldas30}
!
! !REVISION HISTORY:
! 16 May 2025: David Mocko, Initial Specification
!                           (derived from read_merra2.F90)
! 04 Jun 2025: James Geiger, add support read subsets of NLDAS-3 domain
!
! !INTERFACE:
subroutine read_nldas30(n,order,findex,filename,nldasforc,ferror)
! !USES:
   use LDT_logMod
   use LDT_FORC_AttributesMod
   use LDT_coreMod,       only  : LDT_rc,LDT_domain,LDT_masterproc
   use LDT_metforcingMod, only  : LDT_forc
   use nldas30_forcingMod, only : nldas30_struc
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
   use netcdf
#endif

   implicit none
! !ARGUMENTS:
   integer, intent(in)          :: n
   integer, intent(in)          :: order
   integer, intent(in)          :: findex
   character(len=*), intent(in) :: filename
   real, intent(inout)          :: nldasforc(LDT_rc%met_nf(findex),24, &
      LDT_rc%lnc(n)*LDT_rc%lnr(n))
   integer, intent(out)         :: ferror

! !DESCRIPTION:
!  For the given time, reads parameters from NLDAS-3 data,
!  transforms into 8 LDT forcing parameters and interpolates
!  to the LDT domain. \newline
!
! NLDAS-3 FORCING VARIABLES (unless noted, fields are 1-hr): \newline
!  1. T 2m    Temperature interpolated to 2 metres [$K$] \newline
!  2. q 2m    Instantaneous specific humidity interpolated to 2 metres[$kg/kg$] \newline
!  3. swdn    Downward shortwave flux at the ground [$W/m^2$] \newline
!  4. lwdn    Downward longwave radiation at the ground [$W/m^2$] \newline
!  5. u 10m   Instantaneous zonal wind interpolated to 10 metres [$m/s$] \newline
!  6. v 10m   Instantaneous meridional wind interpolated to 10 metres[$m/s$] \newline
!  7. ps      Instantaneous Surface Pressure [$Pa$] \newline
!  8. preacc  Total precipitation [$mm/s$] \newline
!
!  The arguments are:
!  \begin{description}
!  \item[order]
!    flag indicating which data to be read (order=1, read the previous
!    1 hourly instance, order=2, read the next 1 hourly instance)
!  \item[n]
!    index of the nest
!  \item[name]
!    name of the daily NLDAS-3 forcing file
!  \item[tscount]
!    time step count
!  \item[ferror]
!    return error code (0 indicates success)
!  \end{description}
!
!  The routines invoked are:
!  \begin{description}
!  \item[bilinear\_interp](\ref{bilinear_interp}) \newline
!    spatially interpolate the forcing data using bilinear interpolation
!  \item[conserv\_interp](\ref{conserv_interp}) \newline
!    spatially interpolate the forcing data using conservative interpolation
!  \item[neighbor\_interp](\ref{neighbor_interp}) \newline
!    spatially interpolate the forcing data using nearest neighbor interpolation
!  \item[upscaleByAveraging](\ref{upscaleByAveraging}) \newline
!    upscales scalar data from a finer grid to a coarser grid 
!  \end{description}
!EOP

   integer   :: ftn
   integer   :: tmpId,qId,uId,vId,psId,rainfId,swdnId,lwdnId
   logical   :: file_exists
   integer   :: mo

   real      :: tair(nldas30_struc(n)%nc, nldas30_struc(n)%nr,24)
   real      :: qair(nldas30_struc(n)%nc, nldas30_struc(n)%nr,24)
   real      :: uwind(nldas30_struc(n)%nc, nldas30_struc(n)%nr,24)
   real      :: vwind(nldas30_struc(n)%nc, nldas30_struc(n)%nr,24)
   real      :: ps(nldas30_struc(n)%nc, nldas30_struc(n)%nr,24)
   real      :: rainf(nldas30_struc(n)%nc, nldas30_struc(n)%nr,24)
   real      :: swdn(nldas30_struc(n)%nc, nldas30_struc(n)%nr,24)
   real      :: lwdn(nldas30_struc(n)%nc, nldas30_struc(n)%nr,24)

   integer :: sC, sR, cC, cR
!_______________________________________________________________________

#if (defined USE_NETCDF3)
   write(LDT_logunit,*) "[ERR] NLDAS-3 reader requires NetCDF4"
   call LDT_endrun
#endif

#if (defined USE_NETCDF4)
   ferror = 0
   mo = LDT_rc%lnc(n)*LDT_rc%lnr(n) !nldas30_struc(n)%nc*nldas30_struc(n)%nr !nldas30_struc(n)%bb%NLON*nldas30_struc(n)%bb%NLAT

   ! Read NLDAS-3 fields
   inquire(file=filename,exist=file_exists)
   if (file_exists) then
      write(LDT_logunit,*) '[INFO] Reading NLDAS-3 file (bookend,',     &
         order,' ... ',trim(filename)
      call LDT_verify(nf90_open(path=trim(filename),mode=NF90_NOWRITE,  &
         ncid=ftn),'nf90_open failed in read_nldas30')

!      sC = nldas30_struc(n)%bb%i_llon
!      sR = nldas30_struc(n)%bb%i_llat
!      cC = nldas30_struc(n)%bb%NLON
!      cR = nldas30_struc(n)%bb%NLAT

      call LDT_verify(nf90_inq_varid(ftn,'Tair',tmpId), &
         'nf90_inq_varid failed for Tair in read_nldas30')
      call LDT_verify(nf90_inq_varid(ftn,'Qair',qId), &
         'nf90_inq_varid failed for Qair in read_nldas30')
      call LDT_verify(nf90_inq_varid(ftn,'PSurf',psId), &
         'nf90_inq_varid failed for PSurf in read_nldas30')
      call LDT_verify(nf90_inq_varid(ftn,'LWdown',lwdnId), &
         'nf90_inq_varid failed for LWdown in read_nldas30')
      call LDT_verify(nf90_inq_varid(ftn,'SWdown',swdnId), &
         'nf90_inq_varid failed for SWdown in read_nldas30')
      call LDT_verify(nf90_inq_varid(ftn,'Wind_E',uId), &
         'nf90_inq_varid failed for Wind_E in read_nldas30')
      call LDT_verify(nf90_inq_varid(ftn,'Wind_N',vId), &
         'nf90_inq_varid failed for Wind_N in read_nldas30')
      call LDT_verify(nf90_inq_varid(ftn,'Rainf',rainfId), &
         'nf90_inq_varid failed for Rainf in read_nldas30')

      call LDT_verify(nf90_get_var(ftn, tmpId, tair), &
         'nf90_get_var failed for tair in read_nldas30')
      call LDT_verify(nf90_get_var(ftn, qId, qair), &
         'nf90_get_var failed for qair in read_nldas30')
      call LDT_verify(nf90_get_var(ftn, psId, ps), &
         'nf90_get_var failed for ps in read_nldas30')
      call LDT_verify(nf90_get_var(ftn, lwdnId, lwdn), &
         'nf90_get_var failed for lwdn in read_nldas30')
      call LDT_verify(nf90_get_var(ftn, swdnId, swdn), &
         'nf90_get_var failed for swdn in read_nldas30')
      call LDT_verify(nf90_get_var(ftn, uId, uwind), &
         'nf90_get_var failed for uwind in read_nldas30')
      call LDT_verify(nf90_get_var(ftn, vId, vwind), &
         'nf90_get_var failed for vwind in read_nldas30')
      call LDT_verify(nf90_get_var(ftn, rainfId, rainf), &
         'nf90_get_var failed for rainf in read_nldas30')

!      call LDT_verify(nf90_get_var(ftn, tmpId, tair, &
!         start=(/sC,sR,1/), count=(/cC,cR,24/)), &
!         'nf90_get_var failed for tair in read_nldas30')
!      call LDT_verify(nf90_get_var(ftn, qId, qair, &
!         start=(/sC,sR,1/), count=(/cC,cR,24/)), &
!         'nf90_get_var failed for qair in read_nldas30')
!      call LDT_verify(nf90_get_var(ftn, psId, ps, &
!         start=(/sC,sR,1/), count=(/cC,cR,24/)), &
!         'nf90_get_var failed for ps in read_nldas30')
!      call LDT_verify(nf90_get_var(ftn, lwdnId, lwdn, &
!         start=(/sC,sR,1/), count=(/cC,cR,24/)), &
!         'nf90_get_var failed for lwdn in read_nldas30')
!      call LDT_verify(nf90_get_var(ftn, swdnId, swdn, &
!         start=(/sC,sR,1/), count=(/cC,cR,24/)), &
!         'nf90_get_var failed for swdn in read_nldas30')
!      call LDT_verify(nf90_get_var(ftn, uId, uwind, &
!         start=(/sC,sR,1/), count=(/cC,cR,24/)), &
!         'nf90_get_var failed for uwind in read_nldas30')
!      call LDT_verify(nf90_get_var(ftn, vId, vwind, &
!         start=(/sC,sR,1/), count=(/cC,cR,24/)), &
!         'nf90_get_var failed for vwind in read_nldas30')
!      call LDT_verify(nf90_get_var(ftn, rainfId, rainf, &
!         start=(/sC,sR,1/), count=(/cC,cR,24/)), &
!         'nf90_get_var failed for rainf in read_nldas30')

      call LDT_verify(nf90_close(ftn),'failed to close in read_nldas30')

      call interp_nldas30_var(n,findex,tair,  1,.false.,nldasforc)
      call interp_nldas30_var(n,findex,qair,  2,.false.,nldasforc)
      call interp_nldas30_var(n,findex,ps,    7,.false.,nldasforc)
      call interp_nldas30_var(n,findex,lwdn,  4,.false.,nldasforc)
      call interp_nldas30_var(n,findex,swdn,  3,.false.,nldasforc)
      call interp_nldas30_var(n,findex,uwind, 5,.false.,nldasforc)
      call interp_nldas30_var(n,findex,vwind, 6,.false.,nldasforc)
      call interp_nldas30_var(n,findex,rainf, 8, .true.,nldasforc)
   else
      write(LDT_logunit,*) '[ERR] ',trim(filename)//' does not exist'
      call LDT_endrun
   endif
#endif
end subroutine read_nldas30

!BOP
!
! !ROUTINE: interp_nldas30_var
! \label{interp_nldas30_var}
!
! !INTERFACE:
subroutine interp_nldas30_var(n,findex,input_var,var_index, &
   pcp_flag,nldasforc)

! !USES:
   use LDT_coreMod
   use LDT_logMod
   use LDT_spatialDownscalingMod
   use nldas30_forcingMod, only : nldas30_struc
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
   use netcdf
#endif
   implicit none

! !ARGUMENTS:
   integer, intent(in)    :: n
   integer, intent(in)    :: findex
   real,    intent(in)    :: input_var(nldas30_struc(n)%nc, nldas30_struc(n)%nr,24)
   integer, intent(in)    :: var_index
   logical, intent(in)    :: pcp_flag
   real,    intent(inout) :: nldasforc(LDT_rc%met_nf(findex), &
                                      24, LDT_rc%lnc(n)*LDT_rc%lnr(n))
! !DESCRIPTION:
!  This subroutine spatially interpolates a NLDAS-3 field
!  to the LDT running domain
!
!EOP
   integer   :: t,c,r,k,iret
   integer   :: doy
   integer   :: ftn
   integer   :: pcp1Id,pcp2Id,pcp3Id,pcp4Id,pcp5Id,pcp6Id
   real      :: f (nldas30_struc(n)%nc*nldas30_struc(n)%nr)
   logical*1 :: lb(nldas30_struc(n)%nc*nldas30_struc(n)%nr)
   logical*1 :: lo(LDT_rc%lnc(n)*LDT_rc%lnr(n))
   integer   :: input_size
   logical   :: scal_read_flag
! _____________________________________________________________

!   input_size = nldas30_struc(n)%bb%NLON*nldas30_struc(n)%bb%NLAT

   do t = 1,24
      lb = .true.
      do r = 1, nldas30_struc(n)%nr !bb%NLAT
         do c = 1, nldas30_struc(n)%nc !bb%NLON
            k = c+(r-1)*nldas30_struc(n)%nc  !bb%NLON
            f(k) = input_var(c,r,t)
            if (f(k).eq.-9999.0) then
               f(k)  = LDT_rc%udef
               lb(k) = .false.
            endif
         enddo
      enddo

      if (pcp_flag.and. &
         trim(LDT_rc%met_gridtransform(findex)).eq."budget-bilinear") then
         call conserv_interp(LDT_rc%gridDesc(n,:),lb,f,lo,            &
            nldasforc(var_index,t,:),                                 &
            nldas30_struc(n)%mi,LDT_rc%lnc(n)*LDT_rc%lnr(n),                   &
            LDT_domain(n)%lat,LDT_domain(n)%lon,                      &
            nldas30_struc(n)%w112,nldas30_struc(n)%w122,              &
            nldas30_struc(n)%w212,nldas30_struc(n)%w222,              &
            nldas30_struc(n)%n112,nldas30_struc(n)%n122,              &
            nldas30_struc(n)%n212,nldas30_struc(n)%n222,              &
            LDT_rc%udef,iret)
      elseif ((trim(LDT_rc%met_gridtransform(findex)).eq."bilinear").or. &
              (trim(LDT_rc%met_gridtransform(findex)).eq."budget-bilinear")) then
            call bilinear_interp(LDT_rc%gridDesc(n,:),lb,f,lo,           &
               nldasforc(var_index,t,:),                                 &
               nldas30_struc(n)%mi,LDT_rc%lnc(n)*LDT_rc%lnr(n),                   &
               LDT_domain(n)%lat,LDT_domain(n)%lon,                      &
               nldas30_struc(n)%w111,nldas30_struc(n)%w121,              &
               nldas30_struc(n)%w211,nldas30_struc(n)%w221,              &
               nldas30_struc(n)%n111,nldas30_struc(n)%n121,              &
               nldas30_struc(n)%n211,nldas30_struc(n)%n221,              &
               LDT_rc%udef,iret)
      elseif (trim(LDT_rc%met_gridtransform(findex)).eq."neighbor") then
         call neighbor_interp(LDT_rc%gridDesc(n,:),lb,f,lo,           &
            nldasforc(var_index,t,:),nldas30_struc(n)%mi,                      &
            LDT_rc%lnc(n)*LDT_rc%lnr(n),                              &
            LDT_domain(n)%lat,LDT_domain(n)%lon,                      &
            nldas30_struc(n)%n113,LDT_rc%udef,iret)
      elseif (trim(LDT_rc%met_gridtransform(findex)).eq."average") then
         call upscaleByAveraging(nldas30_struc(n)%mi, LDT_rc%lnc(n)*LDT_rc%lnr(n), &
            LDT_rc%udef, nldas30_struc(n)%n111, &
            lb, f, lo, nldasforc(var_index,t,:))
      else
         write(LDT_logunit,*) '[ERR] Spatial interpolation option '// &
            trim(LDT_rc%met_gridtransform(findex))// &
            ' not supported for NLDAS-3'
         call LDT_endrun
      endif
   enddo

end subroutine interp_nldas30_var
