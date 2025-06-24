!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: read_simSnowGlobePMWobs
!  \label{read_simSnowGlobePMWobs}
!
! !REVISION HISTORY:
!  21Jun2006: Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine read_simSnowGlobePMWobs(n, k, OBS_State, OBS_Pert_State) 
  ! !USES:

  use ESMF
  use LIS_mpiMod
  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use LIS_dataAssimMod
  use LIS_pluginIndices
  use LIS_DAobservationsMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use simSnowGlobePMWobs_module, only : SnowGlobePMW_struc
  
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif
  
  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  integer, intent(in) :: k
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!  
!  reads the synthetic SWE observations produced from a LIS control run. 
!  
!  The arguments are: 
!  \begin{description}
!  \item[n]    index of the nest
!  \item[OBS\_State] observations state
!  \end{description}
!
!EOP
  type(ESMF_Field)    :: sweField

  real,    pointer    :: obsl(:)
  real                :: sweobs_out(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real                :: sweobs(SnowGlobePMW_struc(n)%nc,SnowGlobePMW_struc(n)%nr)
  integer             :: gid(LIS_rc%obs_ngrid(k))
  integer             :: assimflag(LIS_rc%obs_ngrid(k))
  real, allocatable       :: dummy(:)

  character(len=LIS_CONST_PATH_LEN) :: sweobsdir
  character(len=LIS_CONST_PATH_LEN) :: name

  logical             :: data_update
  logical                :: data_upd_flag(LIS_npes)
  logical                :: data_upd_flag_local
  logical                :: data_upd

  logical*1, allocatable :: swe_data_b(:)
  logical*1              :: sweobs_b_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  logical             :: file_exists
  integer*2, allocatable :: var(:,:)
  real, allocatable      :: swe1d(:)

  logical             :: readflag
  integer             :: sweid, status
  integer             :: fnd
  integer             :: c,r
  integer             :: ftn,t
  integer             :: iret

  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       sweobsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(OBS_State,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  call simSnowGlobePMWSWE_filename(name,sweobsdir,&
       LIS_rc%yr,LIS_rc%mo,LIS_rc%da,LIS_rc%hr,LIS_rc%mn)
!  print *,'sweobsdir=',sweobsdir
!  print *,'name=',name

  inquire(file=name,exist=file_exists)

  if(file_exists) then 
     readflag = .true.
!     print *,'reading SimObs file' 
  else 
     readflag = .false.
  endif

  if (readflag) then 
     allocate(dummy(LIS_rc%obs_ngrid(k)))
     allocate(swe_data_b(SnowGlobePMW_struc(n)%nc*SnowGlobePMW_struc(n)%nr))
     allocate(var(SnowGlobePMW_struc(n)%nc,SnowGlobePMW_struc(n)%nr))
     allocate(swe1d(SnowGlobePMW_struc(n)%nc*SnowGlobePMW_struc(n)%nr))

     call ESMF_StateGet(OBS_State,"Observation01",sweField,&
          rc=status)
     call LIS_verify(status)
     call ESMF_FieldGet(sweField,localDE=0, farrayPtr=obsl,rc=status)
     call LIS_verify(status)
  
     write(LIS_logunit,*)  'Reading syn data ',trim(name)
     
     call ESMF_StateGet(OBS_State,"Observation01",sweField,&
          rc=status)
     call LIS_verify(status)
     call ESMF_FieldGet(sweField,localDE=0, farrayPtr=obsl,rc=status)
     call LIS_verify(status)
     
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
     call LIS_verify(nf90_open(path=trim(name),mode=NF90_NOWRITE,ncid=ftn),&
          'Error opening file '//trim(name))
     call LIS_verify(nf90_inq_varid(ftn,'SWE_tavg',sweid),&
          'Error nf90_inq_varid: SWE_tavg')

     call LIS_verify(nf90_get_var(ftn,sweid,sweobs, &
          count=(/SnowGlobePMW_struc(n)%nc,SnowGlobePMW_struc(n)%nr,1/),&
          start=(/1,1,1/)),&
          'Error in nf90_get_var')
 !    print *,'within block'
     call LIS_verify(nf90_close(ftn))


!--------------------------------------------------------------------------
! Interpolate to the observation grid
!-------------------------------------------------------------------------- 
        swe1d = LIS_rc%udef
        swe_data_b = .false.

        do r=1, SnowGlobePMW_struc(n)%nr
           do c=1, SnowGlobePMW_struc(n)%nc
              if(sweobs(c,SnowGlobePMW_struc(n)%nr-r+1).ge.0) then
                 swe1d(c+(r-1)*SnowGlobePMW_struc(n)%nc) = &
                      sweobs(c,SnowGlobePMW_struc(n)%nr-r+1)
                 swe_data_b(c+(r-1)*SnowGlobePMW_struc(n)%nc)=.true.
              endif
           enddo
        enddo
        deallocate(var)

!        call upscaleByAveraging(&
!             SnowGlobePMW_struc(n)%nc*SnowGlobePMW_struc(n)%nr,&
!             LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k),&
!             LIS_rc%udef,&
!             SnowGlobePMW_struc(n)%n11,&
!             swe_data_b,swe1d,&
!             sweobs_b_ip,sweobs_out)

        call bilinear_interp(LIS_rc%obs_gridDesc(k,:),&
             swe_data_b, swe1d, sweobs_b_ip, sweobs_out, &
             SnowGlobePMW_struc(n)%nc*SnowGlobePMW_struc(n)%nr, &
             LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k), &
             SnowGlobePMW_struc(n)%rlat,SnowGlobePMW_struc(n)%rlon,&
             SnowGlobePMW_struc(n)%w11,SnowGlobePMW_struc(n)%w12,&
             SnowGlobePMW_struc(n)%w21,SnowGlobePMW_struc(n)%w22,&
             SnowGlobePMW_struc(n)%n11,SnowGlobePMW_struc(n)%n12,&
             SnowGlobePMW_struc(n)%n21,SnowGlobePMW_struc(n)%n22,LIS_rc%udef,iret)

        deallocate(swe1d)
        deallocate(swe_data_b)

!end MLW edits
     do r =1,LIS_rc%obs_lnr(k)
        do c =1,LIS_rc%obs_lnc(k)
           if (LIS_obs_domain(n,k)%gindex(c,r) .ne. -1)then
              obsl(LIS_obs_domain(n,k)%gindex(c,r)) = &
                   sweobs_out(c+(r-1)*LIS_rc%obs_lnc(k))
           end if
        end do
     end do
     
#endif
     !-------------------------------------------------------------------------
     !  Apply LSM based QC and screening of observations
     !-------------------------------------------------------------------------
     call lsmdaqcobsstate(trim(LIS_rc%lsm)//"+" &
          //trim(LIS_SnowGlobePMWSWEId)//char(0), n, k, OBS_state)

     if(SnowGlobePMW_struc(n)%conditionalScreen.eq.1) then 
        call applySnowGlobeConditionalDAflags(n,k,OBS_State)
     endif
     call LIS_checkForValidObs(n, k, obsl, fnd, sweobs)     

     if (fnd .eq. 0) then
        data_upd_flag_local = .false.
     else
        data_upd_flag_local = .true.
     endif
     
#if (defined SPMD)
     call MPI_ALLGATHER(data_upd_flag_local, 1, &
          MPI_LOGICAL, data_upd_flag(:), &
          1, MPI_LOGICAL, LIS_mpi_comm, status)
     data_upd = any(data_upd_flag)
#else
     data_upd = data_upd_flag_local
#endif

     if (data_upd) then          
        
        do t=1,LIS_rc%obs_ngrid(k)
           gid(t) = t
           if(obsl(t).ne.-9999.0) then 
              assimflag(t) = 1 
           else
              assimflag(t) = 0 
           endif
        enddo
        
        call ESMF_AttributeSet(OBS_State,"Data Update Status",&
             .true., rc=status)
        call LIS_verify(status)
        
        if (LIS_rc%obs_ngrid(k) .gt. 0) then
           call ESMF_AttributeSet(sweField,"Grid Number",&
                gid, itemCount=LIS_rc%obs_ngrid(k),rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeSet(sweField,"Assimilation Flag",&
                assimflag,itemCount=LIS_rc%obs_ngrid(k),rc=status)
           call LIS_verify(status)
        endif
        
     else
        call ESMF_AttributeSet(OBS_State,"Data Update Status",&
             .false., rc=status)
        call LIS_verify(status)
        return
     end if
  
     do t=1,LIS_rc%obs_ngrid(k)
        if(obsl(t).ne.-9999.0) then 
           if(obsl(t).lt.0.0) obsl(t) = 0.0
!        if(obsl(t).gt.200.0 ) obsl(t) = 200.0
        endif
     enddo
     deallocate(dummy)
  else
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status)
  end if

end subroutine read_simSnowGlobePMWobs

subroutine simSnowGlobePMWSWE_filename(name, ndir, yr, mo,da,hr,mn)

  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  implicit none
  character(len=LIS_CONST_PATH_LEN) :: name
  integer           :: yr, mo, da, hr,mn
  character (len=*) :: ndir
  character (len=4) :: fyr
  character (len=2) :: fmo,fda,fhr,fmn
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
  write(unit=fhr, fmt='(i2.2)') hr
  write(unit=fmn, fmt='(i2.2)') mn  

  name = trim(ndir)//'/'//trim(fyr)//trim(fmo)//'/SimObs_'//&
       trim(fyr)//trim(fmo)//trim(fda)//trim(fhr)//&
       trim(fmn)//'.nc'  

end subroutine simSnowGlobePMWSWE_filename

subroutine applySnowGlobeConditionalDAflags(n,k,OBS_State)

  use ESMF
  use netcdf
  use LIS_coreMod
  use LIS_timeMgrMod
  use LIS_logMod
  use LIS_DAobservationsMod
  
  implicit none

  integer             :: n
  integer             :: k
  type(ESMF_State)    :: OBS_State

  integer                 :: syr,smo,sda,shr,smn,sss
  integer                 :: eyr,emo,eda,ehr,emn,ess
  integer                 :: yr,mo,da,hr,mn,ss
  integer                 :: status
  type(ESMF_Time)         :: startTime
  type(ESMF_Time)         :: currTime
  type(ESMF_Time)         :: stopTime
  type(ESMF_TimeInterval) :: tw,ts
  logical                 :: file_exists,togo
  integer                 :: ftn,c,r
  character*100           :: innovfile,cdate1,cdate2
  integer                 :: ninnvId
  integer                 :: tw_sec
  real                    :: innovmask(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  real                    :: tmask(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  type(ESMF_Field)        :: sweField
  real, pointer           :: obsl(:)


  call ESMF_StateGet(OBS_State,"Observation01",sweField,&
       rc=status)
  call LIS_verify(status)
  call ESMF_FieldGet(sweField,localDE=0, farrayPtr=obsl,rc=status)
  call LIS_verify(status)

  call ESMF_TimeIntervalSet(ts,s=nint(LIS_rc%ts),rc=status)
  
  eyr = LIS_rc%yr
  emo = LIS_rc%mo
  eda = LIS_rc%da
  ehr = LIS_rc%hr
  emn = LIS_rc%mn
  ess = LIS_rc%ss
  
  call ESMF_TimeSet(stopTime, yy = eyr, &
       mm = emo, &
       dd = eda, &
       h  = ehr, &
       m  = emn, &
       s  = ess, &
       calendar = LIS_calendar, &
       rc = status)

  tw_sec = nint(7*86400.0)
  call ESMF_TimeIntervalSet(tw,s=tw_sec,rc=status)

  startTime = stopTime - tw

  call ESMF_TimeGet(startTime, yy = syr, &
       mm = smo, &
       dd = sda, &
       h  = shr, &
       m  = smn, &
       s  = sss, &
       calendar = LIS_calendar, &
       rc = status)

  togo = .true.
  currTime = startTime

  innovmask = 0.0
  do while (togo)

     call ESMF_TimeGet(currTime, yy = yr, &
          mm = mo, &
          dd = da, &
          h  = hr, &
          m  = mn, &
          s  = ss, &
          calendar = LIS_calendar, &
          rc = status)

     !
     !Caveats: 
     ! assuming that the files are in the 3-level hierachy format
     ! also hardcoding to look for the second assimilation instance
     ! (snowglobe DA)
     ! assuming that observation space and model space are the same
     !
     write(unit=cdate1, fmt='(i4.4, i2.2)') &
          yr, mo
     write(unit=cdate2, fmt='(i4.4, i2.2, i2.2, i2.2, i2.2)') &
          yr, mo,da,hr,mn     

     innovfile = trim(LIS_rc%odir)//'/EnKF/'//trim(cdate1)//&
          '/LIS_DA_EnKF_'//trim(cdate2)//'_innov.a02.d01.nc'

     inquire(file=trim(innovfile),exist=file_exists)
     
     if(file_exists) then
        write(LIS_logunit,*) 'Reading ',trim(innovfile)
        call LIS_verify(nf90_open(path=trim(innovfile),&
             mode=nf90_nowrite,&
             ncid = ftn))
        call LIS_verify(nf90_inq_varid(ftn,"ninnov_02",&
             ninnvId))
        call LIS_verify(nf90_get_var(ftn,ninnvId,&
             tmask,&
             start=(/LIS_ews_obs_halo_ind(n,LIS_localPet+1),&
             LIS_nss_obs_halo_ind(n,LIS_localPet+1)/),&
             count =(/LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k)/)))
        call LIS_verify(nf90_close(ftn))

        do r=1,LIS_rc%obs_lnr(k)
           do c=1,LIS_rc%obs_lnc(k)
              if(tmask(c,r).ne.-9999.0) then
                 innovmask(c,r) = 1.0
              endif
           enddo
        enddo
     endif

     currTime = currTime + ts
     if(currTime.gt.stopTime) then
        togo = .false.
     endif
  enddo

  do r =1,LIS_rc%obs_lnr(k)
     do c =1,LIS_rc%obs_lnc(k)
        if (LIS_obs_domain(n,k)%gindex(c,r) .ne. -1)then
           if(innovmask(c,r).eq.1) then 
              obsl(LIS_obs_domain(n,k)%gindex(c,r)) = -9999.0
           end if
        endif
     end do
  end do
  
  
end subroutine applySnowGlobeConditionalDAflags

