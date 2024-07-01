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
! !ROUTINE: read_AMSRcnnSnow
! \label{read_AMSRcnnSnow}
!
! !REVISION HISTORY:
!  08 Jun 2022: Sujay Kumar; Initial version
!
! !INTERFACE: 
subroutine read_AMSRcnnSnow(n, k, OBS_State, OBS_Pert_State)
! !USES: 
  use ESMF
  use LIS_mpiMod
  use LIS_coreMod
  use LIS_logMod
  use LIS_timeMgrMod
  use LIS_dataAssimMod
  use map_utils
  use LIS_pluginIndices
  use LIS_DAobservationsMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use AMSRcnnSnowMod, only : AMSRcnnSnow_struc
#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
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
!  reads the AMSRcnnSnow observations and prepares them for DA
!  
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest
!  \item[OBS\_State] observations state
!  \end{description}
!
!EOP
  integer                :: ftn,status
  character(len=LIS_CONST_PATH_LEN) :: sndobsdir
  character(len=LIS_CONST_PATH_LEN) :: fname
  logical                :: alarmCheck
  integer                :: t,c,r,i,j,p,jj
  integer                :: sndId,ierr
  real,          pointer :: obsl(:)
  type(ESMF_Field)       :: sndfield, pertField
  logical*1              :: sndobs_b_ip(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  integer                :: gid(LIS_rc%obs_ngrid(k))
  integer                :: assimflag(LIS_rc%obs_ngrid(k))
  logical                :: data_update
  logical                :: data_upd_flag(LIS_npes)
  logical                :: data_upd_flag_local
  logical                :: data_upd
  logical*1, allocatable :: snd_data_b(:)
  real                   :: sndobs(LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k))
  real, allocatable      :: var(:,:)
  real, allocatable      :: snd1d(:)
  real                   :: snd_current(LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k))
  integer                :: fnd
  real                   :: cornerlat1, cornerlon1
  integer                :: lat_off, lon_off
  logical                :: file_exists
    
  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       sndobsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(OBS_State,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  sndobs = LIS_rc%udef
  data_upd = .false. 
!-------------------------------------------------------------------------
!   Read both ascending and descending passes at 0Z and then store
!   the overpass time as 1.30AM for the descending pass and 1.30PM 
!   for the ascending pass. 
!-------------------------------------------------------------------------
  alarmCheck = LIS_isAlarmRinging(LIS_rc, "AMSRcnnSnow read alarm")
  
  if(alarmCheck.or.AMSRcnnSnow_struc(n)%startMode) then
     AMSRcnnSnow_struc(n)%startMode = .false.
     sndobs = LIS_rc%udef
     
     call create_AMSRcnnSnow_filename(sndobsdir, &
          LIS_rc%yr, LIS_rc%mo, &
          LIS_rc%da, fname)
     
     inquire(file=fname,exist=file_exists)

     if(file_exists.and.AMSRcnnSnow_struc(n)%validDomain) then
        
        write(LIS_logunit,*) '[INFO] Reading ',trim(fname)

        
        allocate(snd_data_b(AMSRcnnSnow_struc(n)%nc*AMSRcnnSnow_struc(n)%nr))
        allocate(var(AMSRcnnSnow_struc(n)%nc,AMSRcnnSnow_struc(n)%nr))
        allocate(snd1d(AMSRcnnSnow_struc(n)%nc*AMSRcnnSnow_struc(n)%nr))

        if(AMSRcnnSnow_struc(n)%datares.eq.10) then         
#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
           ierr = nf90_open(path=fname,mode=NF90_NOWRITE,ncid=ftn)
           call LIS_verify(ierr,'error opening AMSR CNN file')
           
           ierr = nf90_inq_varid(ftn,'SD',sndId)
           call LIS_verify(ierr, 'nf90_inq_varid failed for SDin read_AMSRcnnSnow')
           
           ierr = nf90_get_var(ftn,sndId, var)
           call LIS_verify(ierr, 'nf90_get_var failed for SD')
           
           ierr = nf90_close(ftn)
           call LIS_verify(ierr)
#endif
           
           snd1d = LIS_rc%udef
           snd_data_b = .false.
           
           do r=1, AMSRcnnSnow_struc(n)%nr
              do c=1, AMSRcnnSnow_struc(n)%nc
                 if(var(c,r).lt.0.or.var(c,r).gt.1000) then 
                    snd1d(c+(AMSRcnnSnow_struc(n)%nr-r)*&
                         AMSRcnnSnow_struc(n)%nc) = -9999.0
                 else
                    snd1d(c+(AMSRcnnSnow_struc(n)%nr-r)*&
                         AMSRcnnSnow_struc(n)%nc) = var(c,r)/100.0 !to meters
                    snd_data_b(c+(AMSRcnnSnow_struc(n)%nr-r)*&
                         AMSRcnnSnow_struc(n)%nc) = .true. 
                 endif
              enddo
           enddo
           deallocate(var)
        elseif(AMSRcnnSnow_struc(n)%datares.eq.5) then
#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )

           cornerlat1=AMSRcnnSnow_struc(n)%gridDesci(7)
           cornerlon1=AMSRcnnSnow_struc(n)%gridDesci(5)           
           
           lat_off = nint((cornerlat1+0.0)/0.05)+1
           lon_off = nint((cornerlon1+179.95)/0.05)+1

           lat_off = 1661 - lat_off + 1
           ierr = nf90_open(path=fname,mode=NF90_NOWRITE,ncid=ftn)
           call LIS_verify(ierr,'error opening AMSR CNN file')
           
           ierr = nf90_inq_varid(ftn,'SD',sndId)
           call LIS_verify(ierr, 'nf90_inq_varid failed for SDin read_AMSRcnnSnow')
           
           ierr = nf90_get_var(ftn,sndId, var,&
                start=(/lon_off,lat_off/),&
                count=(/AMSRcnnSnow_struc(n)%nc,AMSRcnnSnow_struc(n)%nr/))
           call LIS_verify(ierr, 'nf90_get_var failed for SD')
           
           ierr = nf90_close(ftn)
           call LIS_verify(ierr)
           
#endif
           snd1d = LIS_rc%udef
           snd_data_b = .false.
           
           do r=1, AMSRcnnSnow_struc(n)%nr
              do c=1, AMSRcnnSnow_struc(n)%nc
                 if(var(c,r).lt.0.or.var(c,r).gt.1000) then 
                    snd1d(c+(AMSRcnnSnow_struc(n)%nr-r)*&
                         AMSRcnnSnow_struc(n)%nc) = -9999.0
                 else
                    snd1d(c+(AMSRcnnSnow_struc(n)%nr-r)*&
                         AMSRcnnSnow_struc(n)%nc) = var(c,r)/100.0 !to meters
                    snd_data_b(c+(AMSRcnnSnow_struc(n)%nr-r)*&
                         AMSRcnnSnow_struc(n)%nc) = .true. 
                 endif
              enddo
           enddo
           deallocate(var)
        endif
        
!--------------------------------------------------------------------------
! Interpolate to the observation grid
!-------------------------------------------------------------------------- 
        call bilinear_interp(LIS_rc%obs_gridDesc(k,:),&
             snd_data_b,snd1d,&
             sndobs_b_ip,sndobs,&
             AMSRcnnSnow_struc(n)%nc*AMSRcnnSnow_struc(n)%nr,&
             LIS_rc%obs_lnc(k)*LIS_rc%obs_lnr(k),&
             AMSRcnnSnow_struc(n)%rlat, &
             AMSRcnnSnow_struc(n)%rlon, &
             AMSRcnnSnow_struc(n)%w11,&
             AMSRcnnSnow_struc(n)%w12,&
             AMSRcnnSnow_struc(n)%w21,&
             AMSRcnnSnow_struc(n)%w22,&
             AMSRcnnSnow_struc(n)%n11,&                
             AMSRcnnSnow_struc(n)%n12,&
             AMSRcnnSnow_struc(n)%n21,&                
             AMSRcnnSnow_struc(n)%n22,&
             LIS_rc%udef,status)


!        print*, 'here'
!        open(100,file='test.bin',form='unformatted')
!        write(100) sndobs
!        close(100)
!        stop
           
        deallocate(snd1d)
        deallocate(snd_data_b)

     endif

     call ESMF_StateGet(OBS_State,"Observation01",sndfield,&
             rc=status)
     call LIS_verify(status, 'Error: StateGet Observation01')
     
     call ESMF_FieldGet(sndfield,localDE=0,farrayPtr=obsl,rc=status)
     call LIS_verify(status, 'Error: FieldGet')
     
     fnd = 0 


     do r=1,LIS_rc%obs_lnr(k)
        do c=1,LIS_rc%obs_lnc(k)
           if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
              
              if(sndobs(c+(r-1)*LIS_rc%obs_lnc(k)).ge.0) then 
                 fnd = 1
                 exit
              endif
           endif
        enddo
     enddo

     obsl = LIS_rc%udef 
     if(fnd.ne.0) then 
        do r=1, LIS_rc%obs_lnr(k)
           do c=1, LIS_rc%obs_lnc(k)
              if(LIS_obs_domain(n,k)%gindex(c,r).ne.-1) then 
                 obsl(LIS_obs_domain(n,k)%gindex(c,r))=&
                      sndobs(c+(r-1)*LIS_rc%obs_lnc(k))
              endif
           enddo
        enddo
     endif
     
!-------------------------------------------------------------------------
!  Apply LSM based quality control and screening of observations
!-------------------------------------------------------------------------
 
     call lsmdaqcobsstate(trim(LIS_rc%lsm)//"+"&
          //trim(LIS_AMSRcnnSnowobsId)//char(0),n, k, OBS_state)
     
     snd_current = LIS_rc%udef
     call LIS_checkForValidObs(n,k,obsl,fnd,sndobs)

!     if(LIS_localPet.eq.3) then 
!        print*, 'here', LIS_rc%obs_lnc(k),LIS_rc%obs_lnr(k)
!        open(100,file='test.bin',form='unformatted')
!        write(100) sndobs
!        close(100)
!        stop
!     endif
     
     if(fnd.eq.0) then 
        data_upd_flag_local = .false. 
     else
        data_upd_flag_local = .true. 
     endif
     
#if (defined SPMD)
     call MPI_ALLGATHER(data_upd_flag_local,1, &
          MPI_LOGICAL, data_upd_flag(:),&
          1, MPI_LOGICAL, LIS_mpi_comm, status)
    data_upd = any(data_upd_flag)
#else
    data_upd = data_upd_flag_local
#endif
     
     if(data_upd) then
        do t=1,LIS_rc%obs_ngrid(k)
           gid(t) = t
           if(obsl(t).ne.-9999.0) then 
              assimflag(t) = 1
              if(obsl(t).gt.10.0) then
                  !print *,'WARNING: OBSL > 10.0m; OBSL=',obsl(t)
                  assimflag(t) = 0 !Do not assimilate values > 10m
              endif
           else
              assimflag(t) = 0
           endif
        enddo
        
        call ESMF_AttributeSet(OBS_State,"Data Update Status",&
             .true. , rc=status)
        call LIS_verify(status)
        
        if(LIS_rc%obs_ngrid(k).gt.0) then 
           call ESMF_AttributeSet(sndField,"Grid Number",&
                gid,itemCount=LIS_rc%obs_ngrid(k),rc=status)
           call LIS_verify(status)
           
           call ESMF_AttributeSet(sndField,"Assimilation Flag",&
                assimflag,itemCount=LIS_rc%obs_ngrid(k),rc=status)
           call LIS_verify(status)
           
        endif

     else
        call ESMF_AttributeSet(OBS_State,"Data Update Status",&
             .false., rc=status)
        call LIS_verify(status)     
     endif
  else
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status)     
  endif
end subroutine read_AMSRcnnSnow


!BOP
! !ROUTINE: create_AMSRcnnSnow_filename
! \label{create_AMSRcnnSnow_filename}
! 
! !INTERFACE: 
subroutine create_AMSRcnnSnow_filename(ndir, yr, mo,da, filename)
! !USES:   

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  integer           :: yr, mo, da
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates the AMSRcnnSnow data filename
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the AMSRcnnSnow directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
!  \item[filename] Generated AMSRcnnSnow filename
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
 
  filename = trim(ndir)//'/'//fyr//&
       '/AMSR_CNN_NH_'//trim(fyr)//'-'//&
       trim(fmo)//'-'//trim(fda)//'.nc'
  
end subroutine create_AMSRcnnSnow_filename


