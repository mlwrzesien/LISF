!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: LIS_microMetCorrection
! \label{LIS_microMetCorrection}
!
! !REVISION HISTORY:
!
!  21 Dec 2004: Sujay Kumar; Initial Specification
!  21 Jan 2021: Kristi Arsenault; Update to the code, based on SnowModel
!
! !INTERFACE:
subroutine LIS_microMetCorrection(nest, modelelev, LIS_FORC_Base_State)
! !USES:
  use LIS_logMod
  use ESMF
  use LIS_mpiMod
  use LIS_constantsMod
  use LIS_coreMod
  use LIS_FORC_AttributesMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: nest
  real                :: modelelev(LIS_rc%ngrid(nest))
  type(ESMF_State)    :: LIS_FORC_Base_State
!
! !DESCRIPTION:
!  Uses the MicroMet topographic downscaling methods to correct
!   temperature, pressure, humidity, shortwave
!   and longwave radiation, and wind field redistribution,
!   based on:
!
!  Ref:
!  Liston, G. E., & Elder, K. (2006). A Meteorological Distribution System for High-Resolution 
!   Terrestrial Modeling (MicroMet), Journal of Hydrometeorology, 7(2), 217-234. 
!   https://journals.ametsoc.org/view/journals/hydr/7/2/jhm486_1.xml
!
! The arguments are:  
! \begin{description}
!  \item[nest]
!   index of the nest 
!  \item [modelelev]
!   forcing elevation
! \end{description}
!EOP

   integer :: t, k, index, row, col
   real    :: force_tmp, tcforce   ! Input, corrected temperature
   real    :: force_hum, hcforce   ! Input, corrected spec humidity 
   real    :: force_prs, pcforce   ! Input, corrected pressure
   real    :: force_lwd, lcforce   ! Input, corrected LW down rad
   real    :: force_swd, scforce   ! Input, corrected SW down rad
   real    :: force_pcp, pptcforce ! Input, corrected precip
   real    :: force_u,   uwcforce  ! Input, corrected U-wind comp.
   real    :: force_v,   vwcforce  ! Input, corrected U-wind comp.

   real    :: elevdiff, elevdiff2
   real    :: Td,tbar,e,Tdcforce

   integer :: mbefore, mafter, weight, xlat, zone
   real    :: T_lapse_rate, Td_lapse_rate
   real    :: precip_lapse_rate, precip_lapse_rate_m, alfa
   real    :: theta, omegas, thetad, wind, Ww
   real    :: sigma_c, zenith, dec,tau,mu,cosi, psi_dir,psi_dif
   real    :: lhour,czenith
   real    :: ea, epsilona,epsilonb
   real    :: mee, mfe, ee, fe, ratio
   real    :: femiss,emiss
   real    :: elev700, T700, Td700, Rh700, es
   real    :: windspd, winddir
   real    :: deg2rad, rad2deg
   real    :: wslope_max, wind_slope, wslopemax_glb, curve_max, curvemax_glb
   real    :: slope, aspect, curvature
   real    :: dirdiff, xmult
   real    :: deltax, deltay, deltaxy, inc, topo

!    B1=17.502 ==> over water, following Buck (1981) reference
!    B2=22.452 ==> over ice
   real, parameter :: A=611.21, B1=17.502, B2=22.452, C=240.97
   real, parameter :: grav = 9.81, rdry = 287., Sstar=1367.0
   real, parameter :: sigma = 5.67E-8
   real, parameter :: Tf = 273.16
   integer, parameter :: bb=2016

   integer, parameter :: months = 12
   integer lastday(months)
   data lastday/31,28,31,30,31,30,31,31,30,31,30,31/

   ! MicroMet monthly lapse rates (Liston and Kelly, 2006):
   real lapse_rate(months)
   real lapse_rate_nohem(months)
   real lapse_rate_sohem(months)
   data lapse_rate_nohem /0.0044,0.0059,0.0071,0.0078,0.0081,0.0082,&
                          0.0081,0.0081,0.0077,0.0068,0.0055,0.0047/
   data lapse_rate_sohem /0.0081,0.0081,0.0077,0.0068,0.0055,0.0047,&
                          0.0044,0.0059,0.0071,0.0078,0.0081,0.0082/
   ! Atmospheric vapor pressure coeffs are in units of km-1
   !  Following Liston and Kelly (2006)
   real am(months)
   real am_nohem(months)
   real am_sohem(months)
   data am_nohem /0.41,0.42,0.40,0.39,0.38,0.36,&
                  0.33,0.33,0.36,0.37,0.40,0.40/
   data am_sohem /0.33,0.33,0.36,0.37,0.40,0.40,&
                  0.41,0.42,0.40,0.39,0.38,0.36/

   ! The precipitation lapse rate units are in km-1.
   real prec_lapse_rate(months)
   real precip_lapse_rate_nohem(months)
   real precip_lapse_rate_sohem(months)
   data precip_lapse_rate_nohem /0.35,0.35,0.35,0.30,0.25,0.20,&
  &                              0.20,0.20,0.20,0.25,0.30,0.35/
   data precip_lapse_rate_sohem /0.20,0.20,0.20,0.25,0.30,0.35,&
  &                              0.35,0.35,0.35,0.30,0.25,0.20/

   ! Avoid problems of zero (low) winds (for example, turbulence
   !   theory, log wind profile, etc., says that we must have some
   !   wind.  Thus, some equations blow up when the wind speed gets
   !   very small).  This number defines the value that any wind speed
   !   below this gets set to. From snowmodel.par file ...
   real, parameter :: windspd_min = 0.1

! The curvature and wind_slope values range between -0.5 and +0.5.
!   Valid slopewt and curvewt values are between 0 and 1, with
!   values of 0.5 giving approximately equal weight to slope and
!   curvature. Glen Liston suggests that slopewt and curvewt be set such
!   that slopewt + curvewt = 1.0.  This will limit the total
!   wind weight to between 0.5 and 1.5 (but this is not required).
   real, parameter :: slopewt = 0.58
   real, parameter :: curvewt = 0.42
   real :: windwt   ! <-- 2D in MicroMet! and dependent on wind direction

! The curvature is used as part of the wind model.  Define a length
!   scale that the curvature calculation will be performed on.  This
!   has units of meters, and should be approximately one-half the
!   wavelength of the topographic features within the domain.
   real, parameter :: curve_len_scale = 300.0

   ! ESMF metforcing states:
   integer            :: status, ierr
   type(ESMF_Field)   :: tmpField, q2Field,  lwdField, psurfField
   type(ESMF_Field)   :: swdField, pcpField, uField, vField
!    type(ESMF_Field)   :: cpcpField,snowfField
   real, pointer      :: tmp(:),q2(:),lwd(:),psurf(:)
   real, pointer      :: swd(:),pcp(:),uwind(:),vwind(:)

! .............................................................

    ! ESMF calls for metforcing fields of interest
    ! Temperature
    call ESMF_StateGet(LIS_FORC_Base_State,&
         trim(LIS_FORC_Tair%varname(1)),tmpField,&
         rc=status)
    call LIS_verify(status)
 
    call ESMF_FieldGet(tmpField,localDE=0, farrayPtr=tmp,rc=status)
    call LIS_verify(status)

    ! Surface pressure
    call ESMF_StateGet(LIS_FORC_Base_State,&
         trim(LIS_FORC_Psurf%varname(1)),psurfField,&
         rc=status)
    call LIS_verify(status)

    call ESMF_FieldGet(psurfField,localDE=0, farrayPtr=psurf,rc=status)
    call LIS_verify(status)

    ! Spec humidity
    call ESMF_StateGet(LIS_FORC_Base_State,&
         trim(LIS_FORC_Qair%varname(1)),q2Field,&
         rc=status)
    call LIS_verify(status)

    call ESMF_FieldGet(q2Field,localDE=0, farrayPtr=q2,rc=status)
    call LIS_verify(status)

    ! LW down radiation
    call ESMF_StateGet(LIS_FORC_Base_State,&
         trim(LIS_FORC_LWdown%varname(1)),lwdField,&
         rc=status)
    call LIS_verify(status)

    call ESMF_FieldGet(lwdField,localDE=0, farrayPtr=lwd,rc=status)
    call LIS_verify(status)

    ! SW down radiation
    call ESMF_StateGet(LIS_FORC_Base_State,&
         trim(LIS_FORC_SWdown%varname(1)),swdField,&
         rc=status)
    call LIS_verify(status)

    call ESMF_FieldGet(swdField,localDE=0, farrayPtr=swd,rc=status)
    call LIS_verify(status)

    ! Total precip (or rainfall)
    call ESMF_StateGet(LIS_FORC_Base_State,&
         trim(LIS_FORC_Rainf%varname(1)),pcpField,&
         rc=status)
    call LIS_verify(status)

    call ESMF_FieldGet(pcpField,localDE=0, farrayPtr=pcp,rc=status)
    call LIS_verify(status)

    ! U-wind component:
    call ESMF_StateGet(LIS_FORC_Base_State,&
         trim(LIS_FORC_Wind_E%varname(1)),uField,&
         rc=status)
    call LIS_verify(status)

    call ESMF_FieldGet(uField,localDE=0, farrayPtr=uwind,rc=status)
    call LIS_verify(status)

    ! V-wind component:
    call ESMF_StateGet(LIS_FORC_Base_State,&
         trim(LIS_FORC_Wind_N%varname(1)),vField,&
         rc=status)
    call LIS_verify(status)

    call ESMF_FieldGet(vField,localDE=0, farrayPtr=vwind,rc=status)
    call LIS_verify(status)

! .............................................................

!      lat = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lat
!      lon = LIS_domain(n)%grid(LIS_domain(n)%gindex(col, row))%lon

   ! From MicroMet code (get_lapse_rates routine)
   ! Air, dewpoint temperature, and precipitation.
   xlat = LIS_domain(nest)%grid(LIS_domain(nest)%gindex(1,1))%lat
   do k=1,months
      if( xlat.lt.0.0 ) then   ! Southern Hemisphere
        lapse_rate(k) = lapse_rate_sohem(k)              ! Temp
        am(k) = am_sohem(k)                              ! Td
        prec_lapse_rate(k) = precip_lapse_rate_sohem(k)  ! Precip
      else
        lapse_rate(k) = lapse_rate_nohem(k)
        am(k) = am_nohem(k)
        prec_lapse_rate(k) = precip_lapse_rate_nohem(k)
      endif
   enddo

   ! Find the month before and after the day in question.
   if (LIS_rc%da.le.15) then
      mbefore = LIS_rc%mo - 1
      if (mbefore.eq.0) mbefore = 12
      mafter = LIS_rc%mo 
      weight = (real(lastday(mbefore)) - 15. + real(LIS_rc%da)) / &
              &  real(lastday(mbefore))
   else
      mbefore = LIS_rc%mo 
      mafter = LIS_rc%mo + 1
      if (mafter.eq.13) mafter = 1
      weight = (real(LIS_rc%da) - 15.) / real(lastday(mbefore))
   endif

   ! Define the temperature lapse rate (Kelvin/1000m).
   T_lapse_rate = (- (weight * lapse_rate(mafter) + &
     &  (1. - weight) * lapse_rate(mbefore))) 

   ! Define the dew-point temperature lapse rate (deg C/m).
   Td_lapse_rate = (- ((weight * am(mafter) + &
     &  (1. - weight) * am(mbefore)) * C)) / (B1 * 1000.0)  ! MicroMet

   ! Define the precipitation lapse rate (km-1).
   precip_lapse_rate = weight * prec_lapse_rate(mafter) + &
     &  (1. - weight) * prec_lapse_rate(mbefore)

   ! For topographic wind distribution calculation
   deg2rad = LIS_CONST_PI / 180.0
   rad2deg = 180.0 / LIS_CONST_PI

   ! Loop over model-space tiles:
   do t=1,LIS_rc%ntiles(nest)

      if(tmp(t).gt.0) then
        force_tmp = tmp(t)
        force_prs = psurf(t)
        force_hum = q2(t)
        force_pcp = pcp(t)
        force_lwd = lwd(t)
!        force_swd = swd(t)
!        force_u   = uwind(t)
!        force_v   = vwind(t)

        ! Calculate elev diff between hi-res elev and forcing elev
        index = LIS_domain(nest)%tile(t)%index
        elevdiff = LIS_domain(nest)%tile(t)%elev-&
                       modelelev(index)

        ! Perform lapse-rate correction on temperature:
        tcforce = force_tmp+(T_lapse_rate*elevdiff)

        ! Pressure - elevation change (from Cosgrove et. al (2003):
        tbar=(force_tmp+tcforce)/2.
        pcforce=force_prs/(exp((LIS_CONST_G*elevdiff)/(rdry*tbar)))

        ! Perform lapse-rate correction on humidity fields:
        ! Compute water vapor pressure from spec humidity and pressure
        e = force_hum*force_prs / (0.622+0.378*force_hum)   

        Td = C*log((e/A))/(B1-log((e/A))) + Tf      ! matches MicroMet equation
        tdcforce = Td + (Td_lapse_rate*elevdiff)   

        ! Convert back to spec hum from Td (MicroMet)
        e = A*exp((B1*(tdcforce-Tf))/(C+(tdcforce-Tf)))  
        hcforce = (0.622*e) / (pcforce - (0.378*e))

        ! Use a precipitation "lapse rate", or adjust the
        !  the precipitation on the actual elevation grid.  
        !  The reason the interpolated station elevations are 
        !  used as the topographic reference surface (instead 
        !  of something like sea level), is because precip adjustment
        !  precipitation adjustment factor is a non-linear 
        !  function of elevation difference.
        !
        ! The adjustment factor that is used comes from: Thornton, P. E.,
        !   S. W. Running, and M. A. White, 1997: Generating surfaces of
        !   daily meteorological variables over large regions of complex
        !   terrain.  J. Hydrology, 190, 214-251.

        ! Convert the precipitation "lapse rate" (km-1) to m-1.
        precip_lapse_rate_m = precip_lapse_rate / 1000.0

        ! Apply Glen Liston's original precipitation increase with
        !  elevation scheme. Also, we could implement (in the future) 
        !  Ward van Pelt's scheme (see van Pelt et al. 2016).
        ! Don't let the elevation difference be greater than some number
        !  (like 1800 meters gives a factor of 4.4).  If it is too large
        !   you get huge precipitation adjustment, a divide by zero, or
        !   even negative adjustments for high elevations).

        elevdiff2 = min(elevdiff,1800.0)   
        alfa = precip_lapse_rate_m * elevdiff2
        pptcforce = force_pcp * (1.0 + alfa)/(1.0 - alfa)
        pptcforce = max(0.0,pptcforce)   ! Ensure no negative precip values

        ! .............
        ! Longwave radiation -- Elevation change (source??): -- OLD LIS CODE
!        ea = e * 0.01
!        epsilona = 0.70 + (5.95e-5 * ea * exp(1500 / tcforce))
!        epsilonb = -0.792 + (3.161 * epsilona) - (1.573 * epsilona * epsilona)
!        force_lwd = epsilonb*5.67E-8*(tcforce)**4

        ! From Cogsrove et al. (2003) -- For now!
        !  Will consider including the actual MicroMet routines here for LWdown ...
        !  requires air temp and rel humidity ...
        !  Will have to figure out a more "generic" way to handle the forest-canopy
        !   vegetation type values for LAI ... (have to see if the different land
        !   cover maps have enough overlap between the types ...
        ee = (force_hum*force_prs)/0.622
        fe = (hcforce*pcforce)/0.622
        mee = ee/100.
        mfe = fe/100.
       !----------------------------------------------------------------------
       ! correct for negative vapor pressure at very low temperatures at
       ! high latitudes
       !----------------------------------------------------------------------
        if (mee .le. 0) mee = 1e-08
        if (mfe .le. 0) mfe = 1e-08
        emiss  = 1.08*(1-exp(-mee**(force_tmp/bb)))
        femiss = 1.08*(1-exp(-mfe**(tcforce/bb)))
        ratio  = (femiss*(tcforce**4))/(emiss*(force_tmp**4))
        lcforce= force_lwd*ratio

        ! Downward SW radiation ...
        !  Will use the existing "slope-aspect" correction version for now, and
        !  will incorporate the MM routines in the future ...
        !  "solar" routine from MicroMet calculates the "local" SWdown, and
        !   only set up for radiation values for <= 3 hour model timestep
        !   or daily timestep.
        !  Again, need to figure out what to do with the vegtype-dependent forest-canopy 
        !  code and LAI fractions.  


        ! Wind fields ...
        ! WILL NEED TO RETURN TO THIS AFTER ADDING CURVATURE PARAMETER FIELD TO TO LDT
        ! Calculate topo fields for wind and SW down corrections:
        ! Convert domain dx, dy into meters:
        if( LIS_domain(nest)%dx <= 1.) then  ! Temp. solution for lat/lon dx,dy
          deltax = LIS_domain(nest)%dx*100000.
          deltay = LIS_domain(nest)%dy*100000.
        else
          deltax = LIS_domain(nest)%dx*1000.
          deltay = LIS_domain(nest)%dy*1000.
        endif

        ! Assign local slope, aspect, curvature for winds and radiation calculations:
        slope  = LIS_domain(nest)%tile(t)%slope/deg2rad
        aspect = LIS_domain(nest)%tile(t)%aspect/deg2rad
        curvature = 0.1
!        curvature = LIS_domain(nest)%tile(t)%curvature   ! Start providing via LDT??
!        lat = LIS_domain(nest)%grid(index)%lat*LIS_CONST_PI/180.0

        ! Compute the average grid increment.
        deltaxy = 0.5 * (deltax + deltay)

!! WILL RETURN TO DEAL WITH SWrad and WINDS AFTER WE IMPLEMENT CURVATURE IN LDT!
#if 0
        ! Convert the length scale to an appropriate grid increment.
        inc = max(1,nint(curve_len_scale/deltaxy))

        ! Scale the curvature such that the max abs(curvature) has a value
        !  of abs(0.5).  Include a 1 mm curvature in curve_max to prevent
        !  divisions by zero in flat terrain where the curvature is zero.
        curve_max = 0.0 + 0.001
        curve_max = max(curve_max,abs(curvature))

#if (defined SPMD)
        call MPI_Barrier(LIS_MPI_COMM, ierr)
        call MPI_ALLREDUCE(curve_max, curvemax_glb, 1,&
             MPI_REAL, MPI_MAX,&
             LIS_mpi_comm, ierr)
        curve_max = curvemax_glb
#endif
        curvature = curvature / (2.0 * curve_max)

        ! Wind correction from MicroMet:
        !
        ! If desired, impose a wind speed increase with elevation.  Here
        !   the wind_lapse_rate = the wind speed increase per 1-km elevation
        !   gain. The adjustment is patterned after the precipitation-
        !   elevation adjustment.  
        ! DEFAULT VALUE SET IN SNOWMODEL.PAR: wind_lapse_rate == 0
        !if (wind_lapse_rate.ne.0.0) then   ! Brought in as an option from SnowModel
        !  alfa1 = (wind_lapse_rate - 1.0) / (1.0 + wind_lapse_rate)
        ! Convert to m-1.
        !  alfa1 = alfa1 / 1000.0
        !  do j=1,ny
        !   do i=1,nx
        !     delta_topo = topo(i,j) - topo_ref_grid(i,j)
              ! Impose some limits to the adjustment.
        !     delta_topo = min(delta_topo,1800.0)
        !     alfa2 = alfa1 * delta_topo
        !     u_grid(i,j) = u_grid(i,j) * (1.0 + alfa2)/(1.0 - alfa2)
        !     v_grid(i,j) = v_grid(i,j) * (1.0 + alfa2)/(1.0 - alfa2)
        !   enddo
        !  enddo
        !endif

        ! Convert these u and v components to speed and directions.
        ! Some compilers do not allow both u and v to be 0.0 in
        !   the atan2 computation.
       
        if (abs(force_u).lt.1e-10) force_u = 1e-10
        winddir = rad2deg * atan2(force_u,force_v)
        if (winddir.ge.180.0) then
          winddir = winddir - 180.0
        else
          winddir = winddir + 180.0
        endif
        windspd = sqrt(force_u**2 + force_v**2)

        ! Modify the wind speed and direction according to simple
        !  wind-topography relationships. From "topo_mod_winds"  
        !  subroutine in micromet_code.f90

        ! Compute the wind modification factor which is a function of
        !   topography and wind direction following Liston and Sturm (1998).
        !
        ! Compute the slope in the direction of the wind.
        wind_slope = deg2rad * slope * &
     &      cos(deg2rad * (winddir - aspect))

        ! Scale the wind slope such that the max abs(wind slope) has a value
        !   of abs(0.5).  Include a 1 mm slope in slope_max to prevent
        !   divisions by zero in flat terrain where the slope is zero.
        wslope_max = 0.0 + 0.001
        wslope_max = max(wslope_max,abs(wind_slope))

#if (defined SPMD)
        call MPI_Barrier(LIS_MPI_COMM, ierr)
        call MPI_ALLREDUCE(wslope_max, wslopemax_glb, 1,&
                 MPI_REAL, MPI_MAX, LIS_mpi_comm, ierr)
        wslope_max = wslopemax_glb
#endif

        wind_slope = wind_slope / (2.0 * wslope_max)

        ! Calculate the wind speed and direction adjustments.  The
        !  curvature and wind_slope values range between -0.5 and +0.5.
        !  Valid slopewt and curvewt values are between 0 and 1, with
        !  values of 0.5 giving approximately equal weight to slope and
        !  curvature. Glen Liston suggests that slopewt and curvewt be set such
        !  that slopewt + curvewt = 1.0.  This will limit the total
        !  wind weight to between 0.5 and 1.5 (but this is not required).

        ! Compute the wind weighting factor -- dependent on wind direction input.
        windwt = 1.0 + slopewt * wind_slope + curvewt * curvature

        ! Smooth the wind weighting factor to eliminate any sharp speed
        !   increases resulting from the merging of the curve wt and the
        !   slope wt.  Define the number of times this is done to be a
        !   function of the curvature length scale and the grid increment.
        !   The first 0.5 here just means that half of the caclulated
        !   number appears to be about right (additional loops with this
        !   smoother does not change the results much).  If there are
        !   unwanted wave features in the snow distribution, this 0.5
        !   factor can be increased to 1.0 or more, to get rid of these
        !   waves.  Also see "loops_snowd_smoother" in snowtran_code.f90.
        ! xmult = 0.5
        xmult = 1.0     ! Based on MicroMet code
        ! xmult = 1.5

        loops_windwt_smoother = nint(xmult * curve_len_scale / &
     &                          (0.5 * (deltax + deltay)))

        ! NEED TO HAVE IN 2D SPACE ... NC, NR ... for ntiles ...
        !    LIS_domain(nest)%tile(t)%col,LIS_domain(nest)%tile(t)%row
        ! windwt --> is 2D and dependent on winddir
        !
        ! Don't do this smoothing if the domain is arbitrarily small.
!        if (nx.gt.100 .and. ny.gt.100) then
!        if (LIS_rc%lnc(nest).gt.100 .and. LIS_rc%lnr(nest).gt.100) then
!          do k=1,loops_windwt_smoother
!             call smoother9(LIS_rc%lnc(nest),LIS_rc%lnr(nest),windwt)  !<-- windwt(nc,nr)
!          enddo
!        endif

        ! Generate the terrain-modified wind speed.
!        windspd = windwt(nc,nr) * windspd <-- Windwt is 2D and dependent on wind direction ...
        windspd = windwt * windspd

        ! Further modify the wind speed to account for forest canopies.
        !  if (vegtype(i,j).le.5.0) then
        !     nveg = nint(vegtype(i,j))
        ! As with LW and SW radiation calcs -- need to account for veg type
        !  Which can landcover/use map dependent, or even model-dependent ...
        ! May add this option later ...

        ! Modify the wind direction according to Ryan (1977).  Note that it
        !   is critical that "dirdiff" handles the cases where the slope
        !   azimuth and the wind direction are on different sides of the
        !   360-0 line.

        if( aspect.gt.270.0.and. &
     &      winddir.lt.90.0 ) then
          dirdiff = aspect - winddir - 360.0
        elseif( aspect.lt.90.0.and. &
     &          winddir.gt.270.0 ) then
          dirdiff = aspect - winddir + 360.0
        else
          dirdiff = aspect - winddir
        endif
        if (abs(dirdiff).le.90.0) then
          winddir = winddir - 0.5 * &
     &        min(wind_slope*rad2deg,45.0) * &
     &        sin(deg2rad * (2.0 * dirdiff))
          if (winddir.gt.360.0) then
            winddir = winddir - 360.0
          elseif (winddir.lt.0.0) then
            winddir = winddir + 360.0
          endif
        endif

        ! Avoid problems of zero (low) winds (for example, turbulence
        !   theory, log wind profile, etc., says that we must have some
        !   wind. Thus, some equations blow up when the wind speed gets
        !   very small).
        if (windspd.lt.windspd_min) then
          windspd = windspd_min
        endif
        uwcforce = (- windspd) * sin(deg2rad*winddir)
        vwcforce = (- windspd) * cos(deg2rad*winddir)

#endif

        ! ------------------------------------
        ! Convert back to LIS-model tile space:
        tmp(t)   = tcforce
        psurf(t) = pcforce
        q2(t)    = hcforce
        pcp(t)   = pptcforce
        lwd(t)   = lcforce
!        swd(t)   = scforce
!        uwind(t) = uwcforce
!        vwind(t) = vwcforce

      endif
   enddo

   write(LIS_logunit,*) '[INFO] MicroMet forcing topographic downscaling completed'

end subroutine LIS_microMetCorrection
