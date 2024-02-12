module cb24cnn

!=====================================================================
!
! Purpose: Implement cb24cnn of the model state of U,V,T,Q, and/or PS
!          toward specified values from random analyses increments.
!
! Author: Wiliiam Chapman -- wchapman@ucar.edu
!
! Description:
!
!    This module assumes that the user has {U,V,T,Q,PS} values from analyses
!    which have been preprocessed onto the current model grid and adjusted
!    for differences in topography. It is also assumed that these resulting
!    values and are stored in individual files which are indexed with respect
!    to year, month, day, and second of the day. When the model is inbetween
!    the given begining and ending times, a relaxation forcing is added to
!    cb24cnn the model toward the analyses values determined from the forcing
!    option specified. After the model passes the ending analyses time, the
!    forcing discontinues.
!
!    Some beans analyses products can have gaps in the available data, where values
!    are missing for some interval of time. When files are missing, the cb24cnn
!    force is switched off for that interval of time, so we effectively 'coast'
!    thru the gap.
!
!    Currently, the cb24cnn module is set up to accomodate cb24cnn of PS
!    values, however that functionality requires forcing that is applied in
!    the selected dycore and is not yet implemented.
!
!    The cb24cnn of the model toward the analyses data is controlled by
!    the 'cb24cnn_nl' namelist in 'user_nl_cam'; whose variables control the
!    time interval over which cb24cnn is applied, the strength of the cb24cnn
!    tendencies, and its spatial distribution.
!

! Useful modules
  !------------------
  use shr_kind_mod,   only:r8=>SHR_KIND_R8,cs=>SHR_KIND_CS,cl=>SHR_KIND_CL
  use time_manager,   only:timemgr_time_ge,timemgr_time_inc,timemgr_time_inc_minus,get_curr_date,get_curr_calday,get_step_size
  use, intrinsic :: ieee_exceptions, only: ieee_get_halting_mode, ieee_set_halting_mode, ieee_all, ieee_support_halting, ieee_overflow
  use phys_grid   ,   only:scatter_field_to_chunk, gather_chunk_to_field
  use cam_abortutils, only:endrun
  use spmd_utils  ,   only:masterproc, mstrid=>masterprocid, mpicom
  use cam_logfile ,   only:iulog
#ifdef SPMD
  use mpishorthand
#endif

  use forpy_mod,                 only : module_py,list,ndarray,object,tuple
  use forpy_mod,                 only : err_print
  use forpy_mod,                 only : forpy_initialize,get_sys_path,import_py,print_py
  use forpy_mod,                 only : ndarray_create,tuple_create,call_py,cast
  use forpy_mod,                 only : forpy_finalize, dict, dict_create

! Set all Global values and routines to private by default
  ! and then explicitly set their exposure.
  !----------------------------------------------------------
  implicit none
  private

  public:: cb24cnn_Model!,cb24cnn_Model_Regress,cb24cnn_ON
  !public:: cb24cnn_readnl
  public:: cb24cnn_init
  public:: cb24cnn_timestep_init
  public:: cb24cnn_timestep_tend
  public:: cb24cnn_set_tend
  public:: cb24cnn_finalize


  !public:: regress_daml_timestep_tend
  !public:: regress_diurnal_daml_timestep_tend
  !public:: cb24cnn_diurnal_timestep_tend

  !private::cb24cnn_readin_weights_all
  !private::cb24cnn_readin_weights
  !private::cb24cnn_readin_weights_diurnal
  !private::cb24cnn_readin_normdict
  !private::cb24cnn_set_PSprofile
  !private::cb24cnn_set_profile


  !  cb24cnn Parameters
  !--------------------
  logical          :: cb24cnn_Model       =.true.
  logical          :: halting_mode(5)
  type(module_py)  :: pymodule



  integer :: &
    upwp_clubb_idx, &          ! Cloud fraction
    taux_idx,       &
    pblh_idx,       &
    ustarwec_idx,   &
    flntc_idx,      &
    fsntoa_idx


  ! cb24cnn State Arrays
  integer cb24cnn_nlon,cb24cnn_nlat,cb24cnn_ncol,cb24cnn_nlev,cb24cnn_nlevp
  real(r8),allocatable:: cb24cnn_Ustep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: cb24cnn_Vstep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: cb24cnn_Sstep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: cb24cnn_Qstep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: cb24cnn_PSstep(:,:)    !(pcols,begchunk:endchunk)



  !CNN input vars
  real(r8),allocatable:: cb24cnn_Model_V (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: cb24cnn_Model_T (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: cb24cnn_Model_U (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: cb24cnn_Model_Q (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: cb24cnn_Model_W (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: cb24cnn_Model_UPWP (:,:,:)  !(pcols,pverp,begchunk:endchunk)
  real(r8),allocatable:: cb24cnn_Model_TAUX (:,:)    !(pcols,begchunk:endchunk)
  real(r8),allocatable:: cb24cnn_Model_TAUY (:,:)    !(pcols,begchunk:endchunk)
  real(r8),allocatable:: cb24cnn_Model_USTAR (:,:)   !(pcols,begchunk:endchunk)
  real(r8),allocatable:: cb24cnn_Model_TBOT (:,:) !(pcols,begchunk:endchunk)
  real(r8),allocatable:: cb24cnn_Model_PBLH (:,:) !(pcols,begchunk:endchunk)

  !derived type interfaces:

  !> Control structure for Python interface
  type, public :: python_interface ; private
    type(module_py) :: pymodule
    type(list) :: paths
  end type

contains

   !================================================================
  subroutine cb24cnn_init
    !
    ! cb24cnn_INIT: Allocate space and initialize cb24cnn values
    !===============================================================
    use ppgrid        ,only: pver,pverp,pcols,begchunk,endchunk
    use error_messages,only: alloc_err
    use dycore        ,only: dycore_is
    use dyn_grid      ,only: get_horiz_grid_dim_d
    use phys_grid     ,only: get_rlat_p,get_rlon_p,get_ncols_p
    use cam_history   ,only: addfld
    use shr_const_mod ,only: SHR_CONST_PI
    use filenames     ,only: interpret_filename_spec

    implicit none

    integer :: ierror
    type(list) :: my_list
    type(list) :: paths

    !local Values
    !-------------
    integer  istat,lchnk
    integer  hdim1_d,hdim2_d
    character(len=:), allocatable :: return_string
    
    ! Allocate Space for cb24 data arrays
    !-----------------------------------------

    ! Allocate Space for spatial dependence of
    ! cb24cnn Coefs and cb24cnn Forcing.
    !-------------------------------------------
    allocate(cb24cnn_Ustep(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24cnn_init','cb24cnn_Ustep',pcols*pver*((endchunk-begchunk)+1))
    allocate(cb24cnn_Vstep(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24cnn_init','cb24cnn_Vstep',pcols*pver*((endchunk-begchunk)+1))
    allocate(cb24cnn_Sstep(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24cnn_init','cb24cnn_Sstep',pcols*pver*((endchunk-begchunk)+1))
    allocate(cb24cnn_Qstep(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24cnn_init','cb24cnn_Qstep',pcols*pver*((endchunk-begchunk)+1))
    allocate(cb24cnn_PSstep(pcols,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24cnn_init','cb24cnn_PSstep',pcols*((endchunk-begchunk)+1)) !currently not used.

    allocate(cb24cnn_Model_U(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24cnn_init','cb24cnn_Model_U',pcols*pver*((endchunk-begchunk)+1))
    allocate(cb24cnn_Model_V(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24cnn_init','cb24cnn_Model_V',pcols*pver*((endchunk-begchunk)+1))
    allocate(cb24cnn_Model_T(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24cnn_init','cb24cnn_Model_T',pcols*pver*((endchunk-begchunk)+1))
    allocate(cb24cnn_Model_Q(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24cnn_init','cb24cnn_Model_Q',pcols*pver*((endchunk-begchunk)+1))
    allocate(cb24cnn_Model_W(pcols,pverp,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24cnn_init','cb24cnn_Model_W',pcols*pverp*((endchunk-begchunk)+1))
    allocate(cb24cnn_Model_UPWP(pcols,pverp,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24cnn_init','cb24cnn_Model_Q',pcols*pverp*((endchunk-begchunk)+1))

    allocate(cb24cnn_Model_TAUX(pcols,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24cnn_init','cb24cnn_Model_TAUX',pcols*((endchunk-begchunk)+1)) !currently not used.
    allocate(cb24cnn_Model_TAUY(pcols,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24cnn_init','cb24cnn_Model_TAUY',pcols*((endchunk-begchunk)+1)) !currently not used.
    allocate(cb24cnn_Model_USTAR(pcols,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24cnn_init','cb24cnn_Model_USTAR',pcols*((endchunk-begchunk)+1)) !currently not used.
    allocate(cb24cnn_Model_TBOT(pcols,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24cnn_init','cb24cnn_Model_TBOT',pcols*((endchunk-begchunk)+1)) !currently not used.
    allocate(cb24cnn_Model_PBLH(pcols,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24cnn_init','cb24cnn_Model_PBLH',pcols*((endchunk-begchunk)+1)) !currently not used.

    ! Register output fields with the cam history module
     !-----------------------------------------------------
    call addfld( 'cb24cnn_U',(/ 'lev' /),'A','m/s/s'  ,'U DAMLining Tendency')
    call addfld( 'cb24cnn_V',(/ 'lev' /),'A','m/s/s'  ,'V DAMLining Tendency')
    call addfld( 'cb24cnn_T',(/ 'lev' /),'A','K/s'    ,'T DAMLining Tendency')
    call addfld( 'cb24cnn_Q',(/ 'lev' /),'A','kg/kg/s','Q DAMLining Tendency')

    ! Initialize column and level dimensions
     !--------------------------------------------------------
     call get_horiz_grid_dim_d(hdim1_d,hdim2_d)
     cb24cnn_nlon=hdim1_d
     cb24cnn_nlat=hdim2_d
     cb24cnn_ncol=hdim1_d*hdim2_d
     cb24cnn_nlev=pver
     cb24cnn_nlevp=pverp

    ! Initialize the forpy module and import the python script "DAMLcnn.py"
     !-----------------------------------------------------
    if (masterproc) write(iulog,*) 'cb24cnn_init starting'
    if (.not. ieee_support_halting(ieee_overflow)) then
       call endrun("ieee_halting is not supported")
    endif
    call ieee_get_halting_mode(ieee_all, halting_mode)
    print *,__FILE__,__LINE__,halting_mode
    call ieee_set_halting_mode(ieee_all, .false.)

    ierror = forpy_initialize()
    ierror = get_sys_path(paths)
    ierror = paths%append(".") !stash it in the run directory; future NL

    ierror = import_py(pymodule,"DAMLcnn") !future NL
    call ieee_set_halting_mode(ieee_all, halting_mode)
    call paths%destroy

    ! Init forcing as zeros:
    do lchnk=begchunk,endchunk
      cb24cnn_Ustep(:pcols,:pver,lchnk)=0._r8
      cb24cnn_Vstep(:pcols,:pver,lchnk)=0._r8
      cb24cnn_Sstep(:pcols,:pver,lchnk)=0._r8
      cb24cnn_Qstep(:pcols,:pver,lchnk)=0._r8
      cb24cnn_PSstep(:pcols,lchnk)=0._r8
    end do

  end subroutine cb24cnn_init



  subroutine cb24cnn_timestep_init(phys_state,phys_tend,pbuf,cam_in,cam_out)

    !
    ! cb24cnn_TIMESTEP_INIT:
    !                 Check the current time and update Model/Nudging
    !                 arrays when necessary. Toggle the Nudging flag
    !                 when the time is withing the nudging window.
    !===============================================================
    use physconst    ,only: cpair
    use physics_types,only: physics_state,physics_ptend,physics_ptend_init, &
                            physics_state_copy
    use constituents ,only: cnst_get_ind
    use dycore       ,only: dycore_is
    use ppgrid       ,only: pver,pverp,pcols,begchunk,endchunk
    use filenames    ,only: interpret_filename_spec
    use ESMF
    use camsrfexch     ,only: cam_in_t,cam_out_t
    use physics_buffer ,only: pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field, &
                              physics_buffer_desc
    use cam_history  ,only: outfld

    ! Arguments
    !-------------
    type(physics_state), intent(in) :: phys_state
    type(cam_in_t),      intent(in)    :: cam_in
    type(cam_out_t),     intent(in)    :: cam_out
    type(physics_ptend), intent(out):: phys_tend
    type(physics_buffer_desc), pointer :: pbuf(:)
    
    ! Local values
    !----------------
    integer Year,Month,Day,Sec
    integer YMD1,YMD2,YMD,YMD3,YMD4,YMD5,YMD6
    logical Update_Model,Update_Nudge,Sync_Error
    logical After_Beg   ,Before_End
    integer lchnk,ncol,indw,itim_old

    real(r8), pointer, dimension(:,:) :: upwp      ! upwp                               [no idea]
    real(r8), pointer, dimension(:) :: pblh     ! planetary boundary layer height       [m]
    real(r8), pointer, dimension(:) :: taux     ! Surface zonal momentum flux           [m]
    real(r8), pointer, dimension(:) :: tauy     ! Surface meridional momentum flux
    real(r8), pointer, dimension(:) :: ustarwec ! Surface meridional momentum flux

    !Grab Chunk and Columns
    !----------------------
    lchnk = phys_state%lchnk
    ncol  = phys_state%ncol

    itim_old = pbuf_old_tim_idx()

    upwp_clubb_idx     = pbuf_get_index('UPWP')         ! UPWP
    call pbuf_get_field(pbuf, upwp_clubb_idx,    upwp,    start=(/1,1,itim_old/), kount=(/pcols,pverp,1/))

    ustarwec_idx     = pbuf_get_index('ustarwec')         ! Ustar !16|33
    call pbuf_get_field(pbuf, ustarwec_idx,    ustarwec)

    pblh_idx     = pbuf_get_index('pblh')         ! PBL height
    call pbuf_get_field(pbuf, pblh_idx,    pblh)

    !cam_out%tbot
    !cam_in%wsx
    !cam_in%wsy

    ! Load values at Current into the Model arrays Predictors...
    !-----------------------------------------------
    ! look here: you need to update this:
    !/glade/work/wchapman/cesm/perfect_nudging/f.e21.FHIST.f09_f09_mg17_nudge_it_ERA5_1999/SourceMods/src.cam
    !nudging_calc_tend()

    call cnst_get_ind('Q',indw)

    cb24cnn_Model_V(:ncol,:pver,lchnk)=phys_state%v(:ncol,:pver)
    cb24cnn_Model_T(:ncol,:pver,lchnk)=phys_state%t(:ncol,:pver)
    cb24cnn_Model_U(:ncol,:pver,lchnk)=phys_state%u(:ncol,:pver)
    cb24cnn_Model_Q(:ncol,:pver,lchnk)=phys_state%q(:ncol,:pver,indw)
    cb24cnn_Model_W(:ncol,:pver,lchnk)=phys_state%omega(:ncol,:pver)
    cb24cnn_Model_UPWP(:ncol,:pverp,lchnk)=upwp(:ncol,:pverp)
    cb24cnn_Model_TAUX(:ncol,lchnk)=cam_in%wsx(:ncol)
    cb24cnn_Model_TAUY(:ncol,lchnk)=cam_in%wsy(:ncol)
    cb24cnn_Model_TBOT(:ncol,lchnk)=cam_out%tbot(:ncol)
    cb24cnn_Model_PBLH(:ncol,lchnk)=pblh(:ncol)
    cb24cnn_Model_USTAR(:ncol,lchnk)=ustarwec(:ncol)


  end subroutine cb24cnn_timestep_init



  subroutine cb24cnn_timestep_tend()

    use physics_types ,only: physics_state,physics_ptend,physics_ptend_init,physics_state_copy
    use physics_buffer,only: pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field, &
                              physics_buffer_desc
    use time_manager  ,only:timemgr_time_ge,timemgr_time_inc,timemgr_time_inc_minus,get_curr_date,get_curr_calday,get_step_size
    use ppgrid        ,only: pver,pverp,pcols,begchunk,endchunk
    use cam_history   ,only: outfld
    use netcdf

    integer :: ierror
    type(tuple) :: args
    type(object) :: return_value

    ! Arguments
    !-------------
    
    !local values:
    !---------------
    integer lev
    integer c                               
    integer nlon,nlat,plev,istat,lchnk,indw
    integer ncid,varid
    integer ilat,ilon,ilev
    integer londimid, latdimid, levdimid,vardimid, rhvarid
    integer  Year,Month,Day,Sec,Hour,calday
    real(r8) Xcnn(cb24cnn_nlon,cb24cnn_nlat,11)  !future NL
    real(r8) Vanal(cb24cnn_nlon,cb24cnn_nlat,cb24cnn_nlev)
    real(r8) Uanal(cb24cnn_nlon,cb24cnn_nlat,cb24cnn_nlev)
    real(r8) Tanal(cb24cnn_nlon,cb24cnn_nlat,cb24cnn_nlev)
    real(r8) Qanal(cb24cnn_nlon,cb24cnn_nlat,cb24cnn_nlev)
    real(r8) Wanal(cb24cnn_nlon,cb24cnn_nlat,cb24cnn_nlev)
    real(r8) UPWPanal(cb24cnn_nlon,cb24cnn_nlat,cb24cnn_nlevp)
    real(r8) TAUXanal(cb24cnn_nlon,cb24cnn_nlat)
    real(r8) TAUYanal(cb24cnn_nlon,cb24cnn_nlat)
    real(r8) USTARanal(cb24cnn_nlon,cb24cnn_nlat)
    real(r8) TBOTanal(cb24cnn_nlon,cb24cnn_nlat)
    real(r8) PBLHanal(cb24cnn_nlon,cb24cnn_nlat)
    real(r8) Lat_anal(cb24cnn_nlat)
    real(r8) Lon_anal(cb24cnn_nlon)
    real(r8) Xtrans(cb24cnn_nlon,cb24cnn_nlev,cb24cnn_nlat)
    real(r8) Xtransp(cb24cnn_nlon,cb24cnn_nlevp,cb24cnn_nlat)
    real(r8) Xtransf(1,cb24cnn_nlon,cb24cnn_nlat)
    integer  nn,Nindex
    real(r8) PI

    type(ndarray) :: in1_py,in2_py,in3_py,in4_py,in5_py,in6_py,in7_py        !< variables in the form of numpy array
    type(ndarray) :: in8_py,in9_py,in10_py,in11_py,in12_py,in13_py,out_arr   !< variables in the form of numpy array
    character(len=:), allocatable :: return_string !< outputs from Python module
    real(r8), pointer, dimension(:,:,:,:) :: out_for !< outputs from Python module
    
    if (masterproc) write(iulog,*) 'cb24cnn_timestep_tend starting....'

    PI    = 4.D0*DATAN(1.D0)
    !Grab Chunk and Columns
    !----------------------
    
    !---
    !V
    !---
    call gather_chunk_to_field(1,cb24cnn_nlev,1,cb24cnn_nlon,cb24cnn_Model_V,Xtrans)
    if (masterproc) then
        do ilat=1,cb24cnn_nlat
        do ilev=1,cb24cnn_nlev
        do ilon=1,cb24cnn_nlon
           Vanal(ilon,ilat,ilev)=Xtrans(ilon,ilev,ilat)
        end do
        end do
        end do
    end if 

    !---
    !U
    !---
    call gather_chunk_to_field(1,cb24cnn_nlev,1,cb24cnn_nlon,cb24cnn_Model_U,Xtrans)
    if (masterproc) then
        do ilat=1,cb24cnn_nlat
        do ilev=1,cb24cnn_nlev
        do ilon=1,cb24cnn_nlon
           Uanal(ilon,ilat,ilev)=Xtrans(ilon,ilev,ilat)
        end do
        end do
        end do
    end if 

    !---
    !T
    !---
    call gather_chunk_to_field(1,cb24cnn_nlev,1,cb24cnn_nlon,cb24cnn_Model_T,Xtrans)
    if (masterproc) then
        do ilat=1,cb24cnn_nlat
        do ilev=1,cb24cnn_nlev
        do ilon=1,cb24cnn_nlon
           Tanal(ilon,ilat,ilev)=Xtrans(ilon,ilev,ilat)
        end do
        end do
        end do
    end if 
    
    !---
    !Q
    !---
    call gather_chunk_to_field(1,cb24cnn_nlev,1,cb24cnn_nlon,cb24cnn_Model_Q,Xtrans)
    if (masterproc) then
        do ilat=1,cb24cnn_nlat
        do ilev=1,cb24cnn_nlev
        do ilon=1,cb24cnn_nlon
           Qanal(ilon,ilat,ilev)=Xtrans(ilon,ilev,ilat)
        end do
        end do
        end do
    end if 

    !---
    !OMEGA
    !---
    call gather_chunk_to_field(1,cb24cnn_nlev,1,cb24cnn_nlon,cb24cnn_Model_W,Xtrans)
    if (masterproc) then
        do ilat=1,cb24cnn_nlat
        do ilev=1,cb24cnn_nlev
        do ilon=1,cb24cnn_nlon
           Wanal(ilon,ilat,ilev)=Xtrans(ilon,ilev,ilat)
        end do
        end do
        end do
    end if 
    
    !---
    !UPWP
    !---
    call gather_chunk_to_field(1,cb24cnn_nlevp,1,cb24cnn_nlon,cb24cnn_Model_UPWP,Xtransp)
    if (masterproc) then
        do ilat=1,cb24cnn_nlat
        do ilev=1,cb24cnn_nlevp
        do ilon=1,cb24cnn_nlon
           UPWPanal(ilon,ilat,ilev)=Xtransp(ilon,ilev,ilat)
        end do
        end do
        end do
    end if 

    !---
    !TAUX
    !---
    call gather_chunk_to_field(1,1,1,cb24cnn_nlon,cb24cnn_Model_TAUX,Xtransf)
    if (masterproc) then
        do ilat=1,cb24cnn_nlat
        do ilon=1,cb24cnn_nlon
           TAUXanal(ilon,ilat)=Xtransf(1,ilon,ilat)
        end do
        end do
        !if you want to see if you've done good kid. 
        !istat=nf90_create('CESM_dumpTAUX.nc',NF90_CLOBBER,ncid)
        !istat=nf90_def_dim(ncid, "lon",cb24cnn_nlon , londimid)
        !istat=nf90_def_dim(ncid, "lat",cb24cnn_nlat , latdimid)
        !istat=nf90_def_var(ncid, "TAUX", nf90_double, (/ londimid, latdimid /), rhvarid)
        !istat=nf90_enddef(ncid)
        !istat=nf90_put_var(ncid, rhvarid, TAUXanal) 
        !istat=nf90_close(ncid)
        !write(*,*) 'Finished TAUX dumping field.'
    end if

    !---
    !TAUY
    !---
    call gather_chunk_to_field(1,1,1,cb24cnn_nlon,cb24cnn_Model_TAUY,Xtransf)
    if (masterproc) then
        do ilat=1,cb24cnn_nlat
        do ilon=1,cb24cnn_nlon
           TAUYanal(ilon,ilat)=Xtransf(1,ilon,ilat)
        end do
        end do
    end if

    !---
    !USTAR
    !---
    call gather_chunk_to_field(1,1,1,cb24cnn_nlon,cb24cnn_Model_USTAR,Xtransf)
    if (masterproc) then
        do ilat=1,cb24cnn_nlat
        do ilon=1,cb24cnn_nlon
           USTARanal(ilon,ilat)=Xtransf(1,ilon,ilat)
        end do
        end do
    end if

    !---
    !TBOT
    !---
    call gather_chunk_to_field(1,1,1,cb24cnn_nlon,cb24cnn_Model_TBOT,Xtransf)
    if (masterproc) then
        do ilat=1,cb24cnn_nlat
        do ilon=1,cb24cnn_nlon
           TBOTanal(ilon,ilat)=Xtransf(1,ilon,ilat)
        end do
        end do
    end if

    !---
    !PBLH
    !---
    call gather_chunk_to_field(1,1,1,cb24cnn_nlon,cb24cnn_Model_PBLH,Xtransf)
    if (masterproc) then
        do ilat=1,cb24cnn_nlat
        do ilon=1,cb24cnn_nlon
           PBLHanal(ilon,ilat)=Xtransf(1,ilon,ilat)
        end do
        end do
    end if

    
    
    if (masterproc) write(iulog,*) 'finished gather .....'

    !All your python stuff is here:

    if (masterproc) then

        calday = get_curr_calday()
        call get_curr_date(Year,Month,Day,Sec)
        Hour = Sec/3600.0_r8
        
        call ieee_set_halting_mode(ieee_all, .false.)
    
        ierror = tuple_create(args, 13)
        
        ierror = ndarray_create(in1_py, Vanal)
        if (ierror/=0) then; call err_print; endif
        ierror = args%setitem(0,in1_py)
    
        ierror = ndarray_create(in2_py, Tanal)
        if (ierror/=0) then; call err_print; endif
        ierror = args%setitem(1,in2_py)
    
        ierror = ndarray_create(in3_py, Uanal)
        if (ierror/=0) then; call err_print; endif
        ierror = args%setitem(2,in3_py)
    
        ierror = ndarray_create(in4_py, Qanal)
        if (ierror/=0) then; call err_print; endif
        ierror = args%setitem(3,in4_py)
    
        ierror = ndarray_create(in5_py, Wanal)
        if (ierror/=0) then; call err_print; endif
        ierror = args%setitem(4,in5_py)
    
        ierror = ndarray_create(in6_py, TAUXanal)
        if (ierror/=0) then; call err_print; endif
        ierror = args%setitem(5,in6_py)
    
        ierror = ndarray_create(in7_py, TAUYanal)
        if (ierror/=0) then; call err_print; endif
        ierror = args%setitem(6,in7_py)
    
        ierror = ndarray_create(in8_py, USTARanal)
        if (ierror/=0) then; call err_print; endif
        ierror = args%setitem(7,in8_py)
    
        ierror = ndarray_create(in9_py, UPWPanal)
        if (ierror/=0) then; call err_print; endif
        ierror = args%setitem(8,in9_py)
    
        ierror = ndarray_create(in10_py, TBOTanal)
        if (ierror/=0) then; call err_print; endif
        ierror = args%setitem(9,in10_py)
        
        ierror = ndarray_create(in11_py, PBLHanal)
        if (ierror/=0) then; call err_print; endif
        ierror = args%setitem(10,in11_py)
    
        ierror = args%setitem(11, Hour)
        ierror = args%setitem(12, calday)
    
        if (masterproc) write(iulog,*)'... placed ndarray ...'
        ierror = call_py(return_value,pymodule,"DAMLcnn_run", args)
        if (ierror/=0) then; call err_print; endif
        ierror = cast(out_arr, return_value)
        if (ierror/=0) then; call err_print; endif
        ierror = out_arr%get_data(out_for, order='C')
        if (ierror/=0) then; call err_print; endif
    
        write(iulog,*) 'outfld yahoo:',size(out_for, 1),size(out_for, 2),size(out_for, 3),size(out_for, 4)
        write(iulog,*) 'poutfld yahoo:',size(Uanal, 1),size(Uanal, 2),size(Uanal, 3)
        call args%destroy
        call return_value%destroy
    
         
        !istat=nf90_create('CESM_dumpOUTFOR.nc',NF90_CLOBBER,ncid)
        !!istat=nf90_def_dim(ncid, "lon",cb24cnn_nlon , londimid)
        !istat=nf90_def_dim(ncid, "lat",cb24cnn_nlat , latdimid)
        !istat=nf90_def_dim(ncid, "lev",cb24cnn_nlev , levdimid)
        !istat=nf90_def_dim(ncid, "var",4, vardimid)
        !istat=nf90_def_var(ncid, "TAUX", nf90_double, (/ londimid, latdimid, levdimid, vardimid /), rhvarid)
        !istat=nf90_enddef(ncid)
        !istat=nf90_put_var(ncid, rhvarid, out_for) 
        !istat=nf90_close(ncid)
        !write(*,*) 'Finished out_for dumping field.'
    
        if (masterproc) write(iulog,*)'cb24cnn_init ending '
        !here we need to do a scatter: 
        call ieee_set_halting_mode(ieee_all, halting_mode)

    !All your python stuff is here:
    end if
    
    if (masterproc) then
    do ilat=1,cb24cnn_nlat
        do ilev=1,cb24cnn_nlev
            do ilon=1,cb24cnn_nlon
               Xtrans(ilon,ilev,ilat)=out_for(ilon,ilat,ilev,1) !U
            end do
        end do
    end do
    end if
    call scatter_field_to_chunk(1,cb24cnn_nlev,1,cb24cnn_nlon,Xtrans,    &
                                cb24cnn_Ustep(1,1,begchunk))


    if (masterproc) then
    do ilat=1,cb24cnn_nlat
        do ilev=1,cb24cnn_nlev
            do ilon=1,cb24cnn_nlon
               Xtrans(ilon,ilev,ilat)=out_for(ilon,ilat,ilev,2) !V
            end do
        end do
    end do
    end if

    call scatter_field_to_chunk(1,cb24cnn_nlev,1,cb24cnn_nlon,Xtrans,    &
                                cb24cnn_Vstep(1,1,begchunk))

    if (masterproc) then
    do ilat=1,cb24cnn_nlat
        do ilev=1,cb24cnn_nlev
            do ilon=1,cb24cnn_nlon
               Xtrans(ilon,ilev,ilat)=out_for(ilon,ilat,ilev,3) !T
            end do
        end do
    end do
    end if

    call scatter_field_to_chunk(1,cb24cnn_nlev,1,cb24cnn_nlon,Xtrans,    &
                                cb24cnn_Sstep(1,1,begchunk))

    
    if (masterproc) then
    do ilat=1,cb24cnn_nlat
        do ilev=1,cb24cnn_nlev
            do ilon=1,cb24cnn_nlon
               Xtrans(ilon,ilev,ilat)=out_for(ilon,ilat,ilev,4) !Q
            end do
        end do
    end do
    end if

    call scatter_field_to_chunk(1,cb24cnn_nlev,1,cb24cnn_nlon,Xtrans,    &
                                cb24cnn_Qstep(1,1,begchunk))
    
    do c=begchunk,endchunk
        !update phys_tend
        call outfld('cb24cnn_U',cb24cnn_Ustep(:,:,c),pcols,c)
        call outfld('cb24cnn_V',cb24cnn_Vstep(:,:,c),pcols,c)
        call outfld('cb24cnn_T',cb24cnn_Sstep(:,:,c),pcols,c)
        call outfld('cb24cnn_Q',cb24cnn_Qstep(:,:,c),pcols,c)
    end do

  !write(iulog,*) 'done with that stuff, again.....'
   ! End Routine
   !------------
  end subroutine cb24cnn_timestep_tend 

  subroutine cb24cnn_set_tend(phys_state,phys_tend)
  !
   ! CB24cnn_set_teset:
   !                If Nudging is ON, return the Nudging contributions
   !                to forcing using the current contents of the Nudge
   !                arrays. Send output to the cam history module as well.
   !===============================================================
   use physconst    ,only: cpair
   use physics_types,only: physics_state,physics_ptend,physics_ptend_init
   use constituents ,only: cnst_get_ind,pcnst
   use ppgrid       ,only: pver,pcols,begchunk,endchunk
   use cam_history  ,only: outfld

   ! Arguments
   !-------------
   type(physics_state), intent(in) :: phys_state
   type(physics_ptend), intent(out):: phys_tend

   ! Local values
   !--------------------
   integer indw,ncol,lchnk
   logical lq(pcnst)

   call cnst_get_ind('Q',indw)
   lq(:)   =.false.
   lq(indw)=.true.
   call physics_ptend_init(phys_tend,phys_state%psetcols,'nudging',lu=.true.,lv=.true.,ls=.true.,lq=lq)

   lchnk=phys_state%lchnk
   ncol =phys_state%ncol
   phys_tend%u(:ncol,:pver)     =cb24cnn_Ustep(:ncol,:pver,lchnk)
   phys_tend%v(:ncol,:pver)     =cb24cnn_Vstep(:ncol,:pver,lchnk)
   phys_tend%s(:ncol,:pver)     =cb24cnn_Sstep(:ncol,:pver,lchnk)
   phys_tend%q(:ncol,:pver,indw)=cb24cnn_Qstep(:ncol,:pver,lchnk)

  end subroutine cb24cnn_set_tend 
  
  subroutine cb24cnn_finalize()
    call pymodule%destroy
    call forpy_finalize
  end subroutine cb24cnn_finalize
  end module cb24cnn
