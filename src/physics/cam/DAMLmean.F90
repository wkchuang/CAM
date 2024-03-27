module cb24mean

!=====================================================================
!
! Purpose: Implement cb24mean of the model state of U,V,T,Q, and/or PS
!          toward specified values
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
!    cb24mean the model toward the analyses values determined from the forcing
!    option specified. After the model passes the ending analyses time, the
!    forcing discontinues.
!
!    Some beans analyses products can have gaps in the available data, where values
!    are missing for some interval of time. When files are missing, the cb24mean
!    force is switched off for that interval of time, so we effectively 'coast'
!    thru the gap.
!
!    Currently, the cb24mean module is set up to accomodate cb24mean of PS
!    values, however that functionality requires forcing that is applied in
!    the selected dycore and is not yet implemented.
!
!    The cb24mean of the model toward the analyses data is controlled by
!    the 'cb24mean_nl' namelist in 'user_nl_cam'; whose variables control the
!    time interval over which cb24mean is applied, the strength of the cb24mean
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

  public:: cb24mean_Model!,cb24mean_Model_Regress,cb24mean_ON
  !public:: cb24mean_readnl
  public:: cb24mean_init
  public:: cb24mean_timestep_init
  public:: cb24mean_timestep_tend
  public:: cb24mean_set_tend
  public:: cb24mean_finalize

  !  cb24mean Parameters
  !--------------------
  logical          :: cb24mean_Model       =.true.
  logical          :: halting_mode(5)
  type(module_py)  :: pymodule

  integer :: &
    upwp_clubb_idx, &          ! Cloud fraction
    taux_idx,       &
    pblh_idx,       &
    ustarwec_idx,   &
    flntc_idx,      &
    fsntoa_idx

  ! cb24mean State Arrays
  integer cb24mean_nlon,cb24mean_nlat,cb24mean_ncol,cb24mean_nlev,cb24mean_nlevp
  real(r8),allocatable:: cb24mean_Ustep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: cb24mean_Vstep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: cb24mean_Sstep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: cb24mean_Qstep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: cb24mean_PSstep(:,:)    !(pcols,begchunk:endchunk)

  !CNN input vars
  real(r8),allocatable:: cb24mean_Model_V (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: cb24mean_Model_T (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: cb24mean_Model_U (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: cb24mean_Model_Q (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: cb24mean_Model_W (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: cb24mean_Model_UPWP (:,:,:)  !(pcols,pverp,begchunk:endchunk)
  real(r8),allocatable:: cb24mean_Model_TAUX (:,:)    !(pcols,begchunk:endchunk)
  real(r8),allocatable:: cb24mean_Model_TAUY (:,:)    !(pcols,begchunk:endchunk)
  real(r8),allocatable:: cb24mean_Model_USTAR (:,:)   !(pcols,begchunk:endchunk)
  real(r8),allocatable:: cb24mean_Model_TBOT (:,:) !(pcols,begchunk:endchunk)
  real(r8),allocatable:: cb24mean_Model_PBLH (:,:) !(pcols,begchunk:endchunk)

  !derived type interfaces:

  !> Control structure for Python interface
  type, public :: python_interface ; private
    type(module_py) :: pymodule
    type(list) :: paths
  end type

contains

   !================================================================
  subroutine cb24mean_init
    !
    ! cb24mean_INIT: Allocate space and initialize cb24mean values
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
    ! cb24mean Coefs and cb24mean Forcing.
    !-------------------------------------------
    allocate(cb24mean_Ustep(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24mean_init','cb24mean_Ustep',pcols*pver*((endchunk-begchunk)+1))
    allocate(cb24mean_Vstep(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24mean_init','cb24mean_Vstep',pcols*pver*((endchunk-begchunk)+1))
    allocate(cb24mean_Sstep(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24mean_init','cb24mean_Sstep',pcols*pver*((endchunk-begchunk)+1))
    allocate(cb24mean_Qstep(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24mean_init','cb24mean_Qstep',pcols*pver*((endchunk-begchunk)+1))
    allocate(cb24mean_PSstep(pcols,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24mean_init','cb24mean_PSstep',pcols*((endchunk-begchunk)+1)) !currently not used.

    allocate(cb24mean_Model_U(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24mean_init','cb24mean_Model_U',pcols*pver*((endchunk-begchunk)+1))
    allocate(cb24mean_Model_V(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24mean_init','cb24mean_Model_V',pcols*pver*((endchunk-begchunk)+1))
    allocate(cb24mean_Model_T(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24mean_init','cb24mean_Model_T',pcols*pver*((endchunk-begchunk)+1))
    allocate(cb24mean_Model_Q(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24mean_init','cb24mean_Model_Q',pcols*pver*((endchunk-begchunk)+1))
    allocate(cb24mean_Model_W(pcols,pverp,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24mean_init','cb24mean_Model_W',pcols*pverp*((endchunk-begchunk)+1))
    allocate(cb24mean_Model_UPWP(pcols,pverp,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24mean_init','cb24mean_Model_Q',pcols*pverp*((endchunk-begchunk)+1))

    allocate(cb24mean_Model_TAUX(pcols,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24mean_init','cb24mean_Model_TAUX',pcols*((endchunk-begchunk)+1)) !currently not used.
    allocate(cb24mean_Model_TAUY(pcols,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24mean_init','cb24mean_Model_TAUY',pcols*((endchunk-begchunk)+1)) !currently not used.
    allocate(cb24mean_Model_USTAR(pcols,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24mean_init','cb24mean_Model_USTAR',pcols*((endchunk-begchunk)+1)) !currently not used.
    allocate(cb24mean_Model_TBOT(pcols,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24mean_init','cb24mean_Model_TBOT',pcols*((endchunk-begchunk)+1)) !currently not used.
    allocate(cb24mean_Model_PBLH(pcols,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24mean_init','cb24mean_Model_PBLH',pcols*((endchunk-begchunk)+1)) !currently not used.

    ! Register output fields with the cam history module
     !-----------------------------------------------------
    call addfld( 'cb24mean_U',(/ 'lev' /),'A','m/s/s'  ,'U DAMLmean Tendency')
    call addfld( 'cb24mean_V',(/ 'lev' /),'A','m/s/s'  ,'V DAMLmean Tendency')
    call addfld( 'cb24mean_T',(/ 'lev' /),'A','K/s'    ,'T DAMLmean Tendency')
    call addfld( 'cb24mean_Q',(/ 'lev' /),'A','kg/kg/s','Q DAMLmean Tendency')

    ! Initialize column and level dimensions
     !--------------------------------------------------------
     call get_horiz_grid_dim_d(hdim1_d,hdim2_d)
     cb24mean_nlon=hdim1_d
     cb24mean_nlat=hdim2_d
     cb24mean_ncol=hdim1_d*hdim2_d
     cb24mean_nlev=pver
     cb24mean_nlevp=pverp

    ! Initialize the forpy module and import the python script "DAMLcnn.py"
     !-----------------------------------------------------
    if (masterproc) write(iulog,*) 'cb24mean_init starting'
    if (.not. ieee_support_halting(ieee_overflow)) then
       call endrun("ieee_halting is not supported")
    endif
    call ieee_get_halting_mode(ieee_all, halting_mode)
    print *,__FILE__,__LINE__,halting_mode
    call ieee_set_halting_mode(ieee_all, .false.)

    ierror = forpy_initialize()
    ierror = get_sys_path(paths)
    ierror = paths%append(".") !stash it in the run directory; future NL

    ierror = import_py(pymodule,"DAMLmean") !future NL
    call ieee_set_halting_mode(ieee_all, halting_mode)
    call paths%destroy

    ! Init forcing as zeros:
    do lchnk=begchunk,endchunk
      cb24mean_Ustep(:pcols,:pver,lchnk)=0._r8
      cb24mean_Vstep(:pcols,:pver,lchnk)=0._r8
      cb24mean_Sstep(:pcols,:pver,lchnk)=0._r8
      cb24mean_Qstep(:pcols,:pver,lchnk)=0._r8
      cb24mean_PSstep(:pcols,lchnk)=0._r8
    end do

  end subroutine cb24mean_init

  subroutine cb24mean_timestep_init(phys_state,phys_tend,pbuf,cam_in,cam_out)

    !
    ! cb24mean_TIMESTEP_INIT:
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

    cb24mean_Model_V(:ncol,:pver,lchnk)=phys_state%v(:ncol,:pver)
    cb24mean_Model_T(:ncol,:pver,lchnk)=phys_state%t(:ncol,:pver)
    cb24mean_Model_U(:ncol,:pver,lchnk)=phys_state%u(:ncol,:pver)
    cb24mean_Model_Q(:ncol,:pver,lchnk)=phys_state%q(:ncol,:pver,indw)
    cb24mean_Model_W(:ncol,:pver,lchnk)=phys_state%omega(:ncol,:pver)
    cb24mean_Model_UPWP(:ncol,:pverp,lchnk)=upwp(:ncol,:pverp)
    cb24mean_Model_TAUX(:ncol,lchnk)=cam_in%wsx(:ncol)
    cb24mean_Model_TAUY(:ncol,lchnk)=cam_in%wsy(:ncol)
    cb24mean_Model_TBOT(:ncol,lchnk)=cam_out%tbot(:ncol)
    cb24mean_Model_PBLH(:ncol,lchnk)=pblh(:ncol)
    cb24mean_Model_USTAR(:ncol,lchnk)=ustarwec(:ncol)

  end subroutine cb24mean_timestep_init

  subroutine cb24mean_timestep_tend()

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
    integer  Year,Month,Day,Sec,calday
    real(r8) Hour
    
    real(r8) Xtrans(cb24mean_nlon,cb24mean_nlev,cb24mean_nlat)
    integer  nn,Nindex
    real(r8) PI

    type(ndarray) :: out_arr   !< variables in the form of numpy array
    character(len=:), allocatable :: return_string !< outputs from Python module
    real(r8), pointer, dimension(:,:,:,:) :: out_for !< outputs from Python module
    
    if (masterproc) write(iulog,*) 'cb24mean_timestep_tend starting....'    
    if (masterproc) write(iulog,*) 'finished gather .....'

    !All your python stuff is here:

    if (masterproc) then

        calday = get_curr_calday()
        call get_curr_date(Year,Month,Day,Sec)
        Hour = Sec/3600.0_r8
        
        call ieee_set_halting_mode(ieee_all, .false.)
    
        ierror = tuple_create(args, 2)
        ierror = args%setitem(0, Hour)
        ierror = args%setitem(1, calday)
    
        if (masterproc) write(iulog,*)'... placed ndarray ...'
        ierror = call_py(return_value,pymodule,"DAMLmean_run", args)
        if (masterproc) write(iulog,*)'... came back from dead ...'
        if (ierror/=0) then; call err_print; endif
        ierror = cast(out_arr, return_value)
        if (ierror/=0) then; call err_print; endif
        ierror = out_arr%get_data(out_for, order='C')
        !ierror = out_arr%get_data(out_for)
        if (ierror/=0) then; call err_print; endif
    
        write(iulog,*) 'outfld yahoo-mean:',size(out_for, 1),size(out_for, 2),size(out_for, 3),size(out_for, 4)
        call args%destroy
        call return_value%destroy
    
        !istat=nf90_create('CESM_dumpOUTFOR.nc',NF90_CLOBBER,ncid)
        !!istat=nf90_def_dim(ncid, "lon",cb24mean_nlon , londimid)
        !istat=nf90_def_dim(ncid, "lat",cb24mean_nlat , latdimid)
        !istat=nf90_def_dim(ncid, "lev",cb24mean_nlev , levdimid)
        !istat=nf90_def_dim(ncid, "var",4, vardimid)
        !istat=nf90_def_var(ncid, "TAUX", nf90_double, (/ londimid, latdimid, levdimid, vardimid /), rhvarid)
        !istat=nf90_enddef(ncid)
        !istat=nf90_put_var(ncid, rhvarid, out_for) 
        !istat=nf90_close(ncid)
        !write(*,*) 'Finished out_for dumping field.'
    
        if (masterproc) write(iulog,*)'cb24mean_init ending '
        !here we need to do a scatter: 
        call ieee_set_halting_mode(ieee_all, halting_mode)

    !All your python stuff is here:
    end if
    
    if (masterproc) then
    do ilat=1,cb24mean_nlat
        do ilev=1,cb24mean_nlev
            do ilon=1,cb24mean_nlon
               Xtrans(ilon,ilev,ilat)=out_for(ilon,ilat,ilev,1) !U
            end do
        end do
    end do
    end if
    call scatter_field_to_chunk(1,cb24mean_nlev,1,cb24mean_nlon,Xtrans,    &
                                cb24mean_Ustep(1,1,begchunk))


    if (masterproc) then
    do ilat=1,cb24mean_nlat
        do ilev=1,cb24mean_nlev
            do ilon=1,cb24mean_nlon
               Xtrans(ilon,ilev,ilat)=out_for(ilon,ilat,ilev,2) !V
            end do
        end do
    end do
    end if

    call scatter_field_to_chunk(1,cb24mean_nlev,1,cb24mean_nlon,Xtrans,    &
                                cb24mean_Vstep(1,1,begchunk))

    if (masterproc) then
    do ilat=1,cb24mean_nlat
        do ilev=1,cb24mean_nlev
            do ilon=1,cb24mean_nlon
               Xtrans(ilon,ilev,ilat)=out_for(ilon,ilat,ilev,3) !T
            end do
        end do
    end do
    end if

    call scatter_field_to_chunk(1,cb24mean_nlev,1,cb24mean_nlon,Xtrans,    &
                                cb24mean_Sstep(1,1,begchunk))

    
    if (masterproc) then
    do ilat=1,cb24mean_nlat
        do ilev=1,cb24mean_nlev
            do ilon=1,cb24mean_nlon
               Xtrans(ilon,ilev,ilat)=out_for(ilon,ilat,ilev,4) !Q
            end do
        end do
    end do
    end if

    call scatter_field_to_chunk(1,cb24mean_nlev,1,cb24mean_nlon,Xtrans,    &
                                cb24mean_Qstep(1,1,begchunk))
    
    do c=begchunk,endchunk
        !update phys_tend
        call outfld('cb24mean_U',cb24mean_Ustep(:,:,c),pcols,c)
        call outfld('cb24mean_V',cb24mean_Vstep(:,:,c),pcols,c)
        call outfld('cb24mean_T',cb24mean_Sstep(:,:,c),pcols,c)
        call outfld('cb24mean_Q',cb24mean_Qstep(:,:,c),pcols,c)
    end do

    if (masterproc) write(iulog,*)'exiting cb24mean_timestep_tend'

  !write(iulog,*) 'done with that stuff, again.....'
   ! End Routine
   !------------
  end subroutine cb24mean_timestep_tend 

  subroutine cb24mean_set_tend(phys_state,phys_tend)
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
   call physics_ptend_init(phys_tend,phys_state%psetcols,'DAMLmean',lu=.true.,lv=.true.,ls=.true.,lq=lq)

   lchnk=phys_state%lchnk
   ncol =phys_state%ncol
   phys_tend%u(:ncol,:pver)     =cb24mean_Ustep(:ncol,:pver,lchnk)
   phys_tend%v(:ncol,:pver)     =cb24mean_Vstep(:ncol,:pver,lchnk)
   phys_tend%s(:ncol,:pver)     =cb24mean_Sstep(:ncol,:pver,lchnk)
   phys_tend%q(:ncol,:pver,indw)=cb24mean_Qstep(:ncol,:pver,lchnk)

  end subroutine cb24mean_set_tend 
  
  subroutine cb24mean_finalize()
    call pymodule%destroy
    call forpy_finalize
  end subroutine cb24mean_finalize
  end module cb24mean
