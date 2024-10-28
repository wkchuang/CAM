module cb24mjocnn
!=====================================================================
!
! Purpose: Implement cb24mjocnn of the model state of U,V,T,Q, and/or PS
!          toward specified values from random analyses increments.
!
! Author: Wiliiam Chapman -- wchapman@ucar.edu
!
! Description:
!
!    This module assumes that the user has {Nudging Increments} values from analyses
!    which have been preprocessed onto the current model grid and adjusted
!    for differences in topography. It is also assumed that these resulting
!    values and are stored in individual files which are indexed with respect
!    to year, month, day, and second of the day. When the model is inbetween
!    the given begining and ending times, a relaxation forcing is added to
!    cb24mjocnn the model toward the analyses values determined from the forcing
!    option specified. After the model passes the ending analyses time, the
!    forcing discontinues.

!    The cb24mjocnn of the model toward the analyses data is controlled by
!    the 'cb24mjocnn_nl' namelist in 'user_nl_cam'; whose variables control the
!    time interval over which cb24mjocnn is applied, the strength of the cb24mjocnn
!    tendencies, and its spatial distribution.
!
! Useful modules
  !------------------
  use shr_kind_mod,   only:r8=>SHR_KIND_R8,cs=>SHR_KIND_CS,cl=>SHR_KIND_CL,r4=>SHR_KIND_R4
  use time_manager,   only:timemgr_time_ge,timemgr_time_inc,timemgr_time_inc_minus,get_curr_date,get_curr_calday,get_step_size
  use, intrinsic :: ieee_exceptions, only: ieee_get_halting_mode, ieee_set_halting_mode, ieee_all, ieee_support_halting, ieee_overflow
  use phys_grid   ,   only:scatter_field_to_chunk, gather_chunk_to_field
  use cam_abortutils, only:endrun
  use spmd_utils  ,   only:masterproc, mstrid=>masterprocid, mpicom
  use cam_logfile ,   only:iulog
  use cam_history_support, only: add_hist_coord
#ifdef SPMD
  use mpishorthand
#endif

  use ftorch
  
! Set all Global values and routines to private by default
  ! and then explicitly set their exposure.
  !----------------------------------------------------------
  implicit none
  private

  public:: cb24mjocnn_Model 
  !public:: cb24mjocnn_readnl
  public:: cb24mjocnn_init
  public:: cb24mjocnn_timestep_tender
  public:: cb24mjocnn_timestep_tend
  public:: cb24mjocnn_timestep_init
  public:: cb24mjocnn_finalize
  public:: cb24mjocnn_register_cam

  !private::cb24mjocnn_readin_weights_all
  !private::cb24mjocnn_readin_weights
  !private::cb24mjocnn_readin_weights_diurnal
  !private::cb24mjocnn_readin_normdict
  !private::cb24mjocnn_set_PSprofile
  !private::cb24mjocnn_set_profile

  !  cb24mjocnn Parameters
  !--------------------
  logical          :: cb24mjocnn_Model       =.true.
  logical          :: halting_mode(5)
  character(len=cl):: cb24mjocnn_Path
  character(len=cs):: cb24mjocnn_File1
  character(len=cs):: cb24mjocnn_File2
  character(len=cs):: cb24mjocnn_File3
  character(len=cs):: cb24mjocnn_File4
  character(len=cs):: cb24mjocnn_File5
  character(len=cs):: cb24mjocnn_File6
  character(len=cs):: cb24mjocnn_File7
  character(len=cs):: cb24mjocnn_File8

  integer :: &
    upwp_clubb_idx, &          ! Cloud fraction
    taux_idx,       &
    pblh_idx,       &
    ustarwec_idx,   &
    flutc_idx,      &
    fsntoa_idx


  ! cb24mjocnn State Arrays
  integer cb24mjocnn_nlon,cb24mjocnn_nlat,cb24mjocnn_ncol,cb24mjocnn_nlev,cb24mjocnn_nlevp
  integer cb24mjocnn_NumObs
  real(r8),allocatable:: cb24mjocnn_Ustep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: cb24mjocnn_Vstep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: cb24mjocnn_Sstep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: cb24mjocnn_Qstep (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: cb24mjocnn_PSstep(:,:)    !(pcols,begchunk:endchunk)

  !CNN input vars
  real(r8),allocatable:: cb24mjocnn_Probs (:,:,:)  !(nprobs,pver,begchunk:endchunk)
  real(r8),allocatable:: cb24mjocnn_Model_V (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: cb24mjocnn_Model_U (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: cb24mjocnn_Model_V200 (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: cb24mjocnn_Model_U200 (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: cb24mjocnn_Model_V850 (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: cb24mjocnn_Model_U850 (:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8),allocatable:: cb24mjocnn_Model_FLUTC (:,:) !(pcols,begchunk:endchunk)
  
  !derived type interfaces:
  real(r8),allocatable::MJO_Inc_00 (:,:,:,:) !(pcols,pver,begchunk:endchunk,cb24mjocnn_NumObs)
  real(r8),allocatable::MJO_Inc_01 (:,:,:,:) !(pcols,pver,begchunk:endchunk,cb24mjocnn_NumObs)
  real(r8),allocatable::MJO_Inc_02 (:,:,:,:) !(pcols,pver,begchunk:endchunk,cb24mjocnn_NumObs)
  real(r8),allocatable::MJO_Inc_03 (:,:,:,:) !(pcols,pver,begchunk:endchunk,cb24mjocnn_NumObs)
  real(r8),allocatable::MJO_Inc_04 (:,:,:,:) !(pcols,pver,begchunk:endchunk,cb24mjocnn_NumObs)
  real(r8),allocatable::MJO_Inc_05 (:,:,:,:) !(pcols,pver,begchunk:endchunk,cb24mjocnn_NumObs)
  real(r8),allocatable::MJO_Inc_06 (:,:,:,:) !(pcols,pver,begchunk:endchunk,cb24mjocnn_NumObs)
  real(r8),allocatable::MJO_Inc_07 (:,:,:,:) !(pcols,pver,begchunk:endchunk,cb24mjocnn_NumObs)
  real(r8),allocatable::MJO_Inc_08 (:,:,:,:) !(pcols,pver,begchunk:endchunk,cb24mjocnn_NumObs)

  !FTORCH vectors:
  type(torch_model) :: model_ftorch !ftorch
  real(r4) output_ftorch(9) !ftorch
  real(r8),allocatable:: cb24mjocnn_Model_input (:,:,:,:)  !(pcols,pver,begchunk:endchunk)
  real(r8), dimension(:), allocatable, target :: probabilities_ftorch !ftorch
  
  
contains

  subroutine cb24mjocnn_register_cam()
  !-------------------------------------------------------------------------------
    ! Description:
    !   Register the constituents and fields in the physics buffer
    ! Author: Will Fucking Chapman
    !
    !-------------------------------------------------------------------------------
    use ppgrid        ,only: pver,pverp,pcols,begchunk,endchunk
    use error_messages,only: alloc_err
    use dycore        ,only: dycore_is
    use dyn_grid      ,only: get_horiz_grid_dim_d
    use phys_grid     ,only: get_rlat_p,get_rlon_p,get_ncols_p
    use cam_history   ,only: addfld
    use shr_const_mod ,only: SHR_CONST_PI
    use filenames     ,only: interpret_filename_spec


    if(masterproc) write(iulog,*)'WAGO: Register nprob'
    call add_hist_coord('nprob', 9, 'mjo probability index')
    if(masterproc) write(iulog,*)'Done'

  end subroutine cb24mjocnn_register_cam

   !================================================================
  subroutine cb24mjocnn_init
    !
    ! cb24mjocnn_INIT: Allocate space and initialize cb24mjocnn values
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
    integer            :: ierr
    integer, parameter :: n = 9

    !local Values
    !-------------
    integer  istat,lchnk
    integer  hdim1_d,hdim2_d
    character(len=:), allocatable :: return_string
    character(len=:), allocatable :: Torchpath

    real(r4) :: input_vector(n) = [2.0, 1.0, 0.5, 3.0, 2.5, &
                                   2.0, 0.0, 1.0, 2.0]
    real(r4) :: softmax_vector(n)
    real(r4) :: max_val
    integer :: i
    allocate(probabilities_ftorch(9))
    ierr = 42

    
    ! Allocate Space for cb24 data arrays
    !-----------------------------------------

    ! Set Default Namelist values
    !-----------------------------

    cb24mjocnn_Path   = '/glade/work/wchapman/DA_ML/CESML_AI/MJO_Increments/'
    cb24mjocnn_File1  = 'MJO_phase_increment_001.nc'
    cb24mjocnn_File2  = 'MJO_phase_increment_002.nc'
    cb24mjocnn_File3  = 'MJO_phase_increment_003.nc'
    cb24mjocnn_File4  = 'MJO_phase_increment_004.nc'
    cb24mjocnn_File5  = 'MJO_phase_increment_005.nc'
    cb24mjocnn_File6  = 'MJO_phase_increment_006.nc'
    cb24mjocnn_File7  = 'MJO_phase_increment_007.nc'
    cb24mjocnn_File8  = 'MJO_phase_increment_008.nc'

    ! Allocate Space for spatial dependence of
    ! cb24mjocnn Coefs and cb24mjocnn Forcing.
    !-------------------------------------------
    allocate(cb24mjocnn_Ustep(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24mjocnn_init','cb24mjocnn_Ustep',pcols*pver*((endchunk-begchunk)+1))
    allocate(cb24mjocnn_Vstep(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24mjocnn_init','cb24mjocnn_Vstep',pcols*pver*((endchunk-begchunk)+1))
    allocate(cb24mjocnn_Sstep(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24mjocnn_init','cb24mjocnn_Sstep',pcols*pver*((endchunk-begchunk)+1))
    allocate(cb24mjocnn_Qstep(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24mjocnn_init','cb24mjocnn_Qstep',pcols*pver*((endchunk-begchunk)+1))
    allocate(cb24mjocnn_PSstep(pcols,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24mjocnn_init','cb24mjocnn_PSstep',pcols*((endchunk-begchunk)+1)) !currently not used.


    allocate(cb24mjocnn_Model_U(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24mjocnn_init','cb24mjocnn_Model_U',pcols*pver*((endchunk-begchunk)+1))
    allocate(cb24mjocnn_Model_V(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24mjocnn_init','cb24mjocnn_Model_V',pcols*pver*((endchunk-begchunk)+1))
    allocate(cb24mjocnn_Probs(pcols,9,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24mjocnn_init','cb24mjocnn_Model_Probs',pcols*9*((endchunk-begchunk)+1))
    allocate(cb24mjocnn_Model_U200(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24mjocnn_init','cb24mjocnn_Model_U200',pcols*pver*((endchunk-begchunk)+1))
    allocate(cb24mjocnn_Model_V200(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24mjocnn_init','cb24mjocnn_Model_V200',pcols*pver*((endchunk-begchunk)+1))
    allocate(cb24mjocnn_Model_U850(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24mjocnn_init','cb24mjocnn_Model_U850',pcols*pver*((endchunk-begchunk)+1))
    allocate(cb24mjocnn_Model_V850(pcols,pver,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24mjocnn_init','cb24mjocnn_Model_V850',pcols*pver*((endchunk-begchunk)+1))

    allocate(cb24mjocnn_Model_FLUTC(pcols,begchunk:endchunk),stat=istat)
    call alloc_err(istat,'cb24mjocnn_init','cb24mjocnn_Model_FLUTC',pcols*((endchunk-begchunk)+1))


    ! Register output fields with the cam history module
     !-----------------------------------------------------
    call addfld( 'cb24mjocnn_U',(/ 'lev' /),'A','m/s/s'  ,'U cb24mjocnn Tendency')
    call addfld( 'cb24mjocnn_V',(/ 'lev' /),'A','m/s/s'  ,'V cb24mjocnn Tendency')
    call addfld( 'cb24mjocnn_T',(/ 'lev' /),'A','K/s'    ,'T cb24mjocnn Tendency')
    call addfld( 'cb24mjocnn_Q',(/ 'lev' /),'A','kg/kg/s','Q cb24mjocnn Tendency')
    call addfld( 'state_mjo',(/ 'nprob' /),'A','fraction','MJO state Prob')
    ! Initialize column and level dimensions
    !--------------------------------------------------------
    call get_horiz_grid_dim_d(hdim1_d,hdim2_d)
    cb24mjocnn_nlon=hdim1_d
    cb24mjocnn_nlat=hdim2_d
    cb24mjocnn_ncol=hdim1_d*hdim2_d
    cb24mjocnn_nlev=pver
    cb24mjocnn_nlevp=pverp
    cb24mjocnn_NumObs = 2

    ! Broadcast other variables that have changed
    !---------------------------------------------
#ifdef SPMD
    call mpibcast(cb24mjocnn_NumObs        ,            1, mpiint, 0, mpicom)
#endif

    !Load the six MJO arrays... 
    ! Allocate Space for cb24mjocnn MJO Nudging arrays, initialize with 0's
    !---------------------------------------------------------------------
    allocate(MJO_Inc_00(pcols,pver,begchunk:endchunk,cb24mjocnn_NumObs),stat=istat)
    call alloc_err(istat,'cb24mjocnn_init','MJO_Inc_00',pcols*pver*((endchunk-begchunk)+1)*cb24mjocnn_NumObs)
    allocate(MJO_Inc_01(pcols,pver,begchunk:endchunk,cb24mjocnn_NumObs),stat=istat)
    call alloc_err(istat,'cb24mjocnn_init','MJO_Inc_01',pcols*pver*((endchunk-begchunk)+1)*cb24mjocnn_NumObs)
    allocate(MJO_Inc_02(pcols,pver,begchunk:endchunk,cb24mjocnn_NumObs),stat=istat)
    call alloc_err(istat,'cb24mjocnn_init','MJO_Inc_02',pcols*pver*((endchunk-begchunk)+1)*cb24mjocnn_NumObs)
    allocate(MJO_Inc_03(pcols,pver,begchunk:endchunk,cb24mjocnn_NumObs),stat=istat)
    call alloc_err(istat,'cb24mjocnn_init','MJO_Inc_03',pcols*pver*((endchunk-begchunk)+1)*cb24mjocnn_NumObs)
    allocate(MJO_Inc_04(pcols,pver,begchunk:endchunk,cb24mjocnn_NumObs),stat=istat)
    call alloc_err(istat,'cb24mjocnn_init','MJO_Inc_04',pcols*pver*((endchunk-begchunk)+1)*cb24mjocnn_NumObs)
    allocate(MJO_Inc_05(pcols,pver,begchunk:endchunk,cb24mjocnn_NumObs),stat=istat)
    call alloc_err(istat,'cb24mjocnn_init','MJO_Inc_05',pcols*pver*((endchunk-begchunk)+1)*cb24mjocnn_NumObs)
    allocate(MJO_Inc_06(pcols,pver,begchunk:endchunk,cb24mjocnn_NumObs),stat=istat)
    call alloc_err(istat,'cb24mjocnn_init','MJO_Inc_06',pcols*pver*((endchunk-begchunk)+1)*cb24mjocnn_NumObs)
    allocate(MJO_Inc_07(pcols,pver,begchunk:endchunk,cb24mjocnn_NumObs),stat=istat)
    call alloc_err(istat,'cb24mjocnn_init','MJO_Inc_07',pcols*pver*((endchunk-begchunk)+1)*cb24mjocnn_NumObs)
    allocate(MJO_Inc_08(pcols,pver,begchunk:endchunk,cb24mjocnn_NumObs),stat=istat)
    call alloc_err(istat,'cb24mjocnn_init','MJO_Inc_08',pcols*pver*((endchunk-begchunk)+1)*cb24mjocnn_NumObs)
    
    MJO_Inc_00(:pcols,:pver,begchunk:endchunk,:cb24mjocnn_NumObs)=0._r8
    MJO_Inc_01(:pcols,:pver,begchunk:endchunk,:cb24mjocnn_NumObs)=0._r8
    MJO_Inc_02(:pcols,:pver,begchunk:endchunk,:cb24mjocnn_NumObs)=0._r8
    MJO_Inc_03(:pcols,:pver,begchunk:endchunk,:cb24mjocnn_NumObs)=0._r8
    MJO_Inc_04(:pcols,:pver,begchunk:endchunk,:cb24mjocnn_NumObs)=0._r8
    MJO_Inc_05(:pcols,:pver,begchunk:endchunk,:cb24mjocnn_NumObs)=0._r8
    MJO_Inc_06(:pcols,:pver,begchunk:endchunk,:cb24mjocnn_NumObs)=0._r8
    MJO_Inc_07(:pcols,:pver,begchunk:endchunk,:cb24mjocnn_NumObs)=0._r8
    MJO_Inc_08(:pcols,:pver,begchunk:endchunk,:cb24mjocnn_NumObs)=0._r8

    if(dycore_is('UNSTRUCTURED')) then
        call cb24mjocnn_readin_weights_all(trim(cb24mjocnn_Path)//trim(cb24mjocnn_File1),&
                                       trim(cb24mjocnn_Path)//trim(cb24mjocnn_File2),&
                                       trim(cb24mjocnn_Path)//trim(cb24mjocnn_File3),&
                                       trim(cb24mjocnn_Path)//trim(cb24mjocnn_File4),&
                                       trim(cb24mjocnn_Path)//trim(cb24mjocnn_File5),&
                                       trim(cb24mjocnn_Path)//trim(cb24mjocnn_File6),&
                                       trim(cb24mjocnn_Path)//trim(cb24mjocnn_File7),&
                                       trim(cb24mjocnn_Path)//trim(cb24mjocnn_File8))
    elseif(dycore_is('EUL')) then
         call cb24mjocnn_readin_weights_all(trim(cb24mjocnn_Path)//trim(cb24mjocnn_File1),&
                                       trim(cb24mjocnn_Path)//trim(cb24mjocnn_File2),&
                                       trim(cb24mjocnn_Path)//trim(cb24mjocnn_File3),&
                                       trim(cb24mjocnn_Path)//trim(cb24mjocnn_File4),&
                                       trim(cb24mjocnn_Path)//trim(cb24mjocnn_File5),&
                                       trim(cb24mjocnn_Path)//trim(cb24mjocnn_File6),&
                                       trim(cb24mjocnn_Path)//trim(cb24mjocnn_File7),&
                                       trim(cb24mjocnn_Path)//trim(cb24mjocnn_File8))
    else !if(dycore_is('LR')) then
         call cb24mjocnn_readin_weights_all(trim(cb24mjocnn_Path)//trim(cb24mjocnn_File1),&
                                       trim(cb24mjocnn_Path)//trim(cb24mjocnn_File2),&
                                       trim(cb24mjocnn_Path)//trim(cb24mjocnn_File3),&
                                       trim(cb24mjocnn_Path)//trim(cb24mjocnn_File4),&
                                       trim(cb24mjocnn_Path)//trim(cb24mjocnn_File5),&
                                       trim(cb24mjocnn_Path)//trim(cb24mjocnn_File6),&
                                       trim(cb24mjocnn_Path)//trim(cb24mjocnn_File7),&
                                       trim(cb24mjocnn_Path)//trim(cb24mjocnn_File8))
    endif !read mjo weights in
    
    !

    !Load the Pytorch Model
    Torchpath = '/glade/work/wchapman/DA_ML/CESML_AI/MJO_Increments/saved_MJOcnnClassify_SS_model_cpu.pt'
    call torch_model_load(model_ftorch, Torchpath)
    
    ! Init forcing as zeros:
    do lchnk=begchunk,endchunk
      cb24mjocnn_Ustep(:pcols,:pver,lchnk)=0._r8
      cb24mjocnn_Vstep(:pcols,:pver,lchnk)=0._r8
      cb24mjocnn_Sstep(:pcols,:pver,lchnk)=0._r8
      cb24mjocnn_Qstep(:pcols,:pver,lchnk)=0._r8
      cb24mjocnn_PSstep(:pcols,lchnk)=0._r8
    end do
    probabilities_ftorch(:) = 0._r8
  end subroutine cb24mjocnn_init


  subroutine cb24mjocnn_timestep_tender()
    use physics_types ,only: physics_state,physics_ptend,physics_ptend_init,physics_state_copy
    use physics_buffer,only: pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field, &
                              physics_buffer_desc
    use time_manager  ,only:timemgr_time_ge,timemgr_time_inc,timemgr_time_inc_minus,get_curr_date,get_curr_calday,get_step_size
    use ppgrid        ,only: pver,pverp,pcols,begchunk,endchunk
    use cam_history   ,only: outfld
    use netcdf
    
    !+++ Ftorch tensors
    type(torch_tensor), dimension(1) :: in_tensor, out_tensor !ftorch
    !--- Ftorch tensors

    !local values:
    !---------------
    integer lev
    integer c                               
    integer nlon,nlat,plev,istat,lchnk,indw,ii,ierr
    integer ncid,varid
    integer ilat,ilon,ilev
    integer londimid, latdimid, levdimid,vardimid, rhvarid,singdimid
    integer  Year,Month,Day,Sec,Hour
    real(r8) calday
    real(r8) U200m,U200s,V200m,V200s,U850m,U850s,V850m,V850s,FLUTCm,FLUTCs
    real(r8) V200anal(cb24mjocnn_nlon,cb24mjocnn_nlat)
    real(r8) U200anal(cb24mjocnn_nlon,cb24mjocnn_nlat)
    real(r8) V850anal(cb24mjocnn_nlon,cb24mjocnn_nlat)
    real(r8) U850anal(cb24mjocnn_nlon,cb24mjocnn_nlat)
    real(r8) FLUTCanal(cb24mjocnn_nlon,cb24mjocnn_nlat)
    real(r8) model_input_sub(1,6,32,cb24mjocnn_nlon)
    real(r4) model_input_sub_single(1,6,32,cb24mjocnn_nlon)
    real(r8) Lat_anal(cb24mjocnn_nlat)
    real(r8) Lon_anal(cb24mjocnn_nlon)
    real(r8) Xtransf(1,cb24mjocnn_nlon,cb24mjocnn_nlat)
    real(r8) Xtrans(cb24mjocnn_nlon,cb24mjocnn_nlev,cb24mjocnn_nlat)
    integer  nn,Nindex
    real(r8) PI

    !+++ Ftorch vectors
    integer, parameter :: in_dims = 4
    integer, parameter :: n = 9 !for softmax function 
    integer :: in_shape(in_dims) = [1, 6, 32, 288]
    integer :: in_layout(in_dims) = [1,2,3,4]
    integer, parameter :: out_dims = 2
    integer :: out_shape(out_dims) =  [1,9]
    integer :: out_layout(out_dims) = [1,2]
    integer, parameter :: n_inputs = 1
    real(r4), dimension(:,:), allocatable, target :: out_data !ftorch
    real(r4), dimension(:), allocatable, target :: out_data_squeeze !ftorch
    real(r8), dimension(:), allocatable, target :: out_data_squeeze_double !ftorch
    
    allocate(out_data(out_shape(1), out_shape(2)))
    allocate(out_data_squeeze(out_shape(2)))
    allocate(out_data_squeeze_double(out_shape(2)))
    !--- Ftorch vectors

    U200m = 0.83579528
    U200s = 13.20837984
    V200m = -0.40176933
    V200s = 7.4425693
    U850m = -3.5611767
    U850s = 5.63111131
    V850m = 0.3286065
    V850s = 3.2882435
    FLUTCm = 280.6467022
    FLUTCs = 12.20693268

    !---
    !V
    !---
    call gather_chunk_to_field(1,cb24mjocnn_nlev,1,cb24mjocnn_nlon,cb24mjocnn_Model_V,Xtrans)
    if (masterproc) then
        do ilat=1,cb24mjocnn_nlat
        do ilon=1,cb24mjocnn_nlon
           V850anal(ilon,ilat)=Xtrans(ilon,26,ilat)
           V200anal(ilon,ilat)=Xtrans(ilon,15,ilat)
        end do
        end do
    end if 
    !---
    !U
    !---
    call gather_chunk_to_field(1,cb24mjocnn_nlev,1,cb24mjocnn_nlon,cb24mjocnn_Model_U,Xtrans)
    if (masterproc) then
        do ilat=1,cb24mjocnn_nlat
        do ilon=1,cb24mjocnn_nlon
           U850anal(ilon,ilat)=Xtrans(ilon,26,ilat)
           U200anal(ilon,ilat)=Xtrans(ilon,15,ilat)
        end do
        end do
    end if
    !---
    !FLUTC
    !---
    call gather_chunk_to_field(1,1,1,cb24mjocnn_nlon,cb24mjocnn_Model_FLUTC,Xtransf)
    if (masterproc) then
        do ilat=1,cb24mjocnn_nlat
        do ilon=1,cb24mjocnn_nlon
           FLUTCanal(ilon,ilat)=Xtransf(1,ilon,ilat)
        end do
        end do
    end if

    if (masterproc) then
        calday = get_curr_calday()
        calday = 2.0_r8*SIN(2.0_r8*3.14159265359_r8*calday/365.0_r8)
    end if
    
    !if(masterproc) write(iulog,*)'WAGO:', calday
    if (masterproc) then
        do ilat=81,112
        do ilon=1,cb24mjocnn_nlon
           !['FLUTC','U200','V200','U850','V850','DOY']
           model_input_sub(1,2,ilat-80,ilon) = (U200anal(ilon,ilat)-U200m)/U200s
           model_input_sub(1,4,ilat-80,ilon) = (U850anal(ilon,ilat)-U850m)/U850s
           model_input_sub(1,3,ilat-80,ilon) = (V200anal(ilon,ilat)-V200m)/V200s
           model_input_sub(1,5,ilat-80,ilon) = (V850anal(ilon,ilat)-V850m)/V850s
           model_input_sub(1,1,ilat-80,ilon) = (FLUTCanal(ilon,ilat)-FLUTCm)/FLUTCs
           model_input_sub(1,6,ilat-80,ilon) = calday
        end do
        end do
        !if you want to see if you've done good kid. 
        istat=nf90_create('Categorical_Input.nc',NF90_CLOBBER,ncid)
        istat=nf90_def_dim(ncid, "lon" ,288, londimid)
        istat=nf90_def_dim(ncid, "lat" , 32, latdimid )
        istat=nf90_def_dim(ncid, "var" , 6 , vardimid )
        istat=nf90_def_dim(ncid, "sing", 1 , singdimid)
        istat=nf90_def_var(ncid, "indim", nf90_double, (/singdimid, vardimid, latdimid, londimid /), rhvarid)
        istat=nf90_enddef(ncid)
        istat=nf90_put_var(ncid, rhvarid, model_input_sub) 
        istat=nf90_close(ncid)
    end if

    if (masterproc) then
        !+++ Call FTorch prediction and handle data:
        model_input_sub_single = real(model_input_sub, r4)
        !if (masterproc) write(iulog,*) 'Cowabunga!', model_input_sub_single(1,1,:,20)
        !if (masterproc) write(iulog,*) 'Cowabunga!', model_input_sub_single(1,2,:,20)
        !if (masterproc) write(iulog,*) 'Cowabunga!', model_input_sub_single(1,3,:,20)
        !if (masterproc) write(iulog,*) 'Cowabunga!', model_input_sub_single(1,4,:,20)
        !if (masterproc) write(iulog,*) 'Cowabunga!', model_input_sub_single(1,5,:,20)
        !if (masterproc) write(iulog,*) 'Cowabunga!', model_input_sub_single(1,6,:,20)
        
        call torch_tensor_from_array(in_tensor(1),model_input_sub_single, in_layout, torch_kCPU)
        call torch_tensor_from_array(out_tensor(1), out_data, out_layout, torch_kCPU)
        
        if (masterproc) write(iulog,*) 'Calling Ftorch!'!ftorch
        call torch_model_forward(model_ftorch, in_tensor, out_tensor)
        ! Loop through the vector and replace NaNs with zeros
        do ii = 1, 9
            if (isnan(out_data(1,ii))) then
                out_data(1,ii)=0.0
                out_data(1,1) = 1.0
            end if
        end do
    
        do ii = 1,9
            out_data_squeeze(ii) = out_data(1,ii)
        end do
        out_data_squeeze_double = real(out_data_squeeze, r8)
        
        ! Compute the softmax of the input_vector
        call softmax(out_data_squeeze_double, probabilities_ftorch, n)
        !if (masterproc) write(iulog,*) 'TENSOR INFER Probs! .....', probabilities_ftorch(:) !ftorch
        !if (ierr /= 0) call endrun("well this is a catastrophe")
    end if
#ifdef SPMD
        call mpibcast(probabilities_ftorch, 9, mpir8, 0, mpicom)
#endif 
    !--- Call FTorch prediction and handle data:
    
    call torch_delete(in_tensor)
    call torch_delete(out_tensor)
  end subroutine cb24mjocnn_timestep_tender


  subroutine cb24mjocnn_timestep_init(phys_state,phys_tend,pbuf,cam_in,cam_out)

    !
    ! cb24mjocnn_TIMESTEP_INIT:
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

    real(r8), pointer, dimension(:) :: flutc ! Surface meridional momentum flux
    
    !Grab Chunk and Columns
    !----------------------
    lchnk = phys_state%lchnk
    ncol  = phys_state%ncol

    itim_old = pbuf_old_tim_idx()

    flutc_idx = pbuf_get_index('FLUTC')         ! Ustar !16|33
    call pbuf_get_field(pbuf, flutc_idx,    flutc)
    
    ! Load values at Current into the Model arrays Predictors...
    !-----------------------------------------------
    ! look here: you need to update this:
    !/glade/work/wchapman/cesm/perfect_nudging/f.e21.FHIST.f09_f09_mg17_nudge_it_ERA5_1999/SourceMods/src.cam
    !nudging_calc_tend()

    cb24mjocnn_Model_V(:ncol,:pver,lchnk) = phys_state%v(:ncol,:pver)
    cb24mjocnn_Model_U(:ncol,:pver,lchnk) = phys_state%u(:ncol,:pver)
    cb24mjocnn_Model_FLUTC(:ncol,lchnk) = flutc(:ncol)

    !if(masterproc) write(iulog,*)'Flutey Pie:',flutc(:ncol)
    !if(masterproc) write(iulog,*)'Twoutey Pie:', cb24mjocnn_Model_FLUTC(:ncol,lchnk)
  end subroutine cb24mjocnn_timestep_init

  subroutine cb24mjocnn_readin_weights_all(anal_file1,anal_file2,anal_file3,anal_file4,anal_file5,anal_file6,anal_file7,anal_file8)
   !
   ! NUDGING_UPDATE_ANALYSES:
   !                 Open the given analyses data file, read in
   !                 U,V,T,Q, and PS values and then distribute
   !                 the values to all of the chunks.
   !===============================================================
   use ppgrid          ,only: pcols,pver,begchunk,endchunk
   use cam_pio_utils   ,only: cam_pio_openfile
   use pio             ,only: PIO_BCAST_ERROR,PIO_INTERNAL_ERROR,pio_inq_dimlen,pio_inq_dimid
   use pio             ,only: pio_inq_att, pio_inq_dimid, pio_inq_dimlen, pio_inq_dimname, pio_inq_vardimid, pio_inq_varid,pio_get_var
   use pio             ,only: pio_closefile,pio_seterrorhandling,file_desc_t
   use ncdio_atm       ,only: infld
   use cam_grid_support,only: cam_grid_id,cam_grid_get_dim_names,DLEN=>max_hcoordname_len
   use ioFileMod       ,only: getfil
   ! Arguments
   !-------------
   character(len=*),intent(in):: anal_file1
   character(len=*),intent(in):: anal_file2
   character(len=*),intent(in):: anal_file3
   character(len=*),intent(in):: anal_file4
   character(len=*),intent(in):: anal_file5
   character(len=*),intent(in):: anal_file6
   character(len=*),intent(in):: anal_file7
   character(len=*),intent(in):: anal_file8

   ! Local values
   !-------------
   type(file_desc_t)  :: fileID1,fileID2,fileID3,fileID4,fileID5,fileID6,fileID7,fileID8
   integer            :: nn,Nindex,pvar,c
   integer            :: dimid,dimid2,dimid3
   logical            :: VARflag
   integer            :: grid_id
   integer            :: ierr
   integer            :: err_handling1,err_handling2,err_handling3,err_handling4,err_handling5,err_handling6,err_handling7,err_handling8
   
   
   character(len=256) :: locfn1,locfn2,locfn3,locfn4,locfn5,locfn6,locfn7,locfn8
   character(len=256) :: filespec1,filespec2,filespec3,filespec4,filespec5,filespec6,filespec7,filespec8
   integer nlon,nlat,nlev,npredvar,npredvar2,npredvar3

   real(r8),allocatable:: Tmp3D(:,:,:)
   
   character(len=*), parameter :: prefix = 'nudging_update_analyses: '
   character(len=*), parameter :: sub = 'DAMLdiurnal_weights_init'

   ! Open the file and get the fileID, and define the predvar length
   !-------------------------------------
   filespec1 = trim( anal_file1 )
   call getfil( filespec1, locfn1, 0 )
   call cam_pio_openfile(fileID1,trim(locfn1), 0)
   
   call pio_seterrorhandling(fileID1,PIO_BCAST_ERROR,oldmethod=err_handling1)
   if(masterproc) write(iulog,*)'PIO_OPEN: file=',trim(anal_file1)
   
   ierr = pio_inq_dimid( fileID1, 'lat', dimid )
   ierr = pio_inq_dimlen( fileID1, dimid, nlat )
   
   ierr = pio_inq_dimid( fileID1, 'lon', dimid )
   ierr = pio_inq_dimlen( fileID1, dimid, nlon )
   
   ierr = pio_inq_dimid( fileID1, 'lev', dimid )
   ierr = pio_inq_dimlen( fileID1, dimid, nlev )
   
   allocate(Tmp3D(pcols,pver,begchunk:endchunk))
   
   if(masterproc) write(iulog,*)'pcols:',pcols,'pver:',pver   
   
   ! Get the Variables
   !-------------------------------------
   
   if(masterproc) write(iulog,*)'writing weights'       
   call infld('Nudge_U',fileID1,'lat','lev','lon', &
              1,pcols,1,pver,begchunk,endchunk,Tmp3D, &
              VARflag,gridname='physgrid',timelevel=1 )
              
   if (.not. VARflag) call endrun(sub//': ERROR: Nudge_U not found in anomaly file')
             
   MJO_Inc_01(:,:,begchunk:endchunk,1) = Tmp3D(:,:,begchunk:endchunk)

   if(masterproc) write(iulog,*)'writing weights'       
   call infld('Nudge_V',fileID1,'lat','lev','lon', &
              1,pcols,1,pver,begchunk,endchunk,Tmp3D, &
              VARflag,gridname='physgrid',timelevel=1 )
              
   if (.not. VARflag) call endrun(sub//': ERROR: Nudge_V not found in anomaly file')
             
   MJO_Inc_01(:,:,begchunk:endchunk,2) = Tmp3D(:,:,begchunk:endchunk)
   
   !-------------------------------------
   ! XXXXXXXXXXXXxxd111111
   ! Open the file and get the fileID, and define the predvar length
   !-------------------------------------
   filespec2 = trim( anal_file2 )
   call getfil( filespec2, locfn2, 0 )
   call cam_pio_openfile(fileID2,trim(locfn2), 0)
   
   call pio_seterrorhandling(fileID2,PIO_BCAST_ERROR,oldmethod=err_handling2)
   if(masterproc) write(iulog,*)'PIO_OPEN: file=',trim(anal_file2)
   
   ierr = pio_inq_dimid( fileID2, 'lat', dimid )
   ierr = pio_inq_dimlen( fileID2, dimid, nlat )
   
   ierr = pio_inq_dimid( fileID2, 'lon', dimid )
   ierr = pio_inq_dimlen( fileID2, dimid, nlon )
   
   ierr = pio_inq_dimid( fileID2, 'lev', dimid )
   ierr = pio_inq_dimlen( fileID2, dimid, nlev )
   
   !allocate(Tmp3D(pcols,pver,begchunk:endchunk))
   
   if(masterproc) write(iulog,*)'pcols:',pcols,'pver:',pver  
   ! Get the Variables
   !-------------------------------------
   
   if(masterproc) write(iulog,*)'writing weights'       
   call infld('Nudge_U',fileID2,'lat','lev','lon', &
              1,pcols,1,pver,begchunk,endchunk,Tmp3D, &
              VARflag,gridname='physgrid',timelevel=1 )
              
   if (.not. VARflag) call endrun(sub//': ERROR: Nudge_U not found in anomaly file')
             
   MJO_Inc_02(:,:,begchunk:endchunk,1) = Tmp3D(:,:,begchunk:endchunk)

   if(masterproc) write(iulog,*)'writing weights'       
   call infld('Nudge_V',fileID2,'lat','lev','lon', &
              1,pcols,1,pver,begchunk,endchunk,Tmp3D, &
              VARflag,gridname='physgrid',timelevel=1 )
              
   if (.not. VARflag) call endrun(sub//': ERROR: Nudge_V not found in anomaly file')
             
   MJO_Inc_02(:,:,begchunk:endchunk,2) = Tmp3D(:,:,begchunk:endchunk)


   !-------------------------------------
   ! XXXXXXXXXXXXxxd111111
   ! Open the file and get the fileID, and define the predvar length
   !-------------------------------------
   filespec3 = trim( anal_file3 )
   call getfil( filespec3, locfn3, 0 )
   call cam_pio_openfile(fileID3,trim(locfn3), 0)
   
   call pio_seterrorhandling(fileID3,PIO_BCAST_ERROR,oldmethod=err_handling3)
   if(masterproc) write(iulog,*)'PIO_OPEN: file=',trim(anal_file3)
   
   ierr = pio_inq_dimid( fileID3, 'lat', dimid )
   ierr = pio_inq_dimlen( fileID3, dimid, nlat )
   
   ierr = pio_inq_dimid( fileID3, 'lon', dimid )
   ierr = pio_inq_dimlen( fileID3, dimid, nlon )
   
   ierr = pio_inq_dimid( fileID3, 'lev', dimid )
   ierr = pio_inq_dimlen( fileID3, dimid, nlev )
   
   !allocate(Tmp3D(pcols,pver,begchunk:endchunk))
   
   if(masterproc) write(iulog,*)'pcols:',pcols,'pver:',pver   
   ! Get the Variables
   !-------------------------------------
   
   if(masterproc) write(iulog,*)'writing weights'       
   call infld('Nudge_U',fileID3,'lat','lev','lon', &
              1,pcols,1,pver,begchunk,endchunk,Tmp3D, &
              VARflag,gridname='physgrid',timelevel=1 )
              
   if (.not. VARflag) call endrun(sub//': ERROR: Nudge_U not found in anomaly file')
             
   MJO_Inc_03(:,:,begchunk:endchunk,1) = Tmp3D(:,:,begchunk:endchunk)

   if(masterproc) write(iulog,*)'writing weights'       
   call infld('Nudge_V',fileID3,'lat','lev','lon', &
              1,pcols,1,pver,begchunk,endchunk,Tmp3D, &
              VARflag,gridname='physgrid',timelevel=1 )
              
   if (.not. VARflag) call endrun(sub//': ERROR: Nudge_V not found in anomaly file')
             
   MJO_Inc_03(:,:,begchunk:endchunk,2) = Tmp3D(:,:,begchunk:endchunk)


   !-------------------------------------
   ! XXXXXXXXXXXXxxd111111
   ! Open the file and get the fileID, and define the predvar length
   !-------------------------------------
   filespec4 = trim( anal_file4 )
   call getfil( filespec4, locfn4, 0 )
   call cam_pio_openfile(fileID4,trim(locfn4), 0)
   
   call pio_seterrorhandling(fileID4,PIO_BCAST_ERROR,oldmethod=err_handling4)
   if(masterproc) write(iulog,*)'PIO_OPEN: file=',trim(anal_file4)
   
   ierr = pio_inq_dimid( fileID4, 'lat', dimid )
   ierr = pio_inq_dimlen( fileID4, dimid, nlat )
   
   ierr = pio_inq_dimid( fileID4, 'lon', dimid )
   ierr = pio_inq_dimlen( fileID4, dimid, nlon )
   
   ierr = pio_inq_dimid( fileID4, 'lev', dimid )
   ierr = pio_inq_dimlen( fileID4, dimid, nlev )
   
   !allocate(Tmp3D(pcols,pver,begchunk:endchunk))
   
   if(masterproc) write(iulog,*)'pcols:',pcols,'pver:',pver  
   ! Get the Variables
   !-------------------------------------
   
   if(masterproc) write(iulog,*)'writing weights'       
   call infld('Nudge_U',fileID4,'lat','lev','lon', &
              1,pcols,1,pver,begchunk,endchunk,Tmp3D, &
              VARflag,gridname='physgrid',timelevel=1 )
              
   if (.not. VARflag) call endrun(sub//': ERROR: Nudge_U not found in anomaly file')
             
   MJO_Inc_04(:,:,begchunk:endchunk,1) = Tmp3D(:,:,begchunk:endchunk)

   if(masterproc) write(iulog,*)'writing weights'       
   call infld('Nudge_V',fileID4,'lat','lev','lon', &
              1,pcols,1,pver,begchunk,endchunk,Tmp3D, &
              VARflag,gridname='physgrid',timelevel=1 )
              
   if (.not. VARflag) call endrun(sub//': ERROR: Nudge_V not found in anomaly file')
             
   MJO_Inc_04(:,:,begchunk:endchunk,2) = Tmp3D(:,:,begchunk:endchunk)


   !-------------------------------------
   ! XXXXXXXXXXXXxxd111111
   ! Open the file and get the fileID, and define the predvar length
   !-------------------------------------
   filespec5 = trim( anal_file5 )
   call getfil( filespec5, locfn5, 0 )
   call cam_pio_openfile(fileID5,trim(locfn5), 0)
   
   call pio_seterrorhandling(fileID5,PIO_BCAST_ERROR,oldmethod=err_handling5)
   if(masterproc) write(iulog,*)'PIO_OPEN: file=',trim(anal_file5)
   
   ierr = pio_inq_dimid( fileID5, 'lat', dimid )
   ierr = pio_inq_dimlen( fileID5, dimid, nlat )
   
   ierr = pio_inq_dimid( fileID5, 'lon', dimid )
   ierr = pio_inq_dimlen( fileID5, dimid, nlon )
   
   ierr = pio_inq_dimid( fileID5, 'lev', dimid )
   ierr = pio_inq_dimlen( fileID5, dimid, nlev )
   
   !allocate(Tmp3D(pcols,pver,begchunk:endchunk))
   
   if(masterproc) write(iulog,*)'pcols:',pcols,'pver:',pver  
   ! Get the Variables
   !-------------------------------------
   
   if(masterproc) write(iulog,*)'writing weights'       
   call infld('Nudge_U',fileID5,'lat','lev','lon', &
              1,pcols,1,pver,begchunk,endchunk,Tmp3D, &
              VARflag,gridname='physgrid',timelevel=1 )
              
   if (.not. VARflag) call endrun(sub//': ERROR: Nudge_U not found in anomaly file')
             
   MJO_Inc_05(:,:,begchunk:endchunk,1) = Tmp3D(:,:,begchunk:endchunk)

   if(masterproc) write(iulog,*)'writing weights'       
   call infld('Nudge_V',fileID5,'lat','lev','lon', &
              1,pcols,1,pver,begchunk,endchunk,Tmp3D, &
              VARflag,gridname='physgrid',timelevel=1 )
              
   if (.not. VARflag) call endrun(sub//': ERROR: Nudge_V not found in anomaly file')
             
   MJO_Inc_05(:,:,begchunk:endchunk,2) = Tmp3D(:,:,begchunk:endchunk)

   !-------------------------------------
   ! XXXXXXXXXXXXxxd111111
   ! Open the file and get the fileID, and define the predvar length
   !-------------------------------------
   filespec6 = trim( anal_file6 )
   call getfil( filespec6, locfn6, 0 )
   call cam_pio_openfile(fileID6,trim(locfn6), 0)
   
   call pio_seterrorhandling(fileID6,PIO_BCAST_ERROR,oldmethod=err_handling6)
   if(masterproc) write(iulog,*)'PIO_OPEN: file=',trim(anal_file6)
   
   ierr = pio_inq_dimid( fileID6, 'lat', dimid )
   ierr = pio_inq_dimlen( fileID6, dimid, nlat )
   
   ierr = pio_inq_dimid( fileID6, 'lon', dimid )
   ierr = pio_inq_dimlen( fileID6, dimid, nlon )
   
   ierr = pio_inq_dimid( fileID6, 'lev', dimid )
   ierr = pio_inq_dimlen( fileID6, dimid, nlev )
   
   !allocate(Tmp3D(pcols,pver,begchunk:endchunk))
   
   if(masterproc) write(iulog,*)'pcols:',pcols,'pver:',pver 
   ! Get the Variables
   !-------------------------------------
   
   if(masterproc) write(iulog,*)'writing weights'       
   call infld('Nudge_U',fileID6,'lat','lev','lon', &
              1,pcols,1,pver,begchunk,endchunk,Tmp3D, &
              VARflag,gridname='physgrid',timelevel=1 )
              
   if (.not. VARflag) call endrun(sub//': ERROR: Nudge_U not found in anomaly file')
             
   MJO_Inc_06(:,:,begchunk:endchunk,1) = Tmp3D(:,:,begchunk:endchunk)

   if(masterproc) write(iulog,*)'writing weights'       
   call infld('Nudge_V',fileID6,'lat','lev','lon', &
              1,pcols,1,pver,begchunk,endchunk,Tmp3D, &
              VARflag,gridname='physgrid',timelevel=1 )
              
   if (.not. VARflag) call endrun(sub//': ERROR: Nudge_V not found in anomaly file')
             
   MJO_Inc_06(:,:,begchunk:endchunk,2) = Tmp3D(:,:,begchunk:endchunk)


   !-------------------------------------
   ! XXXXXXXXXXXXxxd111111
   ! Open the file and get the fileID, and define the predvar length
   !-------------------------------------
   filespec7 = trim( anal_file7 )
   call getfil( filespec7, locfn7, 0 )
   call cam_pio_openfile(fileID7,trim(locfn7), 0)
   
   call pio_seterrorhandling(fileID7,PIO_BCAST_ERROR,oldmethod=err_handling7)
   if(masterproc) write(iulog,*)'PIO_OPEN: file=',trim(anal_file7)
   
   ierr = pio_inq_dimid( fileID7, 'lat', dimid )
   ierr = pio_inq_dimlen( fileID7, dimid, nlat )
   
   ierr = pio_inq_dimid( fileID7, 'lon', dimid )
   ierr = pio_inq_dimlen( fileID7, dimid, nlon )
   
   ierr = pio_inq_dimid( fileID7, 'lev', dimid )
   ierr = pio_inq_dimlen( fileID7, dimid, nlev )
   
   !allocate(Tmp3D(pcols,pver,begchunk:endchunk))
   
   if(masterproc) write(iulog,*)'pcols:',pcols,'pver:',pver 
   ! Get the Variables
   !-------------------------------------
   
   if(masterproc) write(iulog,*)'writing weights'       
   call infld('Nudge_U',fileID7,'lat','lev','lon', &
              1,pcols,1,pver,begchunk,endchunk,Tmp3D, &
              VARflag,gridname='physgrid',timelevel=1 )
              
   if (.not. VARflag) call endrun(sub//': ERROR: Nudge_U not found in anomaly file')
             
   MJO_Inc_07(:,:,begchunk:endchunk,1) = Tmp3D(:,:,begchunk:endchunk)

   if(masterproc) write(iulog,*)'writing weights'       
   call infld('Nudge_V',fileID7,'lat','lev','lon', &
              1,pcols,1,pver,begchunk,endchunk,Tmp3D, &
              VARflag,gridname='physgrid',timelevel=1 )
              
   if (.not. VARflag) call endrun(sub//': ERROR: Nudge_V not found in anomaly file')
             
   MJO_Inc_07(:,:,begchunk:endchunk,2) = Tmp3D(:,:,begchunk:endchunk)


   !-------------------------------------
   ! XXXXXXXXXXXXxxd111111
   ! Open the file and get the fileID, and define the predvar length
   !-------------------------------------
   filespec8 = trim( anal_file8 )
   call getfil( filespec8, locfn8, 0 )
   call cam_pio_openfile(fileID8,trim(locfn8), 0)
   
   call pio_seterrorhandling(fileID8,PIO_BCAST_ERROR,oldmethod=err_handling8)
   if(masterproc) write(iulog,*)'PIO_OPEN: file=',trim(anal_file8)
   
   ierr = pio_inq_dimid( fileID8, 'lat', dimid )
   ierr = pio_inq_dimlen( fileID8, dimid, nlat )
   
   ierr = pio_inq_dimid( fileID8, 'lon', dimid )
   ierr = pio_inq_dimlen( fileID8, dimid, nlon )
   
   ierr = pio_inq_dimid( fileID8, 'lev', dimid )
   ierr = pio_inq_dimlen( fileID8, dimid, nlev )
   
   !allocate(Tmp3D(pcols,pver,begchunk:endchunk))
   
   if(masterproc) write(iulog,*)'pcols:',pcols,'pver:',pver  
   ! Get the Variables
   !-------------------------------------
   
   if(masterproc) write(iulog,*)'writing weights'       
   call infld('Nudge_U',fileID8,'lat','lev','lon', &
              1,pcols,1,pver,begchunk,endchunk,Tmp3D, &
              VARflag,gridname='physgrid',timelevel=1 )
              
   if (.not. VARflag) call endrun(sub//': ERROR: Nudge_U not found in anomaly file')
             
   MJO_Inc_08(:,:,begchunk:endchunk,1) = Tmp3D(:,:,begchunk:endchunk)

   if(masterproc) write(iulog,*)'writing weights'       
   call infld('Nudge_V',fileID8,'lat','lev','lon', &
              1,pcols,1,pver,begchunk,endchunk,Tmp3D, &
              VARflag,gridname='physgrid',timelevel=1 )
              
   if (.not. VARflag) call endrun(sub//': ERROR: Nudge_V not found in anomaly file')
             
   MJO_Inc_08(:,:,begchunk:endchunk,2) = Tmp3D(:,:,begchunk:endchunk) 
   
   ! Restore old error handling
   !----------------------------
   call pio_seterrorhandling(fileID1,err_handling1)
   call pio_seterrorhandling(fileID2,err_handling2)
   call pio_seterrorhandling(fileID3,err_handling3)
   call pio_seterrorhandling(fileID4,err_handling4)
   call pio_seterrorhandling(fileID5,err_handling5)
   call pio_seterrorhandling(fileID6,err_handling6)
   call pio_seterrorhandling(fileID7,err_handling7)
   call pio_seterrorhandling(fileID8,err_handling8)

   ! Close the analyses file
   !-----------------------
   deallocate(Tmp3D)
   call pio_closefile(fileID1)
   call pio_closefile(fileID2)
   call pio_closefile(fileID3)
   call pio_closefile(fileID4)
   call pio_closefile(fileID5)
   call pio_closefile(fileID6)
   call pio_closefile(fileID7)
   call pio_closefile(fileID8)

   ! End Routine
   !------------

  end subroutine cb24mjocnn_readin_weights_all

    ! Function to compute the softmax of a vector
    subroutine softmax(vector, softmax_vector, n)
        integer :: n
        real(r8), intent(in) :: vector(n)
        real(r8), intent(out) :: softmax_vector(n)
        real(r8) :: sum_exp, max_val
        integer :: i
    
        ! Numerical stability improvement by subtracting the max value
        max_val = maxval(vector)
        sum_exp = 0.0
    
        ! Compute the sum of the exponentials
        do i = 1, n
            sum_exp = sum_exp + exp(vector(i) - max_val)
        end do
    
        ! Compute the softmax values
        do i = 1, n
            softmax_vector(i) = exp(vector(i) - max_val) / sum_exp
        end do
    end subroutine softmax

  !================================================================
  subroutine cb24mjocnn_timestep_tend(phys_state,phys_tend)
   !
   ! NUDGING_TIMESTEP_TEND:
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
   integer indw,ncol,lchnk,i,k
   logical lq(pcnst)
   real(r8), dimension(:), allocatable :: probabilities_ftorch_double !ftorch
   allocate(probabilities_ftorch_double(9))

   call cnst_get_ind('Q',indw)
   lq(:)   =.false.
   lq(indw)=.true.
   call physics_ptend_init(phys_tend,phys_state%psetcols,'nudging',lu=.true.,lv=.true.,ls=.true.,lq=lq)

   if (masterproc) write(iulog,*) 'Tstep Tend Probs!:', probabilities_ftorch(:) !ftorch
   if (masterproc) write(iulog,*) 'Tstep Tend Probs Eight!:', probabilities_ftorch(9) !ftorch

   probabilities_ftorch_double = real(probabilities_ftorch, r8)

   if (masterproc) write(iulog,*) 'Tstep Tend Probs_double!:', probabilities_ftorch_double(:) !ftorch
   if (masterproc) write(iulog,*) 'Tstep Tend Probs Eight_double!:', probabilities_ftorch_double(9) !ftorch

   lchnk=phys_state%lchnk
   ncol =phys_state%ncol
   
   do i=1,ncol
      do k=1,pver
          phys_tend%u(i,k) = MJO_Inc_00(i,k,lchnk,1)*probabilities_ftorch(1) + &
                             MJO_Inc_01(i,k,lchnk,1)*probabilities_ftorch(2) + &
                             MJO_Inc_02(i,k,lchnk,1)*probabilities_ftorch(3) + &
                             MJO_Inc_03(i,k,lchnk,1)*probabilities_ftorch(4) + &
                             MJO_Inc_04(i,k,lchnk,1)*probabilities_ftorch(5) + &
                             MJO_Inc_05(i,k,lchnk,1)*probabilities_ftorch(6) + &
                             MJO_Inc_06(i,k,lchnk,1)*probabilities_ftorch(7) + &
                             MJO_Inc_07(i,k,lchnk,1)*probabilities_ftorch(8) + &
                             MJO_Inc_08(i,k,lchnk,1)*probabilities_ftorch(9)
                              
          phys_tend%v(i,k) = MJO_Inc_00(i,k,lchnk,2)*probabilities_ftorch(1) + &
                             MJO_Inc_01(i,k,lchnk,2)*probabilities_ftorch(2) + &
                             MJO_Inc_02(i,k,lchnk,2)*probabilities_ftorch(3) + &
                             MJO_Inc_03(i,k,lchnk,2)*probabilities_ftorch(4) + &
                             MJO_Inc_04(i,k,lchnk,2)*probabilities_ftorch(5) + &
                             MJO_Inc_05(i,k,lchnk,2)*probabilities_ftorch(6) + &
                             MJO_Inc_06(i,k,lchnk,2)*probabilities_ftorch(7) + &
                             MJO_Inc_07(i,k,lchnk,2)*probabilities_ftorch(8) + &
                             MJO_Inc_08(i,k,lchnk,2)*probabilities_ftorch(9)
      end do
   end do 

   do i=1,ncol
      do k=1,9
          cb24mjocnn_Probs(i,k,lchnk) = probabilities_ftorch_double(k)
      end do
   end do   

   call outfld( 'cb24mjocnn_U',phys_tend%u                ,pcols,lchnk)
   call outfld( 'cb24mjocnn_V',phys_tend%v                ,pcols,lchnk)
   call outfld( 'cb24mjocnn_T',phys_tend%s/cpair          ,pcols,lchnk) !MUST EDIT TO MESS WITH T->DSE [WEC]
   call outfld( 'cb24mjocnn_Q',phys_tend%q(1,1,indw)      ,pcols,lchnk)
   call outfld( 'state_mjo'   ,cb24mjocnn_Probs,pcols,lchnk)
     
   ! End Routine
   !------------
   return
  end subroutine ! nudging_timestep_tend

  subroutine cb24mjocnn_finalize()
    ! Cleanup
    call torch_delete(model_ftorch)
  end subroutine cb24mjocnn_finalize

  end module cb24mjocnn
