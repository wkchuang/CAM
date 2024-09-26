#define MMF_NN_EMULATOR
module physpkg
  !-----------------------------------------------------------------------
  ! Purpose:
  !
  ! Provides the interface to CAM physics package
  !
  ! Revision history:
  ! Aug  2005,  E. B. Kluzek,  Creation of module from physpkg subroutine
  ! 2005-10-17  B. Eaton       Add contents of inti.F90 to phys_init().  Add
  !                            initialization of grid info in phys_state.
  ! Nov 2010    A. Gettelman   Put micro/macro physics into separate routines
  ! Aug 2024    W. Chuang      Add in SPCAM emulator adapted from E3SM
  !   For the emulator: 
  !   Each parameterization should be implemented with this sequence of calls:
  !    1)  call the physics routine to calculate tendencies
  !    2)  call physics_update_ml() to apply tendencies from ptend to the state
  !    3)  call check_energy_chng() to ensure that the energy and water 
  !        changes match the boundary fluxes
  !-----------------------------------------------------------------------

  use cam_abortutils,   only: endrun
  use shr_kind_mod,     only: i8 => shr_kind_i8, r8 => shr_kind_r8
  use shr_sys_mod,      only: shr_sys_irtc, shr_sys_flush
  use spmd_utils,       only: masterproc, iam
  use physconst,        only: latvap, latice, rh2o
  use physics_types,    only: physics_state, physics_tend, physics_state_set_grid, &
        physics_ptend, &
        physics_tend_init, physics_ptend_dealloc, &
        physics_update, &
        physics_type_alloc, &
        physics_state_alloc, physics_state_dealloc, physics_tend_alloc, physics_tend_dealloc, &
        physics_state_copy
  use physics_update_mod, only: physics_update_ml, physics_update_init, hist_vars, nvars_prtrb_hist, get_var
  use phys_grid,        only: get_ncols_p, print_cost_p, update_cost_p
  use phys_gmean,       only: gmean_mass
  use ppgrid,           only: begchunk, endchunk, pcols, pver, pverp, psubcols
  use constituents,     only: pcnst, cnst_name, cnst_get_ind, cnst_name, setup_moist_indices
  use camsrfexch,       only: cam_out_t, cam_in_t

  use cam_control_mod,  only: ideal_phys, adiabatic
  use phys_control,     only: phys_do_flux_avg, phys_getopts, waccmx_is
  use scamMod,          only: single_column, scm_crm_mode
  use flux_avg,         only: flux_avg_init
  use perf_mod
  use cam_logfile,     only: iulog
  use camsrfexch,      only: cam_export

  use modal_aero_calcsize,    only: modal_aero_calcsize_init, modal_aero_calcsize_diag, modal_aero_calcsize_reg
  use modal_aero_wateruptake, only: modal_aero_wateruptake_init, modal_aero_wateruptake_dr, modal_aero_wateruptake_reg

  implicit none
  private
  save

  ! Public methods
  public phys_register ! was initindx  - register physics methods
  public phys_init   ! Public initialization method
  public phys_run1   ! First phase of the public run method
  public phys_run2   ! Second phase of the public run method
  public phys_final  ! Public finalization method
#ifdef MMF_NN_EMULATOR
  public mmf_nn_emulator_driver ! MMF_NN_EMULATOR NN-emulation driver
#endif
  ! Private module data

  ! Physics package options
  character(len=16) :: shallow_scheme
  character(len=16) :: macrop_scheme
  character(len=16) :: microp_scheme
  integer           :: cld_macmic_num_steps    ! Number of macro/micro substeps
  logical           :: do_clubb_sgs
  logical           :: use_subcol_microp   ! if true, use subcolumns in microphysics
  logical           :: state_debug_checks  ! Debug physics_state.
  logical           :: clim_modal_aero     ! climate controled by prognostic or prescribed modal aerosols
  logical           :: prog_modal_aero     ! Prognostic modal aerosols present

  !  Physics buffer index
  integer ::  teout_idx          = 0

  integer ::  landm_idx          = 0
  integer ::  sgh_idx            = 0
  integer ::  sgh30_idx          = 0

  integer ::  qini_idx           = 0
  integer ::  cldliqini_idx      = 0
  integer ::  cldiceini_idx      = 0

  integer ::  prec_str_idx       = 0
  integer ::  snow_str_idx       = 0
  integer ::  prec_sed_idx       = 0
  integer ::  snow_sed_idx       = 0
  integer ::  prec_pcw_idx       = 0
  integer ::  snow_pcw_idx       = 0
  integer ::  prec_dp_idx        = 0
  integer ::  snow_dp_idx        = 0
  integer ::  prec_sh_idx        = 0
  integer ::  snow_sh_idx        = 0
  integer ::  dlfzm_idx          = 0     ! detrained convective cloud water mixing ratio.

!=======================================================================
contains
!=======================================================================

  subroutine phys_register
    !-----------------------------------------------------------------------
    !
    ! Purpose: Register constituents and physics buffer fields.
    !
    ! Author:    CSM Contact: M. Vertenstein, Aug. 1997
    !            B.A. Boville, Oct 2001
    !            A. Gettelman, Nov 2010 - put micro/macro physics into separate routines
    !
    !-----------------------------------------------------------------------
    use cam_abortutils,     only: endrun
    use physics_buffer,     only: pbuf_init_time
    use physics_buffer,     only: pbuf_add_field, dtype_r8, pbuf_register_subcol
    use shr_kind_mod,       only: r8 => shr_kind_r8
    use spmd_utils,         only: masterproc
    use constituents,       only: pcnst, cnst_add, cnst_chk_dim, cnst_name

    use cam_control_mod,    only: moist_physics
    use chemistry,          only: chem_register
    use cloud_fraction,     only: cldfrc_register
    use rk_stratiform,      only: rk_stratiform_register
    use microp_driver,      only: microp_driver_register
    use microp_aero,        only: microp_aero_register
    use macrop_driver,      only: macrop_driver_register
    use clubb_intr,         only: clubb_register_cam
    use conv_water,         only: conv_water_register
    use physconst,          only: mwdry, cpair, mwh2o, cpwv
    use tracers,            only: tracers_register
    use check_energy,       only: check_energy_register
    use carma_intr,         only: carma_register
    use cam3_aero_data,     only: cam3_aero_data_on, cam3_aero_data_register
    use cam3_ozone_data,    only: cam3_ozone_data_on, cam3_ozone_data_register
    use ghg_data,           only: ghg_data_register
    use vertical_diffusion, only: vd_register
    use convect_deep,       only: convect_deep_register
    use convect_shallow,    only: convect_shallow_register
    use radiation,          only: radiation_register
    use co2_cycle,          only: co2_register
    use flux_avg,           only: flux_avg_register
    use iondrag,            only: iondrag_register
    use waccmx_phys_intr,   only: waccmx_phys_ion_elec_temp_reg
    use string_utils,       only: to_lower
    use prescribed_ozone,   only: prescribed_ozone_register
    use prescribed_volcaero,only: prescribed_volcaero_register
    use prescribed_strataero,only: prescribed_strataero_register
    use prescribed_aero,    only: prescribed_aero_register
    use prescribed_ghg,     only: prescribed_ghg_register
    use sslt_rebin,         only: sslt_rebin_register
    use aoa_tracers,        only: aoa_tracers_register
    use aircraft_emit,      only: aircraft_emit_register
    use cam_diagnostics,    only: diag_register
    use cloud_diagnostics,  only: cloud_diagnostics_register
    use cospsimulator_intr, only: cospsimulator_intr_register
    use rad_constituents,   only: rad_cnst_get_info ! Added to query if it is a modal aero sim or not
    use subcol,             only: subcol_register
    use subcol_utils,       only: is_subcol_on
    use dyn_comp,           only: dyn_register
    use spcam_drivers,      only: spcam_register
    use offline_driver,     only: offline_driver_reg
    use upper_bc,           only: ubc_fixed_conc
    use cb24mjocnn,         only: cb24mjocnn_register_cam

    !---------------------------Local variables-----------------------------
    !
    integer  :: m        ! loop index
    integer  :: mm       ! constituent index
    integer  :: nmodes
    logical  :: has_fixed_ubc ! for upper bndy cond
    !-----------------------------------------------------------------------

    ! Get physics options
    call phys_getopts(shallow_scheme_out       = shallow_scheme, &
                      macrop_scheme_out        = macrop_scheme,   &
                      microp_scheme_out        = microp_scheme,   &
                      cld_macmic_num_steps_out = cld_macmic_num_steps, &
                      do_clubb_sgs_out         = do_clubb_sgs,     &
                      use_subcol_microp_out    = use_subcol_microp, &
                      state_debug_checks_out   = state_debug_checks)

    ! Initialize dyn_time_lvls
    call pbuf_init_time()

    ! Register the subcol scheme
    call subcol_register()

    ! Register water vapor.
    ! ***** N.B. ***** This must be the first call to cnst_add so that
    !                  water vapor is constituent 1.
    has_fixed_ubc = ubc_fixed_conc('Q') ! .false.
    if (moist_physics) then
       call cnst_add('Q', mwh2o, cpwv, 1.E-12_r8, mm, fixed_ubc=has_fixed_ubc, &
            longname='Specific humidity', readiv=.true., is_convtran1=.true.)
    else
       call cnst_add('Q', mwh2o, cpwv, 0.0_r8, mm, fixed_ubc=has_fixed_ubc, &
            longname='Specific humidity', readiv=.false., is_convtran1=.true.)
    end if

    ! Topography file fields.
    call pbuf_add_field('LANDM',     'global',  dtype_r8, (/pcols/),      landm_idx)
    call pbuf_add_field('SGH',       'global',  dtype_r8, (/pcols/),      sgh_idx)
    call pbuf_add_field('SGH30',     'global',  dtype_r8, (/pcols/),      sgh30_idx)

    ! Fields for physics package diagnostics
    call pbuf_add_field('QINI',      'physpkg', dtype_r8, (/pcols,pver/), qini_idx)
    call pbuf_add_field('CLDLIQINI', 'physpkg', dtype_r8, (/pcols,pver/), cldliqini_idx)
    call pbuf_add_field('CLDICEINI', 'physpkg', dtype_r8, (/pcols,pver/), cldiceini_idx)

    ! check energy package
    call check_energy_register

    ! If using a simple physics option (e.g., held_suarez, adiabatic),
    ! the normal CAM physics parameterizations are not called.
    if (moist_physics) then

       ! register fluxes for saving across time
       if (phys_do_flux_avg()) call flux_avg_register()

       call cldfrc_register()

       ! cloud water
       if( microp_scheme == 'RK' ) then
          call rk_stratiform_register()
       elseif( microp_scheme == 'MG' ) then
          if (.not. do_clubb_sgs) call macrop_driver_register()
          call microp_aero_register()
          call microp_driver_register()
       end if

       ! Register CLUBB_SGS here
       if (do_clubb_sgs) call clubb_register_cam()


       call pbuf_add_field('PREC_STR',  'physpkg',dtype_r8,(/pcols/),prec_str_idx)
       call pbuf_add_field('SNOW_STR',  'physpkg',dtype_r8,(/pcols/),snow_str_idx)
       call pbuf_add_field('PREC_PCW',  'physpkg',dtype_r8,(/pcols/),prec_pcw_idx)
       call pbuf_add_field('SNOW_PCW',  'physpkg',dtype_r8,(/pcols/),snow_pcw_idx)
       call pbuf_add_field('PREC_SED',  'physpkg',dtype_r8,(/pcols/),prec_sed_idx)
       call pbuf_add_field('SNOW_SED',  'physpkg',dtype_r8,(/pcols/),snow_sed_idx)
       if (is_subcol_on()) then
         call pbuf_register_subcol('PREC_STR', 'phys_register', prec_str_idx)
         call pbuf_register_subcol('SNOW_STR', 'phys_register', snow_str_idx)
         call pbuf_register_subcol('PREC_PCW', 'phys_register', prec_pcw_idx)
         call pbuf_register_subcol('SNOW_PCW', 'phys_register', snow_pcw_idx)
         call pbuf_register_subcol('PREC_SED', 'phys_register', prec_sed_idx)
         call pbuf_register_subcol('SNOW_SED', 'phys_register', snow_sed_idx)
       end if

    ! Who should add FRACIS?
    ! -- It does not seem that aero_intr should add it since FRACIS is used in convection
    !     even if there are no prognostic aerosols ... so do it here for now
       call pbuf_add_field('FRACIS','physpkg',dtype_r8,(/pcols,pver,pcnst/),m)

       call conv_water_register()

       ! Determine whether its a 'modal' aerosol simulation  or not
       call rad_cnst_get_info(0, nmodes=nmodes)
       clim_modal_aero = (nmodes > 0)

       if (clim_modal_aero) then
          call modal_aero_calcsize_reg()
          call modal_aero_wateruptake_reg()
       endif

       ! register chemical constituents including aerosols ...
       call chem_register()

       ! co2 constituents
       call co2_register()

       ! register data model ozone with pbuf
       if (cam3_ozone_data_on) then
          call cam3_ozone_data_register()
       end if
       call prescribed_volcaero_register()
       call prescribed_strataero_register()
       call prescribed_ozone_register()
       call prescribed_aero_register()
       call prescribed_ghg_register()
       call sslt_rebin_register

       ! CAM3 prescribed aerosols
       if (cam3_aero_data_on) then
          call cam3_aero_data_register()
       end if

       ! register various data model gasses with pbuf
       call ghg_data_register()

       ! carma microphysics
       !
       call carma_register()

       ! Register iondrag variables with pbuf
       call iondrag_register()

       ! Register ionosphere variables with pbuf if mode set to ionosphere
       if( waccmx_is('ionosphere') ) then
          call waccmx_phys_ion_elec_temp_reg()
       endif

       call aircraft_emit_register()

       ! deep convection
       call convect_deep_register

       !  shallow convection
       call convect_shallow_register


       call spcam_register

       ! radiation
       call radiation_register
       call cloud_diagnostics_register

       ! COSP
       call cospsimulator_intr_register

       ! vertical diffusion
       call vd_register()
    else
       ! held_suarez/adiabatic physics option should be in simple_physics
       call endrun('phys_register: moist_physics configuration error')
    end if

    ! Register diagnostics PBUF
    call diag_register()

    ! Register age of air tracers
    call aoa_tracers_register()

    ! Register test tracers
    call tracers_register()

    call dyn_register()

    ! All tracers registered, check that the dimensions are correct
    call cnst_chk_dim()

    ! ***NOTE*** No registering constituents after the call to cnst_chk_dim.

    call offline_driver_reg()
    !call cb24mjocnn_register_cam()

  end subroutine phys_register



  !=======================================================================

  subroutine phys_inidat( cam_out, pbuf2d )
    use cam_abortutils,      only: endrun

    use physics_buffer,      only: pbuf_get_index, pbuf_get_field, physics_buffer_desc, pbuf_set_field, dyn_time_lvls


    use cam_initfiles,       only: initial_file_get_id, topo_file_get_id
    use cam_grid_support,    only: cam_grid_check, cam_grid_id
    use cam_grid_support,    only: cam_grid_get_dim_names
    use pio,                 only: file_desc_t
    use ncdio_atm,           only: infld
    use dycore,              only: dycore_is
    use polar_avg,           only: polar_average
    use short_lived_species, only: initialize_short_lived_species
    use cam_control_mod,     only: aqua_planet
    use waccmx_phys_intr,    only: waccmx_phys_ion_elec_temp_inidat

    type(cam_out_t),     intent(inout) :: cam_out(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    integer :: lchnk, m, n, i, k, ncol
    type(file_desc_t), pointer :: fh_ini, fh_topo
    character(len=8) :: fieldname
    real(r8), pointer :: tptr(:,:), tptr_2(:,:), tptr3d(:,:,:), tptr3d_2(:,:,:)
    real(r8), pointer :: qpert(:,:)

    character(len=11) :: subname='phys_inidat' ! subroutine name
    integer :: tpert_idx, qpert_idx, pblh_idx

    logical :: found=.false., found2=.false.
    integer :: ierr
    character(len=8) :: dim1name, dim2name
    integer :: ixcldice, ixcldliq
    integer                   :: grid_id  ! grid ID for data mapping

    nullify(tptr,tptr_2,tptr3d,tptr3d_2)

    fh_ini  => initial_file_get_id()
    fh_topo => topo_file_get_id()

    !   dynamics variables are handled in dyn_init - here we read variables needed for physics
    !   but not dynamics

    grid_id = cam_grid_id('physgrid')
    if (.not. cam_grid_check(grid_id)) then
      call endrun(trim(subname)//': Internal error, no "physgrid" grid')
    end if
    call cam_grid_get_dim_names(grid_id, dim1name, dim2name)

    allocate(tptr(1:pcols,begchunk:endchunk))

    if (associated(fh_topo) .and. .not. aqua_planet) then
      call infld('SGH', fh_topo, dim1name, dim2name, 1, pcols, begchunk, endchunk, &
           tptr, found, gridname='physgrid')
      if(.not. found) call endrun('ERROR: SGH not found on topo file')

      call pbuf_set_field(pbuf2d, sgh_idx, tptr)

      allocate(tptr_2(1:pcols,begchunk:endchunk))
      call infld('SGH30', fh_topo, dim1name, dim2name, 1, pcols, begchunk, endchunk, &
           tptr_2, found, gridname='physgrid')
      if(found) then
        call pbuf_set_field(pbuf2d, sgh30_idx, tptr_2)
      else
        if (masterproc) write(iulog,*) 'Warning: Error reading SGH30 from topo file.'
        if (masterproc) write(iulog,*) 'The field SGH30 will be filled using data from SGH.'
        call pbuf_set_field(pbuf2d, sgh30_idx, tptr)
      end if

      deallocate(tptr_2)

      call infld('LANDM_COSLAT', fh_topo, dim1name, dim2name, 1, pcols, begchunk, endchunk, &
           tptr, found, gridname='physgrid')

      if(.not.found) call endrun(' ERROR: LANDM_COSLAT not found on topo dataset.')

      call pbuf_set_field(pbuf2d, landm_idx, tptr)

    else
      call pbuf_set_field(pbuf2d, sgh_idx, 0._r8)
      call pbuf_set_field(pbuf2d, sgh30_idx, 0._r8)
      call pbuf_set_field(pbuf2d, landm_idx, 0._r8)
    end if

    call infld('PBLH', fh_ini, dim1name, dim2name, 1, pcols, begchunk, endchunk, &
         tptr(:,:), found, gridname='physgrid')
    if(.not. found) then
       tptr(:,:) = 0._r8
       if (masterproc) write(iulog,*) 'PBLH initialized to 0.'
    end if
    pblh_idx = pbuf_get_index('pblh')

    call pbuf_set_field(pbuf2d, pblh_idx, tptr)

    call infld('TPERT', fh_ini, dim1name, dim2name, 1, pcols, begchunk, endchunk, &
         tptr(:,:), found, gridname='physgrid')
    if(.not. found) then
       tptr(:,:) = 0._r8
       if (masterproc) write(iulog,*) 'TPERT initialized to 0.'
    end if
    tpert_idx = pbuf_get_index( 'tpert')
    call pbuf_set_field(pbuf2d, tpert_idx, tptr)

    fieldname='QPERT'
    qpert_idx = pbuf_get_index( 'qpert',ierr)
    if (qpert_idx > 0) then
       call infld(fieldname, fh_ini, dim1name, dim2name, 1, pcols, begchunk, endchunk, &
            tptr, found, gridname='physgrid')
       if(.not. found) then
          tptr=0_r8
          if (masterproc) write(iulog,*) trim(fieldname), ' initialized to 0.'
       end if

       allocate(tptr3d_2(pcols,pcnst,begchunk:endchunk))
       tptr3d_2 = 0_r8
       tptr3d_2(:,1,:) = tptr(:,:)

       call pbuf_set_field(pbuf2d, qpert_idx, tptr3d_2)
       deallocate(tptr3d_2)
    end if

    fieldname='CUSH'
    m = pbuf_get_index('cush', ierr)
    if (m > 0) then
       call infld(fieldname, fh_ini, dim1name, dim2name, 1, pcols, begchunk, endchunk, &
            tptr, found, gridname='physgrid')
       if(.not.found) then
          if(masterproc) write(iulog,*) trim(fieldname), ' initialized to 1000.'
          tptr=1000._r8
       end if
       do n=1,dyn_time_lvls
          call pbuf_set_field(pbuf2d, m, tptr, start=(/1,n/), kount=(/pcols,1/))
       end do
       deallocate(tptr)
    end if

    !
    ! 3-D fields
    !

    allocate(tptr3d(pcols,pver,begchunk:endchunk))

    fieldname='CLOUD'
    m = pbuf_get_index('CLD')
    call infld(fieldname, fh_ini, dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, &
         tptr3d, found, gridname='physgrid')
    if(found) then
       do n = 1, dyn_time_lvls
          call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
       end do
    else
       call pbuf_set_field(pbuf2d, m, 0._r8)
       if (masterproc) write(iulog,*) trim(fieldname), ' initialized to 0.'
    end if

    fieldname='QCWAT'
    m = pbuf_get_index(fieldname,ierr)
    if (m > 0) then
       call infld(fieldname, fh_ini, dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, &
            tptr3d, found, gridname='physgrid')
       if(.not. found) then
          call infld('Q',fh_ini,dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, &
               tptr3d, found, gridname='physgrid')
          if (found) then
             if (masterproc) write(iulog,*) trim(fieldname), ' initialized with Q'
             if(dycore_is('LR')) call polar_average(pver, tptr3d)
          else
             if (masterproc) write(iulog,*) trim(fieldname), ' initialized to huge()'
             tptr3d = huge(1.0_r8)
          end if
       end if
       do n = 1, dyn_time_lvls
          call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
       end do
    end if

    fieldname = 'ICCWAT'
    m = pbuf_get_index(fieldname, ierr)
    if (m > 0) then
       call infld(fieldname, fh_ini, dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, &
          tptr3d, found, gridname='physgrid')
       if(found) then
          do n = 1, dyn_time_lvls
             call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
          end do
       else
          call cnst_get_ind('CLDICE', ixcldice)
          call infld('CLDICE',fh_ini,dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, &
             tptr3d, found, gridname='physgrid')
          if(found) then
             do n = 1, dyn_time_lvls
                call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
             end do
          else
             call pbuf_set_field(pbuf2d, m, 0._r8)
          end if
          if (masterproc) then
             if (found) then
                write(iulog,*) trim(fieldname), ' initialized with CLDICE'
             else
                write(iulog,*) trim(fieldname), ' initialized to 0.0'
             end if
          end if
       end if
    end if

    fieldname = 'LCWAT'
    m = pbuf_get_index(fieldname,ierr)
    if (m > 0) then
       call infld(fieldname, fh_ini, dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, &
            tptr3d, found, gridname='physgrid')
       if(found) then
          do n = 1, dyn_time_lvls
             call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
          end do
       else
          allocate(tptr3d_2(pcols,pver,begchunk:endchunk))
          call cnst_get_ind('CLDICE', ixcldice)
          call cnst_get_ind('CLDLIQ', ixcldliq)
          call infld('CLDICE',fh_ini,dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, &
               tptr3d, found, gridname='physgrid')
          call infld('CLDLIQ',fh_ini,dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, &
               tptr3d_2, found2, gridname='physgrid')
          if(found .and. found2) then
             do lchnk = begchunk, endchunk
                ncol = get_ncols_p(lchnk)
                tptr3d(:ncol,:,lchnk)=tptr3d(:ncol,:,lchnk)+tptr3d_2(:ncol,:,lchnk)
             end do
             if (masterproc) write(iulog,*) trim(fieldname), ' initialized with CLDICE + CLDLIQ'
          else if (found) then ! Data already loaded in tptr3d
             if (masterproc) write(iulog,*) trim(fieldname), ' initialized with CLDICE only'
          else if (found2) then
             tptr3d(:,:,:)=tptr3d_2(:,:,:)
             if (masterproc) write(iulog,*) trim(fieldname), ' initialized with CLDLIQ only'
          end if

          if (found .or. found2) then
             do n = 1, dyn_time_lvls
                call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
             end do
             if(dycore_is('LR')) call polar_average(pver, tptr3d)
          else
             call pbuf_set_field(pbuf2d, m, 0._r8)
             if (masterproc)  write(iulog,*) trim(fieldname), ' initialized to 0.0'
          end if
          deallocate(tptr3d_2)
       end if
    end if

    deallocate(tptr3d)
    allocate(tptr3d(pcols,pver,begchunk:endchunk))

    fieldname = 'TCWAT'
    m = pbuf_get_index(fieldname,ierr)
    if (m > 0) then
       call infld(fieldname, fh_ini, dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, &
            tptr3d, found, gridname='physgrid')
       if(.not.found) then
          call infld('T', fh_ini, dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, &
               tptr3d, found, gridname='physgrid')
          if (found) then
             if(dycore_is('LR')) call polar_average(pver, tptr3d)
             if (masterproc) write(iulog,*) trim(fieldname), ' initialized with T'
          else
             if (masterproc) write(iulog,*) trim(fieldname), ' initialized to huge()'
             tptr3d = huge(1._r8)
          end if
       end if
       do n = 1, dyn_time_lvls
          call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
       end do
    end if

    deallocate(tptr3d)
    allocate(tptr3d(pcols,pverp,begchunk:endchunk))

    fieldname = 'TKE'
    m = pbuf_get_index( 'tke')
    call infld(fieldname, fh_ini, dim1name, 'ilev', dim2name, 1, pcols, 1, pverp, begchunk, endchunk, &
         tptr3d, found, gridname='physgrid')
    if (found) then
       call pbuf_set_field(pbuf2d, m, tptr3d)
    else
       call pbuf_set_field(pbuf2d, m, 0.01_r8)
       if (masterproc) write(iulog,*) trim(fieldname), ' initialized to 0.01'
    end if


    fieldname = 'KVM'
    m = pbuf_get_index('kvm')
    call infld(fieldname, fh_ini, dim1name, 'ilev', dim2name, 1, pcols, 1, pverp, begchunk, endchunk, &
         tptr3d, found, gridname='physgrid')
    if (found) then
       call pbuf_set_field(pbuf2d, m, tptr3d)
    else
       call pbuf_set_field(pbuf2d, m, 0._r8)
       if (masterproc) write(iulog,*) trim(fieldname), ' initialized to 0.'
    end if


    fieldname = 'KVH'
    m = pbuf_get_index('kvh')
    call infld(fieldname, fh_ini, dim1name, 'ilev', dim2name, 1, pcols, 1, pverp, begchunk, endchunk, &
         tptr3d, found, gridname='physgrid')
    if (found) then
       call pbuf_set_field(pbuf2d, m, tptr3d)
    else
       call pbuf_set_field(pbuf2d, m, 0._r8)
       if (masterproc) write(iulog,*) trim(fieldname), ' initialized to 0.'
    end if

    deallocate(tptr3d)
    allocate(tptr3d(pcols,pver,begchunk:endchunk))

    fieldname = 'CONCLD'
    m = pbuf_get_index('CONCLD',ierr)
    if (m > 0) then
       call infld(fieldname, fh_ini, dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, &
            tptr3d, found, gridname='physgrid')
       if(found) then
          do n = 1, dyn_time_lvls
             call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
          end do
       else
          call pbuf_set_field(pbuf2d, m, 0._r8)
          if (masterproc) write(iulog,*) trim(fieldname), ' initialized to 0.'
       end if

       deallocate (tptr3d)
    end if

    call initialize_short_lived_species(fh_ini, pbuf2d)

    !---------------------------------------------------------------------------------
    !  If needed, get ion and electron temperature fields from initial condition file
    !---------------------------------------------------------------------------------

    call waccmx_phys_ion_elec_temp_inidat(fh_ini,pbuf2d)

  end subroutine phys_inidat


  subroutine phys_init( phys_state, phys_tend, pbuf2d, cam_out )

#ifdef MMF_NN_EMULATOR
   use mmf_nn_emulator,    only: init_neural_net
#endif
    !-----------------------------------------------------------------------
    !
    ! Initialization of physics package.
    !
    !-----------------------------------------------------------------------

    use physics_buffer,     only: physics_buffer_desc, pbuf_initialize, pbuf_get_index
    use physconst,          only: rair, cpair, gravit, stebol, tmelt, &
                                  latvap, latice, rh2o, rhoh2o, pstd, zvir, &
                                  karman, rhodair, physconst_init
    use ref_pres,           only: pref_edge, pref_mid

    use carma_intr,         only: carma_init
    use cam_control_mod,    only: initial_run
    use check_energy,       only: check_energy_init
    use chemistry,          only: chem_init
    use prescribed_ozone,   only: prescribed_ozone_init
    use prescribed_ghg,     only: prescribed_ghg_init
    use prescribed_aero,    only: prescribed_aero_init
    use aerodep_flx,        only: aerodep_flx_init
    use aircraft_emit,      only: aircraft_emit_init
    use prescribed_volcaero,only: prescribed_volcaero_init
    use prescribed_strataero,only: prescribed_strataero_init
    use cloud_fraction,     only: cldfrc_init
    use cldfrc2m,           only: cldfrc2m_init
    use co2_cycle,          only: co2_init, co2_transport
    use convect_deep,       only: convect_deep_init
    use convect_shallow,    only: convect_shallow_init
    use cam_diagnostics,    only: diag_init
    use gw_drag,            only: gw_init
    use cam3_aero_data,     only: cam3_aero_data_on, cam3_aero_data_init
    use cam3_ozone_data,    only: cam3_ozone_data_on, cam3_ozone_data_init
    use radheat,            only: radheat_init
    use radiation,          only: radiation_init
    use cloud_diagnostics,  only: cloud_diagnostics_init
    use rk_stratiform,      only: rk_stratiform_init
    use wv_saturation,      only: wv_sat_init
    use microp_driver,      only: microp_driver_init
    use microp_aero,        only: microp_aero_init
    use macrop_driver,      only: macrop_driver_init
    use conv_water,         only: conv_water_init
    use spcam_drivers,      only: spcam_init
    use tracers,            only: tracers_init
    use aoa_tracers,        only: aoa_tracers_init
    use rayleigh_friction,  only: rayleigh_friction_init
    use pbl_utils,          only: pbl_utils_init
    use vertical_diffusion, only: vertical_diffusion_init
    use phys_debug_util,    only: phys_debug_init
    use phys_debug,         only: phys_debug_state_init
    use rad_constituents,   only: rad_cnst_init
    use aer_rad_props,      only: aer_rad_props_init
    use subcol,             only: subcol_init
    use qbo,                only: qbo_init
    use qneg_module,        only: qneg_init
    use iondrag,            only: iondrag_init, do_waccm_ions
#if ( defined OFFLINE_DYN )
    use metdata,            only: metdata_phys_init
#endif
    use epp_ionization,     only: epp_ionization_init, epp_ionization_active
    use waccmx_phys_intr,   only: waccmx_phys_ion_elec_temp_init  ! Initialization of ionosphere module (WACCM-X)
    use waccmx_phys_intr,   only: waccmx_phys_mspd_init   ! Initialization of major species diffusion module (WACCM-X)
    use clubb_intr,         only: clubb_ini_cam
    use sslt_rebin,         only: sslt_rebin_init
    use tropopause,         only: tropopause_init
    use solar_data,         only: solar_data_init
    use dadadj_cam,         only: dadadj_init
    use cam_abortutils,     only: endrun
    use nudging,            only: Nudge_Model, nudging_init
    use cb24mjocnn,         only: cb24mjocnn_init
    use cb24cnn,            only: cb24cnn_init

    ! Input/output arguments
    type(physics_state), pointer       :: phys_state(:)
    type(physics_tend ), pointer       :: phys_tend(:)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    type(cam_out_t),intent(inout)      :: cam_out(begchunk:endchunk)

    ! local variables
    integer :: lchnk
    integer :: ierr

    !-----------------------------------------------------------------------

    call physics_type_alloc(phys_state, phys_tend, begchunk, endchunk, pcols)

    do lchnk = begchunk, endchunk
       call physics_state_set_grid(lchnk, phys_state(lchnk))
    end do

    !-------------------------------------------------------------------------------------------
    ! Initialize any variables in physconst which are not temporally and/or spatially constant
    !-------------------------------------------------------------------------------------------
    call physconst_init()

    ! Initialize debugging a physics column
    call phys_debug_init()

    call pbuf_initialize(pbuf2d)

    ! Initialize subcol scheme
    call subcol_init(pbuf2d)

    ! diag_init makes addfld calls for dynamics fields that are output from
    ! the physics decomposition
    call diag_init(pbuf2d)

    call check_energy_init()

    call tracers_init()

    ! age of air tracers
    call aoa_tracers_init()

    teout_idx = pbuf_get_index( 'TEOUT')

    ! adiabatic or ideal physics should be only used if in simple_physics
    if (adiabatic .or. ideal_phys) then
      if (adiabatic) then
        call endrun('phys_init: adiabatic configuration error')
      else
        call endrun('phys_init: ideal_phys configuration error')
      end if
    end if

    if (initial_run) then
       call phys_inidat(cam_out, pbuf2d)
    end if

    ! wv_saturation is relatively independent of everything else and
    ! low level, so init it early. Must at least do this before radiation.
    call wv_sat_init

    ! CAM3 prescribed aerosols
    if (cam3_aero_data_on) call cam3_aero_data_init(phys_state)

    ! Initialize rad constituents and their properties
    call rad_cnst_init()
    call aer_rad_props_init()

    ! initialize carma
    call carma_init()

    ! solar irradiance data modules
    call solar_data_init()

    ! Prognostic chemistry.
    call chem_init(phys_state,pbuf2d)

    ! Prescribed tracers
    call prescribed_ozone_init()
    call prescribed_ghg_init()
    call prescribed_aero_init()
    call aerodep_flx_init()
    call aircraft_emit_init()
    call prescribed_volcaero_init()
    call prescribed_strataero_init()

    ! co2 cycle
    if (co2_transport()) then
       call co2_init()
    end if

    ! CAM3 prescribed ozone
    if (cam3_ozone_data_on) call cam3_ozone_data_init(phys_state)

    call gw_init()

    call rayleigh_friction_init()

    call pbl_utils_init(gravit, karman, cpair, rair, zvir)
    call vertical_diffusion_init(pbuf2d)

    if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then
       call waccmx_phys_mspd_init ()
       ! Initialization of ionosphere module if mode set to ionosphere
       if( waccmx_is('ionosphere') ) then
          call waccmx_phys_ion_elec_temp_init(pbuf2d)
       endif
    endif

    call radiation_init(pbuf2d)

    call cloud_diagnostics_init()

    call radheat_init(pref_mid)

    call convect_shallow_init(pref_edge, pbuf2d)

    call cldfrc_init()
    call cldfrc2m_init()

    call convect_deep_init(pref_edge)

    if( microp_scheme == 'RK' ) then
       call rk_stratiform_init()
    elseif( microp_scheme == 'MG' ) then
       if (.not. do_clubb_sgs) call macrop_driver_init(pbuf2d)
       call microp_aero_init()
       call microp_driver_init(pbuf2d)
       call conv_water_init
    elseif( microp_scheme == 'SPCAM_m2005') then
       call conv_water_init
    end if


    ! initiate CLUBB within CAM
    if (do_clubb_sgs) call clubb_ini_cam(pbuf2d)

    call spcam_init(pbuf2d)

    call qbo_init

    call iondrag_init(pref_mid)
    ! Geomagnetic module -- after iondrag_init
    if (epp_ionization_active) then
      call epp_ionization_init()
    endif

#if ( defined OFFLINE_DYN )
    call metdata_phys_init()
#endif
    call sslt_rebin_init()
    call tropopause_init()
    call dadadj_init()

    prec_dp_idx  = pbuf_get_index('PREC_DP')
    snow_dp_idx  = pbuf_get_index('SNOW_DP')
    prec_sh_idx  = pbuf_get_index('PREC_SH')
    snow_sh_idx  = pbuf_get_index('SNOW_SH')

    dlfzm_idx = pbuf_get_index('DLFZM', ierr)

    call phys_getopts(prog_modal_aero_out=prog_modal_aero)

    ! Initialize Nudging Parameters
    !--------------------------------
    if(Nudge_Model) call nudging_init
    

    if (clim_modal_aero) then

       ! If climate calculations are affected by prescribed modal aerosols, the
       ! the initialization routine for the dry mode radius calculation is called
       ! here.  For prognostic MAM the initialization is called from
       ! modal_aero_initialize
       if (.not. prog_modal_aero) then
          call modal_aero_calcsize_init(pbuf2d)
       endif

       call modal_aero_wateruptake_init(pbuf2d)

    end if

    ! Initialize qneg3 and qneg4
    call qneg_init()
    !call cb24mjocnn_init
    call cb24cnn_init

#ifdef MMF_NN_EMULATOR
  call init_neural_net()
#endif

  end subroutine phys_init

!===================================================================================================
!===================================================================================================
#ifdef MMF_NN_EMULATOR
  subroutine mmf_nn_emulator_driver(phys_state, phys_state_aphys1, phys_state_mmf, ztodt, phys_tend, pbuf2d,  cam_in, cam_out, cam_out_mmf)
    !-----------------------------------------------------------------------------
    ! Purpose: mmf_nn_emulator driver
    !-----------------------------------------------------------------------------
    use mmf_nn_emulator,          only: cb_partial_coupling, cb_partial_coupling_vars, cb_spinup_step, cb_do_ramp, cb_ramp_linear_steps, &
                                        cb_ramp_option, cb_ramp_factor, cb_ramp_step_0steps, cb_ramp_step_1steps
    use physics_buffer,   only: physics_buffer_desc, pbuf_get_chunk, &
                                pbuf_allocate, pbuf_get_index, pbuf_get_field
    use time_manager,     only: get_nstep, get_step_size, & 
                                is_first_step,  is_first_restart_step, &
                                get_curr_calday
    use cam_diagnostics,  only: diag_allocate, &
                                diag_mmf_nn_emulator_debug
    use radiation,        only: iradsw, use_rad_dt_cosz 
    use radconstants,     only: nswbands, get_ref_solar_band_irrad
    use rad_solar_var,    only: get_variability
    use orbit,            only: zenith
    use shr_orb_mod,      only: shr_orb_decl
    use phys_grid,        only: get_rlat_all_p, get_rlon_all_p
    use cam_control_mod,  only: lambm0, obliqr, eccen, mvelpp
    use tropopause,       only: tropopause_output
    use camsrfexch,       only: cam_export
    use cam_diagnostics,  only: diag_export
    use geopotential,        only: geopotential_t
    use cam_history_support, only: pflds
    use physconst,        only: cpair, zvir, rair, gravit
    use string_utils,    only: to_lower
    use cam_diagnostics,        only: diag_phys_writeout, diag_conv
    use cloud_diagnostics,      only: cloud_diagnostics_calc

    !-----------------------------------------------------------------------------
    ! Interface arguments
    !-----------------------------------------------------------------------------
    real(r8), intent(in) :: ztodt            ! physics time step unless nstep=0
    type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
    type(physics_state), intent(in),    dimension(begchunk:endchunk) :: phys_state_aphys1
    type(physics_state), intent(inout),    dimension(begchunk:endchunk) :: phys_state_mmf
    type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend
    type(physics_buffer_desc), pointer, dimension(:,:)               :: pbuf2d
    type(cam_in_t),                     dimension(begchunk:endchunk) :: cam_in
    type(cam_out_t),                    dimension(begchunk:endchunk) :: cam_out
    type(cam_out_t),      intent(out),      dimension(begchunk:endchunk) :: cam_out_mmf

    !-----------------------------------------------------------------------------
    ! Local Variables
    !-----------------------------------------------------------------------------
    integer :: lchnk                             ! chunk identifier
    integer :: ncol                              ! number of columns
    integer :: nstep                             ! current timestep number
    type(physics_buffer_desc), pointer :: phys_buffer_chunk(:)

    ! stuff for calling crm_physics_tend
    logical           :: use_ECPP
    type(physics_ptend), dimension(begchunk:endchunk) :: ptend ! indivdual parameterization tendencies

    integer  :: itim_old, cldo_idx, cld_idx   ! pbuf indices
    real(r8), pointer, dimension(:,:) :: cld  ! cloud fraction
    real(r8), pointer, dimension(:,:) :: cldo ! old cloud fraction
    !-----------------------------------------------------------------------------

    !-----------------------------------------------------------------------------
    ! Local Variables (for MMF_NN_EMULATOR)
    !-----------------------------------------------------------------------------
    integer :: c, i, j, k 
    integer, save :: nstep0
    integer       :: nstep_NN, dtime
    logical :: do_mmf_nn_emulator_inference = .false.

    ! - for partial coupling - !
    ! type(physics_state), dimension(begchunk:endchunk)  :: phys_state_nn
    ! type(physics_state), dimension(begchunk:endchunk)  :: phys_state_mmf_backup
    ! type(physics_tend ), dimension(begchunk:endchunk)  :: phys_tend_nn

    type(physics_state), allocatable, dimension(:)  :: phys_state_nn
    type(physics_state), allocatable, dimension(:)  :: phys_state_mmf_backup
    type(physics_tend ), allocatable, dimension(:)  :: phys_tend_nn
    type(cam_out_t),     dimension(begchunk:endchunk)  :: cam_out_nn
    integer :: ixcldice, ixcldliq
    integer :: prec_dp_idx, snow_dp_idx
    real(r8), dimension(:), pointer              :: prec_dp    , snow_dp
    real(r8), dimension(pcols,begchunk:endchunk) :: prec_dp_nn , snow_dp_nn, &
                                                  prec_dp_mmf, snow_dp_mmf
    logical  :: do_geopotential = .false.
    real(r8) :: zvirv_loc(pcols,pver), rairv_loc(pcols,pver)  
    real(r8) :: ramp_ratio 
    ! - !

    real(r8) :: calday       ! current calendar day
    real(r8) :: clat(pcols)  ! current latitudes(radians)
    real(r8) :: clon(pcols)  ! current longitudes(radians)
    real(r8), dimension(pcols,begchunk:endchunk) :: coszrs  ! Cosine solar zenith angle
    real(r8), dimension(pcols,begchunk:endchunk) :: solin   ! Insolation

    real(r8) :: sfac(1:nswbands)  ! time varying scaling factors due to Solar Spectral Irrad at 1 A.U. per band
    real(r8) :: solar_band_irrad(1:nswbands) ! rrtmg-assumed solar irradiance in each sw band
    real(r8) :: dt_avg = 0.0_r8   ! time step to use for the shr_orb_cosz calculation, if use_rad_dt_cosz set to true
    real(r8) :: delta    ! Solar declination angle  in radians
    real(r8) :: eccf     ! Earth orbit eccentricity factor
    integer :: ierr=0
    !-----------------------------------------------------------------------------
    ! phys_run1 opening
    ! - phys_timestep_init advances ghg gases,
    ! - need to advance solar insolation (for NN)
    !-----------------------------------------------------------------------------


    allocate(phys_state_nn(begchunk:endchunk), stat=ierr)
    if (ierr /= 0) then
       ! Handle allocation error
       write(iulog,*) 'Error allocating phys_state_nn error = ',ierr
    end if

    allocate(phys_state_mmf_backup(begchunk:endchunk), stat=ierr)
    if (ierr /= 0) then
       ! Handle allocation error
       write(iulog,*) 'Error allocating phys_state_mmf_backup error = ',ierr
    end if

    allocate(phys_tend_nn(begchunk:endchunk), stat=ierr)
    if (ierr /= 0) then
       ! Handle allocation error
       write(iulog,*) 'Error allocating phys_tend_nn error = ',ierr
    end if

    do lchnk=begchunk,endchunk
       call physics_state_alloc(phys_state_nn(lchnk),lchnk,pcols)
       call physics_state_alloc(phys_state_mmf_backup(lchnk),lchnk,pcols)
       ! call physics_tend_alloc(phys_tend_nn(lchnk),lchnk,pcols)
    end do

    do lchnk=begchunk,endchunk
       call physics_tend_alloc(phys_tend_nn(lchnk),phys_state_nn(lchnk)%psetcols)
    end do

    nstep = get_nstep()
    dtime = get_step_size()

    call pbuf_allocate(pbuf2d, 'physpkg')
    call diag_allocate()

    ! Advance time information
    call t_startf ('phys_timestep_init')
    call phys_timestep_init( phys_state, cam_out, pbuf2d)
    call t_stopf ('phys_timestep_init')

    ! Calculate  COSZRS and SOLIN
    call get_ref_solar_band_irrad( solar_band_irrad ) ! this can move to init subroutine
    call get_variability(sfac)                        ! "
    do lchnk=begchunk,endchunk
       ncol = phys_state(lchnk)%ncol
       calday = get_curr_calday(-dtime) ! get current calendar day with a negative offset to match the time in mli and CRM physics
       ! coszrs
       call get_rlat_all_p(lchnk, ncol, clat)
       call get_rlon_all_p(lchnk, ncol, clon)
       if (use_rad_dt_cosz)  then
         dtime  = get_step_size()
         dt_avg = iradsw*dtime
       end if
       call zenith(calday, clat, clon, coszrs(:,lchnk), ncol, dt_avg)
       ! solin
       call shr_orb_decl(calday  ,eccen     ,mvelpp  ,lambm0  ,obliqr  , &
                         delta   ,eccf      )
       solin(:,lchnk) = sum(sfac(:)*solar_band_irrad(:)) * eccf * coszrs(:,lchnk)
    end do
    ! [TO-DO] Check solin and coszrs from this calculation vs. pbuf_XXX

    prec_dp_idx = pbuf_get_index('PREC_DP', errcode=i) ! Query physics buffer index
    snow_dp_idx = pbuf_get_index('SNOW_DP', errcode=i)

    !-----------------------------------------------------------------------------
    ! phys_run1 main
    !-----------------------------------------------------------------------------

    ! Call init subroutine for neural networks
    ! (loading neural network weights and normalization factors)
    if (is_first_step() .or. is_first_restart_step()) then
       ! call init_neural_net()
       nstep0 = nstep
    end if

    ! Determine if MMF spin-up perioid is over
    ! (currently spin up time is set at 86400 sec ~ 1 day)
    ! [TO-DO] create a namelist variable for mmf spin-up time
    ! nstep_NN = 86400 / get_step_size()
    nstep_NN = cb_spinup_step

    if (nstep-nstep0 .eq. nstep_NN) then
       do_mmf_nn_emulator_inference = .true.
       if (masterproc) then
          write(iulog,*) '---------------------------------------'
          write(iulog,*) '[MMF_NN_EMULATOR] NN coupling starts'
          write(iulog,*) '---------------------------------------'
       end if
    end if

#ifdef MMF_NN_EMULATORDEBUG
  if (masterproc) then
     write (iulog,*) '[MMF_NN_EMULATORDEBUG] nstep - nstep0, nstep_NN, do_mmf_nn_emulator = ', nstep - nstep0, nstep_NN, do_mmf_nn_emulator_inference
  endif
#endif

  !Save original values of subroutine arguments
    if (do_mmf_nn_emulator_inference .and. cb_partial_coupling) then
       do lchnk = begchunk, endchunk
          ! since phys_state_mmf_backup etc is just allocated but have not been initialized (empty), doing this copy won't lead to memory leak
          phys_state_nn(lchnk) = phys_state(lchnk) 
          phys_state_mmf_backup(lchnk) = phys_state_mmf(lchnk)
          phys_tend_nn(lchnk)  = phys_tend(lchnk) 
          cam_out_nn(lchnk)    = cam_out(lchnk) 
       end do
    end if

    ! Run phys_run1 physics
    if (.not. do_mmf_nn_emulator_inference) then  ! MMFspin-up
       call phys_run1(phys_state, ztodt, phys_tend, pbuf2d,  cam_in, cam_out)
       do lchnk = begchunk, endchunk
          call physics_state_dealloc(phys_state_mmf(lchnk)) ! to prevent memory leak
          call physics_state_copy(phys_state(lchnk), phys_state_mmf(lchnk))
          cam_out_mmf(lchnk)    = cam_out(lchnk)
       end do

    else  ! NN inference
       if (cb_partial_coupling) then ! NN partial coupling

          call phys_run1   (phys_state,    ztodt, phys_tend,    pbuf2d, cam_in, cam_out)
          ! store mmf calculation of prec_dp and snow_dp
          do lchnk = begchunk, endchunk
             phys_buffer_chunk => pbuf_get_chunk(pbuf2d, lchnk)
             call pbuf_get_field(phys_buffer_chunk, prec_dp_idx, prec_dp)
             call pbuf_get_field(phys_buffer_chunk, snow_dp_idx, snow_dp)
             prec_dp_mmf(:,lchnk) = prec_dp(:) 
             snow_dp_mmf(:,lchnk) = snow_dp(:)
          end do

          do lchnk = begchunk, endchunk
             ! update phys_state_mmf but cannot overwrite its the mmf adv and phy history
             call physics_state_dealloc(phys_state_mmf(lchnk))
             call physics_state_copy(phys_state(lchnk), phys_state_mmf(lchnk))
             phys_state_mmf(lchnk)%t_adv(:,:,:) = phys_state_mmf_backup(lchnk)%t_adv(:,:,:)
             phys_state_mmf(lchnk)%u_adv(:,:,:) = phys_state_mmf_backup(lchnk)%u_adv(:,:,:)
             phys_state_mmf(lchnk)%t_phy(:,:,:) = phys_state_mmf_backup(lchnk)%t_phy(:,:,:)
             phys_state_mmf(lchnk)%u_phy(:,:,:) = phys_state_mmf_backup(lchnk)%u_phy(:,:,:)
             phys_state_mmf(lchnk)%q_adv(:,:,:,:) = phys_state_mmf_backup(lchnk)%q_adv(:,:,:,:)
             phys_state_mmf(lchnk)%q_phy(:,:,:,:) = phys_state_mmf_backup(lchnk)%q_phy(:,:,:,:)
             cam_out_mmf(lchnk) = cam_out(lchnk)
          end do

          call phys_run1_NN(phys_state_nn, phys_state_aphys1, ztodt, phys_tend_nn, pbuf2d, cam_in, cam_out_nn,&
                            solin, coszrs)
          ! store nn calculation of prec_dp and snow_dp
          do lchnk = begchunk, endchunk
             phys_buffer_chunk => pbuf_get_chunk(pbuf2d, lchnk)
             call pbuf_get_field(phys_buffer_chunk, prec_dp_idx, prec_dp)
             call pbuf_get_field(phys_buffer_chunk, snow_dp_idx, snow_dp)
             prec_dp_nn(:,lchnk) = prec_dp(:)
             snow_dp_nn(:,lchnk) = snow_dp(:)
             prec_dp(:) = prec_dp_mmf(:,lchnk) ! restored to mmf calculation
             snow_dp(:) = snow_dp_mmf(:,lchnk) ! (prep for cb_partial_coupling)
          end do

       else ! NN full coupling
          call phys_run1_NN(phys_state, phys_state_aphys1, ztodt, phys_tend, pbuf2d,  cam_in, cam_out,&
                            solin, coszrs)
          do lchnk = begchunk, endchunk
            ! in fully nn coupling case (no partial coupling), phys_state_mmf is just synced with phys_state but won't be used
            call physics_state_dealloc(phys_state_mmf(lchnk))
            call physics_state_copy(phys_state(lchnk), phys_state_mmf(lchnk))
            cam_out_mmf(lchnk)    = cam_out(lchnk)
          end do
       end if ! (cb_partial_coupling)
    end if ! (.not. do_mmf_nn_emulator_inference)

    ! Partial coupling
    ! NN calculations overide MMF calculations for any variables included in 'cb_partial_coupling_vars'
    ! e.g., [ 'ptend_t','ptend_q0001','ptend_q0002','ptend_q0003', 'ptend_u', 'ptend_v',
    !         'cam_out_NETSW', 'cam_out_FLWDS', 'cam_out_PRECSC', 'cam_out_PRECC',
    !         'cam_out_SOLS', 'cam_out_SOLL', 'cam_out_SOLSD', 'cam_out_SOLLD'           ]
    if (do_mmf_nn_emulator_inference .and. cb_partial_coupling) then
  
       if (cb_do_ramp) then
         write (iulog,*) 'MMF_NN_EMULATOR partial coupling: cb_ramp_option = ', trim(cb_ramp_option)
  
         select case (to_lower(trim(cb_ramp_option)))
           case('constant')
             ramp_ratio = cb_ramp_factor
           case('linear')
          
             if (nstep-nstep0-nstep_NN .le. cb_ramp_linear_steps) then
               ramp_ratio = cb_ramp_factor * (nstep-nstep0-nstep_NN)*1.0/(cb_ramp_linear_steps*1.0)
             else
               ramp_ratio = cb_ramp_factor
             end if
           case('step')
          
             if (mod(nstep-nstep0-nstep_NN, (cb_ramp_step_0steps + cb_ramp_step_1steps)) .le. cb_ramp_step_1steps) then
               ramp_ratio = cb_ramp_factor
             else
               ramp_ratio = 0.0
             end if
         end select
       else
         ramp_ratio = 1.0
       end if
    
      if (cb_do_ramp) then
        write (iulog,*) 'MMF_NN_EMULATOR partial coupling: cb_do_ramp is on'
        write (iulog,*) 'MMF_NN_EMULATOR partial coupling: 1 is fully NN, ramp_ratio = ', ramp_ratio
        write (iulog,*) 'MMF_NN_EMULATOR partial coupling: cb_ramp_option = ', trim(cb_ramp_option)
      else
        write (iulog,*) 'MMF_NN_EMULATOR partial coupling: cb_do_ramp is off'
        write (iulog,*) 'MMF_NN_EMULATOR partial coupling: 1 is fully NN, ramp_ratio = ', ramp_ratio
      end if

      call cnst_get_ind('CLDICE', ixcldice)
      call cnst_get_ind('CLDLIQ', ixcldliq)
      do c = begchunk, endchunk
         k = 1
         do while (k < pflds  .and. cb_partial_coupling_vars(k) /= ' ')
            if (trim(cb_partial_coupling_vars(k)) == 'ptend_t') then
               phys_state(c)%t(:,:)   = phys_state_nn(c)%t(:,:)*ramp_ratio + phys_state(c)%t(:,:)*(1.0_r8-ramp_ratio)
               phys_tend(c)%dtdt(:,:) = phys_tend_nn(c)%dtdt(:,:)*ramp_ratio + phys_tend(c)%dtdt(:,:)*(1.0_r8-ramp_ratio)
               do_geopotential = .true.
               if (nstep-nstep0 .eq. nstep_NN .and. masterproc) then
                  write (iulog,*) 'MMF_NN_EMULATOR partial coupling: ', trim(cb_partial_coupling_vars(k))
               endif
            else if (trim(cb_partial_coupling_vars(k)) == 'ptend_q0001') then
               phys_state(c)%q(:,:,1) = phys_state_nn(c)%q(:,:,1)*ramp_ratio + phys_state(c)%q(:,:,1)*(1.0_r8-ramp_ratio) 
                do_geopotential = .true.
               if (nstep-nstep0 .eq. nstep_NN .and. masterproc) then
                  write (iulog,*) 'MMF_NN_EMULATOR partial coupling: ', trim(cb_partial_coupling_vars(k))
              endif
           else if (trim(cb_partial_coupling_vars(k)) == 'ptend_q0002') then
                phys_state(c)%q(:,:,ixcldliq) = phys_state_nn(c)%q(:,:,ixcldliq)*ramp_ratio + phys_state(c)%q(:,:,ixcldliq)*(1.0_r8-ramp_ratio)
               if (nstep-nstep0 .eq. nstep_NN .and. masterproc) then
                  write (iulog,*) 'MMF_NN_EMULATOR partial coupling: ', trim(cb_partial_coupling_vars(k))
               endif
            else if (trim(cb_partial_coupling_vars(k)) == 'ptend_q0003') then
               phys_state(c)%q(:,:,ixcldice) = phys_state_nn(c)%q(:,:,ixcldice)*ramp_ratio + phys_state(c)%q(:,:,ixcldice)*(1.0_r8-ramp_ratio)
               if (nstep-nstep0 .eq. nstep_NN .and. masterproc) then
                  write (iulog,*) 'MMF_NN_EMULATOR partial coupling: ', trim(cb_partial_coupling_vars(k))
                endif
            else if (trim(cb_partial_coupling_vars(k)) == 'ptend_u') then
               phys_state(c)%u(:,:)   = phys_state_nn(c)%u(:,:)*ramp_ratio + phys_state(c)%u(:,:)*(1.0_r8-ramp_ratio)
               phys_tend(c)%dudt(:,:) = phys_tend_nn(c)%dudt(:,:)*ramp_ratio + phys_tend(c)%dudt(:,:)*(1.0_r8-ramp_ratio)
               if (nstep-nstep0 .eq. nstep_NN .and. masterproc) then
                  write (iulog,*) 'MMF_NN_EMULATOR partial coupling: ', trim(cb_partial_coupling_vars(k))
               endif
            else if (trim(cb_partial_coupling_vars(k)) == 'ptend_v') then
               phys_state(c)%v(:,:)   = phys_state_nn(c)%v(:,:)*ramp_ratio + phys_state(c)%v(:,:)*(1.0_r8-ramp_ratio)
              phys_tend(c)%dvdt(:,:) = phys_tend_nn(c)%dvdt(:,:)*ramp_ratio + phys_tend(c)%dvdt(:,:)*(1.0_r8-ramp_ratio)
               if (nstep-nstep0 .eq. nstep_NN .and. masterproc) then
                  write (iulog,*) 'MMF_NN_EMULATOR partial coupling: ', trim(cb_partial_coupling_vars(k))
              endif
           else if (trim(cb_partial_coupling_vars(k)) == 'cam_out_NETSW') then
               cam_out(c)%netsw(:) = cam_out_nn(c)%netsw(:)*ramp_ratio + cam_out(c)%netsw(:)*(1.0_r8-ramp_ratio)
               if (nstep-nstep0 .eq. nstep_NN .and. masterproc) then
                  write (iulog,*) 'MMF_NN_EMULATOR partial coupling: ', trim(cb_partial_coupling_vars(k))
               endif
           else if (trim(cb_partial_coupling_vars(k)) == 'cam_out_FLWDS') then
              cam_out(c)%flwds(:) = cam_out_nn(c)%flwds(:)*ramp_ratio + cam_out(c)%flwds(:)*(1.0_r8-ramp_ratio)
              if (nstep-nstep0 .eq. nstep_NN .and. masterproc) then
                 write (iulog,*) 'MMF_NN_EMULATOR partial coupling: ', trim(cb_partial_coupling_vars(k))
              endif
           else if (trim(cb_partial_coupling_vars(k)) == 'cam_out_SOLS') then
              cam_out(c)%sols(:) = cam_out_nn(c)%sols(:)*ramp_ratio + cam_out(c)%sols(:)*(1.0_r8-ramp_ratio)
              if (nstep-nstep0 .eq. nstep_NN .and. masterproc) then
                 write (iulog,*) 'MMF_NN_EMULATOR partial coupling: ', trim(cb_partial_coupling_vars(k))
              endif
           else if (trim(cb_partial_coupling_vars(k)) == 'cam_out_SOLL') then
              cam_out(c)%soll(:) = cam_out_nn(c)%soll(:)*ramp_ratio + cam_out(c)%soll(:)*(1.0_r8-ramp_ratio)
              if (nstep-nstep0 .eq. nstep_NN .and. masterproc) then
                 write (iulog,*) 'MMF_NN_EMULATOR partial coupling: ', trim(cb_partial_coupling_vars(k))
              endif
           else if (trim(cb_partial_coupling_vars(k)) == 'cam_out_SOLSD') then
                cam_out(c)%solsd(:) = cam_out_nn(c)%solsd(:)*ramp_ratio + cam_out(c)%solsd(:)*(1.0_r8-ramp_ratio)
                if (nstep-nstep0 .eq. nstep_NN .and. masterproc) then
                   write (iulog,*) 'MMF_NN_EMULATOR partial coupling: ', trim(cb_partial_coupling_vars(k))
                endif
           else if (trim(cb_partial_coupling_vars(k)) == 'cam_out_SOLLD') then
              cam_out(c)%solld(:) = cam_out_nn(c)%solld(:)*ramp_ratio + cam_out(c)%solld(:)*(1.0_r8-ramp_ratio)
                if (nstep-nstep0 .eq. nstep_NN .and. masterproc) then
                 write (iulog,*) 'MMF_NN_EMULATOR partial coupling: ', trim(cb_partial_coupling_vars(k))
                endif
             else if (trim(cb_partial_coupling_vars(k)) == 'cam_out_PRECSC') then
              phys_buffer_chunk => pbuf_get_chunk(pbuf2d, c)
              call pbuf_get_field(phys_buffer_chunk, snow_dp_idx, snow_dp)
                snow_dp(:) = snow_dp_nn(:,c)*ramp_ratio + snow_dp(:)*(1.0_r8-ramp_ratio) 
               if (nstep-nstep0 .eq. nstep_NN .and. masterproc) then
                 write (iulog,*) 'MMF_NN_EMULATOR partial coupling: ', trim(cb_partial_coupling_vars(k))
              endif
           else if (trim(cb_partial_coupling_vars(k)) == 'cam_out_PRECC') then
              phys_buffer_chunk => pbuf_get_chunk(pbuf2d, c)
              call pbuf_get_field(phys_buffer_chunk, prec_dp_idx, prec_dp)
              prec_dp(:) = prec_dp_nn(:,c)*ramp_ratio + prec_dp(:)*(1.0_r8-ramp_ratio) 
              if (nstep-nstep0 .eq. nstep_NN .and. masterproc) then
                 write (iulog,*) 'MMF_NN_EMULATOR partial coupling: ', trim(cb_partial_coupling_vars(k))
              endif
           else
                call endrun('[MMF_NN_EMULATOR: cb_partial_coupling] Wrong variables are included in cb_partial_coupling_vars: ' // trim(cb_partial_coupling_vars(k)))
           end if
           k = k+1
         end do ! k

         if (do_geopotential) then
            ncol = phys_state(c)%ncol
            zvirv_loc(:,:) = zvir
            rairv_loc(:,:) = rair
            call geopotential_t  ( &
                 phys_state(c)%lnpint, phys_state(c)%lnpmid,   phys_state(c)%pint, phys_state(c)%pmid, phys_state(c)%pdel, phys_state(c)%rpdel, &
                 phys_state(c)%t     , phys_state(c)%q(:,:,1), rairv_loc(:,:),  gravit,     zvirv_loc(:,:), &
                 phys_state(c)%zi    , phys_state(c)%zm      , ncol)
            ! update dry static energy for use in next process
            do j = 1, pver
               phys_state(c)%s(:ncol,j) = phys_state(c)%t(:ncol,j)*cpair &
                                          + gravit*phys_state(c)%zm(:ncol,j) + phys_state(c)%phis(:ncol)
            end do ! j
         end if ! (do_geopotential)

      end do ! c
    end if ! (cb_partial coupling)

    ! copy from the tphysbc2 to here. make sure the outputted history file is consistent with the partial coupling
    do lchnk=begchunk, endchunk
      phys_buffer_chunk => pbuf_get_chunk(pbuf2d, lchnk)
      call t_startf('bc_history_write')
      call diag_phys_writeout(phys_state(lchnk), cam_out(lchnk)%psl)
      call diag_conv(phys_state(lchnk), ztodt, phys_buffer_chunk)
      call t_stopf('bc_history_write')

      !-----------------------------------------------------------------------------
      ! Write cloud diagnostics on history file
      !-----------------------------------------------------------------------------
      call t_startf('bc_cld_diag_history_write')
      call cloud_diagnostics_calc(phys_state(lchnk), phys_buffer_chunk)
      call t_stopf('bc_cld_diag_history_write')
    end do 
    !-----------------------------------------------------------------------------
    ! phys_run1 closing
    ! - tphysbc2 diagnostic (including cam_export)
    !-----------------------------------------------------------------------------
    do lchnk=begchunk, endchunk
       ! Diagnose the location of the tropopause
       call tropopause_output(phys_state(lchnk))

       ! Save atmospheric fields to force surface models
       phys_buffer_chunk => pbuf_get_chunk(pbuf2d, lchnk)
       call cam_export(phys_state(lchnk), cam_out(lchnk), phys_buffer_chunk)

       ! Write export state to history file
       call diag_export(cam_out(lchnk))
    end do

    do lchnk=begchunk,endchunk
      call physics_state_dealloc(phys_state_nn(lchnk))
      call physics_state_dealloc(phys_state_mmf_backup(lchnk))
      call physics_tend_dealloc(phys_tend_nn(lchnk))
    end do
    deallocate(phys_state_nn)
    deallocate(phys_state_mmf_backup)
    deallocate(phys_tend_nn)

  end subroutine mmf_nn_emulator_driver

  subroutine phys_run1_NN(phys_state, phys_state_aphys1, ztodt, phys_tend, pbuf2d,  cam_in, cam_out, &
                          solin, coszrs)
    !-----------------------------------------------------------------------------
    ! Purpose: First part of atmos physics before updating of surface components
    !-----------------------------------------------------------------------------
    use mmf_nn_emulator,         only: neural_net, &
         cb_partial_coupling, cb_partial_coupling_vars
    use physics_buffer,  only: physics_buffer_desc, pbuf_get_chunk, pbuf_get_field
    use time_manager,    only: get_nstep
    use check_energy,    only: check_energy_gmean
    use flux_avg,        only: flux_avg_init
    !-----------------------------------------------------------------------------
    ! Interface arguments
    !-----------------------------------------------------------------------------
    real(r8), intent(in) :: ztodt            ! physics time step unless nstep=0
    type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
    type(physics_state), intent(in),    dimension(begchunk:endchunk) :: phys_state_aphys1
    type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend

    type(physics_buffer_desc), pointer, dimension(:,:) :: pbuf2d
    type(cam_in_t),                     dimension(begchunk:endchunk) :: cam_in
    type(cam_out_t),                    dimension(begchunk:endchunk) :: cam_out

    real(r8), intent(in), dimension(pcols,begchunk:endchunk) :: coszrs  ! Cosine solar zenith angle
    real(r8), intent(in), dimension(pcols,begchunk:endchunk) :: solin   ! Insolation

    !-----------------------------------------------------------------------------
    ! Local Variables
    !-----------------------------------------------------------------------------
    type(physics_state)                              :: state
    type(physics_buffer_desc),pointer, dimension(:)  :: phys_buffer_chunk
    type(physics_ptend)                              :: ptend 
    integer :: lchnk                             ! chunk identifier
    integer :: nstep                             ! current timestep number
    integer :: ncol
    integer :: ixcldice, ixcldliq            ! constituent indices for cloud liquid and ice water.
    real(r8), pointer, dimension(:,:) :: tini
    real(r8), pointer, dimension(:,:) :: qini
    real(r8), pointer, dimension(:,:) :: cldliqini
    real(r8), pointer, dimension(:,:) :: cldiceini

    nullify(phys_buffer_chunk)
    nullify(tini)
    nullify(qini)
    nullify(cldliqini)
    nullify(cldiceini)

    !-----------------------------------------------------------------------------
    ! phys_run1 opening
    !-----------------------------------------------------------------------------
    nstep = get_nstep()

    ! Set physics tendencies to 0
    do lchnk=begchunk, endchunk
      ncol = phys_state(lchnk)%ncol
      phys_tend(lchnk)%dtdt(:ncol,:pver)  = 0._r8
      phys_tend(lchnk)%dudt(:ncol,:pver)  = 0._r8
      phys_tend(lchnk)%dvdt(:ncol,:pver)  = 0._r8
    end do

    ! The following initialization depends on the import state (cam_in)
    ! being initialized.  This isn't true when cam_init is called, so need
    ! to postpone this initialization to here.
    if (nstep == 0 .and. phys_do_flux_avg()) call flux_avg_init(cam_in,  pbuf2d)

    ! Compute total energy of input state and previous output state
    call t_startf ('chk_en_gmean')
    call check_energy_gmean(phys_state, pbuf2d, ztodt, nstep)
    call t_stopf ('chk_en_gmean')

#ifdef TRACER_CHECK
  call gmean_mass ('before tphysbc DRY', phys_state)
#endif

   ! these initial states will be used in tphysac diagnostics
    do lchnk=begchunk, endchunk
      state = phys_state(lchnk)
      phys_buffer_chunk => pbuf_get_chunk(pbuf2d, lchnk)
      call pbuf_get_field(phys_buffer_chunk, tini_idx, tini)
      call pbuf_get_field(phys_buffer_chunk, qini_idx, qini)
      call pbuf_get_field(phys_buffer_chunk, cldliqini_idx, cldliqini)
      call pbuf_get_field(phys_buffer_chunk, cldiceini_idx, cldiceini)

      ncol = phys_state(lchnk)%ncol
      tini(:ncol,:pver) = state%t(:ncol,:pver)
      call cnst_get_ind('CLDLIQ', ixcldliq)
      call cnst_get_ind('CLDICE', ixcldice)
      qini     (:ncol,:pver) = state%q(:ncol,:pver,       1)
      cldliqini(:ncol,:pver) = state%q(:ncol,:pver,ixcldliq)
      cldiceini(:ncol,:pver) = state%q(:ncol,:pver,ixcldice)
    end do

    !-----------------------------------------------------------------------------
    ! Neural network
    !-----------------------------------------------------------------------------
    do lchnk=begchunk, endchunk
      phys_buffer_chunk => pbuf_get_chunk(pbuf2d, lchnk)
      call neural_net (ptend, phys_state(lchnk), phys_state_aphys1(lchnk), &
                      phys_buffer_chunk, &
                      cam_in(lchnk), cam_out(lchnk), &
                      coszrs(:,lchnk), solin(:,lchnk), &
                      ztodt, lchnk)
      call physics_update_ml (phys_state(lchnk), ptend, ztodt, phys_tend(lchnk))
    end do

  !-----------------------------------------------------------------------------
  ! phys_run1 closing
  !-----------------------------------------------------------------------------
#ifdef TRACER_CHECK
  call gmean_mass ('between DRY', phys_state)
#endif

  end subroutine phys_run1_NN
#endif /* MMF_NN_EMULATOR */

  !
  !-----------------------------------------------------------------------
  !

  subroutine phys_run1(phys_state, ztodt, phys_tend, pbuf2d,  cam_in, cam_out)
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! First part of atmospheric physics package before updating of surface models
    !
    !-----------------------------------------------------------------------
    use time_manager,   only: get_nstep
    use cam_diagnostics,only: diag_allocate, diag_physvar_ic
    use check_energy,   only: check_energy_gmean
    use phys_control,   only: phys_getopts
    use spcam_drivers,  only: tphysbc_spcam
    use spmd_utils,     only: mpicom
    use physics_buffer, only: physics_buffer_desc, pbuf_get_chunk, pbuf_allocate
    use check_energy,   only: check_energy_gmean, check_energy_chng
    use cb24mjocnn,     only: cb24mjocnn_timestep_tender
    use cb24cnn,        only: cb24cnn_timestep_tend,cb24cnn_set_tend
    
#if (defined BFB_CAM_SCAM_IOP )
    use cam_history,    only: outfld
#endif
    use cam_abortutils, only: endrun
#if ( defined OFFLINE_DYN )
     use metdata,       only: get_met_srf1
#endif
    !
    ! Input arguments
    !
    real(r8), intent(in) :: ztodt            ! physics time step unless nstep=0
    !
    ! Input/Output arguments
    !
    type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
    type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend

    type(physics_buffer_desc), pointer, dimension(:,:) :: pbuf2d
    type(cam_in_t),                     dimension(begchunk:endchunk) :: cam_in
    type(cam_out_t),                    dimension(begchunk:endchunk) :: cam_out
    !-----------------------------------------------------------------------
    !
    !---------------------------Local workspace-----------------------------
    !
    integer :: c                                 ! indices
    integer :: ncol                              ! number of columns
    integer :: nstep                             ! current timestep number
    logical :: use_spcam
    type(physics_ptend)     :: ptend               ! WEC
    real(r8) :: zero(pcols)                        ! array of zeros
    type(physics_buffer_desc), pointer :: phys_buffer_chunk(:)

    call t_startf ('physpkg_st1')
    nstep = get_nstep()

#if ( defined OFFLINE_DYN )
    !
    ! if offline mode set SNOWH and TS for micro-phys
    !
    call get_met_srf1( cam_in )
#endif

    ! The following initialization depends on the import state (cam_in)
    ! being initialized.  This isn't true when cam_init is called, so need
    ! to postpone this initialization to here.
    if (nstep == 0 .and. phys_do_flux_avg()) call flux_avg_init(cam_in,  pbuf2d)

    ! Compute total energy of input state and previous output state
    call t_startf ('chk_en_gmean')
    call check_energy_gmean(phys_state, pbuf2d, ztodt, nstep)
    call t_stopf ('chk_en_gmean')
#ifndef MMF_NN_EMULATOR
    call t_stopf ('physpkg_st1')

    call t_startf ('physpkg_st1')

    call pbuf_allocate(pbuf2d, 'physpkg')
    call diag_allocate()

    !-----------------------------------------------------------------------
    ! Advance time information
    !-----------------------------------------------------------------------

    call phys_timestep_init(phys_state, cam_in, cam_out, pbuf2d)

    call t_stopf ('physpkg_st1')
#endif /* MMF_NN_EMULATOR */
#ifdef TRACER_CHECK
    call gmean_mass ('before tphysbc DRY', phys_state)
#endif

#ifdef MMF_NN_EMULATOR
  !-----------------------------------------------------------------------------
  ! Physics tendency before coupler - Phase 1
  !-----------------------------------------------------------------------------

    call t_barrierf('sync_bc_physics', mpicom)
    call t_startf ('bc_physics')
    call t_startf ('bc_physics1')

!$OMP PARALLEL DO PRIVATE (C, beg_count, phys_buffer_chunk, end_count, chunk_cost)
    do c=begchunk, endchunk

      beg_count = shr_sys_irtc(irtc_rate)

      phys_buffer_chunk => pbuf_get_chunk(pbuf2d, c)

      ! Output physics terms to IC file
      call t_startf ('diag_physvar_ic')
      call diag_physvar_ic ( c,  phys_buffer_chunk, cam_out(c), cam_in(c) )
      call t_stopf ('diag_physvar_ic')

      call tphysbc1(ztodt, fsns(1,c), fsnt(1,c), flns(1,c), flnt(1,c), &
                    phys_state(c), phys_tend(c), phys_buffer_chunk, &
                    fsds(1,c), sgh(1,c), sgh30(1,c), &
                    cam_out(c), cam_in(c) )

      end_count = shr_sys_irtc(irtc_rate)
      chunk_cost = real( (end_count-beg_count), r8)/real(irtc_rate, r8)
      call update_cost_p(c, chunk_cost)

    end do

    call t_stopf ('bc_physics1')

  !-----------------------------------------------------------------------------
  ! CRM physics
  !-----------------------------------------------------------------------------  

    if (use_ECPP) then
      call crm_ecpp_output_initialize(crm_ecpp_output_chunk,pcols,pver)
      call crm_ecpp_output_initialize(crm_ecpp_output,      ncrms,pver)
    end if

    call t_startf('crm_physics_tend')
    call crm_physics_tend(ztodt, phys_state, phys_tend, ptend, pbuf2d, cam_in, cam_out, &
                          species_class, crm_ecpp_output, &
                          mmf_qchk_prec_dp, mmf_qchk_snow_dp, mmf_rad_flux )
    call t_stopf('crm_physics_tend')

    do c=begchunk, endchunk
      call physics_update(phys_state(c), ptend(c), ztodt, phys_tend(c))
      call check_energy_chng(phys_state(c), phys_tend(c), "crm_tend", nstep, ztodt, zero, &
                             mmf_qchk_prec_dp(c,:), mmf_qchk_snow_dp(c,:), mmf_rad_flux(c,:))
    end do

  !-----------------------------------------------------------------------------
  ! Physics tendency before coupler - Phase 2
  !-----------------------------------------------------------------------------

    call t_barrierf('sync_bc_physics', mpicom)
    call t_startf ('bc_physics2')

    ! Determine column start and end indices for crm_ecpp_output
    ncol_sum = 0
    do c=begchunk, endchunk
      icol_beg(c) = ncol_sum + 1
      icol_end(c) = ncol_sum + phys_state(c)%ncol
      ncol_sum = ncol_sum + phys_state(c)%ncol
    end do

!$OMP PARALLEL DO PRIVATE (C, beg_count, phys_buffer_chunk, end_count, chunk_cost)
    do c=begchunk, endchunk

      beg_count = shr_sys_irtc(irtc_rate)

      ! Output physics terms to IC file
      phys_buffer_chunk => pbuf_get_chunk(pbuf2d, c)

      if (use_ECPP) call crm_ecpp_output_copy( crm_ecpp_output, crm_ecpp_output_chunk, icol_beg(c), icol_end(c))

      call tphysbc2(ztodt, fsns(1,c), fsnt(1,c), flns(1,c), flnt(1,c), &
                    phys_state(c), phys_tend(c), phys_buffer_chunk, &
                    fsds(1,c), sgh(1,c), sgh30(1,c), &
                    cam_out(c), cam_in(c), crm_ecpp_output_chunk )

      end_count = shr_sys_irtc(irtc_rate)
      chunk_cost = real( (end_count-beg_count), r8)/real(irtc_rate, r8)
      call update_cost_p(c, chunk_cost)

    end do

    call t_stopf ('bc_physics2')
    call t_stopf ('bc_physics')

    if (use_ECPP) then
      call crm_ecpp_output_finalize(crm_ecpp_output_chunk)
      call crm_ecpp_output_finalize(crm_ecpp_output)
    end if

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------

#else
    !-----------------------------------------------------------------------
    ! Tendency physics before flux coupler invocation
    !-----------------------------------------------------------------------
    !

#if (defined BFB_CAM_SCAM_IOP )
    do c=begchunk, endchunk
      call outfld('Tg',cam_in(c)%ts,pcols   ,c     )
    end do
#endif /* defined BFB_CAM_SCAM_IOP */

    call t_barrierf('sync_bc_physics', mpicom)
    call t_startf ('bc_physics')
    call t_adj_detailf(+1)

    call phys_getopts( use_spcam_out = use_spcam)

!$OMP PARALLEL DO PRIVATE (C, phys_buffer_chunk)
    do c=begchunk, endchunk
      !
      ! Output physics terms to IC file
      !
      phys_buffer_chunk => pbuf_get_chunk(pbuf2d, c)

      call t_startf ('diag_physvar_ic')
      call diag_physvar_ic ( c,  phys_buffer_chunk, cam_out(c), cam_in(c) )
      call t_stopf ('diag_physvar_ic')

      if (use_spcam) then
        call tphysbc_spcam (ztodt, phys_state(c),     &
             phys_tend(c), phys_buffer_chunk, &
             cam_out(c), cam_in(c) )
      else
        call tphysbc (ztodt, phys_state(c),           &
             phys_tend(c), phys_buffer_chunk, &
             cam_out(c), cam_in(c) )
      end if

    end do

    !++WEC Call here so python is only ever called once per timestep 
    if (masterproc) write(iulog,*)'nintendo past cb24mjocnn_timestep_init '
    !call cb24mjocnn_timestep_tender() 
    !++WEC Call here so python is only ever called once per timestep 
    if (masterproc) write(iulog,*)'nintendo past cb24cnn_timestep_init '
    call cb24cnn_timestep_tend() 
    do c=begchunk, endchunk
        call cb24cnn_set_tend(phys_state(c),ptend)
        call physics_update(phys_state(c),ptend,ztodt,phys_tend(c)) 
        call check_energy_chng(phys_state(c),phys_tend(c), "cb24cnn", nstep, ztodt, zero, zero, zero, zero)
    end do
    !--WEC

    call t_adj_detailf(-1)
    call t_stopf ('bc_physics')

    ! Don't call the rest in CRM mode
    if(single_column.and.scm_crm_mode) return
#endif /* ifndef MMF_NN_EMULATOR */
#ifdef TRACER_CHECK
    call gmean_mass ('between DRY', phys_state)
#endif

  end subroutine phys_run1

  !
  !-----------------------------------------------------------------------
  !

  subroutine phys_run2(phys_state, ztodt, phys_tend, pbuf2d,  cam_out, &
       cam_in )
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! Second part of atmospheric physics package after updating of surface models
    !
    !-----------------------------------------------------------------------
    use physics_buffer,  only: physics_buffer_desc, pbuf_get_chunk, pbuf_deallocate, pbuf_update_tim_idx
    use mo_lightning,    only: lightning_no_prod
    use cam_diagnostics, only: diag_deallocate, diag_surf
    use physconst,       only: stebol, latvap
    use carma_intr,      only: carma_accumulate_stats
    use spmd_utils,      only: mpicom
    use iop_forcing,     only: scam_use_iop_srf
#if ( defined OFFLINE_DYN )
    use metdata,         only: get_met_srf2
#endif
    !
    ! Input arguments
    !
    real(r8), intent(in) :: ztodt                       ! physics time step unless nstep=0
    !
    ! Input/Output arguments
    !
    type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
    type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend
    type(physics_buffer_desc),pointer,  dimension(:,:)               :: pbuf2d

    type(cam_out_t),     intent(inout), dimension(begchunk:endchunk) :: cam_out
    type(cam_in_t),      intent(inout), dimension(begchunk:endchunk) :: cam_in
    !
    !-----------------------------------------------------------------------
    !---------------------------Local workspace-----------------------------
    !
    integer :: c                                 ! chunk index
    integer :: ncol                              ! number of columns
    type(physics_buffer_desc),pointer, dimension(:)     :: phys_buffer_chunk
    !
    ! If exit condition just return
    !

    if(single_column.and.scm_crm_mode) then
       call diag_deallocate()
       return
    end if
    !-----------------------------------------------------------------------
    ! if using IOP values for surface fluxes overwrite here after surface components run
    !-----------------------------------------------------------------------
    if (single_column) call scam_use_iop_srf(cam_in)

    !-----------------------------------------------------------------------
    ! Tendency physics after coupler
    ! Not necessary at terminal timestep.
    !-----------------------------------------------------------------------
    !
#if ( defined OFFLINE_DYN )
    !
    ! if offline mode set SHFLX QFLX TAUX TAUY for vert diffusion
    !
    call get_met_srf2( cam_in )
#endif
    ! Set lightning production of NO
    call t_startf ('lightning_no_prod')
    call lightning_no_prod( phys_state, pbuf2d,  cam_in )
    call t_stopf ('lightning_no_prod')

    call t_barrierf('sync_ac_physics', mpicom)
    call t_startf ('ac_physics')
    call t_adj_detailf(+1)

!$OMP PARALLEL DO PRIVATE (C, NCOL, phys_buffer_chunk)

    do c=begchunk,endchunk
       ncol = get_ncols_p(c)
       phys_buffer_chunk => pbuf_get_chunk(pbuf2d, c)
       !
       ! surface diagnostics for history files
       !
       call t_startf('diag_surf')
       call diag_surf(cam_in(c), cam_out(c), phys_state(c), phys_buffer_chunk)
       call t_stopf('diag_surf')

       call tphysac(ztodt, cam_in(c),  &
            cam_out(c),                              &
            phys_state(c), phys_tend(c), phys_buffer_chunk)
    end do                    ! Chunk loop

    call t_adj_detailf(-1)
    call t_stopf('ac_physics')

#ifdef TRACER_CHECK
    call gmean_mass ('after tphysac FV:WET)', phys_state)
#endif

    call t_startf ('carma_accumulate_stats')
    call carma_accumulate_stats()
    call t_stopf ('carma_accumulate_stats')

    call t_startf ('physpkg_st2')
    call pbuf_deallocate(pbuf2d, 'physpkg')

    call pbuf_update_tim_idx()
    call diag_deallocate()
    call t_stopf ('physpkg_st2')

  end subroutine phys_run2

  !
  !-----------------------------------------------------------------------
  !

  subroutine phys_final( phys_state, phys_tend, pbuf2d )
    use physics_buffer, only : physics_buffer_desc, pbuf_deallocate
    use chemistry, only : chem_final
    use carma_intr, only : carma_final
    use wv_saturation, only : wv_sat_final
    use cb24mjocnn,    only : cb24mjocnn_finalize
    use cb24cnn,       only : cb24cnn_finalize
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! Finalization of physics package
    !
    !-----------------------------------------------------------------------
    ! Input/output arguments
    type(physics_state), pointer :: phys_state(:)
    type(physics_tend ), pointer :: phys_tend(:)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    if(associated(pbuf2d)) then
       call pbuf_deallocate(pbuf2d,'global')
       deallocate(pbuf2d)
    end if
    deallocate(phys_state)
    deallocate(phys_tend)
    call chem_final
    call carma_final
    call wv_sat_final
    !call cb24mjocnn_finalize
    call cb24cnn_finalize
  end subroutine phys_final


  subroutine tphysac (ztodt,   cam_in,  &
       cam_out,  state,   tend,    pbuf)
    !-----------------------------------------------------------------------
    !
    ! Tendency physics after coupling to land, sea, and ice models.
    !
    ! Computes the following:
    !
    !   o Aerosol Emission at Surface
    !   o Source-Sink for Advected Tracers
    !   o Symmetric Turbulence Scheme - Vertical Diffusion
    !   o Rayleigh Friction
    !   o Dry Deposition of Aerosol
    !   o Enforce Charge Neutrality ( Only for WACCM )
    !   o Gravity Wave Drag
    !   o QBO Relaxation ( Only for WACCM )
    !   o Ion Drag ( Only for WACCM )
    !   o Scale Dry Mass Energy
    !-----------------------------------------------------------------------
    use physics_buffer, only: physics_buffer_desc, pbuf_set_field, pbuf_get_index, pbuf_get_field, pbuf_old_tim_idx
    use shr_kind_mod,       only: r8 => shr_kind_r8
    use chemistry,          only: chem_is_active, chem_timestep_tend, chem_emissions
    use cam_diagnostics,    only: diag_phys_tend_writeout
    use gw_drag,            only: gw_tend
    use vertical_diffusion, only: vertical_diffusion_tend
    use rayleigh_friction,  only: rayleigh_friction_tend
    use constituents,       only: cnst_get_ind
    use physics_types,      only: physics_state, physics_tend, physics_ptend, physics_update,    &
         physics_dme_adjust, set_dry_to_wet, physics_state_check
    use waccmx_phys_intr,   only: waccmx_phys_mspd_tend  ! WACCM-X major diffusion
    use waccmx_phys_intr,   only: waccmx_phys_ion_elec_temp_tend ! WACCM-X
    use aoa_tracers,        only: aoa_tracers_timestep_tend
    use physconst,          only: rhoh2o, latvap,latice
    use aero_model,         only: aero_model_drydep
    use carma_intr,         only: carma_emission_tend, carma_timestep_tend
    use carma_flags_mod,    only: carma_do_aerosol, carma_do_emission
    use check_energy,       only: check_energy_chng, calc_te_and_aam_budgets
    use check_energy,       only: check_tracers_data, check_tracers_init, check_tracers_chng
    use time_manager,       only: get_nstep
    use cam_abortutils,     only: endrun
    use dycore,             only: dycore_is
    use cam_control_mod,    only: aqua_planet
    use mo_gas_phase_chemdr,only: map2chm
    use clybry_fam,         only: clybry_fam_set
    use charge_neutrality,  only: charge_balance
    use qbo,                only: qbo_relax
    use iondrag,            only: iondrag_calc, do_waccm_ions
    use perf_mod
    use flux_avg,           only: flux_avg_run
    use unicon_cam,         only: unicon_cam_org_diags
    use cam_history,        only: hist_fld_active
    use qneg_module,        only: qneg4
    use co2_cycle,          only: co2_cycle_set_ptend
    use nudging,            only: Nudge_Model,Nudge_ON,nudging_timestep_tend

    !
    ! Arguments
    !
    real(r8), intent(in) :: ztodt                  ! Two times model timestep (2 delta-t)

    type(cam_in_t),      intent(inout) :: cam_in
    type(cam_out_t),     intent(inout) :: cam_out
    type(physics_state), intent(inout) :: state
    type(physics_tend ), intent(inout) :: tend
    type(physics_buffer_desc), pointer :: pbuf(:)


    type(check_tracers_data):: tracerint             ! tracer mass integrals and cummulative boundary fluxes

    !
    !---------------------------Local workspace-----------------------------
    !
    type(physics_ptend)     :: ptend               ! indivdual parameterization tendencies

    integer  :: nstep                              ! current timestep number
    real(r8) :: zero(pcols)                        ! array of zeros

    integer :: lchnk                                ! chunk identifier
    integer :: ncol                                 ! number of atmospheric columns
    integer i,k,m                 ! Longitude, level indices
    integer :: yr, mon, day, tod       ! components of a date
    integer :: ixcldice, ixcldliq      ! constituent indices for cloud liquid and ice water.

    logical :: labort                            ! abort flag

    real(r8) tvm(pcols,pver)           ! virtual temperature
    real(r8) prect(pcols)              ! total precipitation
    real(r8) surfric(pcols)            ! surface friction velocity
    real(r8) obklen(pcols)             ! Obukhov length
    real(r8) :: fh2o(pcols)            ! h2o flux to balance source from methane chemistry
    real(r8) :: flx_heat(pcols)        ! Heat flux for check_energy_chng.
    real(r8) :: tmp_q     (pcols,pver) ! tmp space
    real(r8) :: tmp_cldliq(pcols,pver) ! tmp space
    real(r8) :: tmp_cldice(pcols,pver) ! tmp space
    real(r8) :: tmp_trac  (pcols,pver,pcnst) ! tmp space
    real(r8) :: tmp_pdel  (pcols,pver) ! tmp space
    real(r8) :: tmp_ps    (pcols)      ! tmp space

    ! physics buffer fields for total energy and mass adjustment
    integer itim_old, ifld

    real(r8), pointer, dimension(:,:) :: cld
    real(r8), pointer, dimension(:,:) :: qini
    real(r8), pointer, dimension(:,:) :: cldliqini
    real(r8), pointer, dimension(:,:) :: cldiceini
    real(r8), pointer, dimension(:,:) :: dtcore
    real(r8), pointer, dimension(:,:) :: ast     ! relative humidity cloud fraction

    !-----------------------------------------------------------------------
    lchnk = state%lchnk
    ncol  = state%ncol

    nstep = get_nstep()

    ! Adjust the surface fluxes to reduce instabilities in near sfc layer
    if (phys_do_flux_avg()) then
       call flux_avg_run(state, cam_in,  pbuf, nstep, ztodt)
    endif

    ! Validate the physics state.
    if (state_debug_checks) &
         call physics_state_check(state, name="before tphysac")

    call t_startf('tphysac_init')
    ! Associate pointers with physics buffer fields
    itim_old = pbuf_old_tim_idx()


    ifld = pbuf_get_index('DTCORE')
    call pbuf_get_field(pbuf, ifld, dtcore, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

    call pbuf_get_field(pbuf, qini_idx, qini)
    call pbuf_get_field(pbuf, cldliqini_idx, cldliqini)
    call pbuf_get_field(pbuf, cldiceini_idx, cldiceini)

    ifld = pbuf_get_index('CLD')
    call pbuf_get_field(pbuf, ifld, cld, start=(/1,1,itim_old/),kount=(/pcols,pver,1/))

    ifld = pbuf_get_index('AST')
    call pbuf_get_field(pbuf, ifld, ast, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

    !
    ! accumulate fluxes into net flux array for spectral dycores
    ! jrm Include latent heat of fusion for snow
    !
    do i=1,ncol
       tend%flx_net(i) = tend%flx_net(i) + cam_in%shf(i) + (cam_out%precc(i) &
            + cam_out%precl(i))*latvap*rhoh2o &
            + (cam_out%precsc(i) + cam_out%precsl(i))*latice*rhoh2o
    end do

    ! emissions of aerosols and gas-phase chemistry constituents at surface
    call chem_emissions( state, cam_in )

    if (carma_do_emission) then
       ! carma emissions
       call carma_emission_tend (state, ptend, cam_in, ztodt)
       call physics_update(state, ptend, ztodt, tend)
    end if

    ! get nstep and zero array for energy checker
    zero = 0._r8
    nstep = get_nstep()
    call check_tracers_init(state, tracerint)

    ! Check if latent heat flux exceeds the total moisture content of the
    ! lowest model layer, thereby creating negative moisture.

    call qneg4('TPHYSAC', lchnk, ncol, ztodt ,                                &
         state%q(1,pver,1), state%rpdel(1,pver),                              &
         cam_in%shf, cam_in%lhf, cam_in%cflx)

    call t_stopf('tphysac_init')
    !===================================================
    ! Source/sink terms for advected tracers.
    !===================================================
    call t_startf('adv_tracer_src_snk')
    ! Test tracers

    call aoa_tracers_timestep_tend(state, ptend, cam_in%cflx, cam_in%landfrac, ztodt)
    call physics_update(state, ptend, ztodt, tend)
    call check_tracers_chng(state, tracerint, "aoa_tracers_timestep_tend", nstep, ztodt,   &
         cam_in%cflx)

    call co2_cycle_set_ptend(state, pbuf, ptend)
    call physics_update(state, ptend, ztodt, tend)

    !===================================================
    ! Chemistry and MAM calculation
    ! MAM core aerosol conversion process is performed in the below 'chem_timestep_tend'.
    ! In addition, surface flux of aerosol species other than 'dust' and 'sea salt', and
    ! elevated emission of aerosol species are treated in 'chem_timestep_tend' before
    ! Gas chemistry and MAM core aerosol conversion.
    ! Note that surface flux is not added into the atmosphere, but elevated emission is
    ! added into the atmosphere as tendency.
    !===================================================
    if (chem_is_active()) then
       call chem_timestep_tend(state, ptend, cam_in, cam_out, ztodt, &
            pbuf,  fh2o=fh2o)

       call physics_update(state, ptend, ztodt, tend)
       call check_energy_chng(state, tend, "chem", nstep, ztodt, fh2o, zero, zero, zero)
       call check_tracers_chng(state, tracerint, "chem_timestep_tend", nstep, ztodt, &
            cam_in%cflx)
    end if
    call t_stopf('adv_tracer_src_snk')

    !===================================================
    ! Vertical diffusion/pbl calculation
    ! Call vertical diffusion code (pbl, free atmosphere and molecular)
    !===================================================

    call t_startf('vertical_diffusion_tend')
    call vertical_diffusion_tend (ztodt ,state , cam_in, &
         surfric  ,obklen   ,ptend    ,ast    ,pbuf )

   !------------------------------------------
   ! Call major diffusion for extended model
   !------------------------------------------
    if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then
       call waccmx_phys_mspd_tend (ztodt    ,state    ,ptend)
    endif

    call physics_update(state, ptend, ztodt, tend)

    call t_stopf ('vertical_diffusion_tend')

    !===================================================
    ! Rayleigh friction calculation
    !===================================================
    call t_startf('rayleigh_friction')
    call rayleigh_friction_tend( ztodt, state, ptend)
    call physics_update(state, ptend, ztodt, tend)
    call t_stopf('rayleigh_friction')

    if (do_clubb_sgs) then
      call check_energy_chng(state, tend, "vdiff", nstep, ztodt, zero, zero, zero, zero)
    else
      call check_energy_chng(state, tend, "vdiff", nstep, ztodt, cam_in%cflx(:,1), zero, &
           zero, cam_in%shf)
    endif

    call check_tracers_chng(state, tracerint, "vdiff", nstep, ztodt, cam_in%cflx)

    !  aerosol dry deposition processes
    call t_startf('aero_drydep')
    call aero_model_drydep( state, pbuf, obklen, surfric, cam_in, ztodt, cam_out, ptend )
    call physics_update(state, ptend, ztodt, tend)
    call t_stopf('aero_drydep')

   ! CARMA microphysics
   !
   ! NOTE: This does both the timestep_tend for CARMA aerosols as well as doing the dry
   ! deposition for CARMA aerosols. It needs to follow vertical_diffusion_tend, so that
   ! obklen and surfric have been calculated. It needs to follow aero_model_drydep, so
   ! that cam_out%xxxdryxxx fields have already been set for CAM aerosols and cam_out
   ! can be added to for CARMA aerosols.
   if (carma_do_aerosol) then
     call t_startf('carma_timestep_tend')
     call carma_timestep_tend(state, cam_in, cam_out, ptend, ztodt, pbuf, obklen=obklen, ustar=surfric)
     call physics_update(state, ptend, ztodt, tend)

     call check_energy_chng(state, tend, "carma_tend", nstep, ztodt, zero, zero, zero, zero)
     call t_stopf('carma_timestep_tend')
   end if


    !---------------------------------------------------------------------------------
    !   ... enforce charge neutrality
    !---------------------------------------------------------------------------------
    call charge_balance(state, pbuf)

    !===================================================
    ! Gravity wave drag
    !===================================================
    call t_startf('gw_tend')

    call gw_tend(state, pbuf, ztodt, ptend, cam_in, flx_heat)

    call physics_update(state, ptend, ztodt, tend)
    ! Check energy integrals
    call check_energy_chng(state, tend, "gwdrag", nstep, ztodt, zero, &
         zero, zero, flx_heat)
    call t_stopf('gw_tend')

    ! QBO relaxation
    call qbo_relax(state, pbuf, ptend)
    call physics_update(state, ptend, ztodt, tend)
    ! Check energy integrals
    call check_energy_chng(state, tend, "qborelax", nstep, ztodt, zero, zero, zero, zero)

    ! Ion drag calculation
    call t_startf ( 'iondrag' )

    if ( do_waccm_ions ) then
       call iondrag_calc( lchnk, ncol, state, ptend, pbuf,  ztodt )
    else
       call iondrag_calc( lchnk, ncol, state, ptend)
    endif
    !----------------------------------------------------------------------------
    ! Call ionosphere routines for extended model if mode is set to ionosphere
    !----------------------------------------------------------------------------
    if( waccmx_is('ionosphere') ) then
       call waccmx_phys_ion_elec_temp_tend(state, ptend, pbuf, ztodt)
    endif

    call physics_update(state, ptend, ztodt, tend)
    call calc_te_and_aam_budgets(state, 'pAP')

    !---------------------------------------------------------------------------------
    ! Enforce charge neutrality after O+ change from ionos_tend
    !---------------------------------------------------------------------------------
    if( waccmx_is('ionosphere') ) then
       call charge_balance(state, pbuf)
    endif

    ! Check energy integrals
    call check_energy_chng(state, tend, "iondrag", nstep, ztodt, zero, zero, zero, zero)

    call t_stopf  ( 'iondrag' )

    ! Update Nudging values, if needed
    !----------------------------------
    !if((Nudge_Model).and.(Nudge_ON)) then
    !  call nudging_timestep_tend(state,ptend)
    !  call physics_update(state,ptend,ztodt,tend)
    !  call check_energy_chng(state, tend, "nudging", nstep, ztodt, zero, zero, zero, zero)
    !endif

    !-------------- Energy budget checks vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

    ! Save total energy for global fixer in next timestep (FV and SE dycores)
    call pbuf_set_field(pbuf, teout_idx, state%te_cur, (/1,itim_old/),(/pcols,1/))

    if (shallow_scheme .eq. 'UNICON') then

       ! ------------------------------------------------------------------------
       ! Insert the organization-related heterogeneities computed inside the
       ! UNICON into the tracer arrays here before performing advection.
       ! This is necessary to prevent any modifications of organization-related
       ! heterogeneities by non convection-advection process, such as
       ! dry and wet deposition of aerosols, MAM, etc.
       ! Again, note that only UNICON and advection schemes are allowed to
       ! changes to organization at this stage, although we can include the
       ! effects of other physical processes in future.
       ! ------------------------------------------------------------------------

       call unicon_cam_org_diags(state, pbuf)

    end if
    !
    ! FV: convert dry-type mixing ratios to moist here because physics_dme_adjust
    !     assumes moist. This is done in p_d_coupling for other dynamics. Bundy, Feb 2004.
    if ( dycore_is('LR')) call set_dry_to_wet(state)    ! Physics had dry, dynamics wants moist

    ! Scale dry mass and energy (does nothing if dycore is EUL or SLD)
    call cnst_get_ind('CLDLIQ', ixcldliq)
    call cnst_get_ind('CLDICE', ixcldice)
    tmp_q     (:ncol,:pver) = state%q(:ncol,:pver,1)
    tmp_cldliq(:ncol,:pver) = state%q(:ncol,:pver,ixcldliq)
    tmp_cldice(:ncol,:pver) = state%q(:ncol,:pver,ixcldice)
    ! For not 'FV', physics_dme_adjust is called for energy diagnostic purposes only.  So, save off tracers
    if (.not.dycore_is('FV').and.&
         (hist_fld_active('SE_pAM').or.hist_fld_active('KE_pAM').or.hist_fld_active('WV_pAM').or.&
         hist_fld_active('WL_pAM').or.hist_fld_active('WI_pAM'))) then
      tmp_trac(:ncol,:pver,:pcnst) = state%q(:ncol,:pver,:pcnst)
      tmp_pdel(:ncol,:pver)        = state%pdel(:ncol,:pver)
      tmp_ps(:ncol)                = state%ps(:ncol)
      !
      ! pint, lnpint,rpdel are altered by dme_adjust but not used for tendencies in dynamics of SE
      ! we do not reset them to pre-dme_adjust values
      !
      if (dycore_is('SE')) call set_dry_to_wet(state)
      call physics_dme_adjust(state, tend, qini, ztodt)
      call calc_te_and_aam_budgets(state, 'pAM')
      ! Restore pre-"physics_dme_adjust" tracers
      state%q(:ncol,:pver,:pcnst) = tmp_trac(:ncol,:pver,:pcnst)
      state%pdel(:ncol,:pver)     = tmp_pdel(:ncol,:pver)
      state%ps(:ncol)             = tmp_ps(:ncol)
    end if

    if (dycore_is('LR')) then
      call physics_dme_adjust(state, tend, qini, ztodt)
      call calc_te_and_aam_budgets(state, 'pAM')
    endif

    !!!   REMOVE THIS CALL, SINCE ONLY Q IS BEING ADJUSTED. WON'T BALANCE ENERGY. TE IS SAVED BEFORE THIS
    !!!   call check_energy_chng(state, tend, "drymass", nstep, ztodt, zero, zero, zero, zero)

    ! store T in buffer for use in computing dynamics T-tendency in next timestep
    do k = 1,pver
       dtcore(:ncol,k) = state%t(:ncol,k)
    end do

    !-------------- Energy budget checks ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    if (aqua_planet) then
       labort = .false.
       do i=1,ncol
          if (cam_in%ocnfrac(i) /= 1._r8) labort = .true.
       end do
       if (labort) then
          call endrun ('TPHYSAC error:  grid contains non-ocean point')
       endif
    endif

    call diag_phys_tend_writeout (state, pbuf,  tend, ztodt, tmp_q, tmp_cldliq, tmp_cldice, &
         qini, cldliqini, cldiceini)

    call clybry_fam_set( ncol, lchnk, map2chm, state%q, pbuf )

  end subroutine tphysac

  subroutine tphysbc (ztodt, state,  &
       tend,    pbuf,              &
       cam_out, cam_in )
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    ! Evaluate and apply physical processes that are calculated BEFORE
    ! coupling to land, sea, and ice models.
    !
    ! Processes currently included are:
    !
    !  o Resetting Negative Tracers to Positive
    !  o Global Mean Total Energy Fixer
    !  o Dry Adjustment
    !  o Asymmetric Turbulence Scheme : Deep Convection & Shallow Convection
    !  o Stratiform Macro-Microphysics
    !  o Wet Scavenging of Aerosol
    !  o Radiation
    !
    ! Method:
    !
    ! Each parameterization should be implemented with this sequence of calls:
    !  1)  Call physics interface
    !  2)  Check energy
    !  3)  Call physics_update
    ! See Interface to Column Physics and Chemistry Packages
    !   http://www.ccsm.ucar.edu/models/atm-cam/docs/phys-interface/index.html
    !
    !-----------------------------------------------------------------------

    use physics_buffer,  only: physics_buffer_desc, pbuf_get_field
    use physics_buffer,  only: pbuf_get_index, pbuf_old_tim_idx
    use physics_buffer,  only: col_type_subcol, dyn_time_lvls
    use shr_kind_mod,    only: r8 => shr_kind_r8

    use dadadj_cam,      only: dadadj_tend
    use rk_stratiform,   only: rk_stratiform_tend
    use microp_driver,   only: microp_driver_tend
    use microp_aero,     only: microp_aero_run
    use macrop_driver,   only: macrop_driver_tend
    use physics_types,   only: physics_state, physics_tend, physics_ptend, &
         physics_update, physics_ptend_init, physics_ptend_sum, &
         physics_state_check, physics_ptend_scale
    use cam_diagnostics, only: diag_conv_tend_ini, diag_phys_writeout, diag_conv, diag_export, diag_state_b4_phys_write
    use cam_history,     only: outfld
    use physconst,       only: cpair, latvap
    use constituents,    only: pcnst, qmin, cnst_get_ind
    use convect_deep,    only: convect_deep_tend, convect_deep_tend_2, deep_scheme_does_scav_trans
    use time_manager,    only: is_first_step, get_nstep
    use convect_shallow, only: convect_shallow_tend
    use check_energy,    only: check_energy_chng, check_energy_fix, check_energy_timestep_init
    use check_energy,    only: check_tracers_data, check_tracers_init, check_tracers_chng
    use check_energy,    only: calc_te_and_aam_budgets
    use dycore,          only: dycore_is
    use aero_model,      only: aero_model_wetdep
    use carma_intr,      only: carma_wetdep_tend, carma_timestep_tend
    use carma_flags_mod, only: carma_do_detrain, carma_do_cldice, carma_do_cldliq,  carma_do_wetdep
    use radiation,       only: radiation_tend
    use cloud_diagnostics, only: cloud_diagnostics_calc
    use perf_mod
    use mo_gas_phase_chemdr,only: map2chm
    use clybry_fam,         only: clybry_fam_adj
    use clubb_intr,      only: clubb_tend_cam
    use sslt_rebin,      only: sslt_rebin_adv
    use tropopause,      only: tropopause_output
    use cam_abortutils,  only: endrun
    use subcol,          only: subcol_gen, subcol_ptend_avg
    use subcol_utils,    only: subcol_ptend_copy, is_subcol_on
    use qneg_module,     only: qneg3
    use nudging,         only: Nudge_Model,Nudge_ON,nudging_timestep_tend !++WEC
    use cb24mjocnn,      only: cb24mjocnn_timestep_tend,cb24mjocnn_timestep_init,cb24mjocnn_Model !++WEC
    use cb24cnn,         only: cb24cnn_timestep_init !++WEC

    ! Arguments

    real(r8), intent(in) :: ztodt                          ! 2 delta t (model time increment)

    type(physics_state), intent(inout) :: state
    type(physics_tend ), intent(inout) :: tend
    type(physics_buffer_desc), pointer :: pbuf(:)

    type(cam_out_t),     intent(inout) :: cam_out
    type(cam_in_t),      intent(in)    :: cam_in


    !
    !---------------------------Local workspace-----------------------------
    !

    type(physics_ptend)   :: ptend            ! indivdual parameterization tendencies
    type(physics_state)   :: state_sc         ! state for sub-columns
    type(physics_ptend)   :: ptend_sc         ! ptend for sub-columns
    type(physics_ptend)   :: ptend_aero       ! ptend for microp_aero
    type(physics_ptend)   :: ptend_aero_sc    ! ptend for microp_aero on sub-columns
    type(physics_tend)    :: tend_sc          ! tend for sub-columns

    integer :: nstep                          ! current timestep number

    real(r8) :: net_flx(pcols)

    real(r8) :: zdu(pcols,pver)               ! detraining mass flux from deep convection
    real(r8) :: cmfmc(pcols,pverp)            ! Convective mass flux--m sub c

    real(r8) cmfcme(pcols,pver)                ! cmf condensation - evaporation

    real(r8) dlf(pcols,pver)                   ! Detraining cld H20 from shallow + deep convections
    real(r8) dlf2(pcols,pver)                  ! Detraining cld H20 from shallow convections
    real(r8) pflx(pcols,pverp)                 ! Conv rain flux thru out btm of lev

    integer lchnk                              ! chunk identifier
    integer ncol                               ! number of atmospheric columns

    integer :: i                               ! column indicex
    integer :: ixcldice, ixcldliq              ! constituent indices for cloud liquid and ice water.
    ! for macro/micro co-substepping
    integer :: macmic_it                       ! iteration variables
    real(r8) :: cld_macmic_ztodt               ! modified timestep
    ! physics buffer fields to compute tendencies for stratiform package
    integer itim_old, ifld
    real(r8), pointer, dimension(:,:) :: cld        ! cloud fraction


    ! physics buffer fields for total energy and mass adjustment
    real(r8), pointer, dimension(:  ) :: teout
    real(r8), pointer, dimension(:,:) :: qini
    real(r8), pointer, dimension(:,:) :: cldliqini
    real(r8), pointer, dimension(:,:) :: cldiceini
    real(r8), pointer, dimension(:,:) :: dtcore

    real(r8), pointer, dimension(:,:,:) :: fracis  ! fraction of transported species that are insoluble

    real(r8), pointer :: dlfzm(:,:)                ! ZM detrained convective cloud water mixing ratio.

    ! convective precipitation variables
    real(r8),pointer :: prec_dp(:)                ! total precipitation from ZM convection
    real(r8),pointer :: snow_dp(:)                ! snow from ZM convection
    real(r8),pointer :: prec_sh(:)                ! total precipitation from Hack convection
    real(r8),pointer :: snow_sh(:)                ! snow from Hack convection

    ! carma precipitation variables
    real(r8) :: prec_sed_carma(pcols)          ! total precip from cloud sedimentation (CARMA)
    real(r8) :: snow_sed_carma(pcols)          ! snow from cloud ice sedimentation (CARMA)

    ! stratiform precipitation variables
    real(r8),pointer :: prec_str(:)    ! sfc flux of precip from stratiform (m/s)
    real(r8),pointer :: snow_str(:)     ! sfc flux of snow from stratiform   (m/s)
    real(r8),pointer :: prec_str_sc(:)  ! sfc flux of precip from stratiform (m/s) -- for subcolumns
    real(r8),pointer :: snow_str_sc(:)  ! sfc flux of snow from stratiform   (m/s) -- for subcolumns
    real(r8),pointer :: prec_pcw(:)     ! total precip from prognostic cloud scheme
    real(r8),pointer :: snow_pcw(:)     ! snow from prognostic cloud scheme
    real(r8),pointer :: prec_sed(:)     ! total precip from cloud sedimentation
    real(r8),pointer :: snow_sed(:)     ! snow from cloud ice sedimentation

    ! Local copies for substepping
    real(r8) :: prec_pcw_macmic(pcols)
    real(r8) :: snow_pcw_macmic(pcols)
    real(r8) :: prec_sed_macmic(pcols)
    real(r8) :: snow_sed_macmic(pcols)

    ! energy checking variables
    real(r8) :: zero(pcols)                    ! array of zeros
    real(r8) :: zero_sc(pcols*psubcols)        ! array of zeros
    real(r8) :: rliq(pcols)                    ! vertical integral of liquid not yet in q(ixcldliq)
    real(r8) :: rice(pcols)                    ! vertical integral of ice not yet in q(ixcldice)
    real(r8) :: rliq2(pcols)                   ! vertical integral of liquid from shallow scheme
    real(r8) :: det_s  (pcols)                 ! vertical integral of detrained static energy from ice
    real(r8) :: det_ice(pcols)                 ! vertical integral of detrained ice
    real(r8) :: flx_cnd(pcols)
    real(r8) :: flx_heat(pcols)
    type(check_tracers_data):: tracerint             ! energy integrals and cummulative boundary fluxes
    real(r8) :: zero_tracers(pcols,pcnst)

    logical   :: lq(pcnst)
    !-----------------------------------------------------------------------

    call t_startf('bc_init')

    zero = 0._r8
    zero_tracers(:,:) = 0._r8
    zero_sc(:) = 0._r8

    lchnk = state%lchnk
    ncol  = state%ncol

    nstep = get_nstep()

    ! Associate pointers with physics buffer fields
    itim_old = pbuf_old_tim_idx()
    ifld = pbuf_get_index('CLD')
    call pbuf_get_field(pbuf, ifld, cld, (/1,1,itim_old/),(/pcols,pver,1/))

    call pbuf_get_field(pbuf, teout_idx, teout, (/1,itim_old/), (/pcols,1/))

    call pbuf_get_field(pbuf, qini_idx, qini)
    call pbuf_get_field(pbuf, cldliqini_idx, cldliqini)
    call pbuf_get_field(pbuf, cldiceini_idx, cldiceini)

    ifld   =  pbuf_get_index('DTCORE')
    call pbuf_get_field(pbuf, ifld, dtcore, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

    ifld    = pbuf_get_index('FRACIS')
    call pbuf_get_field(pbuf, ifld, fracis, start=(/1,1,1/), kount=(/pcols, pver, pcnst/)  )
    fracis (:ncol,:,1:pcnst) = 1._r8

    ! Set physics tendencies to 0
    tend %dTdt(:ncol,:pver)  = 0._r8
    tend %dudt(:ncol,:pver)  = 0._r8
    tend %dvdt(:ncol,:pver)  = 0._r8

    ! Verify state coming from the dynamics
    if (state_debug_checks) &
         call physics_state_check(state, name="before tphysbc (dycore?)")

    call clybry_fam_adj( ncol, lchnk, map2chm, state%q, pbuf )

    ! Since clybry_fam_adj operates directly on the tracers, and has no
    ! physics_update call, re-run qneg3.

    call qneg3('TPHYSBCc',lchnk  ,ncol    ,pcols   ,pver    , &
         1, pcnst, qmin  ,state%q )

    ! Validate output of clybry_fam_adj.
    if (state_debug_checks) &
         call physics_state_check(state, name="clybry_fam_adj")

    !
    ! Dump out "before physics" state
    !
    call diag_state_b4_phys_write (state)

    ! compute mass integrals of input tracers state
    call check_tracers_init(state, tracerint)

    call t_stopf('bc_init')

    !===================================================
    ! Global mean total energy fixer
    !===================================================
    call t_startf('energy_fixer')

    call calc_te_and_aam_budgets(state, 'pBF')
    if (dycore_is('LR') .or. dycore_is('SE'))  then
       call check_energy_fix(state, ptend, nstep, flx_heat)
       call physics_update(state, ptend, ztodt, tend)
       call check_energy_chng(state, tend, "chkengyfix", nstep, ztodt, zero, zero, zero, flx_heat)
       call outfld( 'EFIX', flx_heat    , pcols, lchnk   )
    end if
    call calc_te_and_aam_budgets(state, 'pBP')
    ! Save state for convective tendency calculations.
    call diag_conv_tend_ini(state, pbuf)

    call cnst_get_ind('CLDLIQ', ixcldliq)
    call cnst_get_ind('CLDICE', ixcldice)
    qini     (:ncol,:pver) = state%q(:ncol,:pver,       1)
    cldliqini(:ncol,:pver) = state%q(:ncol,:pver,ixcldliq)
    cldiceini(:ncol,:pver) = state%q(:ncol,:pver,ixcldice)

    call outfld('TEOUT', teout       , pcols, lchnk   )
    call outfld('TEINP', state%te_ini, pcols, lchnk   )
    call outfld('TEFIX', state%te_cur, pcols, lchnk   )

    ! T tendency due to dynamics
    if( nstep > dyn_time_lvls-1 ) then
       dtcore(:ncol,:pver) = (state%t(:ncol,:pver) - dtcore(:ncol,:pver))/ztodt
       call outfld( 'DTCORE', dtcore, pcols, lchnk )
    end if

    call t_stopf('energy_fixer')
    !
    !===================================================
    ! Dry adjustment
    !===================================================
    call t_startf('dry_adjustment')

    call dadadj_tend(ztodt, state, ptend)

    call physics_update(state, ptend, ztodt, tend)

    call t_stopf('dry_adjustment')

    !===================================================
    ! Moist convection
    !===================================================
    call t_startf('moist_convection')

    call t_startf ('convect_deep_tend')

    call convect_deep_tend(  &
         cmfmc,      cmfcme,             &
         pflx,    zdu,       &
         rliq,    rice,      &
         ztodt,   &
         state,   ptend, cam_in%landfrac, pbuf)

    call physics_update(state, ptend, ztodt, tend)

    call t_stopf('convect_deep_tend')

    call pbuf_get_field(pbuf, prec_dp_idx, prec_dp )
    call pbuf_get_field(pbuf, snow_dp_idx, snow_dp )
    call pbuf_get_field(pbuf, prec_sh_idx, prec_sh )
    call pbuf_get_field(pbuf, snow_sh_idx, snow_sh )
    call pbuf_get_field(pbuf, prec_str_idx, prec_str )
    call pbuf_get_field(pbuf, snow_str_idx, snow_str )
    call pbuf_get_field(pbuf, prec_sed_idx, prec_sed )
    call pbuf_get_field(pbuf, snow_sed_idx, snow_sed )
    call pbuf_get_field(pbuf, prec_pcw_idx, prec_pcw )
    call pbuf_get_field(pbuf, snow_pcw_idx, snow_pcw )

    if (use_subcol_microp) then
      call pbuf_get_field(pbuf, prec_str_idx, prec_str_sc, col_type=col_type_subcol)
      call pbuf_get_field(pbuf, snow_str_idx, snow_str_sc, col_type=col_type_subcol)
    end if

    ! Check energy integrals, including "reserved liquid"
    flx_cnd(:ncol) = prec_dp(:ncol) + rliq(:ncol)
    snow_dp(:ncol) = snow_dp(:ncol) + rice(:ncol)
    call check_energy_chng(state, tend, "convect_deep", nstep, ztodt, zero, flx_cnd, snow_dp, zero)
    snow_dp(:ncol) = snow_dp(:ncol) - rice(:ncol)

    !
    ! Call Hack (1994) convection scheme to deal with shallow/mid-level convection
    !
    call t_startf ('convect_shallow_tend')

    if (dlfzm_idx > 0) then
       call pbuf_get_field(pbuf, dlfzm_idx, dlfzm)
       dlf(:ncol,:) = dlfzm(:ncol,:)
    else
       dlf(:,:) = 0._r8
    end if

    call convect_shallow_tend (ztodt   , cmfmc, &
         dlf        , dlf2   ,  rliq   , rliq2, &
         state      , ptend  ,  pbuf, cam_in)
    call t_stopf ('convect_shallow_tend')

    call physics_update(state, ptend, ztodt, tend)

    flx_cnd(:ncol) = prec_sh(:ncol) + rliq2(:ncol)
    call check_energy_chng(state, tend, "convect_shallow", nstep, ztodt, zero, flx_cnd, snow_sh, zero)

    call check_tracers_chng(state, tracerint, "convect_shallow", nstep, ztodt, zero_tracers)

    call t_stopf('moist_convection')

    ! Rebin the 4-bin version of sea salt into bins for coarse and accumulation
    ! modes that correspond to the available optics data.  This is only necessary
    ! for CAM-RT.  But it's done here so that the microphysics code which is called
    ! from the stratiform interface has access to the same aerosols as the radiation
    ! code.
    call sslt_rebin_adv(pbuf,  state)

    !===================================================
    ! Calculate tendencies from CARMA bin microphysics.
    !===================================================
    !
    ! If CARMA is doing detrainment, then on output, rliq no longer represents water reserved
    ! for detrainment, but instead represents potential snow fall. The mass and number of the
    ! snow are stored in the physics buffer and will be incorporated by the MG microphysics.
    !
    ! Currently CARMA cloud microphysics is only supported with the MG microphysics.
    call t_startf('carma_timestep_tend')

    if (carma_do_cldice .or. carma_do_cldliq) then
       call carma_timestep_tend(state, cam_in, cam_out, ptend, ztodt, pbuf, dlf=dlf, rliq=rliq, &
            prec_str=prec_str, snow_str=snow_str, prec_sed=prec_sed_carma, snow_sed=snow_sed_carma)
       call physics_update(state, ptend, ztodt, tend)

       ! Before the detrainment, the reserved condensate is all liquid, but if CARMA is doing
       ! detrainment, then the reserved condensate is snow.
       if (carma_do_detrain) then
          call check_energy_chng(state, tend, "carma_tend", nstep, ztodt, zero, prec_str+rliq, snow_str+rliq, zero)
       else
          call check_energy_chng(state, tend, "carma_tend", nstep, ztodt, zero, prec_str, snow_str, zero)
       end if
    end if

    call t_stopf('carma_timestep_tend')

    if( microp_scheme == 'RK' ) then

       !===================================================
       ! Calculate stratiform tendency (sedimentation, detrain, cloud fraction and microphysics )
       !===================================================
       call t_startf('rk_stratiform_tend')

       call rk_stratiform_tend(state, ptend, pbuf, ztodt, &
            cam_in%icefrac, cam_in%landfrac, cam_in%ocnfrac, &
            cam_in%snowhland, & ! sediment
            dlf, dlf2, & ! detrain
            rliq  , & ! check energy after detrain
            cmfmc,  &
            cam_in%ts,      cam_in%sst,        zdu)

       call physics_update(state, ptend, ztodt, tend)
       call check_energy_chng(state, tend, "cldwat_tend", nstep, ztodt, zero, prec_str, snow_str, zero)

       call t_stopf('rk_stratiform_tend')

    elseif( microp_scheme == 'MG' ) then
       ! Start co-substepping of macrophysics and microphysics
       cld_macmic_ztodt = ztodt/cld_macmic_num_steps

       ! Clear precip fields that should accumulate.
       prec_sed_macmic = 0._r8
       snow_sed_macmic = 0._r8
       prec_pcw_macmic = 0._r8
       snow_pcw_macmic = 0._r8

       do macmic_it = 1, cld_macmic_num_steps

          !===================================================
          ! Calculate macrophysical tendency (sedimentation, detrain, cloud fraction)
          !===================================================

          call t_startf('macrop_tend')

          ! don't call Park macrophysics if CLUBB is called
          if (macrop_scheme .ne. 'CLUBB_SGS') then

             call macrop_driver_tend( &
                  state,           ptend,          cld_macmic_ztodt, &
                  cam_in%landfrac, cam_in%ocnfrac, cam_in%snowhland, & ! sediment
                  dlf,             dlf2,                             & ! detrain
                  cmfmc,                                             &
                  cam_in%ts,       cam_in%sst,     zdu,              &
                  pbuf,            det_s,          det_ice)

             ! Since we "added" the reserved liquid back in this routine, we need
             ! to account for it in the energy checker
             flx_cnd(:ncol) = -1._r8*rliq(:ncol)
             flx_heat(:ncol) = det_s(:ncol)

             ! Unfortunately, physics_update does not know what time period
             ! "tend" is supposed to cover, and therefore can't update it
             ! with substeps correctly. For now, work around this by scaling
             ! ptend down by the number of substeps, then applying it for
             ! the full time (ztodt).
             call physics_ptend_scale(ptend, 1._r8/cld_macmic_num_steps, ncol)
             call physics_update(state, ptend, ztodt, tend)
             call check_energy_chng(state, tend, "macrop_tend", nstep, ztodt, &
                  zero, flx_cnd(:ncol)/cld_macmic_num_steps, &
                  det_ice(:ncol)/cld_macmic_num_steps, &
                  flx_heat(:ncol)/cld_macmic_num_steps)

          else ! Calculate CLUBB macrophysics

             ! =====================================================
             !    CLUBB call (PBL, shallow convection, macrophysics)
             ! =====================================================

             call clubb_tend_cam(state, ptend, pbuf, cld_macmic_ztodt,&
                cmfmc, cam_in, macmic_it, cld_macmic_num_steps, &
                dlf, det_s, det_ice)

             ! Since we "added" the reserved liquid back in this routine, we need
             ! to account for it in the energy checker
             flx_cnd(:ncol) = -1._r8*rliq(:ncol)
             flx_heat(:ncol) = cam_in%shf(:ncol) + det_s(:ncol)

             ! Unfortunately, physics_update does not know what time period
             ! "tend" is supposed to cover, and therefore can't update it
             ! with substeps correctly. For now, work around this by scaling
             ! ptend down by the number of substeps, then applying it for
             ! the full time (ztodt).
             call physics_ptend_scale(ptend, 1._r8/cld_macmic_num_steps, ncol)

             ! Update physics tendencies and copy state to state_eq, because that is
             ! input for microphysics
             call physics_update(state, ptend, ztodt, tend)

             ! Use actual qflux (not lhf/latvap) for consistency with surface fluxes and revised code
             call check_energy_chng(state, tend, "clubb_tend", nstep, ztodt, &
                cam_in%cflx(:ncol,1)/cld_macmic_num_steps, &
                flx_cnd(:ncol)/cld_macmic_num_steps, &
                det_ice(:ncol)/cld_macmic_num_steps, &
                flx_heat(:ncol)/cld_macmic_num_steps)

          endif

          call t_stopf('macrop_tend')

          !===================================================
          ! Calculate cloud microphysics
          !===================================================

          if (is_subcol_on()) then
             ! Allocate sub-column structures.
             call physics_state_alloc(state_sc, lchnk, psubcols*pcols)
             call physics_tend_alloc(tend_sc, psubcols*pcols)

             ! Generate sub-columns using the requested scheme
             call subcol_gen(state, tend, state_sc, tend_sc, pbuf)

             !Initialize check energy for subcolumns
             call check_energy_timestep_init(state_sc, tend_sc, pbuf, col_type_subcol)
          end if

          call t_startf('microp_aero_run')
          call microp_aero_run(state, ptend_aero, cld_macmic_ztodt, pbuf)
          call t_stopf('microp_aero_run')

          call t_startf('microp_tend')

          if (use_subcol_microp) then
             call microp_driver_tend(state_sc, ptend_sc, cld_macmic_ztodt, pbuf)

             ! Average the sub-column ptend for use in gridded update - will not contain ptend_aero
             call subcol_ptend_avg(ptend_sc, state_sc%ngrdcol, lchnk, ptend)

             ! Copy ptend_aero field to one dimensioned by sub-columns before summing with ptend
             call subcol_ptend_copy(ptend_aero, state_sc, ptend_aero_sc)
             call physics_ptend_sum(ptend_aero_sc, ptend_sc, state_sc%ncol)
             call physics_ptend_dealloc(ptend_aero_sc)

             ! Have to scale and apply for full timestep to get tend right
             ! (see above note for macrophysics).
             call physics_ptend_scale(ptend_sc, 1._r8/cld_macmic_num_steps, ncol)

             call physics_update (state_sc, ptend_sc, ztodt, tend_sc)
             call check_energy_chng(state_sc, tend_sc, "microp_tend_subcol", &
                  nstep, ztodt, zero_sc, &
                  prec_str_sc(:state_sc%ncol)/cld_macmic_num_steps, &
                  snow_str_sc(:state_sc%ncol)/cld_macmic_num_steps, zero_sc)

             call physics_state_dealloc(state_sc)
             call physics_tend_dealloc(tend_sc)
             call physics_ptend_dealloc(ptend_sc)
          else
             call microp_driver_tend(state, ptend, cld_macmic_ztodt, pbuf)
          end if
          ! combine aero and micro tendencies for the grid
          call physics_ptend_sum(ptend_aero, ptend, ncol)
          call physics_ptend_dealloc(ptend_aero)

          ! Have to scale and apply for full timestep to get tend right
          ! (see above note for macrophysics).
          call physics_ptend_scale(ptend, 1._r8/cld_macmic_num_steps, ncol)

          call physics_update (state, ptend, ztodt, tend)
          call check_energy_chng(state, tend, "microp_tend", nstep, ztodt, &
               zero, prec_str(:ncol)/cld_macmic_num_steps, &
               snow_str(:ncol)/cld_macmic_num_steps, zero)

          call t_stopf('microp_tend')
          prec_sed_macmic(:ncol) = prec_sed_macmic(:ncol) + prec_sed(:ncol)
          snow_sed_macmic(:ncol) = snow_sed_macmic(:ncol) + snow_sed(:ncol)
          prec_pcw_macmic(:ncol) = prec_pcw_macmic(:ncol) + prec_pcw(:ncol)
          snow_pcw_macmic(:ncol) = snow_pcw_macmic(:ncol) + snow_pcw(:ncol)

       end do ! end substepping over macrophysics/microphysics

       prec_sed(:ncol) = prec_sed_macmic(:ncol)/cld_macmic_num_steps
       snow_sed(:ncol) = snow_sed_macmic(:ncol)/cld_macmic_num_steps
       prec_pcw(:ncol) = prec_pcw_macmic(:ncol)/cld_macmic_num_steps
       snow_pcw(:ncol) = snow_pcw_macmic(:ncol)/cld_macmic_num_steps
       prec_str(:ncol) = prec_pcw(:ncol) + prec_sed(:ncol)
       snow_str(:ncol) = snow_pcw(:ncol) + snow_sed(:ncol)

    endif

    ! Add the precipitation from CARMA to the precipitation from stratiform.
    if (carma_do_cldice .or. carma_do_cldliq) then
       prec_sed(:ncol) = prec_sed(:ncol) + prec_sed_carma(:ncol)
       snow_sed(:ncol) = snow_sed(:ncol) + snow_sed_carma(:ncol)
    end if

    if ( .not. deep_scheme_does_scav_trans() ) then

       ! -------------------------------------------------------------------------------
       ! 1. Wet Scavenging of Aerosols by Convective and Stratiform Precipitation.
       ! 2. Convective Transport of Non-Water Aerosol Species.
       !
       !  . Aerosol wet chemistry determines scavenging fractions, and transformations
       !  . Then do convective transport of all trace species except qv,ql,qi.
       !  . We needed to do the scavenging first to determine the interstitial fraction.
       !  . When UNICON is used as unified convection, we should still perform
       !    wet scavenging but not 'convect_deep_tend2'.
       ! -------------------------------------------------------------------------------

       call t_startf('bc_aerosols')
       if (clim_modal_aero .and. .not. prog_modal_aero) then
          call modal_aero_calcsize_diag(state, pbuf)
          call modal_aero_wateruptake_dr(state, pbuf)
       endif
       call aero_model_wetdep( state, ztodt, dlf, cam_out, ptend, pbuf)
       call physics_update(state, ptend, ztodt, tend)


       if (carma_do_wetdep) then
          ! CARMA wet deposition
          !
          ! NOTE: It needs to follow aero_model_wetdep, so that cam_out%xxxwetxxx
          ! fields have already been set for CAM aerosols and cam_out can be added
          ! to for CARMA aerosols.
          call t_startf ('carma_wetdep_tend')
          call carma_wetdep_tend(state, ptend, ztodt, pbuf, dlf, cam_out)
          call physics_update(state, ptend, ztodt, tend)
          call t_stopf ('carma_wetdep_tend')
       end if

       call t_startf ('convect_deep_tend2')
       call convect_deep_tend_2( state,   ptend,  ztodt,  pbuf )
       call physics_update(state, ptend, ztodt, tend)
       call t_stopf ('convect_deep_tend2')

       ! check tracer integrals
       call check_tracers_chng(state, tracerint, "cmfmca", nstep, ztodt,  zero_tracers)

       call t_stopf('bc_aerosols')

   endif
   
    !++WEC 
    ! Update Nudging values, if needed ++WEC
    !----------------------------------
    if((Nudge_Model).and.(Nudge_ON)) then
      call nudging_timestep_tend(state,ptend)
      call physics_update(state,ptend,ztodt,tend)
      call check_energy_chng(state, tend, "nudging", nstep, ztodt, zero, zero, zero, zero)
    endif

    

    !if((cb24mjocnn_Model).and.(cb24mjocnn_Model)) then
    !  call t_startf('cb24mjocnn_init')
    !  call cb24mjocnn_timestep_init(state,ptend,pbuf,cam_in,cam_out)
    !  if (masterproc) write(iulog,*)'nintendo past cb24mjocnn_timestep_init '
    !  call t_stopf('cb24mjocnn_init')
    !  call cb24mjocnn_timestep_tend(state,ptend)
    !  call physics_update(state,ptend,ztodt,tend)
    !  call check_energy_chng(state, tend, "CNNmjo", nstep, ztodt, zero, zero, zero, zero)
    !endif
    !===================================================
    ! Moist physical parameteriztions complete:
    ! send dynamical variables, and derived variables to history file
    !===================================================

    call t_startf('cb24cnn_init')
    call cb24cnn_timestep_init(state,ptend,pbuf,cam_in,cam_out)
    if (masterproc) write(iulog,*)'nintendo past cb24cnn_timestep_init '
    !call cb24cnn_timestep_tend(state,ptend,pbuf,cam_in,cam_out) !WEC please turn off later...
    call t_stopf('cb24cnn_init')
    

    call t_startf('bc_history_write')
    call diag_phys_writeout(state, pbuf)
    call diag_conv(state, ztodt, pbuf)

    call t_stopf('bc_history_write')

    !===================================================
    ! Write cloud diagnostics on history file
    !===================================================

    call t_startf('bc_cld_diag_history_write')

    call cloud_diagnostics_calc(state, pbuf)

    call t_stopf('bc_cld_diag_history_write')

    !===================================================
    ! Radiation computations
    !===================================================
    call t_startf('radiation')


    call radiation_tend( &
       state, ptend, pbuf, cam_out, cam_in, net_flx)

    ! Set net flux used by spectral dycores
    do i=1,ncol
       tend%flx_net(i) = net_flx(i)
    end do
    call physics_update(state, ptend, ztodt, tend)
    call check_energy_chng(state, tend, "radheat", nstep, ztodt, zero, zero, zero, net_flx)

    call t_stopf('radiation')

    ! Diagnose the location of the tropopause and its location to the history file(s).
    call t_startf('tropopause')
    call tropopause_output(state)
    call t_stopf('tropopause')

    ! Save atmospheric fields to force surface models
    call t_startf('cam_export')
    call cam_export (state,cam_out,pbuf)
    call t_stopf('cam_export')

    ! Write export state to history file
    call t_startf('diag_export')
    call diag_export(cam_out)
    call t_stopf('diag_export')

  end subroutine tphysbc

  subroutine tphysbc1(ztodt, fsns, fsnt, flns, flnt, &
                      state, tend, pbuf, fsds, &
                      sgh, sgh30, cam_out, cam_in )
    !----------------------------------------------------------------------------- 
    ! Purpose: Evaluate physics processes BEFORE coupling to sfc components
    !          Phase 1 - energy fixer and dry adjustment
    !
    ! Pass surface fields for separate surface flux calculations
    ! Dump appropriate fields to history file.
    !-----------------------------------------------------------------------------
    use physics_buffer,         only: physics_buffer_desc, pbuf_get_field
    use physics_buffer,         only: pbuf_get_index, pbuf_old_tim_idx
    use physics_buffer,         only: dyn_time_lvls, pbuf_set_field
    use physics_types,          only: physics_ptend_init, physics_ptend_sum, &
                                    physics_state_check, physics_ptend_scale
    use cam_diagnostics,        only: diag_conv_tend_ini, diag_state_b4_phys_write
    use cam_history,            only: outfld, fieldname_len
    use physconst,              only: cpair, latvap, rga
    use constituents,           only: pcnst, qmin, cnst_get_ind
    use time_manager,           only: get_nstep
    use check_energy,           only: check_energy_chng, check_energy_fix, & 
                                      check_water, check_qflx
    use mo_gas_phase_chemdr,    only: map2chm
    use clybry_fam,             only: clybry_fam_adj
    use output_aerocom_aie,     only: do_aerocom_ind3
    use phys_control,           only: use_qqflx_fixer, use_mass_borrower
    use crm_physics,            only: crm_physics_tend, crm_surface_flux_bypass_tend
    use cloud_diagnostics,      only: cloud_diagnostics_calc

    implicit none
    !-----------------------------------------------------------------------------
    ! Interface Arguments
    !-----------------------------------------------------------------------------
    real(r8),            intent(in   ) :: ztodt         ! 2 delta t (model time increment)
    real(r8),            intent(inout) :: fsns(pcols)   ! Surface solar absorbed flux
    real(r8),            intent(inout) :: fsnt(pcols)   ! Net column abs solar flux at model top
    real(r8),            intent(inout) :: flns(pcols)   ! Srf longwave cooling (up-down) flux
    real(r8),            intent(inout) :: flnt(pcols)   ! Net outgoing lw flux at model top
    real(r8),            intent(inout) :: fsds(pcols)   ! Surface solar down flux
    real(r8),            intent(in   ) :: sgh(pcols)    ! Std. deviation of orography
    real(r8),            intent(in   ) :: sgh30(pcols)  ! Std. deviation of 30 s orography for tms
    type(physics_state), intent(inout) :: state
    type(physics_tend ), intent(inout) :: tend
    type(physics_buffer_desc), pointer :: pbuf(:)
    type(cam_out_t),     intent(inout) :: cam_out
    type(cam_in_t),      intent(in)    :: cam_in
    !-----------------------------------------------------------------------------
    ! Local variables
    !-----------------------------------------------------------------------------
    type(physics_ptend)   :: ptend            ! indivdual parameterization tendencies
    type(physics_state)   :: state_alt        ! alt state for CRM input
    integer  :: nstep                         ! current timestep number
    real(r8) :: net_flx(pcols)
    real(r8) :: rtdt                          ! 1./ztodt
    integer  :: lchnk                         ! chunk identifier
    integer  :: ncol                          ! number of atmospheric columns
    integer  :: ierr
    integer  :: i,k,m,ihist                   ! Longitude, level, constituent indices
    integer  :: ixcldice, ixcldliq            ! constituent indices for cloud liquid and ice water.

    ! physics buffer indices
    integer itim_old, ifld

    ! physics buffer fields for total energy and mass adjustment
    real(r8), pointer, dimension(:  ) :: teout
    real(r8), pointer, dimension(:,:) :: tini
    real(r8), pointer, dimension(:,:) :: qini
    real(r8), pointer, dimension(:,:) :: cldliqini
    real(r8), pointer, dimension(:,:) :: cldiceini
    real(r8), pointer, dimension(:,:) :: dtcore
    real(r8), pointer, dimension(:,:,:) :: fracis  ! fraction of transported species that are insoluble

    ! energy checking variables
    real(r8) :: zero(pcols)                    ! array of zeros
    real(r8) :: flx_heat(pcols)

    logical   :: lq(pcnst)

    character(len=fieldname_len)   :: varname, vsuffix

    real(r8) :: ftem(pcols,pver)         ! tmp space
    real(r8), pointer, dimension(:) :: static_ener_ac_2d ! Vertically integrated static energy
    real(r8), pointer, dimension(:) :: water_vap_ac_2d   ! Vertically integrated water vapor
    real(r8) :: CIDiff(pcols)            ! Difference in vertically integrated static energy

    logical :: l_bc_energy_fix, l_dry_adj

    call phys_getopts( l_bc_energy_fix_out    = l_bc_energy_fix    &
                      ,l_dry_adj_out          = l_dry_adj          )

    !-----------------------------------------------------------------------------
    ! Initialize stuff
    !-----------------------------------------------------------------------------
    call t_startf('bc_init')

    zero = 0._r8
    lchnk = state%lchnk
    ncol  = state%ncol
    rtdt = 1._r8/ztodt
    nstep = get_nstep()

    if (pergro_test_active) then 
      !call outfld calls
      do ihist = 1 , nvars_prtrb_hist
        vsuffix  = trim(adjustl(hist_vars(ihist)))
        varname  = trim(adjustl(vsuffix))//'_topphysbc1' ! form variable name
        call outfld( trim(adjustl(varname)),get_var(state,vsuffix), pcols , lchnk )
      enddo
    endif

    call pbuf_get_field(pbuf, pbuf_get_index('static_ener_ac'), static_ener_ac_2d )
    call pbuf_get_field(pbuf, pbuf_get_index('water_vap_ac'), water_vap_ac_2d )

    ! Integrate and compute the difference
    ! CIDiff = difference of column integrated values
    if( nstep == 0 ) then
      CIDiff(:ncol) = 0.0_r8
      call outfld('DTENDTH', CIDiff, pcols, lchnk )
      call outfld('DTENDTQ', CIDiff, pcols, lchnk )
    else
      ! MSE first
      ftem(:ncol,:) = (state%s(:ncol,:) + latvap*state%q(:ncol,:,1)) * state%pdel(:ncol,:)*rga
      do k=2,pver
        ftem(:ncol,1) = ftem(:ncol,1) + ftem(:ncol,k)
      end do
      CIDiff(:ncol) = (ftem(:ncol,1) - static_ener_ac_2d(:ncol))*rtdt

      call outfld('DTENDTH', CIDiff, pcols, lchnk )
      ! Water vapor second
      ftem(:ncol,:) = state%q(:ncol,:,1)*state%pdel(:ncol,:)*rga
      do k=2,pver
        ftem(:ncol,1) = ftem(:ncol,1) + ftem(:ncol,k)
      end do
      CIDiff(:ncol) = (ftem(:ncol,1) - water_vap_ac_2d(:ncol))*rtdt

      call outfld('DTENDTQ', CIDiff, pcols, lchnk )
    end if

    ! Associate pointers with physics buffer fields
    itim_old = pbuf_old_tim_idx()
    call pbuf_get_field(pbuf, teout_idx, teout, (/1,itim_old/), (/pcols,1/))
    call pbuf_get_field(pbuf, tini_idx, tini)
    call pbuf_get_field(pbuf, qini_idx, qini)
    call pbuf_get_field(pbuf, cldliqini_idx, cldliqini)
    call pbuf_get_field(pbuf, cldiceini_idx, cldiceini)
    call pbuf_get_field(pbuf, pbuf_get_index('DTCORE'), dtcore, (/1,1,itim_old/), (/pcols,pver,1/) )
    call pbuf_get_field(pbuf, pbuf_get_index('FRACIS'), fracis, (/1,1,1/),        (/pcols, pver, pcnst/)  )
    fracis (:ncol,:,1:pcnst) = 1._r8

    ! Set physics tendencies to 0
    tend%dTdt(:ncol,:pver)  = 0._r8
    tend%dudt(:ncol,:pver)  = 0._r8
    tend%dvdt(:ncol,:pver)  = 0._r8

    !-----------------------------------------------------------------------------
    ! Mass checks and fixers
    !-----------------------------------------------------------------------------

    call check_qflx (state, tend, "PHYBC01", nstep, ztodt, cam_in%cflx(:,1))
    call check_water(state, tend, "PHYBC01", nstep, ztodt)

    ! make sure tracers are all positive - if use_mass_borrower then just print diagnostics
    call qneg3('TPHYSBCb', lchnk, ncol, pcols, pver, 1, pcnst, qmin, state%q, .not.use_mass_borrower )

    if(use_mass_borrower) then 
      ! tracer borrower for mass conservation 
      do m = 1, pcnst 
        call massborrow("PHYBC01",lchnk,ncol,state%psetcols,m,m,qmin(m),state%q(1,1,m),state%pdel)
      end do
    end if 

    ! Validate state coming from the dynamics.
    if (state_debug_checks) call physics_state_check(state, name="before tphysbc (dycore?)")

    ! Adjust chemistry for conservation issues
    call clybry_fam_adj( ncol, lchnk, map2chm, state%q, pbuf )

    ! Validate output of clybry_fam_adj
    if (state_debug_checks) call physics_state_check(state, name="clybry_fam_adj")

    ! make sure tracers are all positive, again - if use_mass_borrower then just print diagnostics
    call qneg3('TPHYSBCc',lchnk  ,ncol, pcols, pver, 1, pcnst, qmin, state%q, .not.use_mass_borrower )

    if(use_mass_borrower) then
      ! tracer borrower for mass conservation 
      do m = 1, pcnst
        call massborrow("PHYBC02",lchnk,ncol,state%psetcols,m,m,qmin(m),state%q(1,1,m),state%pdel)
      end do
    end if

    call check_water(state, tend, "PHYBC02", nstep, ztodt)

    ! Dump out "before physics" state
    call diag_state_b4_phys_write (state)

    call t_stopf('bc_init')

    !-----------------------------------------------------------------------------
    ! Global mean total energy fixer
    !-----------------------------------------------------------------------------
    if (l_bc_energy_fix) then

      call t_startf('energy_fixer')

      tini(:ncol,:pver) = state%t(:ncol,:pver)

      call check_energy_fix(state, ptend, nstep, flx_heat)
      call physics_update(state, ptend, ztodt, tend)
      call check_energy_chng(state, tend, "chkengyfix", nstep, ztodt, zero, zero, zero, flx_heat)

      ! Save state for convective tendency calculations.
      call diag_conv_tend_ini(state, pbuf)

      call cnst_get_ind('CLDLIQ', ixcldliq)
      call cnst_get_ind('CLDICE', ixcldice)

      qini     (:ncol,:pver) = state%q(:ncol,:pver,       1)
      cldliqini(:ncol,:pver) = state%q(:ncol,:pver,ixcldliq)
      cldiceini(:ncol,:pver) = state%q(:ncol,:pver,ixcldice)

      call outfld('TEOUT', teout       , pcols, lchnk   )
      call outfld('TEINP', state%te_ini, pcols, lchnk   )
      call outfld('TEFIX', state%te_cur, pcols, lchnk   )

      ! set and output the dse change due to dynpkg
      if( nstep > dyn_time_lvls-1 ) then
        do k = 1,pver
          dtcore(:ncol,k) = (tini(:ncol,k) - dtcore(:ncol,k))/(ztodt) + tend%dTdt(:ncol,k)
        end do
        call outfld( 'DTCORE', dtcore, pcols, lchnk )
      end if

      call t_stopf('energy_fixer')

    end if

    !-----------------------------------------------------------------------------
    ! Dry adjustment
    !-----------------------------------------------------------------------------
    if (l_dry_adj) then

      call t_startf('dry_adjustment')

      ! Copy state info for input to dadadj
      ! This is a kludge so dadadj doesn't have to be reformulated for DSE
      ! This code block is not a good example of interfacing a parameterization
      lq(:) = .FALSE.
      lq(1) = .TRUE.
      call physics_ptend_init(ptend, state%psetcols, 'dadadj', ls=.true., lq=lq)
      ptend%s(:ncol,:pver)   = state%t(:ncol,:pver)
      ptend%q(:ncol,:pver,1) = state%q(:ncol,:pver,1)

      call dadadj (lchnk, ncol, state%pmid, state%pint, state%pdel, ptend%s, ptend%q(1,1,1))

      ptend%s(:ncol,:)   = (ptend%s(:ncol,:)   - state%t(:ncol,:)  )/ztodt * cpair
      ptend%q(:ncol,:,1) = (ptend%q(:ncol,:,1) - state%q(:ncol,:,1))/ztodt

      call physics_update(state, ptend, ztodt, tend)
      call t_stopf('dry_adjustment')

    end if

    !-----------------------------------------------------------------------------
    ! MMF surface flux bypass
    !-----------------------------------------------------------------------------
#if defined( MMF_FLUX_BYPASS )

  ! Check if LHF exceeds the total moisture content of the lowest layer
  call qneg4('TPHYSBC ', lchnk, ncol, ztodt, &
              state%q(1,pver,1), state%rpdel(1,pver), &
              cam_in%shf, cam_in%lhf, cam_in%cflx )

  call crm_surface_flux_bypass_tend(state, cam_in, ptend)
  call physics_update(state, ptend, ztodt, tend)  
  call check_energy_chng(state, tend, "crm_tend", nstep, ztodt,  &
                         cam_in%shf(:), zero, zero, cam_in%cflx(:,1)) 
#endif

  end subroutine tphysbc1

!===================================================================================================
!===================================================================================================

  subroutine tphysbc2(ztodt, fsns, fsnt, flns, flnt, &
                      state, tend, pbuf, fsds, &
                      sgh, sgh30, cam_out, cam_in, crm_ecpp_output )
    !----------------------------------------------------------------------------- 
    ! Purpose: Evaluate physics processes BEFORE coupling to sfc components
    !          Phase 2 - aerosols, radiation, and diagnostics
    !
    ! Pass surface fields for separate surface flux calculations
    ! Dump appropriate fields to history file.
    !-----------------------------------------------------------------------------
    use physics_buffer,         only: physics_buffer_desc, pbuf_get_field
    use physics_buffer,         only: pbuf_get_index, pbuf_old_tim_idx
    use physics_buffer,         only: dyn_time_lvls
    use physics_types,          only: physics_ptend_init, physics_ptend_sum, &
                        physics_state_check, physics_ptend_scale
    use cam_diagnostics,        only: diag_conv_tend_ini, diag_phys_writeout, diag_conv, &
                        diag_export, diag_state_b4_phys_write
    use cam_history,            only: outfld, fieldname_len
    use physconst,              only: cpair, latvap, gravit, rga
    use constituents,           only: pcnst, qmin, cnst_get_ind
    use time_manager,           only: get_nstep
    use check_energy,           only: check_energy_chng, & 
                        check_tracers_data, check_tracers_init, &
                        check_tracers_chng, check_tracers_fini
    use aero_model,             only: aero_model_wetdep
    use modal_aero_calcsize,    only: modal_aero_calcsize_sub
    use modal_aero_wateruptake, only: modal_aero_wateruptake_dr
    use radiation,              only: radiation_tend
    use tropopause,             only: tropopause_output
    use output_aerocom_aie,     only: do_aerocom_ind3, cloud_top_aerocom
    use cloud_diagnostics,      only: cloud_diagnostics_calc
    use crm_ecpp_output_module, only: crm_ecpp_output_type
    use camsrfexch,             only: cam_export
#if defined( ECPP )
    use module_ecpp_ppdriver2, only: parampollu_driver2
    use module_data_ecpp1,     only: dtstep_pp_input
    use crmclouds_camaerosols, only: crmclouds_mixnuc_tend
#endif
    implicit none
    !-----------------------------------------------------------------------------
    ! Interface Arguments
    !-----------------------------------------------------------------------------
    real(r8),                  intent(in   ) :: ztodt         ! 2 delta t (model time increment)
    real(r8),                  intent(inout) :: fsns(pcols)   ! Surface solar absorbed flux
    real(r8),                  intent(inout) :: fsnt(pcols)   ! Net column abs solar flux at model top
    real(r8),                  intent(inout) :: flns(pcols)   ! Srf longwave cooling (up-down) flux
    real(r8),                  intent(inout) :: flnt(pcols)   ! Net outgoing lw flux at model top
    real(r8),                  intent(inout) :: fsds(pcols)   ! Surface solar down flux
    real(r8),                  intent(in   ) :: sgh(pcols)    ! Std. deviation of orography
    real(r8),                  intent(in   ) :: sgh30(pcols)  ! Std. deviation of 30 s orography for tms
    type(physics_state),       intent(inout) :: state
    type(physics_tend ),       intent(inout) :: tend
    type(physics_buffer_desc), pointer       :: pbuf(:)
    type(cam_out_t),           intent(inout) :: cam_out
    type(cam_in_t),            intent(in   ) :: cam_in
    type(crm_ecpp_output_type),intent(inout) :: crm_ecpp_output   ! CRM output data for ECPP calculations
    !-----------------------------------------------------------------------------
    ! Local variables
    !-----------------------------------------------------------------------------
    type(physics_ptend)   :: ptend            ! indivdual parameterization tendencies
    type(physics_state)   :: state_alt        ! alt state for CRM input
    integer  :: nstep                         ! current timestep number
    real(r8) :: net_flx(pcols)
    real(r8) :: cmfmc2(pcols,pverp)           ! Moist convection cloud mass flux
    real(r8) :: dlf(pcols,pver)               ! Detraining cld H20 from shallow + deep convections
    real(r8) :: dlf2(pcols,pver)              ! Detraining cld H20 from shallow convections
    integer  :: lchnk                         ! chunk identifier
    integer  :: ncol                          ! number of atmospheric columns
    integer  :: ierr
    integer  :: i,k,m,ihist                   ! Longitude, level, constituent indices

    real(r8) :: sh_e_ed_ratio(pcols,pver)      ! shallow conv [ent/(ent+det)] ratio 

    ! pbuf fields
    integer itim_old
    real(r8), pointer, dimension(:,:) :: cld  ! cloud fraction
    real(r8), pointer, dimension(:,:) :: cldo ! old cloud fraction

    ! energy checking variables
    real(r8) :: zero(pcols)                    ! array of zeros
    type(check_tracers_data):: tracerint       ! energy integrals and cummulative boundary fluxes
    real(r8) :: zero_tracers(pcols,pcnst)

    logical   :: lq(pcnst)

    character(len=fieldname_len)   :: varname, vsuffix

    !BSINGH - these were moved from zm_conv_intr because they are  used by aero_model_wetdep 
    real(r8), dimension(pcols,pver) :: mu, eu, du, md, ed, dp
    real(r8):: dsubcld(pcols) ! wg layer thickness in mbs (between upper/lower interface).
    integer :: jt(pcols)      ! wg layer thickness in mbs between lcl and maxi.    
    integer :: maxg(pcols)    ! wg top  level index of deep cumulus convection.
    integer :: ideep(pcols)   ! wg gathered values of maxi.
    integer :: lengath        ! w holds position of gathered points vs longitude index

    logical :: l_tracer_aero, l_rad

    logical           :: use_ECPP
    real(r8), pointer :: mmf_clear_rh(:,:) ! CRM clear air relative humidity used for aerosol water uptake

    ! ECPP variables
    real(r8),pointer,dimension(:)   :: pblh              ! PBL height (for ECPP)
    real(r8),pointer,dimension(:,:) :: acldy_cen_tbeg    ! cloud fraction
    real(r8)                        :: dtstep_pp         ! ECPP time step (seconds)
    integer                         :: necpp             ! number of GCM time steps in which ECPP is called once

    call phys_getopts( l_tracer_aero_out      = l_tracer_aero      &
                      ,l_rad_out              = l_rad              &
                      ,use_ECPP_out           = use_ECPP           )

    !-----------------------------------------------------------------------------
    ! Initialize stuff
    !-----------------------------------------------------------------------------
    call t_startf('bc_init')

    zero = 0._r8
    zero_tracers(:,:) = 0._r8
    lchnk = state%lchnk
    ncol  = state%ncol
    nstep = get_nstep()

    if (pergro_test_active) then 
      !call outfld calls
      do ihist = 1 , nvars_prtrb_hist
        vsuffix  = trim(adjustl(hist_vars(ihist)))
        varname  = trim(adjustl(vsuffix))//'_topphysbc2' ! form variable name
        call outfld( trim(adjustl(varname)),get_var(state,vsuffix), pcols , lchnk )
      enddo
    endif

    call pbuf_get_field(pbuf, mmf_clear_rh_idx, mmf_clear_rh )

    ! compute mass integrals of input tracers state
    call check_tracers_init(state, tracerint)

    !-----------------------------------------------------------------------------
    ! Modal aerosol wet radius for radiative calculation
    !-----------------------------------------------------------------------------
#if defined( ECPP ) && defined(MODAL_AERO)
    ! temporarily turn on all lq, so it is allocated
    lq(:) = .true.
    call physics_ptend_init(ptend, state%psetcols, 'crm - modal_aero_wateruptake_dr', lq=lq)

    call t_startf('modal_aero_mmf')

    ! set all ptend%lq to false as they will be set in modal_aero_calcsize_sub
    ptend%lq(:) = .false.
    call modal_aero_calcsize_sub (state, ztodt, pbuf, ptend)
    call modal_aero_wateruptake_dr(state, pbuf, clear_rh_in=mmf_clear_rh)

    ! ECPP handles aerosol wet deposition, so tendency from wet depostion is 
    ! not updated in mz_aero_wet_intr (mz_aerosols_intr.F90), but tendencies
    ! from other parts of crmclouds_aerosol_wet_intr() are still updated here.
    call physics_update (state, ptend, ztodt, tend)
    call t_stopf('modal_aero_mmf')

    call check_energy_chng(state, tend, "modal_aero_mmf", nstep, ztodt, &
             zero, zero, zero, zero)

#endif /* ECPP and MODAL_AERO */
!-----------------------------------------------------------------------------
! ECPP - Explicit-Cloud Parameterized-Pollutant
!-----------------------------------------------------------------------------
#if defined( ECPP )
    if (use_ECPP) then

      call pbuf_get_field(pbuf, pbuf_get_index('pblh'), pblh)
      call pbuf_get_field(pbuf, pbuf_get_index('ACLDY_CEN'), acldy_cen_tbeg)

      dtstep_pp = dtstep_pp_input
      necpp = dtstep_pp/ztodt

      if (nstep.ne.0 .and. mod(nstep, necpp).eq.0) then

        ! aerosol tendency from droplet activation and mixing
        ! cldo and cldn are set to be the same in crmclouds_mixnuc_tend,
        ! So only turbulence mixing is done here.
        call t_startf('crmclouds_mixnuc')
        call crmclouds_mixnuc_tend(state, ptend, dtstep_pp,           &
                         cam_in%cflx, pblh, pbuf,           &
                         crm_ecpp_output%wwqui_cen,         &
                         crm_ecpp_output%wwqui_cloudy_cen,  &
                         crm_ecpp_output%wwqui_bnd,         &
                         crm_ecpp_output%wwqui_cloudy_bnd,  &
                         species_class)
        call physics_update(state, ptend, dtstep_pp, tend)
        call t_stopf('crmclouds_mixnuc')

        ! ECPP interface
        call t_startf('ecpp')
        call parampollu_driver2(state, ptend, pbuf, dtstep_pp, dtstep_pp,   &
                      crm_ecpp_output%acen,       crm_ecpp_output%abnd,         &
                      crm_ecpp_output%acen_tf,    crm_ecpp_output%abnd_tf,      &
                      crm_ecpp_output%massflxbnd, crm_ecpp_output%rhcen,        &
                      crm_ecpp_output%qcloudcen,  crm_ecpp_output%qlsinkcen,    &
                      crm_ecpp_output%precrcen,   crm_ecpp_output%precsolidcen, &
                      acldy_cen_tbeg)
        call physics_update(state, ptend, dtstep_pp, tend)
        call t_stopf ('ecpp')

      end if ! nstep.ne.0 .and. mod(nstep, necpp).eq.0

    end if ! use_ECPP
#endif /* ECPP */
    !-----------------------------------------------------------------------------
    ! save old CRM cloud fraction - w/o CRM, this is done in cldwat2m.F90
    !-----------------------------------------------------------------------------
    itim_old = pbuf_old_tim_idx()
    call pbuf_get_field(pbuf,pbuf_get_index('CLDO'),cldo,start=(/1,1,itim_old/),kount=(/pcols,pver,1/))
    call pbuf_get_field(pbuf,pbuf_get_index('CLD') ,cld ,start=(/1,1,itim_old/),kount=(/pcols,pver,1/))

    cldo(1:ncol,1:pver) = cld(1:ncol,1:pver)

    !-----------------------------------------------------------------------------
    ! Aerosol stuff
    !-----------------------------------------------------------------------------
    if (l_tracer_aero) then
      if (use_ECPP) then
        ! With MMF + ECPP we can skip the conventional aerosol routines
      else
        ! Aerosol wet chemistry determines scavenging and transformations.
        ! This is followed by convective transport of all trace species except
        ! water vapor and condensate. Scavenging needs to occur prior to
        ! transport in order to determine interstitial fraction.

        ! Without ECPP we should be using prescribed aerosols, so we only
        ! need to consider the wet deposition and water uptake for radiation

        ! Aerosol wet removal (including aerosol water uptake)
        call t_startf('aero_model_wetdep')
        call aero_model_wetdep( ztodt, dlf, dlf2, cmfmc2, state,  & ! inputs
                                sh_e_ed_ratio, mu, md, du, eu, ed, dp, dsubcld,    &
                                jt, maxg, ideep, lengath, species_class,           &
                                cam_out, pbuf, ptend,                              & ! outputs
                                clear_rh=mmf_clear_rh) ! clear air relative humidity for water uptake
        call physics_update(state, ptend, ztodt, tend)
        call t_stopf('aero_model_wetdep')

        ! check tracer integrals
        call check_tracers_chng(state, tracerint, "aero_model_wetdep", nstep, ztodt,  zero_tracers)

      end if
    end if ! l_tracer_aero

#ifndef MMF_NN_EMULATOR
    ! in mmf_nn_emulator_driver, this block is included, so need to skip when MMF_NN_EMULATOR is defined
    ! even not executing this block, it won't influence the radiation_tend, which won't use cam_out%psl
    !-----------------------------------------------------------------------------
    ! Moist physical parameteriztions complete: 
    ! send dynamical variables, and derived variables to history file
    !-----------------------------------------------------------------------------
    call t_startf('bc_history_write')
    call diag_phys_writeout(state, cam_out%psl)
    call diag_conv(state, ztodt, pbuf)
    call t_stopf('bc_history_write')

    !-----------------------------------------------------------------------------
    ! Write cloud diagnostics on history file
    !-----------------------------------------------------------------------------
    call t_startf('bc_cld_diag_history_write')
    call cloud_diagnostics_calc(state, pbuf)
    call t_stopf('bc_cld_diag_history_write')
#endif

    !-----------------------------------------------------------------------------
    ! Radiation computations
    !-----------------------------------------------------------------------------
    if (l_rad) then
      call t_startf('radiation')
      call radiation_tend(state,ptend, pbuf, cam_out, cam_in, &
              cam_in%landfrac, cam_in%icefrac, cam_in%snowhland, &
              fsns, fsnt, flns, flnt, fsds, &
              net_flx, is_cmip6_volc, ztodt, clear_rh=mmf_clear_rh)
      call t_stopf('radiation')

      ! We don't need to call physics_update or check_energy_chng for the MMF
      ! because the radiative tendency is added within the call to crm_physics_tend
      ! seems that it will update the crm_rad etc in pbuf, this will update the state in the next call to crm_physics_tend and crm_module

    end if ! l_rad

    !-----------------------------------------------------------------------------
    ! Diagnostics
    !-----------------------------------------------------------------------------

    call t_startf('tphysbc_diagnostics')

    if(do_aerocom_ind3) call cloud_top_aerocom(state, pbuf) 

#ifndef MMF_NN_EMULATOR
    ! in mmf_nn_emulator_driver, this block is included, so need to skip when MMF_NN_EMULATOR is defined
    ! Diagnose the location of the tropopause
    call tropopause_output(state)

    ! Save atmospheric fields to force surface models
    call cam_export(state,cam_out,pbuf)

    ! Write export state to history file
    call diag_export(cam_out)
#endif

    call check_tracers_fini(tracerint)

    call t_stopf('tphysbc_diagnostics')

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
  end subroutine tphysbc2

subroutine phys_timestep_init(phys_state, cam_in, cam_out, pbuf2d)
!-----------------------------------------------------------------------------------
!
! Purpose: The place for parameterizations to call per timestep initializations.
!          Generally this is used to update time interpolated fields from boundary
!          datasets.
!
!-----------------------------------------------------------------------------------
  use shr_kind_mod,        only: r8 => shr_kind_r8
  use chemistry,           only: chem_timestep_init
  use chem_surfvals,       only: chem_surfvals_set
  use physics_types,       only: physics_state
  use physics_buffer,      only: physics_buffer_desc
  use carma_intr,          only: carma_timestep_init
  use ghg_data,            only: ghg_data_timestep_init
  use cam3_aero_data,      only: cam3_aero_data_on, cam3_aero_data_timestep_init
  use cam3_ozone_data,     only: cam3_ozone_data_on, cam3_ozone_data_timestep_init
  use aoa_tracers,         only: aoa_tracers_timestep_init
  use vertical_diffusion,  only: vertical_diffusion_ts_init
  use radheat,             only: radheat_timestep_init
  use solar_data,          only: solar_data_advance
  use qbo,                 only: qbo_timestep_init
  use iondrag,             only: do_waccm_ions, iondrag_timestep_init
  use perf_mod

  use prescribed_ozone,    only: prescribed_ozone_adv
  use prescribed_ghg,      only: prescribed_ghg_adv
  use prescribed_aero,     only: prescribed_aero_adv
  use aerodep_flx,         only: aerodep_flx_adv
  use aircraft_emit,       only: aircraft_emit_adv
  use prescribed_volcaero, only: prescribed_volcaero_adv
  use prescribed_strataero,only: prescribed_strataero_adv
  use mo_apex,             only: mo_apex_init
  use epp_ionization,      only: epp_ionization_active
  use iop_forcing,         only: scam_use_iop_srf
  use nudging,             only: Nudge_Model, nudging_timestep_init

  implicit none

  type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
  type(cam_in_t),      intent(inout), dimension(begchunk:endchunk) :: cam_in
  type(cam_out_t),     intent(inout), dimension(begchunk:endchunk) :: cam_out

  type(physics_buffer_desc), pointer                 :: pbuf2d(:,:)

  !-----------------------------------------------------------------------------

  if (single_column) call scam_use_iop_srf(cam_in)

  ! update geomagnetic coordinates
  if (epp_ionization_active .or. do_waccm_ions) then
     call mo_apex_init(phys_state)
  endif

  ! Chemistry surface values
  call chem_surfvals_set()

  ! Solar irradiance
  call solar_data_advance()

  ! Time interpolate for chemistry.
  call chem_timestep_init(phys_state, pbuf2d)

  ! Prescribed tracers
  call prescribed_ozone_adv(phys_state, pbuf2d)
  call prescribed_ghg_adv(phys_state, pbuf2d)
  call prescribed_aero_adv(phys_state, pbuf2d)
  call aircraft_emit_adv(phys_state, pbuf2d)
  call prescribed_volcaero_adv(phys_state, pbuf2d)
  call prescribed_strataero_adv(phys_state, pbuf2d)

  ! prescribed aerosol deposition fluxes
  call aerodep_flx_adv(phys_state, pbuf2d, cam_out)

  ! CAM3 prescribed aerosol masses
  if (cam3_aero_data_on) call cam3_aero_data_timestep_init(pbuf2d,  phys_state)

  ! CAM3 prescribed ozone data
  if (cam3_ozone_data_on) call cam3_ozone_data_timestep_init(pbuf2d,  phys_state)

  ! Time interpolate data models of gasses in pbuf2d
  call ghg_data_timestep_init(pbuf2d,  phys_state)

  ! Upper atmosphere radiative processes
  call radheat_timestep_init(phys_state, pbuf2d)

  ! Time interpolate for vertical diffusion upper boundary condition
  call vertical_diffusion_ts_init(pbuf2d, phys_state)

  !----------------------------------------------------------------------
  ! update QBO data for this time step
  !----------------------------------------------------------------------
  call qbo_timestep_init

  call iondrag_timestep_init()

  call carma_timestep_init()

  ! age of air tracers
  call aoa_tracers_timestep_init(phys_state)

  ! Update Nudging values, if needed
  !----------------------------------
  if(Nudge_Model) call nudging_timestep_init(phys_state)

end subroutine phys_timestep_init

end module physpkg
