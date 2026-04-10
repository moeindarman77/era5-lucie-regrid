program main

  use stringtype, only : string

  implicit none 

  integer :: year_i, month_i, start_year,end_year,start_month,end_month

  ! Base directory for inputs/, outputs/, aux/ subdirectories.
  ! Change this to your repo clone location.
  character(len=256) :: base_dir = '/glade/derecho/scratch/mdarman/sst_data/era5_data_download/gen_lucie_3d'

  character(len=3) :: file_end='.nc'
  character(len=7) :: file_begin = 'era_5_m'
  character(len=10) :: p_file = '_surface_p'
  character(len=8) :: temp2m_file = '_2m_temp'
  character(len=9) :: vort_div_file = '_div_vort'
  character(len=5) :: full_file = '_full'
  character(len=10) :: regrid_file = '_regridded'
  character(len=22) :: regrid_file_test = '_regridded_speedyptest'
  character(len=2) :: q_file ='_q'
  character(len=2) :: mid_file='_y'
  character(len=1) :: month_1
  character(len=2) :: month_2
  character(len=4) :: year
  character(len=:), allocatable :: file_path 
  character(len=:), allocatable :: file_path_tiered
  character(len=:), allocatable :: full_file_name
  character(len=:), allocatable :: p_file_name 
  character(len=:), allocatable :: regrid_file_name
  character(len=:), allocatable :: specific_h_file_name
  character(len=:), allocatable :: temp2m_file_name
  character(len=:), allocatable :: vort_div_file_name
  character(len=:), allocatable :: toa_tisr_name
  character(len=:), allocatable :: mpi_final_file_name
  character(len=:), allocatable :: format_month
  character(len=:), allocatable :: month

  type(string), allocatable :: regrid_files(:)

  allocate(regrid_files(12))

  start_year = 2022
  end_year = 2022
  start_month = 7
  end_month = 7
  
  do year_i=start_year,end_year
     do month_i=start_month,end_month
        if(month_i >= 10) then
          format_month = '(I2)'
          
          write(year,'(I4)') year_i
          write(month_2,'(I2)') month_i
          
          month = month_2

          file_path = trim(base_dir)//'/inputs/'//year//'/'

          full_file_name = file_path//file_begin//month//mid_file//year//full_file//file_end
          p_file_name = file_path//file_begin//month//mid_file//year//p_file//file_end
          regrid_file_name = trim(base_dir)//'/outputs/'//year//'/'//file_begin//month//mid_file//year//regrid_file//'_spectral_mpi'//file_end
          specific_h_file_name = file_path//file_begin//month//mid_file//year//q_file//file_end
          temp2m_file_name = file_path//file_begin//month//mid_file//year//temp2m_file//file_end

          regrid_files(month_i)%str = regrid_file_name
          print *, month_i, regrid_files(month_i)%str
        else
          format_month = '(I1)'

          write(year,'(I4)') year_i
          write(month_1,'(I1)') month_i

          month = month_1

          file_path = trim(base_dir)//'/inputs/'//year//'/'

          full_file_name = file_path//file_begin//month//mid_file//year//full_file//file_end
          p_file_name = file_path//file_begin//month//mid_file//year//p_file//file_end
          regrid_file_name = trim(base_dir)//'/outputs/'//year//'/'//file_begin//month//mid_file//year//regrid_file//'_spectral_mpi'//file_end
          specific_h_file_name = file_path//file_begin//month//mid_file//year//q_file//file_end
          temp2m_file_name = file_path//file_begin//month//mid_file//year//temp2m_file//file_end

          regrid_files(month_i)%str = regrid_file_name
          print *, month_i, regrid_files(month_i)%str
        endif
       
        call month_year_regrid(full_file_name,p_file_name,specific_h_file_name,temp2m_file_name,regrid_file_name,base_dir)
     enddo
     print *, shape(regrid_files)
     !call year_regrid_tisr(toa_tisr_name,regrid_files)
  enddo
end program 

subroutine month_year_regrid(full_file_name,surface_p_file,q_file_name,temp2m_file_name,regrid_file,base_dir)
  use mod_interp2d, only : bilinear_interp, bspline_interface, era_p_level_to_speedy_sigma_level, speedy_nlvls, era_p_level_to_speedy_p_level
  use mod_utilities, only : dp, sp, gridx, gridy, speedy_lat, flip, hydrostatic_factor, g0, e_constant
  use mod_io, only : read_netcdf_4d_sp,read_netcdf_4d_dp_era, read_netcdf_1d_dp, write_era_5_regridded, read_netcdf_3d_dp_era,read_netcdf_3d_dp, read_netcdf_2d_dp,nc_write_2d_dp,add_var_to_file_4d
  use stringtype, only : string


  implicit none

  character(len=*), intent(in) :: full_file_name,surface_p_file,regrid_file,q_file_name,temp2m_file_name
  character(len=*), intent(in) :: base_dir

  real(kind=dp), allocatable :: var(:,:,:,:), var3d(:,:,:)
  real(kind=dp), allocatable :: era_lats(:),era_lons(:),speedy_lons(:),flipped_era_lats(:),levels(:),times_var(:)
  real(kind=dp), allocatable :: regrid_var(:,:,:,:,:), regrid_surfacep(:,:,:)
  real(kind=dp), allocatable :: target_vertical_indices(:,:,:), averaged_speedy_surfacep(:,:)
  real(kind=dp), allocatable :: era_orography(:,:), regridded_era_orography(:,:), speedy_orography(:,:)
  real(kind=dp), allocatable :: temp3d(:,:,:)
  real(kind=dp), allocatable :: era_regrid_2mtemp(:,:,:)
  real(kind=dp), allocatable :: hydrostatic_factor_var(:,:)

  integer :: height_levels, counter, t, z, times, numvar, num_era_vars, era_x, era_y

  type(string), allocatable :: era_vars(:)

  !Get year average speedy log(surfacep)
  call read_netcdf_3d_dp('average_logp',trim(base_dir)//'/aux/year_average_logp.nc',temp3d)
  averaged_speedy_surfacep = temp3d(:,:,1)
  averaged_speedy_surfacep = e_constant**averaged_speedy_surfacep
  averaged_speedy_surfacep = averaged_speedy_surfacep*1000.0_dp
  deallocate(temp3d)
  !print *, shape(averaged_speedy_logp)

  !Get ERA orography data in meters
  call read_netcdf_3d_dp_era('z',trim(base_dir)//'/aux/era_orography_jan01_1980.nc',temp3d)
  era_orography = temp3d(:,:,1)
  deallocate(temp3d)
  era_orography = era_orography/g0

  call read_netcdf_2d_dp('orography',trim(base_dir)//'/aux/speedy_orography.nc',speedy_orography)

  num_era_vars = 4
  allocate(era_vars(num_era_vars))

  era_vars(1)%str = 't'
  era_vars(2)%str = 'u'
  era_vars(3)%str = 'v'
  era_vars(4)%str = 'q'

  !Get dimensions of vars for allocation 
  print *, 'full_file_name',full_file_name,'regrid',regrid_file
  call read_netcdf_1d_dp('latitude',full_file_name,era_lats)
  call read_netcdf_1d_dp('longitude',full_file_name,era_lons)
  call read_netcdf_1d_dp('level',full_file_name,levels)
  call read_netcdf_1d_dp('time',full_file_name,times_var)

  !Era lats are sorted from biggest to smallest
  !needs to be flipped for everything to work
  allocate(flipped_era_lats,mold=era_lats)

  call flip(era_lats,flipped_era_lats)

  era_lats = flipped_era_lats
  deallocate(flipped_era_lats)

  !Get speedy lons
  allocate(speedy_lons(gridx))
  speedy_lons = (/(real(counter)*3.75,counter=0,95)/)

  !Get size of each dimension
  era_x = size(era_lons,1) 
  era_y = size(era_lats,1)
  height_levels = size(levels,1)
  times = size(times_var,1)

  allocate(regrid_var(num_era_vars,gridx,gridy,speedy_nlvls,times))
  allocate(regrid_surfacep(gridx,gridy,times))

  !TODO NEED to impose a hydrostatic balance 
  !should take place after horizontal regridded 
  !then do the vertical interploation 
  !allocate(regridded_era_orography(gridx,gridy))
  !call bspline_interface(era_lons,era_lats,era_orography,speedy_lons,speedy_lat,regridded_era_orography(:,gridy:1:-1))

  !call read_netcdf_3d_dp_era('tisr',temp2m_file_name,var3d)
  !allocate(era_regrid_2mtemp(gridx,gridy,times))
  !do t=1,times
  !   call bspline_interface(era_lons,era_lats,var3d(:,:,t),speedy_lons,speedy_lat,era_regrid_2mtemp(:,gridy:1:-1,t))
  !enddo
  !stop 
  !deallocate(var3d)
  
  !allocate(hydrostatic_factor_var(gridx,gridy)) 
  call read_netcdf_3d_dp_era('sp',surface_p_file,var3d)
  print *, 'doing surface p interp'
  do t=1,times
     call bspline_interface(era_lons,era_lats,var3d(:,:,t),speedy_lons,speedy_lat,regrid_surfacep(:,gridy:1:-1,t))
     !call hydrostatic_factor(regridded_era_orography,speedy_orography,era_regrid_2mtemp(:,:,t),hydrostatic_factor_var)
     !regrid_surfacep(:,:,t) = regrid_surfacep(:,:,t)*hydrostatic_factor_var
  enddo

  do numvar=1,num_era_vars
     print *, era_vars(numvar)%str
     !Changed file structure so for 1991+ the if statement needs to be commented
     !out 
     !if(numvar == 4) then
     !  call read_netcdf_4d_dp_era(era_vars(numvar)%str,q_file_name,var)
     !else
     call read_netcdf_4d_dp_era(era_vars(numvar)%str,full_file_name,var)
     !endif 
     print *, var(16,18,:,1)
     do t=1,times
        print *, 't',t
        allocate(temp3d(gridx,gridy,height_levels))

        do z=1,height_levels
           call bspline_interface(era_lons,era_lats,var(:,:,z,t),speedy_lons,speedy_lat,temp3d(:,gridy:1:-1,z))
           !"gridy:1:-1" is a trick to flip the matrix over the y axis see np.flipud
        enddo
        ! Troy's paper: "compute the value of σ at each horizontal grid point" — this
        ! is per-gridpoint TRUE sigma (P_target = sigma * actual local SP).
        call era_p_level_to_speedy_p_level(averaged_speedy_surfacep,temp3d,levels,regrid_var(numvar,:,:,:,t))
        !call era_p_level_to_speedy_sigma_level(regrid_surfacep(:,:,t),temp3d,levels,regrid_var(numvar,:,:,:,t))
        deallocate(temp3d)
     enddo
     deallocate(var)
  enddo

  regrid_surfacep = log(regrid_surfacep/100000.0_dp)
  call write_era_5_regridded(regrid_file,regrid_var,regrid_surfacep)
end subroutine

subroutine year_regrid_tisr(full_file_name,regrid_files_test)
  use mod_interp2d, only : bilinear_interp, bspline_interface, era_p_level_to_speedy_sigma_level, speedy_nlvls, era_p_level_to_speedy_p_level
  use mod_utilities, only : dp, sp, gridx, gridy, speedy_lat, flip, hydrostatic_factor, g0, e_constant
  use mod_io, only : read_netcdf_4d_sp,read_netcdf_4d_dp_era, read_netcdf_1d_dp, write_era_5_regridded, read_netcdf_3d_dp_era,read_netcdf_3d_dp, read_netcdf_2d_dp,nc_write_2d_dp,add_var_to_file_4d,add_var_to_file_3d,nc_write_3d_dp
  use stringtype, only : string

  implicit none 

  character(len=*), intent(in) :: full_file_name

  type(string), intent(in), dimension(12) :: regrid_files_test

  real(kind=dp), allocatable :: var(:,:,:,:), var3d(:,:,:)
  real(kind=dp), allocatable :: era_lats(:),era_lons(:),speedy_lons(:),flipped_era_lats(:),levels(:),times_var(:)
  real(kind=dp), allocatable :: regrid_var(:,:,:,:,:), regrid_surfacep(:,:,:)
  real(kind=dp), allocatable :: target_vertical_indices(:,:,:), averaged_speedy_surfacep(:,:)
  real(kind=dp), allocatable :: era_orography(:,:), regridded_era_orography(:,:), speedy_orography(:,:)
  real(kind=dp), allocatable :: temp3d(:,:,:)
  real(kind=dp), allocatable :: era_regrid_2mtemp(:,:,:)
  real(kind=dp), allocatable :: hydrostatic_factor_var(:,:)

  integer :: height_levels, counter, t, z, times, numvar, num_era_vars, era_x, era_y
  integer :: hour_counter 

  type(string), allocatable :: era_vars(:)

  num_era_vars = 1
  allocate(era_vars(num_era_vars))
  era_vars(1)%str = 'tisr'

  print *, shape(regrid_files_test)
  !Get dimensions of vars for allocation 
  print *, 'full_file_name',full_file_name
  print *, 'regrid',regrid_files_test(1)%str

  call read_netcdf_1d_dp('lat',full_file_name,era_lats)
  call read_netcdf_1d_dp('lon',full_file_name,era_lons)
  call read_netcdf_1d_dp('time',full_file_name,times_var)

  !Era lats are sorted from biggest to smallest
  !needs to be flipped for everything to work
  allocate(flipped_era_lats,mold=era_lats)

  !call flip(era_lats,flipped_era_lats)

  era_lats = era_lats!flipped_era_lats
  deallocate(flipped_era_lats)

  !Get speedy lons
  allocate(speedy_lons(gridx))
  speedy_lons = (/(real(counter)*3.75,counter=0,95)/)

  !Get size of each dimension
  era_x = size(era_lons,1) 
  era_y = size(era_lats,1)
  times = size(times_var,1)

  allocate(regrid_var(num_era_vars,gridx,gridy,speedy_nlvls,times))
  allocate(regrid_surfacep(gridx,gridy,times))

  !TODO NEED to impose a hydrostatic balance 
  !should take place after horizontal regridded 
  !then do the vertical interploation 
  !allocate(regridded_era_orography(gridx,gridy))
  !call bspline_interface(era_lons,era_lats,era_orography,speedy_lons,speedy_lat,regridded_era_orography(:,gridy:1:-1))

  call read_netcdf_3d_dp_era('tisr',full_file_name,var3d)
  print *, 'here;'
  allocate(era_regrid_2mtemp(gridx,gridy,times))
  do t=1,times
     call bspline_interface(era_lons,era_lats,var3d(:,:,t),speedy_lons,speedy_lat,regrid_surfacep(:,:,t))
  enddo
 
  call nc_write_3d_dp('test.nc','tisr','J m**-2',regrid_surfacep) 
  deallocate(var3d)

  print *, 'length of regrid_surfacep',size(regrid_surfacep,3)
   
  hour_counter = 0 
  do t=1, 12
     print *, 'i,',t,regrid_files_test(t)%str 
     call read_netcdf_3d_dp_era('logp',regrid_files_test(t)%str,var3d)
     times = size(var3d,3)
     print *, times
     call add_var_to_file_3d(regrid_files_test(t)%str,era_vars(1)%str,'J m**-2',regrid_surfacep(:,:,hour_counter+1:times+hour_counter))
     hour_counter = times + hour_counter
  enddo 

  if(hour_counter /= size(regrid_surfacep,3)) then
    print *, 'full_file_name',full_file_name 
    print *, 'the time dimensions are the same lengthe something is wrong'
    stop
  endif  
end subroutine 


  
