module mod_utilities
  implicit none
  
  integer, parameter :: dp=selected_real_kind(14)  
  integer, parameter :: sp = selected_real_kind(6, 37)

  integer, parameter :: gridy=48, gridx=96
  real(kind=dp), parameter :: speedy_lat(gridy) = (/-87.159, -83.479, -79.777,  -76.070, -72.362, -68.652, &
                 -64.942, -61.232, -57.521, -53.810, -50.099, -46.389, -42.678,-38.967, &
                 -35.256, -31.545 , -27.833, -24.122, -20.411, -16.700, -12.989,-9.278, &
                 -5.567, -1.856, 1.856, 5.567, 9.278, 12.989, 16.700, 20.411, &
                  24.122, 27.833, 31.545, 35.256, 38.967, 42.678, 46.389, 50.099, &
                  53.810, 57.521, 61.232, 64.942, 68.652, 72.362, 76.070, 79.777, &
                  83.479, 87.159/)
  real(kind=dp), parameter :: g0 =  9.80665 !average surface gravity m/s**2
  real(kind=dp), parameter :: e_constant = 2.718281828459045235360287471352 !eulers number 
 
  contains
   
  subroutine flip(original_list,flipped_list)
     real(kind=dp), intent(in) :: original_list(:)
     real(kind=dp), intent(inout) :: flipped_list(:)

     integer :: i, list_size

     list_size = size(original_list)
     
     flipped_list = original_list(list_size:1:-1)
     
     return
   end subroutine

  subroutine hydrostatic_factor(era_orography,speedy_orography,era_2m_temp,factor)
    !Hydrostatic  balance operator based on Baek et al 2009
    real(kind=dp), intent(in)  :: era_orography(:,:)
    real(kind=dp), intent(in)  :: speedy_orography(:,:) !Need to make sur
    real(kind=dp), intent(in)  :: era_2m_temp(:,:)
    
    real(kind=dp), intent(out) :: factor(:,:)

    !local vars
    real(kind=dp), parameter :: g0=9.0866 !m/s**2
    real(kind=dp), parameter :: gamma1=6.5e-3 !meter**-1
    real(kind=dp), parameter :: sigma_level1 = 0.95 !speedy sigma level 1 (0.95 of surfacep)
    real(kind=dp), parameter :: R=287 !J kg**-1 K**-1

    real(kind=dp), allocatable :: average_temp(:,:)

    real(kind=dp), allocatable :: delta_z(:,:)

   
    allocate(delta_z,mold=era_orography)
    delta_z = era_orography - speedy_orography 
   
    allocate(average_temp,mold=era_2m_temp)
    average_temp = era_2m_temp + 0.5*gamma1*delta_z

    factor = (1+ ((g0*delta_z)/(R*average_temp)) + (((g0*delta_z)/(R*average_temp))**2)/2.0_dp)
    
    return 
  end subroutine 

end module mod_utilities
