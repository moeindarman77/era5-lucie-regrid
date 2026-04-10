module mod_interp2d
  use mod_utilities, only : dp

  implicit none

  integer, parameter :: speedy_nlvls = 8

  real(kind=dp), parameter :: speedy_sigma_levels(speedy_nlvls) = (/25,95,200,340,510,685,835,950/)/1000.0_dp
  
  contains 

  pure function binary_search_dp(array,val,tol) result(idx) 
    !binary search of sorted array to find the index of array closest to val 
    real(kind=dp), intent(in) :: array(:)
    real(kind=dp), intent(in) :: val
    
    real(kind=dp), intent(in), optional :: tol
    
    integer                   :: idx

    !Local stuff
    real(kind=dp) :: d

    integer :: left, middle, right

    if(present(tol) == .true.) then
      d = tol
    else 
      d = 1e-2
    endif 

    left = 1
    right = size(array) 

    do 
       if(left > right) then
         exit
       endif 

       middle = nint((left+right)/2.0)

       if(abs(array(middle) - val) <= d) then
          idx = middle 
          return
       elseif(array(middle) > val) then
          right = middle -1
       else
          left = middle + 1
       endif 
    end do  
    idx = right
  end function binary_search_dp
  
  pure function binary_search_dp_reverse(array,val,tol) result(idx)
    !binary search of reversed sorted array to find the index of array closest to val
    !e.g 90 ... -89 -90 instead of -90 .... 89 90
    real(kind=dp), intent(in) :: array(:)
    real(kind=dp), intent(in) :: val

    real(kind=dp), intent(in), optional :: tol

    integer                   :: idx

    !Local stuff
    real(kind=dp) :: d

    integer :: left, middle, right

    if(present(tol) == .true.) then
      d = tol
    else
      d = 1e-2
    endif

    left = size(array)
    right = 1

    do
       if(left<right) then
         exit
       endif

       middle = nint((left+right)/2.0)

       if(abs(array(middle) - val) <= d) then
          idx = middle
          return
       elseif(array(middle) > val) then
          right = middle +1
       else
          left = middle - 1
       endif
    end do
    idx = left
  end function binary_search_dp_reverse 

  pure function bilinear_interp(x_grid,y_grid,f,x_target,y_target) result(interp_val)
    
    real(kind=dp), intent(in)  :: x_grid(:)
    real(kind=dp), intent(in)  :: y_grid(:)
    
    real(kind=dp), intent(in)  :: f(:,:)
    
    real(kind=dp), intent(in)  :: x_target, y_target

    real(kind=dp)              :: interp_val

    !local vars 
    real(kind=dp) :: norm
    
    real(kind=dp) :: x1, x2
    real(kind=dp) :: y1, y2
    
    integer       :: i, j

    i = binary_search_dp(x_grid,x_target)
    j = binary_search_dp_reverse(y_grid,y_target) !Because era is stupid and lats are in reverse order 

    x1 = x_grid(i)
    x2 = x_grid(i+1)

    y1 = y_grid(j)
    y2 = y_grid(j-1) !!Because era is stupid and lats are in reverse order
   
    norm = (x2 - x1)*(y2 - y1)
    
    interp_val = (f(i,j)*(x2-x_target)*(y2-y_target) + f(i+1,j)*(x_target-x1)*(y2-y_target) + &
            f(i,j+1)*(x2-x_target)*(y_target-y1) + f(i+1, j+1)*(x_target-x1)*(y_target-y1))/norm
  end function bilinear_interp
  
  function log_interp1d(x1,x2,val1,val2,new_x) result(newval)
    real(kind=dp), intent(in)  :: x1
    real(kind=dp), intent(in)  :: x2

    real(kind=dp), intent(in)  :: val1,val2

    real(kind=dp), intent(in)  :: new_x

    real(kind=dp)              :: newval
    
    newval = val1 + (new_x - x1)*((val2 - val1)/(x2 - x1))
 
  end function log_interp1d

  subroutine bspline_interface1d(oldx,var,newx,newvar)
    use bspline_module
    use bspline_kinds_module, only: wp

    real(kind=dp), intent(in) :: oldx(:)
    real(kind=dp), intent(in) :: var(:)

    real(kind=dp), intent(in) :: newx(:)
    real(kind=dp), intent(inout) :: newvar(:)

    integer                   :: x_size, y_size
    integer                   :: i, j

    !Bspline stuff
    ! Troy's paper says "1D cubic B-spline" but the original code (which produced
    ! the LUCIE training data we're matching against) used kx=3, inherited from
    ! FITPACK convention where k=3 means cubic. In bspline-fortran kx=3 = quadratic.
    ! Empirically kx=3 reproduces the training data best, so we keep it.
    integer,parameter :: kx     = 3    !! x bspline order (FITPACK "cubic" = bspline-fortran kx=3)
    integer,parameter :: idx    = 0    !! [[db2val]] input
    integer           :: nx            !! number of points in x dimension in original grid
    integer           :: nx_new        !! number of points in x dimension for new grid
    integer,parameter :: iknot  = 0    !! automatically select the knots

    real(wp), allocatable            :: tx(:)       !! x knots
    real(wp), allocatable            :: bcoefs(:)
    real(wp), allocatable            :: w0(:)       !! work array for db1val (size 3*kx)
    real(wp) :: val,tru,err,errmax
    integer  :: iflag  !! status flag
    integer  :: inbvx,inbvy,iloy

    logical, parameter  :: extrap=.True.

    x_size = size(newx)

    nx = size(oldx)

    allocate(tx(nx+kx))
    allocate(bcoefs(nx))
    allocate(w0(3*kx))

    call db1ink(oldx,nx,var,kx,iknot,tx,bcoefs,iflag)

    inbvx = 1
    do i=1,x_size
          call db1val(newx(i),idx,tx,nx,kx,bcoefs,val,iflag,inbvx,w0,extrap)
          if (iflag/=0) print *,'error calling db1val',get_status_message(iflag)
          if(iflag/=0) stop
          newvar(i) = val
    enddo
    return
  end subroutine 

  subroutine bspline_interface(oldx,oldy,var,newx,newy,newvar)
    !Regrids data using a nth order bspline 
    !input oldx,oldy,and var(oldx,oldy) and the target new x and y
    !output var(newx,newy)

    use bspline_module
    use bspline_kinds_module, only: wp

    real(kind=dp), intent(in) :: oldx(:), oldy(:)
    real(kind=dp), intent(inout) :: var(:,:)

    real(kind=dp), intent(in) :: newx(:), newy(:)
    real(kind=dp), intent(inout) :: newvar(:,:)

    integer                   :: x_size, y_size
    integer                   :: i, j

    !Bspline stuff
    ! Troy's paper says "2D quadratic B-spline" but the original code (which
    ! produced the LUCIE training data) used kx=ky=2, inherited from FITPACK
    ! convention where k=2 means quadratic. In bspline-fortran kx=2 = linear.
    ! Empirically kx=ky=2 reproduces logp at machine precision, so we keep it.
    integer,parameter :: kx     = 2    !! x bspline order (FITPACK "quadratic" = bspline-fortran kx=2)
    integer,parameter :: ky     = 2    !! y bspline order (FITPACK "quadratic" = bspline-fortran kx=2)
    integer,parameter :: idx    = 0    !! [[db2val]] input
    integer,parameter :: idy    = 0    !! [[db2val]] input
    integer           :: nx            !! number of points in x dimension in original grid
    integer           :: ny            !! number of points in y dimension in original grid
    integer           :: nx_new        !! number of points in x dimension for new grid
    integer           :: ny_new        !! number of points in y dimension for new grid
    integer,parameter :: iknot  = 0    !! automatically select the knots

    real(wp), allocatable            :: tx(:)       !! x knots
    real(wp), allocatable            :: ty(:)       !! y knots
    real(wp), allocatable            :: bcoef(:,:)  !! b-spline coefficients (separate from input fcn)
    real(wp), allocatable            :: w1(:)       !! work array for db2val (size ky)
    real(wp), allocatable            :: w0(:)       !! work array for db2val (size 3*max(kx,ky))
    real(wp) :: val,tru,err,errmax
    integer  :: iflag  !! status flag
    integer  :: inbvx,inbvy,iloy

    logical, parameter :: extrap=.True.

    x_size = size(newx)
    y_size = size(newy)

    nx = size(oldx)
    ny = size(oldy)

    allocate(tx(nx+kx))
    allocate(ty(ny+ky))
    ! Use a separate bcoef buffer instead of aliasing `var` (the input fcn) — newer
    ! bspline-fortran versions disallow aliasing fcn and bcoef in db2ink (intent
    ! conflict on `pure` subroutines). The aliased pattern was the original code.
    allocate(bcoef(nx,ny))
    allocate(w1(ky))
    allocate(w0(3*max(kx,ky)))

    call db2ink(oldx,nx,oldy,ny,var,kx,ky,iknot,tx,ty,bcoef,iflag)
    if (iflag/=0) print *,'error calling db2ink',get_status_message(iflag)
    if (iflag/=0) print *, oldy
    inbvx = 1
    inbvy = 1
    iloy  = 1
    do i=1,x_size
       do j=1,y_size
          call db2val(newx(i),newy(j),idx,idy,tx,ty,nx,ny,kx,ky,bcoef,val,iflag,&
               inbvx,inbvy,iloy,w1,w0,extrap)
          if (iflag/=0) print *,'error calling db2val',get_status_message(iflag)
          if(iflag/=0) stop
          newvar(i,j) = val
       enddo
    enddo
    
    return
  end subroutine 
  
  subroutine era_p_level_to_speedy_sigma_level(era_surface_p,era_var,era_levels,regridded_data)
    !Subroutine that interpolates era pressure level data to speedy sigma level
    !data using bsplines 

    real(kind=dp), intent(in) :: era_surface_p(:,:), era_var(:,:,:)
    real(kind=dp), intent(in) :: era_levels(:)
    
    real(kind=dp), intent(inout) :: regridded_data(:,:,:) !x y and sigma

    !local stuff
    real(kind=dp), allocatable :: era_sigma_levels_at_xy(:)

    real(kind=dp)              :: p1,p2
    integer                    :: i, j, k
    integer                    :: era_x, era_y, speedy_level_size
    
    era_x = size(era_surface_p,1)
    era_y = size(era_surface_p,2)
    
    !allocate(regridded_data(era_x,era_y,speedy_nlvls))
    allocate(era_sigma_levels_at_xy(speedy_nlvls))

    do i=1,era_x
       do j=1, era_y
          era_sigma_levels_at_xy = (era_surface_p(i,j)*speedy_sigma_levels)/100.0_dp
          call bspline_interface1d(era_levels,era_var(i,j,:),era_sigma_levels_at_xy,regridded_data(i,j,:))
       enddo 
    enddo  
  end subroutine
 
  subroutine era_p_level_to_speedy_p_level(speedy_surface_p,era_var,era_levels,regridded_data)
    !Subroutine that interpolates era pressure level data to closest model
    !pressure level using bsplines
    !NOTE that speedy_surface_p is a monthly averaged array from jan 1980

    real(kind=dp), intent(in) :: speedy_surface_p(:,:), era_var(:,:,:)
    real(kind=dp), intent(in) :: era_levels(:)

    real(kind=dp), intent(out) :: regridded_data(:,:,:) !x y and sigma

    !local stuff
    real(kind=dp), allocatable :: speedy_p_at_sigma_levels_at_xy(:)

    real(kind=dp)              :: p1,p2
    integer                    :: i, j, k
    integer                    :: speedy_x, speedy_y, speedy_level_size

    speedy_x = size(speedy_surface_p,1)
    speedy_y = size(speedy_surface_p,2)

    !allocate(regridded_data(era_x,era_y,speedy_nlvls))
    allocate(speedy_p_at_sigma_levels_at_xy(speedy_nlvls))
    do i=1,speedy_x
       do j=1, speedy_y
          speedy_p_at_sigma_levels_at_xy = (speedy_surface_p(i,j)*speedy_sigma_levels)
          !Sometimes speedy_p_at_sigma_levels_at_xy is out of levels range !TODO
          call bspline_interface1d(era_levels,era_var(i,j,:),speedy_p_at_sigma_levels_at_xy,regridded_data(i,j,:))
          !if(any(regridded_data(i,j,:) > 350.0)) then
          !  print *, 'speedy_p_at_sigma_levels_at_xy',speedy_p_at_sigma_levels_at_xy
          !  print *, 'era_var(i,j,:)',era_var(i,j,:)
          !  print *, 'regridded_data(i,j,:)',regridded_data(i,j,:)
          !elseif(any(regridded_data(i,j,:) < 150.0)) then
          !  print *, 'speedy_p_at_sigma_levels_at_xy',speedy_p_at_sigma_levels_at_xy
          !  print *, 'era_var(i,j,:)',era_var(i,j,:)
          !  print *, 'regridded_data(i,j,:)',regridded_data(i,j,:)
          !endif 
       enddo
    enddo
  end subroutine 
end module mod_interp2d
