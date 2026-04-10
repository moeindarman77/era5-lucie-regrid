module mod_io
   use netcdf
   use mod_utilities, only : dp, sp

   implicit none
   
   integer :: error  

   contains 
     
     function file_exists(filename) result(result_bool)
       character(len=*), intent(in) :: filename
       logical                      :: result_bool

       !Very basic way to check if a file exists
       inquire(file=trim(filename),exist=result_bool)
       
     end function
    
     subroutine read_netcdf_4d_dp_era(varname,filename,var)
        character(len=*), intent(in) :: varname
        character(len=*), intent(in) :: filename
      
        real(kind=dp), allocatable, intent(out) :: var(:,:,:,:) 
        
        !Parmeter
        integer, parameter :: numofdims = 4 !We can assume its a 4d variable

        !Local netcdf variables
        integer :: ncid
        integer :: varid
 
        integer :: dimids(numofdims), dim_length(numofdims)
    
        integer :: i

        integer                       :: scale_length, offset_length
        real(kind=dp)                 :: scale_val, off_val
        character(len=:), allocatable :: offset_attribute, scale_attribute
         
        scale_attribute = 'scale_factor'
        offset_attribute = 'add_offset' 

        call nc_check(nf90_open(filename, nf90_nowrite, ncid))

        call nc_check(nf90_inq_varid(ncid,varname,varid))

        call nc_check(nf90_inquire_variable(ncid,varid,dimids=dimids))
    
        call nc_check(nf90_get_att(ncid,varid,scale_attribute,scale_val))
     
        call nc_check(nf90_get_att(ncid,varid,offset_attribute,off_val))

        do i=1,numofdims
           call nc_check(nf90_inquire_dimension(ncid,dimids(i),len=dim_length(i)))
        enddo 
      
        allocate(var(dim_length(1),dim_length(2),dim_length(3),dim_length(4)))

        call nc_check(nf90_get_var(ncid,varid,var))
        call nc_check(nf90_close(ncid))

        var = var * scale_val + off_val
        return 
     end subroutine 
    
     subroutine read_netcdf_4d_sp(varname,filename,var)
        character(len=*), intent(in) :: varname
        character(len=*), intent(in) :: filename

        real(kind=sp), allocatable, intent(out) :: var(:,:,:,:)

        !Parmeter
        integer, parameter :: numofdims = 4 !We can assume its a 4d variable

        !Local netcdf variables
        integer :: ncid
        integer :: varid

        integer :: dimids(numofdims), dim_length(numofdims)

        integer :: i

        call nc_check(nf90_open(filename, nf90_nowrite, ncid))

        call nc_check(nf90_inq_varid(ncid,varname,varid))

        call nc_check(nf90_inquire_variable(ncid,varid,dimids=dimids))

        do i=1,numofdims
           call nc_check(nf90_inquire_dimension(ncid,dimids(i),len=dim_length(i)))
        enddo

        allocate(var(dim_length(1),dim_length(2),dim_length(3),dim_length(4)))
        call nc_check(nf90_get_var(ncid,varid,var))
        call nc_check(nf90_close(ncid))
        return
     end subroutine

     subroutine read_netcdf_3d_dp(varname,filename,var)
        character(len=*), intent(in) :: varname
        character(len=*), intent(in) :: filename

        real(kind=dp), allocatable, intent(out) :: var(:,:,:)

        !Parmeter
        integer, parameter :: numofdims = 3 !We can assume its a 3d variable

        !Local netcdf variables
        integer :: ncid
        integer :: varid

        integer :: dimids(numofdims), dim_length(numofdims)

        integer :: i

        call nc_check(nf90_open(filename, nf90_nowrite, ncid))

        call nc_check(nf90_inq_varid(ncid,varname,varid))

        call nc_check(nf90_inquire_variable(ncid,varid,dimids=dimids))

        do i=1,numofdims
           call nc_check(nf90_inquire_dimension(ncid,dimids(i),len=dim_length(i)))
        enddo

        allocate(var(dim_length(1),dim_length(2),dim_length(3)))

        call nc_check(nf90_get_var(ncid,varid,var))
        call nc_check(nf90_close(ncid))
        return
     end subroutine
       
     subroutine read_netcdf_3d_dp_era(varname,filename,var)
        character(len=*), intent(in) :: varname
        character(len=*), intent(in) :: filename

        real(kind=dp), allocatable, intent(out) :: var(:,:,:)

        !Parmeter
        integer, parameter :: numofdims = 3 !We can assume its a 3d variable

        !Local netcdf variables
        integer :: ncid
        integer :: varid

        integer :: dimids(numofdims), dim_length(numofdims)

        integer :: i

        integer                       :: scale_length, offset_length
        real(kind=dp)                 :: scale_val, off_val
        character(len=:), allocatable :: offset_attribute, scale_attribute

        scale_attribute = 'scale_factor'
        offset_attribute = 'add_offset'

        call nc_check(nf90_open(filename, nf90_nowrite, ncid))

        call nc_check(nf90_inq_varid(ncid,varname,varid))

        call nc_check(nf90_inquire_variable(ncid,varid,dimids=dimids))

        if(varname /= 'tisr') then
          call nc_check(nf90_get_att(ncid,varid,scale_attribute,scale_val))

          call nc_check(nf90_get_att(ncid,varid,offset_attribute,off_val))
        endif 

        do i=1,numofdims
           call nc_check(nf90_inquire_dimension(ncid,dimids(i),len=dim_length(i)))
        enddo

        allocate(var(dim_length(1),dim_length(2),dim_length(3)))

        call nc_check(nf90_get_var(ncid,varid,var))
        call nc_check(nf90_close(ncid))

        if(varname == 'tisr') then
           var = var !* scale_val + off_val
        else 
           var = var * scale_val + off_val
        endif 
        return
     end subroutine 

     subroutine read_netcdf_3d_sp(varname,filename,var)
        character(len=*), intent(in) :: varname
        character(len=*), intent(in) :: filename

        real(kind=sp), allocatable, intent(out) :: var(:,:,:)

        !Parmeter
        integer, parameter :: numofdims = 3 !We can assume its a 3d variable

        !Local netcdf variables
        integer :: ncid
        integer :: varid

        integer :: dimids(numofdims), dim_length(numofdims)

        integer :: i

        call nc_check(nf90_open(filename, nf90_nowrite, ncid))

        call nc_check(nf90_inq_varid(ncid,varname,varid))

        call nc_check(nf90_inquire_variable(ncid,varid,dimids=dimids))

        do i=1,numofdims
           call nc_check(nf90_inquire_dimension(ncid,dimids(i),len=dim_length(i)))
        enddo

        allocate(var(dim_length(1),dim_length(2),dim_length(3)))

        call nc_check(nf90_get_var(ncid,varid,var))
        call nc_check(nf90_close(ncid))
        return
     end subroutine

     subroutine read_netcdf_2d_dp(varname,filename,var)
        character(len=*), intent(in) :: varname
        character(len=*), intent(in) :: filename

        real(kind=dp), allocatable, intent(out) :: var(:,:)

        !Parmeter
        integer, parameter :: numofdims = 2 !We can assume its a 2d variable

        !Local netcdf variables
        integer :: ncid
        integer :: varid

        integer :: dimids(numofdims), dim_length(numofdims)

        integer :: i

        call nc_check(nf90_open(filename, nf90_nowrite, ncid))

        call nc_check(nf90_inq_varid(ncid,varname,varid))

        call nc_check(nf90_inquire_variable(ncid,varid,dimids=dimids))

        do i=1,numofdims
           call nc_check(nf90_inquire_dimension(ncid,dimids(i),len=dim_length(i)))
        enddo

        allocate(var(dim_length(1),dim_length(2)))

        call nc_check(nf90_get_var(ncid,varid,var))
        call nc_check(nf90_close(ncid))
        return
     end subroutine 

     subroutine read_netcdf_1d_dp(varname,filename,var)
        character(len=*), intent(in) :: varname
        character(len=*), intent(in) :: filename

        real(kind=dp), allocatable, intent(out) :: var(:)

        !Parmeter
        integer, parameter :: numofdims = 1 !We can assume its a 3d variable

        !Local netcdf variables
        integer :: ncid
        integer :: varid

        integer :: dimids(numofdims), dim_length(numofdims)

        integer :: i

        call nc_check(nf90_open(filename, nf90_nowrite, ncid))

        call nc_check(nf90_inq_varid(ncid,varname,varid))

        call nc_check(nf90_inquire_variable(ncid,varid,dimids=dimids))

        do i=1,numofdims
           call nc_check(nf90_inquire_dimension(ncid,dimids(i),len=dim_length(i)))
        enddo

        allocate(var(dim_length(1)))

        call nc_check(nf90_get_var(ncid,varid,var))
        call nc_check(nf90_close(ncid))
        return
     end subroutine
     
     subroutine write_era_5_regridded(filename,grid5d,grid3d)
       use stringtype, only : string
       use mod_utilities, only : speedy_lat, gridx, gridy

       character(len=*), intent(in) :: filename
       
       real(kind=dp), intent(inout)    :: grid5d(:,:,:,:,:)
       real(kind=dp), intent(in)    :: grid3d(:,:,:)

       integer, parameter :: numdims4d=4, numdims3d=3, numspeedyvars=5
       integer :: dimsx,dimsy,dimsz,dimst

       integer :: file_id, xdim_id, ydim_id, zdim_id, timedim_id
       integer :: array_id, xvar_id, yvar_id
       integer :: start4d(numdims4d),varcount4d(numdims4d),start3d(numdims3d),varcount3d(numdims3d)
       integer :: arrdims4d(numdims4d),arrdims3d(numdims3d)

       integer :: i, counter

       real(kind=dp) :: lon(gridx)

       type(string) :: units(numspeedyvars)
       type(string) :: varname(numspeedyvars)

       units(1)%str = 'Kelvin'
       varname(1)%str = 'Temperature'

       units(2)%str = 'm/s'
       varname(2)%str = 'U-wind'

       units(3)%str = 'm/s'
       varname(3)%str = 'V-wind'

       units(4)%str = 'kg/kg'
       varname(4)%str = 'Specific_Humidity'

       units(5)%str = 'log(surfacepressure)'
       varname(5)%str = 'logp'  

       dimsx = size(grid5d,2)
       dimsy = size(grid5d,3)
       dimsz = size(grid5d,4)
       dimst = size(grid5d,5)

       varcount4d = (/ dimsx, dimsy, dimsz, dimst /)
       start4d = (/ 1, 1, 1, 1 /)

       varcount3d = (/ dimsx, dimsy, dimst /)
       start3d = (/ 1, 1, 1 /)

       lon = (/(real(counter)*3.75,counter=0,95)/)
       
       call nc_check(nf90_create(path=filename,cmode=NF90_CLOBBER,ncid=file_id))

       ! define the dimensions
       call nc_check(nf90_def_dim(file_id, 'Lon', dimsx, xdim_id))
       call nc_check(nf90_def_dim(file_id, 'Lat', dimsy, ydim_id))
       call nc_check(nf90_def_dim(file_id, 'Sigma_Level', dimsz, zdim_id))
       print *,'here'
       call nc_check(nf90_def_dim(file_id, 'Timestep', NF90_UNLIMITED, timedim_id))
 
       !Assign lat and lon ids and units
       call nc_check(nf90_def_var(file_id,'Lon',NF90_REAL,xdim_id,xvar_id))
       call nc_check(nf90_def_var(file_id,'Lat',NF90_REAL,ydim_id,yvar_id))

       call nc_check(nf90_put_att(file_id,xvar_id,"units",'degrees_north'))
       call nc_check(nf90_put_att(file_id,yvar_id,"units",'degrees_east'))
         
       ! now that the dimensions are defined, we can define variables on
       ! them,...
       arrdims4d = (/ xdim_id, ydim_id, zdim_id, timedim_id /)
       arrdims3d = (/ xdim_id, ydim_id, timedim_id /)

       do i=1, numspeedyvars-1
          call nc_check(nf90_def_var(file_id,varname(i)%str,NF90_REAL,arrdims4d,array_id))
          ! ...and assign units to them as an attribute 
          call nc_check(nf90_put_att(file_id, array_id, "units", units(i)%str))

          call nc_check(nf90_enddef(file_id))

          !Write out the values
          call nc_check(nf90_put_var(file_id, array_id, grid5d(i,:,:,:,:),start=start4d, count=varcount4d))

          call nc_check(nf90_redef(file_id))
       enddo
       !Lets do logp
       call nc_check(nf90_def_var(file_id,varname(5)%str,NF90_REAL,arrdims3d,array_id))

       call nc_check(nf90_put_att(file_id, array_id, "units", units(5)%str))

       call nc_check(nf90_enddef(file_id))

       !Write out the values
       call nc_check(nf90_put_var(file_id, array_id, grid3d,start=start3d, count=varcount3d))

       call nc_check(nf90_put_var(file_id, xvar_id, lon))

       call nc_check(nf90_put_var(file_id, yvar_id, speedy_lat))

       ! close; done
       call nc_check(nf90_close(file_id)) 
     end subroutine 

     subroutine nc_write_3d_dp(filename,varname,unitname,var3d)
       use mod_utilities, only : speedy_lat, gridx, gridy

       character(len=*), intent(in) :: filename,varname,unitname

       real(kind=dp), intent(in)    :: var3d(:,:,:)

       integer, parameter :: numdims3d=3
       integer :: dimsx,dimsy,dimsz,dimst

       integer :: file_id, xdim_id, ydim_id, zdim_id, timedim_id
       integer :: array_id, xvar_id, yvar_id
       integer :: start3d(numdims3d),varcount3d(numdims3d)
       integer :: arrdims3d(numdims3d)

       integer :: i, counter

       real(kind=dp) :: lon(gridx)

       dimsx = size(var3d,1)
       dimsy = size(var3d,2)
       dimst = size(var3d,3)


       varcount3d = (/ dimsx, dimsy, dimst /)
       start3d = (/ 1, 1, 1 /)

       lon = (/(real(counter)*3.75,counter=0,95)/)

       call nc_check(nf90_create(path=filename,cmode=NF90_CLOBBER,ncid=file_id))

       ! define the dimensions
       call nc_check(nf90_def_dim(file_id, 'Lon', dimsx, xdim_id))
       call nc_check(nf90_def_dim(file_id, 'Lat', dimsy, ydim_id))
       call nc_check(nf90_def_dim(file_id, 'Timestep', NF90_UNLIMITED,timedim_id))

       !Assign lat and lon ids and units
       call nc_check(nf90_def_var(file_id,'Lon',NF90_REAL,xdim_id,xvar_id))
       call nc_check(nf90_def_var(file_id,'Lat',NF90_REAL,ydim_id,yvar_id))

       call nc_check(nf90_put_att(file_id,xvar_id,"units",'degrees_north'))
       call nc_check(nf90_put_att(file_id,yvar_id,"units",'degrees_east'))
       ! now that the dimensions are defined, we can define variables on
       ! them,...
       arrdims3d = (/ xdim_id, ydim_id, timedim_id /)
        
       call nc_check(nf90_def_var(file_id,varname,NF90_REAL,arrdims3d,array_id))

       call nc_check(nf90_put_att(file_id, array_id, "units", unitname))

       call nc_check(nf90_enddef(file_id))

       !Write out the values
       call nc_check(nf90_put_var(file_id, array_id, var3d,start=start3d, count=varcount3d))

       call nc_check(nf90_put_var(file_id, xvar_id, lon))

       call nc_check(nf90_put_var(file_id, yvar_id, speedy_lat))
       ! close; done
       call nc_check(nf90_close(file_id))
     end subroutine 

     subroutine nc_write_2d_dp(filename,varname,unitname,var2d)
       use mod_utilities, only : speedy_lat, gridx, gridy

       character(len=*), intent(in) :: filename,varname,unitname

       real(kind=dp), intent(in)    :: var2d(:,:)

       integer, parameter :: numdims2d=2
       integer :: dimsx,dimsy,dimsz,dimst

       integer :: file_id, xdim_id, ydim_id
       integer :: array_id, xvar_id, yvar_id
       integer :: start2d(numdims2d),varcount2d(numdims2d)
       integer :: arrdims2d(numdims2d)

       integer :: i, counter

       real(kind=dp) :: lon(gridx)

       dimsx = size(var2d,1)
       dimsy = size(var2d,2)

       varcount2d = (/ dimsx, dimsy/)
       start2d = (/ 1, 1 /)

       lon = (/(real(counter)*3.75,counter=0,95)/)

       call nc_check(nf90_create(path=filename,cmode=NF90_CLOBBER,ncid=file_id))

       ! define the dimensions
       call nc_check(nf90_def_dim(file_id, 'Lon', dimsx, xdim_id))
       call nc_check(nf90_def_dim(file_id, 'Lat', dimsy, ydim_id))

       !Assign lat and lon ids and units
       call nc_check(nf90_def_var(file_id,'Lon',NF90_REAL,xdim_id,xvar_id))
       call nc_check(nf90_def_var(file_id,'Lat',NF90_REAL,ydim_id,yvar_id))

       call nc_check(nf90_put_att(file_id,xvar_id,"units",'degrees_north'))
       call nc_check(nf90_put_att(file_id,yvar_id,"units",'degrees_east'))
       ! now that the dimensions are defined, we can define variables on
       ! them,...
       arrdims2d = (/ xdim_id, ydim_id/)

       call nc_check(nf90_def_var(file_id,varname,NF90_REAL,arrdims2d,array_id))
      
       call nc_check(nf90_put_att(file_id, array_id, "units", unitname))

       call nc_check(nf90_enddef(file_id))

       !Write out the values
       call nc_check(nf90_put_var(file_id, array_id, var2d,start=start2d,count=varcount2d))

       call nc_check(nf90_put_var(file_id, xvar_id, lon))

       call nc_check(nf90_put_var(file_id, yvar_id, speedy_lat))
       ! close; done
       call nc_check(nf90_close(file_id))
     end subroutine

     subroutine add_var_to_file_4d(filename,varname,unitname,var4d)
       character(len=*), intent(in) :: filename,varname,unitname

       real(kind=dp), intent(in)    :: var4d(:,:,:,:)

       integer, parameter :: numdims4d=4
       integer :: dimsx,dimsy,dimsz,dimst

       integer :: file_id, xdim_id, ydim_id, zdim_id, timedim_id
       integer :: array_id, xvar_id, yvar_id
       integer :: start4d(numdims4d),varcount4d(numdims4d)
       integer :: arrdims4d(numdims4d)

       integer :: i, counter

       dimsx = size(var4d,1)
       dimsy = size(var4d,2)
       dimsz = size(var4d,3)
       dimst = size(var4d,4)
   
       start4d = (/1,1,1,1/)
       varcount4d = (/ dimsx, dimsy, dimsz, dimst /)

       call nc_check(nf90_open(filename,nf90_write,file_id))

       call nc_check(nf90_inq_dimid(file_id, 'Lon',xdim_id))
       call nc_check(nf90_inq_dimid(file_id, 'Lat',ydim_id))        
       call nc_check(nf90_inq_dimid(file_id, 'Sigma_Level',zdim_id))
       call nc_check(nf90_inq_dimid(file_id, 'Timestep',timedim_id))
  
       arrdims4d = (/ xdim_id, ydim_id, zdim_id, timedim_id /) 
      
       call nc_check(nf90_redef(file_id)) 
      
       call nc_check(nf90_def_var(file_id,varname,NF90_REAL,arrdims4d,array_id))
       ! ...and assign units to them as an attribute
       call nc_check(nf90_put_att(file_id, array_id, "units", unitname))

       call nc_check(nf90_enddef(file_id))

       !Write out the values
       call nc_check(nf90_put_var(file_id, array_id, var4d,start=start4d, count=varcount4d))

       !call nc_check(nf90_put_var(file_id, xvar_id, lon))

       !call nc_check(nf90_put_var(file_id, yvar_id, speedy_lat))
        call nc_check(nf90_close(file_id))
     end subroutine 

     subroutine add_var_to_file_3d(filename,varname,unitname,var3d)
       character(len=*), intent(in) :: filename,varname,unitname

       real(kind=dp), intent(in)    :: var3d(:,:,:)

       integer, parameter :: numdims3d=3
       integer :: dimsx,dimsy,dimsz,dimst

       integer :: file_id, xdim_id, ydim_id, zdim_id, timedim_id
       integer :: array_id, xvar_id, yvar_id
       integer :: start3d(numdims3d),varcount3d(numdims3d)
       integer :: arrdims3d(numdims3d)

       integer :: i, counter

       dimsx = size(var3d,1)
       dimsy = size(var3d,2)
       dimst = size(var3d,3)

       start3d = (/1,1,1/)
       varcount3d = (/ dimsx, dimsy, dimst /)

       call nc_check(nf90_open(filename,nf90_write,file_id))

       call nc_check(nf90_inq_dimid(file_id, 'Lon',xdim_id))
       call nc_check(nf90_inq_dimid(file_id, 'Lat',ydim_id))
       call nc_check(nf90_inq_dimid(file_id, 'Timestep',timedim_id))

       arrdims3d = (/ xdim_id, ydim_id, timedim_id /)

       call nc_check(nf90_redef(file_id))

       call nc_check(nf90_def_var(file_id,varname,NF90_REAL,arrdims3d,array_id))
       ! ...and assign units to them as an attribute
       call nc_check(nf90_put_att(file_id, array_id, "units", unitname))

       call nc_check(nf90_enddef(file_id))

       !Write out the values
       call nc_check(nf90_put_var(file_id, array_id, var3d,start=start3d, count=varcount3d))

       !call nc_check(nf90_put_var(file_id, xvar_id, lon))

       !call nc_check(nf90_put_var(file_id, yvar_id, speedy_lat))
       call nc_check(nf90_close(file_id))
     end subroutine
 
     subroutine nc_check(status)
         use netcdf
         integer, intent (in) :: status

         if(status /= nf90_noerr) then
            print *, trim(nf90_strerror(status))
            stop "Stopped"
         end if
      end subroutine nc_check 
end module mod_io
