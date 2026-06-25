module mod_thakkar
    use :: mod_types
    use :: mpi
    implicit none
    real(rp),allocatable    :: h_map(:,:)
    real(rp)                :: lx,ly,lz,dx,dy,dz,dxi,dyi,dzi
    integer                 :: nx,ny,nz,nx_data,ny_data
    contains
    subroutine read_thakkar_bin(myid)
        character(len=120)      ::  path_bin
        integer,intent(in)      ::  myid
        integer                 :: iunit,ierr,numarg
        character(len=1024)     :: c_iomsg
        logical                 :: is_io_fallback
        path_bin="../thakkar_roughness/roughness.bin"
        open(newunit=iunit,file=trim(path_bin),status="old",action="read",&
        access="stream", form="unformatted",iostat=ierr,iomsg=c_iomsg)
            if(ierr/=0)then
                if(myid == 0) print*, 'ERROR: opening the bin input file: ', trim(c_iomsg)
                if(myid == 0) print*, 'Aborting...'
                call MPI_FINALIZE(ierr)
                close(iunit)
                error stop
            end if
            ! read .bin file and fill h_map
            read(iunit,iostat=ierr,iomsg=c_iomsg) h_map
            if(ierr /= 0) then
                if(myid == 0) print*, 'ERROR: reading `bin file` namelist: ', trim(c_iomsg)
                if(myid == 0) print*, 'Aborting...'
                call MPI_FINALIZE(ierr)
                close(iunit)
                error stop
            end if            
        
    end subroutine read_thakkar_bin
    subroutine read_thakkar_nfo(myid)
        character(len=120)      ::  path_nfo
        integer,intent(in)      ::  myid
        integer                 :: iunit,ierr,numarg
        character(len=1024)     :: c_iomsg
        logical                 :: is_io_fallback
        namelist/roughnessinfo/&
                            lx,&
                            ly,&
                            lz,&
                            nx,&
                            ny,&
                            nz,&
                            nx_data,&
                            ny_data
        path_nfo="../thakkar_roughness/roughness.nfo"
        open(newunit=iunit,file=trim(path_nfo),status="old",action="read",&
        iostat=ierr,iomsg=c_iomsg)
            if(ierr/=0)then
                if(myid == 0) print*, 'ERROR: opening the roughness input file: ', trim(c_iomsg)
                if(myid == 0) print*, 'Aborting...'
                call MPI_FINALIZE(ierr)
                close(iunit)
                error stop
            end if
            read(iunit,nml=roughnessinfo,iostat=ierr,iomsg=c_iomsg)
            if(ierr /= 0) then
                if(myid == 0) print*, 'ERROR: reading `roughnessinfo` namelist: ', trim(c_iomsg)
                if(myid == 0) print*, 'Aborting...'
                call MPI_FINALIZE(ierr)
                close(iunit)
                error stop
            end if
            rewind(iunit)
            close(iunit)
            !
            dx=real(lx/nx,kind=rp)
            dy=real(ly/ny,kind=rp)
            dz=real(lz/nz,kind=rp)
            dxi=dx**(-1)
            dyi=dy**(-1)
            dzi=dz**(-1)
            allocate(h_map(nx_data,ny_data))
    end subroutine read_thakkar_nfo 
    subroutine correct_hmap(lo,hmap,l,n,dl,dix,diy,diz)
        ! the problem CaNS'velocities dont live in the same place 
        ! as the data of the hmap so wee need to adjust the data according to
        implicit none
        real(rp),intent(inout)                   :: hmap(:,:)
        integer , intent(in   ), dimension(3)    :: lo
        integer, intent(in), dimension(3)       :: n
        real(rp), intent(in), dimension(3)       :: l,dl
        integer,intent(in) ::  dix,diy,diz
        real(rp)    :: h_nx,h_ny
        real(rp)    :: x,y,z,gx,gy,gz
        integer     :: i,j,k,ii,jj,kk
        h_nx=size(hmap,1)
        h_ny=size(hmap,2)
        do k=0,n(3)
            do j=0,n(2)
                do i=0,n(1)
                    ii = lo(1)+i-1
                    jj = lo(2)+j-1
                    kk = lo(3)+k-1
                    ! these coordinates where the velocities of CaNS lives
                    x=(real(ii,rp) -0.5d0+ real(dix,rp)*0.5d0)*dl(1)
                    y=(real(jj,rp) -0.5d0+ real(diy,rp)*0.5d0)*dl(2)
                    z=(real(kk,rp) -0.5d0+ real(diz,rp)*0.5d0)*dl(3)

                    gx=x*h_nx/l(1)
                    gy=y*h_ny/l(2)
                    print*, x,gx,y,gy
                end do 
            end do 
        end do 

    end subroutine correct_hmap
    
end module mod_thakkar