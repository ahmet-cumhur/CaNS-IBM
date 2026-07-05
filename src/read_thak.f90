module mod_thakkar
    use :: mod_types
    use :: mpi
    implicit none
    real(rp),allocatable    :: hmap(:,:)
    real(rp)                :: lx,ly,lz,dx,dy,dz,dxi,dyi,dzi
    integer                 :: nx,ny,nz,nx_hmap,ny_hmap
    contains
    subroutine read_thakkar_bin(myid)
        character(len=120)      ::  path_bin
        integer,intent(in)      ::  myid
        integer                 :: iunit,ierr,numarg
        character(len=1024)     :: c_iomsg
        logical                 :: is_io_fallback
        path_bin="/home/cumhur/codes/CaNS-IBM/thakkar_roughness/roughness.bin"
        open(newunit=iunit,file=trim(path_bin),status="old",action="read",&
        access="stream", form="unformatted",iostat=ierr,iomsg=c_iomsg)
            if(ierr/=0)then
                if(myid == 0) print*, 'ERROR: opening the bin input file: ', trim(c_iomsg)
                if(myid == 0) print*, 'Aborting...'
                call MPI_FINALIZE(ierr)
                close(iunit)
                error stop
            end if
            ! read .bin file and fill hmap
            read(iunit,iostat=ierr,iomsg=c_iomsg) hmap
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
                            nx_hmap,&
                            ny_hmap
        path_nfo="/home/cumhur/codes/CaNS-IBM/thakkar_roughness/roughness.nfo"
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
            allocate(hmap(0:nx_hmap-1,0:ny_hmap-1))
    end subroutine read_thakkar_nfo 
end module mod_thakkar