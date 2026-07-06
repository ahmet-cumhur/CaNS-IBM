module mod_thakkar
    use :: mod_types
    use :: mpi
    implicit none
    real(rp),allocatable    :: hmap_tha(:,:)
    real(rp)                :: lx_tha,ly_tha,lz_tha,dx_tha,dy_tha,dz_tha,dxi_tha,dyi_tha,dzi_tha
    integer                 :: nx_tha,ny_tha,nz_tha,nx_hmap_tha,ny_hmap_tha
    contains
    subroutine read_thakkar_bin(myid)
        character(len=120)      ::  path_bin
        integer,intent(in)      ::  myid
        integer                 :: iunit,ierr,numarg
        character(len=1024)     :: c_iomsg
        logical                 :: is_io_fallback
        path_bin="../src/thakkar_roughness/roughness.bin"
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
            read(iunit,iostat=ierr,iomsg=c_iomsg) hmap_tha
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
                            lx_tha,&
                            ly_tha,&
                            lz_tha,&
                            nx_tha,&
                            ny_tha,&
                            nz_tha,&
                            nx_hmap_tha,&
                            ny_hmap_tha
        path_nfo="../src/thakkar_roughness/roughness.nfo"
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
            dx_tha=real(lx_tha/nx_tha,kind=rp)
            dy_tha=real(ly_tha/ny_tha,kind=rp)
            dz_tha=real(lz_tha/nz_tha,kind=rp)
            dxi_tha=dx_tha**(-1)
            dyi_tha=dy_tha**(-1)
            dzi_tha=dz_tha**(-1)
            allocate(hmap_tha(0:nx_hmap_tha-1,0:ny_hmap_tha-1))
    end subroutine read_thakkar_nfo 
end module mod_thakkar