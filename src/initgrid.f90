! -
!
! SPDX-FileCopyrightText: Pedro Costa and the CaNS contributors
! SPDX-License-Identifier: MIT
!
! -
module mod_initgrid
  use mod_param, only:pi,is_gridpoint_natural_channel
  use mod_types
  implicit none
  private
  public initgrid,save_grid
  contains
  subroutine initgrid(gtype,n,gr,lz,dzc,dzf,zc,zf,is_periodic,&
                      visc,dzp_min_goal,dzp_max_goal,hch,hmax,myid)
    !
    ! initializes the non-uniform grid along z
    !
    implicit none
    integer, parameter :: CLUSTER_TWO_END   = 1, &
                          CLUSTER_ONE_END   = 2, &
                          CLUSTER_ONE_END_R = 3, &
                          CLUSTER_MIDDLE    = 4, &
                          CLUSTER_ROUGH     = 5

    integer , intent(in ) :: gtype,n
    real(rp), intent(in ) :: gr,lz
    real(rp), intent(out), dimension(0:n+1) :: dzc,dzf,zc,zf
    logical , intent(in ) :: is_periodic
    real(rp) :: z0
    integer :: k
    procedure (), pointer :: gridpoint => null()
    real(rp),intent(in),optional            :: visc,dzp_min_goal,dzp_max_goal
    real(rp),intent(in),optional            :: hch,hmax
    real(rp)                                :: re_tau
    logical                                 :: diagnostics=.true.
    integer, intent(in),optional            :: myid
    real(rp)                                :: dzfp,dzcp
    integer                                 :: iunit
    character(len=*), parameter :: fmt_dp = '(*(es24.16e3,1x))', &
                                   fmt_sp = '(*(es15.8e2,1x))'
#if !defined(_SINGLE_PRECISION)
    character(len=*), parameter :: fmt_rp = fmt_dp
#else
    character(len=*), parameter :: fmt_rp = fmt_sp
#endif
    select case(gtype)
    case(CLUSTER_TWO_END)
      gridpoint => gridpoint_cluster_two_end
    case(CLUSTER_ONE_END)
      gridpoint => gridpoint_cluster_one_end
    case(CLUSTER_ONE_END_R)
      gridpoint => gridpoint_cluster_one_end_r
    case(CLUSTER_MIDDLE)
      gridpoint => gridpoint_cluster_middle
    case default
      gridpoint => gridpoint_cluster_two_end
    end select
    !present check
    if(gtype==CLUSTER_ROUGH)then
      if(.not.present(visc).or.&
         .not.present(dzp_min_goal).or.&
         .not.present(dzp_max_goal).or.&
         .not.present(hch).or.&
         .not.present(hmax))then
        error stop "Needed rough parameters arent accesible"
      end if
    endif
    if(gtype==CLUSTER_ROUGH.and.diagnostics)then
      !calc re_tau for diag
      re_tau=0._rp
      !we assume u_tau=1.0 after scaling
      re_tau=hch/visc
    endif
    !
    ! step 1) determine coordinates of cell faces zf
    !
    zf(0) = 0.
    if(.not.is_gridpoint_natural_channel.and.gtype/=CLUSTER_ROUGH) then
      do k=1,n
        z0  = (k-0.)/(1.*n)
        call gridpoint(gr,z0,zf(k))
      end do
    elseif(gtype==CLUSTER_ROUGH)then
      call gridpoint_rough(n,lz,zf,visc,dzp_min_goal,dzp_max_goal,hch,hmax)
    else
      do k=1,n
        call gridpoint_natural(k,n,zf(k))
      end do
    end if
    if(gtype/=CLUSTER_ROUGH)then
      zf(1:n) = zf(1:n)*lz
    endif
    !
    if(abs(gr) < epsilon(1._rp) .and. gtype/=CLUSTER_ROUGH) then
      !
      ! for uniform grids, set constant spacing to avoid round-off errors
      !
      dzf(:) = lz/(1.*n)
      dzc(:) = lz/(1.*n)
    else
      !
      ! step 2) determine grid spacing between faces dzf
      !
      do k=1,n
        dzf(k) = zf(k)-zf(k-1)
      end do
      if(.not.is_periodic) then
        dzf(0  ) = dzf(1)
        dzf(n+1) = dzf(n)
      else
        dzf(0  ) = dzf(n)
        dzf(n+1) = dzf(1)
      end if
      !
      ! step 3) determine grid spacing between centers dzc
      !
      do k=0,n
        dzc(k) = .5*(dzf(k)+dzf(k+1))
      end do
      if(.not.is_periodic) then
        dzc(n+1) = dzc(n)
      else
        dzc(n+1) = dzc(1)
      end if
    end if
    !
    ! step 4) compute coordinates of cell centers zc and faces zf
    !
    zc(0)    = -dzc(0)/2.
    zf(0)    = 0.
    do k=1,n+1
      zc(k) = zc(k-1) + dzc(k-1)
      zf(k) = zf(k-1) + dzf(k)
    end do
    !lets get diagnostics 
    if(diagnostics.eqv..true..and.gtype==CLUSTER_ROUGH.and.myid==0)then
      open(newunit=iunit,file=trim("grid_diagnostic.out"),&
          status="replace", action="write")
      do k=0,n
        dzfp=dzf(k)*re_tau/hch
        dzcp=dzc(k)*re_tau/hch
        write(iunit,fmt_rp) zf(k),zc(k),dzf(k),dzc(k),dzfp,dzcp
      end do
      close(iunit)
    end if
  end subroutine initgrid
  !IBM aware grid stretching with geometric growth
  subroutine gridpoint_rough(nz,lz,zf,visc,dzp_min_goal,dzp_max_goal,hch,hmax)
        implicit none
        integer, intent(in)                         :: nz
        real(rp), intent(in)                        :: lz
        real(rp), intent(out), dimension(0:nz+1)    :: zf
        integer                                     :: k,j,i
        real(rp),intent(in)                         :: visc,hch,hmax
        real(rp),intent(in)                         :: dzp_max_goal,dzp_min_goal
        real(rp)                                    :: re_tau,u_tau,hmax_p
        integer                                     :: nz_fine,nz_tra,nz_coarse_it,nz_fine_tra
        real(rp)                                    :: dz_fine,dzp_coarse_it
        real(rp)                                    :: gr
        real(rp)                                    :: lz_fine,lz_tra,lz_coarse_it
        real(rp)                                    :: err
        logical                                     :: odd
        real(rp)                                    :: lzhp,lzh
        !init
        lzh=real(lz/2,kind=rp)
        gr=1.1_rp
        i=0
        err=0.0
        dzp_coarse_it=0.0
        !lets get the basic number for fine region
        zf(0) = 0.0
        !change if needed
        u_tau=1.0;
        re_tau=u_tau*hch/visc;
        dz_fine=dzp_min_goal*hch/re_tau
        lzhp=lzh*re_tau/hch 
        !now we need to divide the grid into 2 equal parts
        if (modulo(nz,2)/=0)then
            odd=.true.
        else
            odd=.false.
        endif 
        !fine region
        hmax_p=(hmax/hch*re_tau)
        nz_fine=ceiling(hmax_p/dzp_min_goal)
        lz_fine=real(nz_fine,kind=rp)*dzp_min_goal
        do k=1,nz_fine
            zf(k)=zf(k-1)+dz_fine
        end do
        !geometric growth region iterative calc.
        err=(dzp_max_goal-dzp_coarse_it)/dzp_max_goal
        do while(abs(err)>1e-4)
            i=i+1
            gr=gr-0.001*err
            err=(dzp_max_goal-dzp_coarse_it)/dzp_max_goal
            !print*,lz_fine,nz_fine
            nz_tra=floor(log(dzp_max_goal/dzp_min_goal)/log(gr))
            lz_tra=dzp_min_goal*gr*((gr**nz_tra)-1)/(gr-1)
            !print*, lz_tra,nz_tra
            if(odd.eqv..false.)then
                lz_coarse_it=lzhp-(lz_fine+lz_tra)
                nz_coarse_it=(nz/2)-(nz_fine+nz_tra)
                dzp_coarse_it=lz_coarse_it/nz_coarse_it
            else
                lz_coarse_it=lzhp-(lz_fine+lz_tra)
                nz_coarse_it=(nz/2)-(nz_fine+nz_tra)
                dzp_coarse_it=lz_coarse_it/real(nz_coarse_it+0.5,kind=rp)!this is the half cell since the nz is odd
            endif

            !print*, lz_coarse_it,nz_coarse_it,dzp_coarse_it
            !print*, gr
            !print*, i
            if(i>1e6)then
                print*, "gr couldnt be estimated please use other method or different number of cells."
                exit
            endif
        end do 
        !now we now number of cells that we need to apply
        nz_fine_tra=nz_fine+nz_tra
        do k=nz_fine+1,nz_fine_tra
            j=k-nz_fine
            zf(k)=zf(k-1)+dz_fine*gr**j
        end do
        
        do k=nz_fine_tra+1,nz/2
            zf(k)=zf(k-1)+dzp_coarse_it*hch/re_tau
        end do
        do k=0,(nz/2)
            zf(nz-k)=lz-zf(k)
        end do
        zf(nz+1)=zf(nz)+(zf(nz-1)-zf(nz-2))
        !print*,zf
    end subroutine gridpoint_rough
  !
  ! grid stretching functions
  ! see e.g., Fluid Flow Phenomena -- A Numerical Toolkit, by P. Orlandi
  !           Pirozzoli et al. JFM 788, 614–639 (commented)
  !
  subroutine gridpoint_cluster_two_end(alpha,z0,z)
    !
    ! clustered at the two sides
    !
    implicit none
    real(rp), intent(in) :: alpha,z0
    real(rp), intent(out) :: z
    if(alpha > epsilon(0._rp)) then
      z = 0.5*(1.+tanh((z0-0.5)*alpha)/tanh(alpha/2.))
      !z = 0.5*(1.+erf( (z0-0.5)*alpha)/erf( alpha/2.))
    else
      z = z0
    end if
  end subroutine gridpoint_cluster_two_end
  subroutine gridpoint_cluster_one_end(alpha,z0,z)
    !
    ! clustered at the lower side
    !
    implicit none
    real(rp), intent(in) :: alpha,z0
    real(rp), intent(out) :: z
    if(alpha > epsilon(0._rp)) then
      z = 1.0*(1.+tanh((z0-1.0)*alpha)/tanh(alpha/1.))
      !z = 1.0*(1.+erf( (z0-1.0)*alpha)/erf( alpha/1.))
    else
      z = z0
    end if
  end subroutine gridpoint_cluster_one_end
  subroutine gridpoint_cluster_one_end_r(alpha,r0,r)
    !
    ! clustered at the upper side
    !
    implicit none
    real(rp), intent(in ) :: alpha,r0
    real(rp), intent(out) :: r
    if(alpha > epsilon(0._rp)) then
      r = 1._rp-1.0_rp*(1._rp+tanh((1._rp-r0-1.0_rp)*alpha)/tanh(alpha/1._rp))
      !r = 1._rp-1.0_rp*(1._rp+erf( (1._rp-r0-1.0_rp)*alpha)/erf( alpha/1._rp))
    else
      r = r0
    end if
  end subroutine gridpoint_cluster_one_end_r
  subroutine gridpoint_cluster_middle(alpha,z0,z)
    !
    ! clustered in the middle
    !
    implicit none
    real(rp), intent(in) :: alpha,z0
    real(rp), intent(out) :: z
    if(alpha > epsilon(0._rp)) then
      if(     z0 <= 0.5) then
        z = 0.5*(1.-1.+tanh(2.*alpha*(z0-0.))/tanh(alpha))
        !z = 0.5*(1.-1.+erf( 2.*alpha*(z0-0.))/erf( alpha))
      else if(z0  > 0.5) then
        z = 0.5*(1.+1.+tanh(2.*alpha*(z0-1.))/tanh(alpha))
        !z = 0.5*(1.+1.+erf( 2.*alpha*(z0-1.))/erf( alpha))
      end if
    else
      z = z0
    end if
  end subroutine gridpoint_cluster_middle
  subroutine gridpoint_natural(kg,nzg,z,kb_a,alpha_a,c_eta_a,dyp_a)
    !
    ! a physics-based, 'natural' grid stretching function for wall-bounded turbulence
    ! see Pirozzoli & Orlandi, JCP 439 - 110408 (2021)
    !
    ! clustered at the two sides
    !
    implicit none
    real(rp), parameter :: kb_p     = 32._rp,    &
                           alpha_p  = pi/1.5_rp, &
                           c_eta_p  = 0.8_rp,    &
                           dyp_p    = 0.05_rp
    integer , intent(in ) :: kg,nzg
    real(rp), intent(out) :: z
    real(rp), intent(in ), optional :: kb_a,alpha_a,c_eta_a,dyp_a
    real(rp)                        :: kb  ,alpha  ,c_eta  ,dyp
    real(rp) :: retau,n,k
    !
    ! handle input parameters
    !
    kb    = kb_p   ; if(present(kb_a   )) kb    = kb_a
    alpha = alpha_p; if(present(alpha_a)) alpha = alpha_a
    c_eta = c_eta_p; if(present(c_eta_a)) c_eta = c_eta_a
    dyp   = dyp_p  ; if(present(dyp_a  )) dyp   = dyp_a
    !
    ! determine retau
    !
    n = nzg/2._rp
    retau = 1._rp/(1._rp+(n/kb)**2)*(dyp*n+(3._rp/4._rp*alpha*c_eta*n)**(4._rp/3._rp)*(n/kb)**2)
#if 0
    if(kg==1) print*,'Grid targeting Retau = ',retau
#endif
    k = 1._rp*min(kg,(nzg-kg))
    !
    ! dermine z/(2h)
    !
    z = 1._rp/(1._rp+(k/kb)**2)*(dyp*k+(3._rp/4._rp*alpha*c_eta*k)**(4._rp/3._rp)*(k/kb)**2)/(2._rp*retau)
    if( kg > nzg-kg ) z = 1._rp-z
  end subroutine gridpoint_natural
  !
  subroutine save_grid(datadir,fname,ng,zc,zf,dzc,dzf)
    !
    ! saves grid for post-processing
    !
#if defined(_USE_HDF5)
    use hdf5
#endif
    implicit none
    character(len=*), intent(in) :: datadir,fname
    integer , intent(in), dimension(3) :: ng
    real(rp), intent(in), dimension(0:) :: zc,zf,dzc,dzf
    integer :: iunit,k
#if defined(_USE_HDF5)
    integer :: ierr_h5
    integer(HID_T) :: file_id,dset,space
    integer(HSIZE_T), dimension(1) :: dims
#endif
    !
    open(newunit=iunit,file=trim(datadir)//trim(fname)//'.bin',action='write',form='unformatted',access='stream',status='replace')
    write(iunit) dzc(1:ng(3)),dzf(1:ng(3)),zc(1:ng(3)),zf(1:ng(3))
    close(iunit)
    open(newunit=iunit,file=trim(datadir)//trim(fname)//'.out',status='replace')
    do k=0,ng(3)+1
      write(iunit,*) 0.,zf(k),zc(k),dzf(k),dzc(k)
    end do
    close(iunit)
#if defined(_USE_HDF5)
    call h5open_f(ierr_h5)
    call h5fcreate_f(trim(datadir)//trim(fname)//'.h5',H5F_ACC_TRUNC_F,file_id,ierr_h5)
    dims(1) = int(ng(3),HSIZE_T)
    call h5screate_simple_f(1,dims,space,ierr_h5)
    call h5dcreate_f(file_id,'rc',HDF5_REAL_RP(),space,dset,ierr_h5)
    call h5dwrite_f(dset,HDF5_REAL_RP(),zc(1:ng(3)),dims,ierr_h5)
    call h5dclose_f(dset,ierr_h5)
    call h5dcreate_f(file_id,'rf',HDF5_REAL_RP(),space,dset,ierr_h5)
    call h5dwrite_f(dset,HDF5_REAL_RP(),zf(1:ng(3)),dims,ierr_h5)
    call h5dclose_f(dset,ierr_h5)
    call h5dcreate_f(file_id,'drc',HDF5_REAL_RP(),space,dset,ierr_h5)
    call h5dwrite_f(dset,HDF5_REAL_RP(),dzc(1:ng(3)),dims,ierr_h5)
    call h5dclose_f(dset,ierr_h5)
    call h5dcreate_f(file_id,'drf',HDF5_REAL_RP(),space,dset,ierr_h5)
    call h5dwrite_f(dset,HDF5_REAL_RP(),dzf(1:ng(3)),dims,ierr_h5)
    call h5dclose_f(dset,ierr_h5)
    call h5sclose_f(space,ierr_h5)
    call h5fclose_f(file_id,ierr_h5)
    call h5close_f(ierr_h5)
#endif
  end subroutine save_grid
#if defined(_USE_HDF5)
  integer(HID_T) function HDF5_REAL_RP()
    use hdf5
    implicit none
    if(rp == dp) then
      HDF5_REAL_RP = H5T_NATIVE_DOUBLE
    else
      HDF5_REAL_RP = H5T_NATIVE_REAL
    end if
  end function HDF5_REAL_RP
#endif
end module mod_initgrid
