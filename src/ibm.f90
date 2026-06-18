module mod_ibm
    use :: mod_types 
    use :: mod_param
    implicit none
    contains 
    ! we initalize the ibm coef. here
    ! at the main.f90
    ! check if the real location(x,y,z) is in the given body shape 
    logical function isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,x,y,z,n,l)
        implicit none
        logical , intent(in), dimension(0:1,3)   :: ibm_direction
        real(rp), intent(in), dimension(0:1,3)   :: amp_l
        integer , intent(in), dimension(0:1,3)   :: n_wave
        real(rp), intent(in), dimension(0:1,3)   :: l_0
        real(rp), intent(in), dimension(0:1,3)   :: phase_l
        integer , intent(in), dimension(3)       :: n
        real(rp), intent(in), dimension(3)       :: l
        real(rp),intent(in)                      :: x,y,z 
        real(rp),parameter                       :: pi = 3.141592653589793_rp
        real(rp)                                 :: height(0:1,3)
        integer                                  :: side,t
        real(rp)                                 :: xyz(3)
        integer                                  :: i,ii
        xyz = [x,y,z]
        do side = 0,1
            do t = 1,3
                i=modulo(t,3)+1
                ii=modulo(t+1,3)+1
                if(ibm_direction(side,t))then
                    height(side,t)=amp_l(side,i)*0.5_rp*(1._rp+sin(2._rp*pi*real(n_wave(side,i)*xyz(i)/l(i),rp)+phase_l(side,i)))+&
                                   amp_l(side,ii)*0.5_rp*(1._rp+sin(2._rp*pi*real(n_wave(side,ii)*xyz(ii)/l(ii),rp)+phase_l(side,ii)))
                else
                    height(side,t)=0._rp
                endif
            end do 
        end do 
        isInBody=.false.
        do t = 1,3
            if(ibm_direction(0,t))then
                height(0,t)=l_0(0,t)+height(0,t)
                if(xyz(t)<=height(0,t))then
                    isInBody=.true.
                endif
            endif
            if(ibm_direction(1,t))then
                height(1,t)=l_0(1,t)-height(1,t)
                if(xyz(t)>=height(1,t))then
                    isInBody=.true.
                endif
            endif   
        end do 
    end function isInbody
    ! 1st order IBM 
    ! we change the diL depending the velocity mask we are handling 
    ! e.g. we need to apply 1,0,0 for mask_u and we need to apply 0,1,0 for mask_v 
    ! 0,0,1 for mask_w
    subroutine set_ibm_staircase(lo,mask_id,dix,diy,diz,n,l,dl,ibm_direction,amp_l,n_wave,l_0,phase_l)
        implicit none
        real(rp), intent(in   ), dimension(3)       :: l
        real(rp), intent(in   ), dimension(3)       :: dl
        integer , intent(in   ), dimension(3)       :: n
        integer , intent(in   ), dimension(3)       :: lo
        logical,intent(inout),dimension(0:,0:,0:)   :: mask_id
        integer,intent(in)                          :: dix,diy,diz 
        integer                                     :: i,j,k
        integer                                     :: ii,jj,kk
        real(rp)                                    :: x,y,z
        logical , intent(in), dimension(0:1,3)      :: ibm_direction
        real(rp), intent(in), dimension(0:1,3)      :: amp_l
        integer , intent(in), dimension(0:1,3)      :: n_wave
        real(rp), intent(in), dimension(0:1,3)      :: l_0
        real(rp), intent(in), dimension(0:1,3)      :: phase_l
        print*, "***IBM coefficients are deploying***"
        do k = lbound(mask_id,3),ubound(mask_id,3)
            do j = lbound(mask_id,2),ubound(mask_id,2)
                do i = lbound(mask_id,1),ubound(mask_id,1)
                    ii = lo(1)+i-1
                    jj = lo(2)+j-1
                    kk = lo(3)+k-1
                    ! we create the real location of each velocity here
                    x = (real(ii,rp) -0.5d0+ real(dix,rp)*0.5d0)*dl(1)
                    y = (real(jj,rp) -0.5d0+ real(diy,rp)*0.5d0)*dl(2)
                    z = (real(kk,rp) -0.5d0+ real(diz,rp)*0.5d0)*dl(3)
                    if(isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,x,y,z,n,l).eqv..true.)then
                        mask_id(i,j,k) = .true.
                    endif 
                end do 
            end do 
        end do
    end subroutine set_ibm_staircase
    !2nd order scheme--laplacian settings
    subroutine set_ibm_2nd(lo,mask_id,laplacian_id,dix,diy,diz&
        ,n,l,dl,ibm_direction,amp_l,n_wave,l_0,phase_l)
        implicit none
        logical,intent(inout)                       :: mask_id(0:,0:,0:)
        real(rp), intent(in   ), dimension(3)       :: l
        real(rp), intent(in   ), dimension(3)       :: dl
        integer , intent(in   ), dimension(3)       :: n
        integer , intent(in   ), dimension(3)       :: lo
        real(rp),intent(inout),dimension(0:,0:,0:)  :: laplacian_id
        integer,intent(in)                          :: dix,diy,diz 
        integer                                     :: i,j,k
        integer                                     :: ii,jj,kk
        real(rp)                                    :: x,y,z,xp,xm,yp,ym,zp,zm
        logical , intent(in), dimension(0:1,3)      :: ibm_direction
        real(rp), intent(in), dimension(0:1,3)      :: amp_l
        integer , intent(in), dimension(0:1,3)      :: n_wave
        real(rp), intent(in), dimension(0:1,3)      :: l_0
        real(rp), intent(in), dimension(0:1,3)      :: phase_l
        real(rp)                                    :: lambda
        integer                                     :: n_dir
        do k = lbound(mask_id,3),ubound(mask_id,3)
            do j = lbound(mask_id,2),ubound(mask_id,2)
                do i = lbound(mask_id,1),ubound(mask_id,1)
                    ii = lo(1)+i-1
                    jj = lo(2)+j-1
                    kk = lo(3)+k-1
                    ! we create the real location of each velocity here
                    x = (real(ii,rp) -0.5d0+ real(dix,rp)*0.5d0)*dl(1)
                    y = (real(jj,rp) -0.5d0+ real(diy,rp)*0.5d0)*dl(2)
                    z = (real(kk,rp) -0.5d0+ real(diz,rp)*0.5d0)*dl(3)
                    if(isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,x,y,z,n,l).eqv..true.)then
                        mask_id(i,j,k) = .true.
                    endif 
                    xp=x+dl(1);xm=x-dl(1);yp=y+dl(2);ym=y-dl(2);zp=z+dl(3);zm=z-dl(3)
                    do n_dir=1,6
                        select case(n_dir)
                            case(1)
                                if(.not.mask_id(i,j,k).and.isInBody(ibm_direction,amp_l,n_wave,&
                                        l_0,phase_l,xp,y,z,n,l))then
                                            call calc_lambda(x,y,z,xp,1,lambda,ibm_direction,&
                                                            amp_l,n_wave,l_0,phase_l,n,l,dl)
                                            laplacian_id(i,j,k)=laplacian_id(i,j,k)+lambda
                                endif
                            case(2)
                                if(.not.mask_id(i,j,k).and.isInBody(ibm_direction,amp_l,n_wave,&
                                        l_0,phase_l,xm,y,z,n,l))then
                                            call calc_lambda(x,y,z,xm,1,lambda,ibm_direction,&
                                                            amp_l,n_wave,l_0,phase_l,n,l,dl)
                                            laplacian_id(i,j,k)=laplacian_id(i,j,k)+lambda
                                endif
                            case(3)
                                if(.not.mask_id(i,j,k).and.isInBody(ibm_direction,amp_l,n_wave,&
                                        l_0,phase_l,x,yp,z,n,l))then
                                            call calc_lambda(x,y,z,yp,2,lambda,ibm_direction,&
                                                            amp_l,n_wave,l_0,phase_l,n,l,dl)
                                            laplacian_id(i,j,k)=laplacian_id(i,j,k)+lambda
                                endif
                            case(4)
                                if(.not.mask_id(i,j,k).and.isInBody(ibm_direction,amp_l,n_wave,&
                                        l_0,phase_l,x,ym,z,n,l))then
                                            call calc_lambda(x,y,z,ym,2,lambda,ibm_direction,&
                                                            amp_l,n_wave,l_0,phase_l,n,l,dl)
                                            laplacian_id(i,j,k)=laplacian_id(i,j,k)+lambda
                                endif       
                            case(5)
                                if(.not.mask_id(i,j,k).and.isInBody(ibm_direction,amp_l,n_wave,&
                                        l_0,phase_l,x,y,zp,n,l))then
                                            call calc_lambda(x,y,z,zp,3,lambda,ibm_direction,&
                                                            amp_l,n_wave,l_0,phase_l,n,l,dl)
                                            laplacian_id(i,j,k)=laplacian_id(i,j,k)+lambda
                                endif
                            case(6)
                                if(.not.mask_id(i,j,k).and.isInBody(ibm_direction,amp_l,n_wave,&
                                        l_0,phase_l,x,y,zm,n,l))then
                                            call calc_lambda(x,y,z,zm,3,lambda,ibm_direction,&
                                                            amp_l,n_wave,l_0,phase_l,n,l,dl)
                                            laplacian_id(i,j,k)=laplacian_id(i,j,k)+lambda
                                endif
                        end select
                    end do
                end do 
            end do 
        end do
    end subroutine set_ibm_2nd

    subroutine calc_lambda(x,y,z,l_n,case_num,lambda,ibm_direction,amp_l,n_wave,l_0,phase_l,n,l,dl)
        implicit none
        logical , intent(in), dimension(0:1,3)      :: ibm_direction
        real(rp), intent(in), dimension(0:1,3)      :: amp_l
        integer , intent(in), dimension(0:1,3)      :: n_wave
        real(rp), intent(in), dimension(0:1,3)      :: l_0
        real(rp), intent(in), dimension(0:1,3)      :: phase_l
        real(rp), intent(in   ), dimension(3)       :: l
        real(rp), intent(in   ), dimension(3)       :: dl
        integer , intent(in   ), dimension(3)       :: n
        integer,intent(in)                          ::case_num
        real(rp),intent(in)                         :: x,y,z,l_n
        real(rp),intent(out)                        :: lambda
        real(rp)                                    :: l_fluid,l_solid,l_int,l_diff
        real(rp)                                    :: eps
        integer                                     ::  n_iter
        lambda=0._rp
        select case(case_num)
            case(1)
                eps = 1.e-10_rp*dl(1)
            case(2)
                eps = 1.e-10_rp*dl(2)
            case(3)
                eps = 1.e-10_rp*dl(3)
        end select
        select case(case_num)
            case(1)
                l_fluid = x;l_solid = l_n;l_int=0._rp;l_diff=0._rp
            case(2)
                l_fluid = y;l_solid = l_n;l_int=0._rp;l_diff=0._rp
            case(3)
                l_fluid = z;l_solid = l_n;l_int=0._rp;l_diff=0._rp
        end select
        do n_iter=1,60
            l_int=real((l_solid+l_fluid)/2._rp,kind=rp)
            select case(case_num)
                case(1)!x
                    if(isInBody(ibm_direction,amp_l,&
                    n_wave,l_0,phase_l,l_int,y,z,n,l))then
                        l_solid=l_int
                    else
                        l_fluid=l_int
                    endif

                case(2)!y
                    if(isInBody(ibm_direction,amp_l,&
                    n_wave,l_0,phase_l,x,l_int,z,n,l))then
                        l_solid=l_int
                    else
                        l_fluid=l_int
                    endif

                case(3)!z
                    if(isInBody(ibm_direction,amp_l,&
                    n_wave,l_0,phase_l,x,y,l_int,n,l))then
                        l_solid=l_int
                    else
                        l_fluid=l_int
                    endif
            end select
        enddo
        select case(case_num)
            case(1)
                l_diff=abs(x-l_int)
            case(2)
                l_diff=abs(y-l_int)
            case(3)
                l_diff=abs(z-l_int)
        end select
        if(l_diff<=eps)then
            l_diff=eps
        endif
        select case(case_num)
            case(1)
                lambda=real((1._rp/dl(1)**2)*(dl(1)/l_diff-1._rp),kind=rp)
            case(2)
                lambda=real((1._rp/dl(2)**2)*(dl(2)/l_diff-1._rp),kind=rp)
            case(3)
                lambda=real((1._rp/dl(3)**2)*(dl(3)/l_diff-1._rp),kind=rp)
        end select  
    end subroutine calc_lambda

    subroutine apply_ibm_staircase(field,mask_id,dt)
        implicit none
        real(rp),intent(inout),dimension(0:,0:,0:)  :: field
        logical,intent(in),dimension(0:,0:,0:)      :: mask_id
        real(rp),intent(in)                         :: dt
        integer :: i,j,k
        do k = lbound(field,3)+1,ubound(field,3)-1
            do j = lbound(field,2)+1,ubound(field,2)-1
                do i = lbound(field,1)+1,ubound(field,1)-1
                    if (mask_id(i,j,k).eqv..true.)then
                        field(i,j,k) = 0._rp 
                    endif
                end do 
            end do 
        end do
    end subroutine apply_ibm_staircase
end module mod_ibm