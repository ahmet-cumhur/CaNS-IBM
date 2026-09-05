module mod_ibm
    use mod_types 
    use mod_param
    use mod_thakkar
    implicit none
    contains 
    ! we initalize the ibm coef. here
    ! at the main.f90
    ! check if the real location(x,y,z) is in the given body shape 
    logical function isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,x,y,z,n,l,&
                                hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap)
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
        real(rp),intent(in),optional             :: hmap(0:,0:)
        real(rp),intent(in),optional             :: l1_hmap,l2_hmap 
        integer,intent(in),optional              :: n1_hmap,n2_hmap
        real(rp)                                 :: dl_hmap(2)
        integer                                  :: i1(2),i2(2)
        real(rp)                                 :: r1(2),r2(2)
        real(rp)                                 :: w1(2),w2(2)
        logical                                  :: macdonald
        macdonald=.true.
        xyz = [x,y,z]
        do side = 0,1
            do t = 1,3
                i=modulo(t,3)+1
                ii=modulo(t+1,3)+1
                if(ibm_direction(side,t))then
                    ! use wave wall  
                    if(.not. use_hmap)then
                        if(.not.macdonald)then
                        height(side,t)=amp_l(side,i)*0.5_rp*(1._rp+sin(2._rp*pi*&
                                        real(n_wave(side,i)*xyz(i)/l(i),rp)+phase_l(side,i)))+&
                                        amp_l(side,ii)*0.5_rp*(1._rp+sin(2._rp*pi*&
                                        real(n_wave(side,ii)*xyz(ii)/l(ii),rp)+phase_l(side,ii)))
                        else
                            if(side == 0)then
                                height(side,t)=amp_l(side,i)*(1._rp + &
                                                cos(2._rp*pi*real(n_wave(side,i)*xyz(i)/l(i),rp))*&
                                                cos(2._rp*pi*real(n_wave(side,ii)*xyz(ii)/l(ii),rp)))
                            else
                                height(side,t)=amp_l(side,i)*(1._rp - &
                                                cos(2._rp*pi*real(n_wave(side,i)*xyz(i)/l(i),rp))*&
                                                cos(2._rp*pi*real(n_wave(side,ii)*xyz(ii)/l(ii),rp)))
                            endif
                        endif
                    ! use hmap                
                    else
                         if(trim(hmap_mode)=="fit")then
                            call get_hmap_loc(side,n1_hmap,n2_hmap,xyz(i),xyz(ii),l(i),l(ii),dl_hmap,i1,i2,r1,r2,w1,w2)
                            ! Bilinear interpolation 
                            height(side,t)=(hmap(i1(1),i2(2))*(w1(1))+hmap(i1(2),i2(2))*(w1(2)))*(w2(2))+&
                                           (hmap(i1(1),i2(1))*(w1(1))+hmap(i1(2),i2(1))*(w1(2)))*(w2(1))
                        elseif(trim(hmap_mode)=="normal")then
                            ! here we need to change the shape of the hmap otherwise we will use different sized dl_hmap...
                            call get_hmap_loc(side,n1_hmap,n2_hmap,xyz(i),xyz(ii),l1_hmap,l2_hmap,dl_hmap,i1,i2,r1,r2,w1,w2)
                            ! Bilinear interpolation 
                            height(side,t)=(hmap(i1(1),i2(2))*(w1(1))+hmap(i1(2),i2(2))*(w1(2)))*(w2(2))+&
                                           (hmap(i1(1),i2(1))*(w1(1))+hmap(i1(2),i2(1))*(w1(2)))*(w2(1))
                        else
                            print*,"unknown hmap mode entered..."
                            stop
                        endif               
                    endif
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
    subroutine get_hmap_loc(side,n1_hmap,n2_hmap,loc_1,loc_2,&
                            l1_hmap,l2_hmap,dl_hmap,i1,i2,r1,r2,w1,w2)
        integer,intent(in)                              :: n1_hmap,n2_hmap
        real(rp),intent(inout)                          :: loc_1,loc_2
        real(rp),intent(in)                             :: l1_hmap,l2_hmap
        real(rp),intent(out)                            :: dl_hmap(2)
        integer,intent(out)                             :: i1(2),i2(2)
        real(rp),intent(out)                            :: r1(2),r2(2)
        real(rp),intent(out)                            :: w1(2),w2(2)
        integer,intent(in)                              :: side
        !     i1(1),i2(2)------------|------------i1(2),i2(2)
        !          |                                    |
        !          |        loc_1,loc_2--->hmap_val     |
        !          |                                    |
        !     i1(1),i2(1)------------|-------------i1(2),i2(1)
        ! additionally we need to shift upper wall lx_hmap/2,ly_hmap/2
        ! so we add a check
        i1(:)=0;i2(:)=0;
        r1(:)=0;r2(:)=0;
        dl_hmap(:)=0._rp;
        dl_hmap(1)=l1_hmap/n1_hmap
        dl_hmap(2)=l2_hmap/n2_hmap

        select case(side)
            case(0)
                loc_1=loc_1
                loc_2=loc_2
            case(1)
                loc_1=modulo(loc_1+(0.5*n1_hmap*dl_hmap(1)),(n1_hmap*dl_hmap(1)))
                loc_2=modulo(loc_2+(0.5*n2_hmap*dl_hmap(2)),(n2_hmap*dl_hmap(2)))
        end select
        
        i1(1)=floor(real(loc_1/dl_hmap(1),kind=rp))
        i1(2)=i1(1)+1

        r1(1)=i1(1)*dl_hmap(1)
        r1(2)=i1(2)*dl_hmap(1)

        i1(1)=modulo(i1(1),n1_hmap)
        i1(2)=modulo(i1(2),n1_hmap)

        i2(1)=floor(real(loc_2/dl_hmap(2),kind=rp))
        i2(2)=i2(1)+1

        r2(1)=i2(1)*dl_hmap(2)
        r2(2)=i2(2)*dl_hmap(2)

        i2(1)=modulo(i2(1),n2_hmap)
        i2(2)=modulo(i2(2),n2_hmap)

        w1(1)=abs(loc_1-r1(2))/dl_hmap(1)
        w1(2)=abs(loc_1-r1(1))/dl_hmap(1)
        w2(1)=abs(loc_2-r2(2))/dl_hmap(2)
        w2(2)=abs(loc_2-r2(1))/dl_hmap(2)

    end subroutine get_hmap_loc

    ! 1st order IBM 
    ! we change the diL depending the velocity mask we are handling 
    ! e.g. we need to apply 1,0,0 for mask_u and we need to apply 0,1,0 for mask_v 
    ! 0,0,1 for mask_w
    subroutine set_ibm_staircase(lo,mask_id,dix,diy,diz,n,l,dl,ibm_direction,amp_l,n_wave,l_0,phase_l,&
                                hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap,zc,zf)
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
        real(rp),intent(in),optional                :: hmap(0:,0:)
        integer,intent(in),optional                 :: n1_hmap,n2_hmap
        real(rp),intent(in),optional                :: l1_hmap,l2_hmap
        real(rp),intent(in),dimension(0:),optional  :: zc,zf           

        do k = lbound(mask_id,3),ubound(mask_id,3)
            do j = lbound(mask_id,2),ubound(mask_id,2)
                do i = lbound(mask_id,1),ubound(mask_id,1)
                    ii = lo(1)+i-1
                    jj = lo(2)+j-1
                    kk = lo(3)+k-1
                    ! we create the real location of each velocity here
                    x = (real(ii,rp) -0.5d0+ real(dix,rp)*0.5d0)*dl(1)
                    y = (real(jj,rp) -0.5d0+ real(diy,rp)*0.5d0)*dl(2)
                    if(diz/=1)then
                    ! this means we are looking for either u or v so their location is at z center
                    ! we gonna use the senter of zc
                        z = zc(k)
                        ! we use k inestead of kk since kk is the global and k is the local array index
                    else
                    ! else than we are looking for the w which is located on the z face
                        z = zf(k)
                    endif
                    if(isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,x,y,z,n,l,&
                                hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap).eqv..true.)then
                        mask_id(i,j,k) = .true. ! this means we are in the solid
                    endif                    
                end do 
            end do 
        end do
    end subroutine set_ibm_staircase
    !2nd order scheme--laplacian settings
    subroutine set_ibm_2nd(lo,mask_id,laplacian_id,dix,diy,diz&
        ,n,l,dl,ibm_direction,amp_l,n_wave,l_0,phase_l,hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap,&
        zc,zf,dzc,dzf,use_hmap,band_id)
        implicit none
        logical,intent(in)                          :: mask_id(0:,0:,0:)
        real(rp), intent(in   ), dimension(3)       :: l
        real(rp), intent(in   ), dimension(3)       :: dl
        integer , intent(in   ), dimension(3)       :: n
        integer , intent(in   ), dimension(3)       :: lo
        real(rp),intent(inout),dimension(0:,0:,0:)  :: laplacian_id
        integer,intent(in)                          :: dix,diy,diz 
        integer                                     :: i,j,k,ip,im,jp,jm,kp,km
        integer                                     :: ii,jj,kk,c
        real(rp)                                    :: x,y,z,xp,xm,yp,ym,zp,zm
        logical , intent(in), dimension(0:1,3)      :: ibm_direction
        real(rp), intent(in), dimension(0:1,3)      :: amp_l
        integer , intent(in), dimension(0:1,3)      :: n_wave
        real(rp), intent(in), dimension(0:1,3)      :: l_0
        real(rp), intent(in), dimension(0:1,3)      :: phase_l
        real(rp)                                    :: lambda
        integer                                     :: n_dir
        real(rp),intent(in),optional                :: hmap(0:,0:)
        integer,intent(in),optional                 :: n1_hmap,n2_hmap
        real(rp),intent(in),optional                :: l1_hmap,l2_hmap
        real(rp),intent(in),dimension(0:),optional  :: zc,zf
        real(rp),intent(in),dimension(0:),optional  :: dzc,dzf   
        real(rp)                                    :: dzf_l,dzc_l
        logical, intent(in)                         :: use_hmap
        real(rp)                                    :: hmax,hmin
        real(rp)                                    :: dl_int,l_int
        integer                                     :: n_hidden,ncand
        logical                                     :: calc_inBetween
        logical                                     :: ibm_diagnostic
        logical,dimension(0:,0:,0:),intent(inout)   :: band_id
        calc_inBetween=.false.
        ibm_diagnostic=.false.
        if(use_hmap)then
            if (.not.present(hmap)) then
                error stop "use_hmap=T but hmap not present"
            endif
            n_hidden=0;
            ncand=0;
            dl_int=0._rp;
            l_int=0._rp;
            hmax=0._rp;
            hmin=0._rp;
            hmin=minval(hmap);
            hmax=maxval(hmap);
            !hmax=hmax-hmin;
        else
            print*,"some problem between use_hmap and hmap"! doesnt mean problem
            n_hidden=0;
            ncand=0;
            dl_int=0._rp;
            l_int=0._rp;
            hmax=0._rp;
            hmin=0._rp;
        endif
        do k = 1,n(3)
            do j = lbound(mask_id,2),ubound(mask_id,2)
                do i = lbound(mask_id,1),ubound(mask_id,1)
                    ii = lo(1)+i-1
                    jj = lo(2)+j-1
                    kk = lo(3)+k-1
                    dzf_l=0._rp;dzc_l=0._rp
                    xp=x+dl(1);xm=x-dl(1);yp=y+dl(2);ym=y-dl(2);
                    call get_grid_loc(lo,i,j,k,dl,zc,zf,dix,diy,diz,x,y,z,xp,xm,yp,ym,zp,zm)
                    do n_dir=1,6
                        select case(n_dir)
                            case(1)
                                ! xp
                                if(.not.mask_id(i,j,k).and.isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,xp,y,z,n,l,&
                                hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then
                                            call calc_lambda(x,y,z,xp,1,lambda,ibm_direction,amp_l,n_wave,l_0,&
                                                            phase_l,n,l,dl,hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap)
                                            laplacian_id(i,j,k)=laplacian_id(i,j,k)+lambda
                                            band_id(i,j,k) = .true. ! this means on band and fluid ! these are for friction calc.
                                endif
                                !mask_id(i,j,k) already has the velocity body information
                                if(calc_inBetween)then 
                                    ! this means we are close to body region and both cells that we check are fluid now
                                    ! we check if in between part has body!
                                    if(z<hmax.and..not.mask_id(i,j,k).and..not.isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,&
                                                                                xp,y,z,n,l,hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then
                                        ncand=ncand+1
                                        do c=1,14
                                            ! we check 15 times
                                            dl_int=dl(1)*real(c,kind=rp)/15
                                            l_int=x+dl_int
                                            if(isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,l_int,y,z,n,l,&
                                                    hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then
                                                        !call calc_lambda(x,y,z,l_int,1,lambda,ibm_direction,amp_l,n_wave,l_0,&
                                                        !        phase_l,n,l,dl,hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap)
                                                        !laplacian_id(i,j,k)=laplacian_id(i,j,k)+lambda
                                                        n_hidden=n_hidden+1
                                                exit ! we need to exit otherwise we gonna keep adding to laplacian !note#2 not sure what should be done
                                            endif
                                        end do
                                        
                                    endif
                                endif
                            case(2)
                                ! xm
                                if(.not.mask_id(i,j,k).and.isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,xm,y,z,n,l,&
                                hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then
                                            call calc_lambda(x,y,z,xm,1,lambda,ibm_direction,amp_l,n_wave,l_0,&
                                                            phase_l,n,l,dl,hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap)
                                            laplacian_id(i,j,k)=laplacian_id(i,j,k)+lambda
                                            band_id(i,j,k) = .true. ! this means on band
                                endif
                                if(calc_inBetween)then 
                                    if(z<hmax.and..not.mask_id(i,j,k).and..not.isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,&
                                                                                xm,y,z,n,l,hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then
                                        ncand=ncand+1
                                        do c=1,14
                                            ! we check 4 times since 5th time is on the neigbour which it is in fluid
                                            dl_int=dl(1)*real(c,kind=rp)/15
                                            l_int=x-dl_int
                                            if(isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,l_int,y,z,n,l,&
                                                    hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then
                                                        !call calc_lambda(x,y,z,l_int,1,lambda,ibm_direction,amp_l,n_wave,l_0,&
                                                        !        phase_l,n,l,dl,hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap)
                                                        !laplacian_id(i,j,k)=laplacian_id(i,j,k)+lambda
                                                        n_hidden=n_hidden+1
                                                exit 
                                            endif
                                        end do
                                    endif
                                endif
                            case(3)
                                ! yp
                                if(.not.mask_id(i,j,k).and.isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,x,yp,z,n,l,&
                                hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then
                                            call calc_lambda(x,y,z,yp,2,lambda,ibm_direction,amp_l,n_wave,l_0,&
                                                            phase_l,n,l,dl,hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap)
                                            laplacian_id(i,j,k)=laplacian_id(i,j,k)+lambda
                                            band_id(i,j,k) = .true. ! this means on band
                                endif
                                if(calc_inBetween)then 
                                    if(z<hmax.and..not.mask_id(i,j,k).and..not.isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,&
                                                                                x,yp,z,n,l,hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then
                                        ncand=ncand+1
                                        do c=1,14
                                            dl_int=dl(2)*real(c,kind=rp)/15
                                            l_int=y+dl_int
                                            if(isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,x,l_int,z,n,l,&
                                                    hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then
                                                        !call calc_lambda(x,y,z,l_int,2,lambda,ibm_direction,amp_l,n_wave,l_0,&
                                                        !        phase_l,n,l,dl,hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap)
                                                        !laplacian_id(i,j,k)=laplacian_id(i,j,k)+lambda
                                                        n_hidden=n_hidden+1
                                                exit 
                                            endif
                                        end do
                                    endif
                                endif
                            case(4)
                                ! ym
                                if(.not.mask_id(i,j,k).and.isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,x,ym,z,n,l,&
                                hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then
                                            call calc_lambda(x,y,z,ym,2,lambda,ibm_direction,amp_l,n_wave,l_0,&
                                                            phase_l,n,l,dl,hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap)
                                            laplacian_id(i,j,k)=laplacian_id(i,j,k)+lambda
                                            band_id(i,j,k) = .true. ! this means on band
                                endif
                                if(calc_inBetween)then 
                                    if(z<hmax.and..not.mask_id(i,j,k).and..not.isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,&
                                                                                x,ym,z,n,l,hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then
                                        ncand=ncand+1
                                        do c=1,14
                                            dl_int=dl(2)*real(c,kind=rp)/15
                                            l_int=y-dl_int
                                            if(isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,x,l_int,z,n,l,&
                                                    hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then
                                                        !call calc_lambda(x,y,z,l_int,2,lambda,ibm_direction,amp_l,n_wave,l_0,&
                                                        !        phase_l,n,l,dl,hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap)
                                                        !laplacian_id(i,j,k)=laplacian_id(i,j,k)+lambda
                                                        n_hidden=n_hidden+1
                                                exit 
                                            endif
                                        end do
                                    endif
                                endif       
                            case(5)
                                ! zp 
                                if(diz==0)then
                                    dzc_l=dzc(k)
                                    dzf_l=dzf(k)
                                else
                                    dzf_l=dzf(k+1)
                                    dzc_l=dzc(k)
                                endif
                                if(.not.mask_id(i,j,k).and.isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,x,y,zp,n,l,&
                                hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then
                                            call calc_lambda(x,y,z,zp,3,lambda,ibm_direction,amp_l,n_wave,l_0,&
                                                            phase_l,n,l,dl,hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap,dzf_l,dzc_l,diz)
                                            laplacian_id(i,j,k)=laplacian_id(i,j,k)+lambda
                                            band_id(i,j,k) = .true. ! this means on band
                                endif
                            case(6)
                                ! zm
                                if(diz==0)then
                                    dzc_l=dzc(k-1)
                                    dzf_l=dzf(k)
                                else
                                    dzf_l=dzf(k)
                                    dzc_l=dzc(k)
                                endif 
                                if(.not.mask_id(i,j,k).and.isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,x,y,zm,n,l,&
                                hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then
                                            call calc_lambda(x,y,z,zm,3,lambda,ibm_direction,amp_l,n_wave,l_0,&
                                                            phase_l,n,l,dl,hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap,dzf_l,dzc_l,diz)
                                            laplacian_id(i,j,k)=laplacian_id(i,j,k)+lambda
                                            band_id(i,j,k) = .true. ! this means on band
                                endif
                        end select
                    end do
                end do 
            end do 
        end do
        if(ibm_diagnostic)then
            c=0
            print*,"max height of the hmap: ",hmax,hmin
            print*, "possible problematic locations: ",n_hidden,ncand
            print*, "number of band locations: ",count(band_id)
            print*, "number of IBM locations: ",count(mask_id)
            do k=1,n(3)
                do j=1,n(2)
                    do i=1,n(1)
                        if(band_id(i,j,k).and.mask_id(i,j,k))then
                            print*,"some problems with ibm processing"
                            print*,"some points are in the band and solid mask"
                            c=c+1
                        endif
                    enddo
                enddo
            enddo
            print*,"overlapping cells: ", c
        endif
    end subroutine set_ibm_2nd

    subroutine calc_lambda(x,y,z,l_n,case_num,lambda,ibm_direction,amp_l,n_wave,l_0,phase_l,n,l,dl,&
                            hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap,dzf,dzc,diz,wall_loc)
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
        real(rp),intent(in),optional                :: hmap(0:,0:)
        integer,intent(in),optional                 :: n1_hmap,n2_hmap
        real(rp),intent(in),optional                :: l1_hmap,l2_hmap
        real(rp),intent(in),optional                :: dzf,dzc
        integer,intent(in),optional                 :: diz   
        real(rp)                                    :: dz
        real(rp),intent(out),optional,dimension(3)  :: wall_loc
        lambda=0._rp
        dz=0._rp
        if(case_num==3)then
            if(diz==0)then
                ! this means we are looking for u or v
                ! their location is at zc 
                dz=dzc
            else
                ! this means we are looking for w
                ! its location is at zf 
                dz=dzf
            endif
        endif
        select case(case_num)
            case(1)
                eps = 1.e-10_rp*dl(1)
            case(2)
                eps = 1.e-10_rp*dl(2)
            case(3)
                eps = 1.e-10_rp*dz
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
                    if(isInBody(ibm_direction,amp_l,n_wave,l_0,phase_l,l_int,y,z,n,l,&
                                hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then
                        l_solid=l_int
                    else
                        l_fluid=l_int
                    endif

                case(2)!y
                    if(isInBody(ibm_direction,amp_l,n_wave,l_0,phase_l,x,l_int,z,n,l,&
                                hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then
                        l_solid=l_int
                    else
                        l_fluid=l_int
                    endif

                case(3)!z
                    if(isInBody(ibm_direction,amp_l,n_wave,l_0,phase_l,x,y,l_int,n,l,&
                                hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then
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
                ! imporant note here  we need both dzc and
                ! dzf since their gradinet either at cell center or face or vice versa
                !                      !
                lambda=real((1._rp/(dzc*dzf))*((dz)/l_diff-1._rp),kind=rp)
        end select 
        ! return the wall location
        if(present(wall_loc))then
            select case(case_num)
            case(1)
                wall_loc=[l_int,y,z]
            case(2)
                wall_loc=[x,l_int,z]
            case(3)
                wall_loc=[x,y,l_int]
            end select
        endif 
    end subroutine calc_lambda
    subroutine calc_grad_dist(grad_dist_id,lo,ibm_direction,amp_l,n_wave,l_0,phase_l,n,l,dl,zc,zf,&
                                band_id,dix,diy,diz,hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap,dzf,dzc)
        ! routine creates the tangent and normal plane for 3D IBM
        real(rp),dimension(:,:),intent(inout)                   :: grad_dist_id                      
        logical , intent(in), dimension(0:1,3)                  :: ibm_direction
        real(rp), intent(in), dimension(0:1,3)                  :: amp_l
        integer , intent(in), dimension(0:1,3)                  :: n_wave
        real(rp), intent(in), dimension(0:1,3)                  :: l_0
        real(rp), intent(in), dimension(0:1,3)                  :: phase_l
        real(rp), intent(in   ), dimension(3)                   :: l
        integer , intent(in   ), dimension(3)                   :: n
        real(rp),intent(in),optional                            :: hmap(0:,0:)
        integer,intent(in),optional                             :: n1_hmap,n2_hmap
        real(rp),intent(in),optional                            :: l1_hmap,l2_hmap
        real(rp),intent(in),dimension(0:),optional              :: dzf,dzc

        logical, intent(in), dimension(0:,0:,0:)                :: band_id
        integer,intent(in)                                      :: dix,diy,diz
        integer                                                 :: i,j,k,c,wc,nc
        real(rp),dimension(0:2,6,0:26)                          :: wall_loc
        logical,dimension(6,0:26)                               :: wall_loc_log
        integer,intent(in),dimension(3)                         :: lo
        real(rp),intent(in),dimension(0:)                       :: zc,zf
        real(rp), intent(in),dimension(3)                       :: dl
        real(rp)                                                :: x,y,z
        real(rp)                                                :: xp,xm,yp,ym,zp,zm
        real(rp)                                                :: lambda!just to fullfill the calc_lamda
        integer                                                 :: di,dj,dk,in,jn,kn,count
        real(rp),dimension(3)                                   :: wall_loc_real
        real(rp)                                                :: dzf_l,dzc_l
        real(rp),dimension(3,2)                                 :: plane
        real(rp),dimension(3)                                   :: nv,cp,bp
        real(rp)                                                :: angle_plane,out_plane,grad

        ! init vars
        grad_dist_id(:,:) = 0._rp 
        count = 0
        angle_plane=0._rp;out_plane=0._rp;grad=0._rp
        wall_loc_log(:,:)=.false.
        wall_loc_real(:)=0._rp
        wall_loc(:,:,:)=0._rp
        x=0._rp;y=0._rp;z=0._rp;
        dzc_l=0._rp;dzf_l=0._rp
        xp=0._rp;xm=0._rp;yp=0._rp;ym=0._rp;zp=0._rp;zm=0._rp
        ! calc part
        do k=(lbound(band_id,3)+1),(ubound(band_id,3)-1)
            do j=(lbound(band_id,2)+1),(ubound(band_id,2)-1)
                do i=(lbound(band_id,1)+1),(ubound(band_id,1)-1)
                    if(band_id(i,j,k))then ! this means we are in band
                        wall_loc(:,:,:)=0._rp
                        wall_loc_log(:,:)=.false.
                        wc=0
                       
                        do nc=0,26 ! time to check all possible near cells
                            
                            di=modulo(nc,3)-1
                            dj=modulo(nc/3,3)-1
                            dk=modulo(nc/9,3)-1
                            in=i+di;jn=j+dj;kn=k+dk;
                            if(band_id(in,jn,kn))then
                                x=0._rp;y=0._rp;z=0._rp;xp=0._rp;xm=0._rp;yp=0._rp;ym=0._rp;zp=0._rp;zm=0._rp
                                call get_grid_loc(lo,in,jn,kn,dl,zc,zf,dix,diy,diz,x,y,z,xp,xm,yp,ym,zp,zm) 
                                do c=1,6    
                                    select case(c)
                                        case(1)!xp
                                            if(isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,xp,y,z,n,l,&
                                            hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then 

                                                call calc_lambda(x,y,z,xp,1,lambda,ibm_direction,amp_l,n_wave,l_0,phase_l,n,l,dl,&
                                                hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap,dzf_l,dzc_l,diz,wall_loc_real)
                                                wc=wc+1 !wc describes the wall counter
                                                wall_loc(:,c,nc)=wall_loc_real
                                                wall_loc_log(c,nc)=.true.
                                                wall_loc_real(:)=0._rp

                                            endif   
                                        case(2)!i-1
                                            if(isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,xm,y,z,n,l,&
                                            hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then 

                                                call calc_lambda(x,y,z,xm,1,lambda,ibm_direction,amp_l,n_wave,l_0,phase_l,n,l,dl,&
                                                hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap,dzf_l,dzc_l,diz,wall_loc_real)
                                                wc=wc+1 !wc describes the wall counter
                                                wall_loc(:,c,nc)=wall_loc_real
                                                wall_loc_log(c,nc)=.true.
                                                wall_loc_real(:)=0._rp

                                            endif   
                                        case(3)!j+1
                                            if(isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,x,yp,z,n,l,&
                                            hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then 
                                            
                                                call calc_lambda(x,y,z,yp,2,lambda,ibm_direction,amp_l,n_wave,l_0,phase_l,n,l,dl,&
                                                hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap,dzf_l,dzc_l,diz,wall_loc_real)
                                                wc=wc+1 !wc describes the wall counter
                                                wall_loc(:,c,nc)=wall_loc_real
                                                wall_loc_log(c,nc)=.true.
                                                wall_loc_real(:)=0._rp

                                            endif    
                                        case(4)!j-1
                                            if(isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,x,ym,z,n,l,&
                                            hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then 
                                                call calc_lambda(x,y,z,ym,2,lambda,ibm_direction,amp_l,n_wave,l_0,phase_l,n,l,dl,&
                                                hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap,dzf_l,dzc_l,diz,wall_loc_real)
                                                wc=wc+1 !wc describes the wall counter
                                                wall_loc(:,c,nc)=wall_loc_real
                                                wall_loc_log(c,nc)=.true.
                                                wall_loc_real(:)=0._rp
                                            end if    
                                        case(5)!k+1
                                            if(diz==0)then
                                                dzc_l=dzc(kn)
                                                dzf_l=dzf(kn)
                                            else
                                                dzf_l=dzf(kn+1)
                                                dzc_l=dzc(kn)
                                            endif
                                            if(isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,x,y,zp,n,l,&
                                            hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then 
                                                call calc_lambda(x,y,z,zp,3,lambda,ibm_direction,amp_l,n_wave,l_0,phase_l,n,l,dl,&
                                                hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap,dzf_l,dzc_l,diz,wall_loc_real)
                                                wc=wc+1 !wc describes the wall counter
                                                wall_loc(:,c,nc)=wall_loc_real
                                                wall_loc_log(c,nc)=.true.
                                                wall_loc_real(:)=0._rp 
                                            endif  
                                        case(6)!k-1
                                            if(diz==0)then
                                                dzc_l=dzc(kn-1)
                                                dzf_l=dzf(kn)
                                            else
                                                dzf_l=dzf(kn)
                                                dzc_l=dzc(kn)
                                            endif
                                            if(isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,x,y,zm,n,l,&
                                            hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then 
                                                call calc_lambda(x,y,z,zm,3,lambda,ibm_direction,amp_l,n_wave,l_0,phase_l,n,l,dl,&
                                                hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap,dzf_l,dzc_l,diz,wall_loc_real)
                                                wc=wc+1 !wc describes the wall counter
                                                wall_loc(:,c,nc)=wall_loc_real
                                                wall_loc_log(c,nc)=.true.
                                                wall_loc_real(:)=0._rp
                                            endif
                                        case default
                                            print*,"something is off.... at ibm.f90"   
                                    end select
                                end do
                            endif
                        end do
                        if(wc<=3)cycle
                        count=count+1
                        call get_plane(wall_loc,wall_loc_log,plane)
                        !now we have the normal and center now lets get the distance
                        !lets get the location of the band point
                        x=0._rp;y=0._rp;z=0._rp;xp=0._rp;xm=0._rp;yp=0._rp;ym=0._rp;zp=0._rp;zm=0._rp
                        call get_grid_loc(lo,i,j,k,dl,zc,zf,dix,diy,diz,x,y,z,xp,xm,yp,ym,zp,zm)
                        !output of these are the location of the bandpoint
                        bp(1)=x;bp(2)=y;bp(3)=z;
                        cp(1)=bp(1)-plane(1,1);
                        cp(2)=bp(2)-plane(2,1);
                        cp(3)=bp(3)-plane(3,1);
                        !now we have exact distance from the center point of the plane
                        !to the band point! lets get them scalar product
                        nv=plane(:,2)
                        call comp_sca(nv,cp,angle_plane,out_plane)
                        if(out_plane<0._rp)nv=-nv ! we check if the output is neg.
                        ! this means the normal vector is facing inwards, we want outwards facing
                        call comp_sca(nv,cp,angle_plane,out_plane)
                        grad_dist_id(count,1)=plane(1,1)
                        grad_dist_id(count,2)=plane(2,1)
                        grad_dist_id(count,3)=plane(3,1)
                        grad_dist_id(count,4)=nv(1)
                        grad_dist_id(count,5)=nv(2)
                        grad_dist_id(count,6)=nv(3)
                        grad_dist_id(count,7)=i!these are the locations of the velocities
                        grad_dist_id(count,8)=j!these are the locations of the velocities
                        grad_dist_id(count,9)=k!these are the locations of the velocities
                        grad_dist_id(count,10)=out_plane
                        grad_dist_id(count,11)=0._rp
                        grad_dist_id(count,12)=0._rp
                        ! here we only compute the geometrical distances! 
                        ! since we have a stat. wall we dont need to compute the distances
                        ! in each time step!
                    endif
                end do 
            end do 
        end do
    end subroutine calc_grad_dist
    subroutine calc_shear_st(grad_dist_id,vel_id)
        ! we assume the walls arent moving therefore we do the distance calc.
        ! time and compute the gradient here at each wanted time step
        real(rp),dimension(0:,0:,0:),intent(in)                 :: vel_id
        real(rp),dimension(:,:),intent(inout)                   :: grad_dist_id
        integer                                                 :: i,j,k,kk
        real(rp)                                                :: grad,dist_grad,shear_st
        grad=0._rp;i=0;j=0;k=0;
        !$acc parallel loop private(grad,dist_grad,shear_st,kk,i,j,k)
        do kk=lbound(grad_dist_id,1),ubound(grad_dist_id,1)
            grad=0._rp
            ! now we go through the data
            i=grad_dist_id(kk,7)
            j=grad_dist_id(kk,8)
            k=grad_dist_id(kk,9)
            dist_grad=grad_dist_id(kk,10)
            !we have the location of the i,j,k for the velocity
            call comp_grad(dist_grad,vel_id(i,j,k),grad)
            shear_st=visc*grad
            grad_dist_id(kk,11)=grad
            grad_dist_id(kk,12)=shear_st
        end do 
        
    end subroutine calc_shear_st
    subroutine comp_grad(d,vel,grad)
        !$acc routine seq
        real(rp), intent(in)                      :: vel
        real(rp),intent(in)                       :: d
        real(rp),intent(out)                      :: grad
        if(abs(d)>1e-10)then
            grad=vel*d**(-1)
        else
            grad=0._rp
        endif
    end subroutine comp_grad
    subroutine get_plane(wall_loc,wall_loc_log,plane)
        real(rp),intent(in),dimension(0:2,6,0:26)           :: wall_loc
        logical,intent(in),dimension(6,0:26)                :: wall_loc_log
        integer                                             :: o,p,m,n,co,g,f   
        real(rp),dimension(3)                               :: center
        real(rp),dimension(0:2,6,0:26)                      :: vec
        real(rp),dimension(3)                               :: vec1,vec2,vec3,n1,n2,n3,n_av,n_ref,norm_nav
        real(rp),dimension(0:2,6,0:26)                      :: normal
        logical,dimension(6,0:26)                           :: normal_loc
        real(rp)                                            :: angle1,angle2,angle3
        real(rp),intent(out),dimension(3,2)                 :: plane
        real(rp)                                            :: out_dum,out_n,angle_dum
        integer                                             :: c_n

        !first lets create a center
        angle1=90._rp;angle2=90._rp;angle3=90._rp;
        plane(:,:)=0._rp
        normal_loc(:,:)=.false.
        vec(:,:,:) = 0._rp
        vec1(:) = 0._rp;vec2(:) = 0._rp;vec3(:) = 0._rp;
        n1(:) = 0._rp;n2(:) = 0._rp;n3(:) = 0._rp;
        n_av(:) = 0._rp;n_ref(:) = 0._rp;norm_nav(:)=0._rp;
        center(:)=0._rp
        co=0;c_n=0
        do n=0,26
            do m=1,6
                if(wall_loc_log(m,n))then
                    center(:)=center(:)+wall_loc(:,m,n)
                    co=co+1
                endif
            enddo
        enddo
        if(co>1)then
            center(:)=center(:)/co
        endif
        !now lets get the vectors from the center
        do n=0,26
            do m=1,6
                if(wall_loc_log(m,n))then
                    vec(:,m,n)=wall_loc(:,m,n)-center(:)
                endif
            enddo
        enddo
        ! now we have the vectors lets compute their normals
        ! we do it with a cross product
        !print*, center
        do n=0,26
            do m=1,6
                if(.not.(wall_loc_log(m,n)))cycle
                vec1(:)=vec(:,m,n)
                do g=0,26
                    do f=1,6
                        if(wall_loc_log(f,g))then
                            if(g/=n.or.f/=m)then
                                vec2(:)=vec(:,f,g)
                                call comp_cross(vec1,vec2,n1)
                                if(c_n<1)n_ref=n1
                                ! now we have a normal lets check this one for all vectors available
                                ! lets check their direction then if they are positive lets add up 
                                ! if they are neg. then lets change them (so they look outwards)
                                ! after that simply average so we get an averaged normal
                                if(c_n>1)then
                                    ! here we check the new normal w/ with the first normal
                                    call comp_sca(n1,n_ref,angle_dum,out_n)
                                    if(out_n<0._rp)n1=-n1
                                    ! now we have positive defined normal vector
                                endif
                                c_n=c_n+1 
                                n_av=n_av+n1
                            endif
                        endif
                    enddo
                enddo    
            enddo
        enddo
        n_av=n_av/c_n
        call nomrlalize_vec(n_av,norm_nav)
        plane(:,1)=center(:)
        plane(:,2)=norm_nav(:)
        ! now we have our normals lets check each normal w/ each vector 
    end subroutine get_plane
    subroutine comp_cross(vec1,vec2,vec3)
        real(rp),intent(in),dimension(3) :: vec1,vec2
        real(rp),intent(out),dimension(3) :: vec3
        real(rp)                          :: i,j,k
        real(rp)                          :: n
        k=(vec1(1)*vec2(2))-(vec1(2)*vec2(1))
        j=(vec1(3)*vec2(1))-(vec1(1)*vec2(3))
        i=(vec1(2)*vec2(3))-(vec1(3)*vec2(2))
        n=sqrt(i**2+j**2+k**2)
        if(abs(n)<1e-10)then
            vec3=0._rp
        else
            i=i*n**(-1);j=j*n**(-1);k=k*n**(-1)
            vec3(:)=[i,j,k]
        endif
        
    end subroutine comp_cross
    subroutine nomrlalize_vec(vec1,vec2)
        real(rp),intent(in),dimension(3)     :: vec1 
        real(rp),intent(out),dimension(3)    :: vec2
        real(rp)                             :: i,j,k,n
        n=sqrt(vec1(1)**2+vec1(2)**2+vec1(3)**2)
        if(n>epsilon(1._rp))then
            vec2(1)=vec1(1)*n**(-1)
            vec2(2)=vec1(2)*n**(-1)
            vec2(3)=vec1(3)*n**(-1)
        else
            vec2(1)=0._rp;vec2(2)=0._rp;vec2(3)=0._rp
        endif
    end subroutine nomrlalize_vec
    subroutine comp_sca(vec1,vec2,angle,out)
        real(rp),intent(in),dimension(3) :: vec1,vec2
        real(rp),intent(out)             :: angle
        real(rp),intent(out)             :: out
        real(rp)                         :: b1,b2
        out=vec1(1)*vec2(1)+vec1(2)*vec2(2)+vec1(3)*vec2(3)
        b1=sqrt(vec1(1)**2+vec1(2)**2+vec1(3)**2)
        b2=sqrt(vec2(1)**2+vec2(2)**2+vec2(3)**2)
        if(abs(b1)>1e-10.and.abs(b2)>1e-10)then
            angle=out*(b1*b2)**(-1)
            angle=acos(angle)*180/(3.1415) !radians to degree
        else
            !we give out as epsilon since we will use it to divide
            out=epsilon(0._rp)
        endif
    end subroutine comp_sca
    subroutine get_grid_loc(lo,i,j,k,dl,zc,zf,dix,diy,diz,x,y,z,xp,xm,yp,ym,zp,zm)
        !$acc routine seq
        integer,intent(in)                      :: i,j,k
        real(rp),intent(out)                    :: x,y,z,xp,xm,yp,ym,zp,zm
        integer,intent(in)                      :: dix,diy,diz
        integer                                 :: ii,jj,kk
        integer,intent(in),dimension(3)         :: lo
        real(rp),intent(in),dimension(0:),optional  :: zc,zf
        real(rp), intent(in),dimension(3)       :: dl
        x=0;y=0;z=0;xp=0;xm=0;yp=0;ym=0;zp=0;zm=0;    
        ii = lo(1)+i-1
        jj = lo(2)+j-1
        kk = lo(3)+k-1
        x = (real(ii,rp) -0.5d0+ real(dix,rp)*0.5d0)*dl(1)
        y = (real(jj,rp) -0.5d0+ real(diy,rp)*0.5d0)*dl(2)
        xp=x+dl(1);xm=x-dl(1);yp=y+dl(2);ym=y-dl(2);
        if(diz/=1)then
            ! this means we are looking for either u or v so their location is at z center
            ! we gonna use the senter of zc
            z = zc(k)
            zp = zc(k+1)
            zm = zc(k-1)
            ! we use k inestead of kk since kk is the global and k is the local array index
        else
            ! else than we are looking for the w which is located on the z face
            z = zf(k)
            zp = zf(k+1)
            zm = zf(k-1)
        endif                 
    end subroutine get_grid_loc
    ! here we do the wall pressure calc stuff

    subroutine get_wall_pres(p,p_grad,lo,dl,zc,zf,mask_s,band_s)
        real(rp),intent(in),dimension(0:,0:,0:)     :: p                        
        real(rp),dimension(:,:),intent(inout)       :: p_grad
        integer,intent(in),dimension(3)             :: lo
        real(rp),intent(in),dimension(0:)           :: zc,zf
        real(rp), intent(in),dimension(3)           :: dl
        logical,intent(in)                          :: mask_s(0:,0:,0:)
        logical, intent(in), dimension(0:,0:,0:)    :: band_s
        real(rp),dimension(3)                       :: c0,p_normalv
        integer                                     :: kk,dix,diy,diz 
        integer                                     :: i,j,k,m,in,jn,kn,di,dj,dk,count,gen_count,wr
        real(rp)                                    :: dist,x,y,z,xp,xm,yp,ym,zp,zm
        real(rp),dimension(3)                       :: vec                                    
        real(rp)                                    :: pw,b !b is the the pressure gradient
        real(rp)                                    :: s0,s1,s2,t0,t1
        kk=0;dix=0;diy=0;diz=0;m=0;count=0;gen_count=0;wr=0;
        dist=0._rp;vec(:)=0._rp
        ! now we have the pressures that are on the surface but arent in the solid
        ! lets get their planes etc. we can directly call the grad distance
        ! to create a normal and a plane for each pressure band point 
        ! like we did w/ the velocities
        ! this subroutine directly returns us an array where the pressure band locations that we located
        ! have their normal and a center point 

        !$acc parallel loop private(m,c0,p_normalv,i,j,k,s0,s1,s2,t0,t1) &
        !$acc& private(vec,pw,b,dist,x,y,z,xp,xm,yp,ym,zp,zm) &
        !$acc& private(di,dj,dk,in,jn,kn,count,dix,diy,diz)
        do kk=lbound(p_grad,1),ubound(p_grad,1)
            gen_count=gen_count+1
            c0(1)=p_grad(kk,1)
            c0(2)=p_grad(kk,2)
            c0(3)=p_grad(kk,3)
            p_normalv(1)=p_grad(kk,4)
            p_normalv(2)=p_grad(kk,5)
            p_normalv(3)=p_grad(kk,6)
            i=p_grad(kk,7)
            j=p_grad(kk,8)
            k=p_grad(kk,9)
            ! so lets search the around the 
            s0=0._rp;s1=0._rp;s2=0._rp;t0=0._rp;t1=0._rp;count=0._rp;
            do m=0,26
                di=modulo(m,3)-1
                dj=modulo(m/3,3)-1
                dk=modulo(m/9,3)-1
                in=i+di;jn=j+dj;kn=k+dk
                !so we check if the neighbour of the band poin is in fluid
                if(kn==k.and.jn==j.and.in==i)cycle ! we avoid the band point we are on
                if(in<lbound(p,1).or.in>ubound(p,1)) cycle
                if(jn<lbound(p,2).or.jn>ubound(p,2)) cycle
                if(kn<lbound(p,3).or.kn>ubound(p,3)) cycle 
                if(mask_s(in,jn,kn))cycle!mask_s is only true if the point is in body
                if(band_s(in,jn,kn))cycle!if the location is also a band point cycle
                call get_grid_loc(lo,in,jn,kn,dl,zc,zf,dix,diy,diz,x,y,z,xp,xm,yp,ym,zp,zm)
                ! comp vec
                count=count+1
                vec(1)=x-c0(1);vec(2)=y-c0(2);vec(3)=z-c0(3)
                dist=dot_product(vec,p_normalv)
                !now we have the distance to the center
                s0=count
                s2=s2+dist**2
                s1=s1+dist
                t0=t0+p(in,jn,kn)
                t1=t1+dist*p(in,jn,kn)
            end do
            if(s0<3)cycle
            ! we use direct mse from the gatti's source 
            pw=(s2*t0-s1*t1)*(s0*s2-s1**2)**(-1)
            b=(s0*t1-s1*t0)*(s0*s2-s1**2)**(-1)
            p_grad(kk,11)=pw
            p_grad(kk,12)=b
        end do 
    end subroutine get_wall_pres
    subroutine init_fric_cubes(band_s,mask_s,lo,dl,zc,zf,dzc,&
                                        dzf,xg,yg,dix,diy,diz,ibm_direction,amp_l,&
                                        n_wave,l_0,phase_l,n,l,hmap,l1_hmap,l2_hmap,&
                                        n1_hmap,n2_hmap,output)
        logical, intent(in),dimension(0:,0:,0:)             :: band_s
        logical,intent(in)                                  :: mask_s(0:,0:,0:)
        integer,intent(in),dimension(3)                     :: lo
        real(rp), intent(in),dimension(3)                   :: dl
        real(rp),intent(in),dimension(0:)                   :: zc,zf,xg,yg,dzf,dzc
        logical , intent(in), dimension(0:1,3)              :: ibm_direction
        real(rp), intent(in), dimension(0:1,3)              :: amp_l
        integer , intent(in), dimension(0:1,3)              :: n_wave
        real(rp), intent(in), dimension(0:1,3)              :: l_0
        real(rp), intent(in), dimension(0:1,3)              :: phase_l
        integer , intent(in), dimension(3)                  :: n
        real(rp), intent(in), dimension(3)                  :: l
        real(rp),intent(in),optional                        :: hmap(0:,0:)
        real(rp),intent(in),optional                        :: l1_hmap,l2_hmap 
        integer,intent(in),optional                         :: n1_hmap,n2_hmap
        integer,intent(in)                                  :: dix,diy,diz
        real(rp),intent(inout),dimension(0:,0:,0:,:)        :: output
        integer                                             :: i,j,k
        output(:,:,:,:) = 0._rp;
        do k=lbound(band_s,3)+1,ubound(band_s,3)-1
            do j=lbound(band_s,2)+1,ubound(band_s,2)-1
                do i=lbound(band_s,1)+1,ubound(band_s,1)-1
                    if(band_s(i,j,k))then
                        call get_fric_cubes(band_s,i,j,k,mask_s,lo,dl,zc,zf,dzc,dzf,xg,yg,dix,diy,diz,&
                                            ibm_direction,amp_l,n_wave,l_0,phase_l,n,l,hmap,l1_hmap,&
                                            l2_hmap,n1_hmap,n2_hmap,output)
                    endif
                end do 
            end do
        enddo


    end subroutine init_fric_cubes
    subroutine get_fric_cubes(band_s,i,j,k,mask_s,lo,dl,zc,zf,dzc,dzf,xg,yg,dix,diy,diz,ibm_direction,amp_l,n_wave,l_0,&
        phase_l,n,l,hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap,output)
        logical, intent(in),dimension(0:,0:,0:)             :: band_s
        logical,intent(in)                                  :: mask_s(0:,0:,0:)
        integer,intent(in),dimension(3)                     :: lo
        real(rp), intent(in),dimension(3)                   :: dl
        real(rp),intent(in),dimension(0:)                   :: zc,zf,xg,yg,dzf,dzc
        logical , intent(in), dimension(0:1,3)              :: ibm_direction
        real(rp), intent(in), dimension(0:1,3)              :: amp_l
        integer , intent(in), dimension(0:1,3)              :: n_wave
        real(rp), intent(in), dimension(0:1,3)              :: l_0
        real(rp), intent(in), dimension(0:1,3)              :: phase_l
        integer , intent(in), dimension(3)                  :: n
        real(rp), intent(in), dimension(3)                  :: l
        real(rp),intent(in),optional                        :: hmap(0:,0:)
        real(rp),intent(in),optional                        :: l1_hmap,l2_hmap 
        integer,intent(in),optional                         :: n1_hmap,n2_hmap
        integer,intent(in)                                  :: i,j,k,dix,diy,diz
        integer                                             :: ic,jc,kc
        real(rp)                                            :: x,y,z,xp,xm,yp,ym,zp,zm,x_int,y_int,z_int
        real(rp)                                            :: dx_int,dy_int,dz_int
        integer                                             :: co
        real(rp)                                            :: cell_vol,fluid_vol,solid_vol
        logical,dimension(10,10,10)                         :: mask_subcell
        integer                                             :: s_xp,s_xm,s_yp,s_ym,s_zp,s_zm 
        real(rp)                                            :: A_dx,A_dy,A_dz,A_x,A_y,A_z
        real(rp)                                            :: A_xsi,A_ysi,A_zsi
        real(rp)                                            :: x_in,x_out,y_in,y_out,z_in,z_out
        real(rp),dimension(3)                               :: normal
        real(rp),dimension(2,3)                             :: A_inlets
        real(rp)                                            :: A_tot 
        real(rp),intent(inout),dimension(0:,0:,0:,:)          :: output
        logical :: diag,inloop
        !options
        diag=.false.;inloop=.false.;
        s_xp=0;s_xm=0;s_yp=0;s_ym=0;s_zp=0;s_zm=0;
        A_dx=0._rp;A_dy=0._rp;A_dz=0._rp;
        A_x=0._rp;A_y=0._rp;A_z=0._rp;
        A_xsi=0;A_ysi=0;A_zsi=0;
        if(band_s(i,j,k))then ! we have a band point
            if(i==0.or.j==0.or.k==0)return !guard
            mask_subcell(:,:,:)=.false.
            x=0._rp;y=0._rp;z=0._rp;xp=0._rp;
            xm=0._rp;yp=0._rp;ym=0._rp;zp=0._rp;zm=0._rp;
            x_int=0._rp;y_int=0._rp;z_int=0._rp;
            co = 0;
            call get_grid_loc(lo,i,j,k,dl,zc,zf,dix,diy,diz,x,y,z,xp,xm,yp,ym,zp,zm)
            ! now lets check the inside of the band point with smaller cubes
            dz_int=dzf(k);dy_int=(yg(j)-yg(j-1));dx_int=(xg(i)-xg(i-1))!get dx,dy,dz from global
            A_dx=dy_int*dz_int*0.01_rp;A_dy=dx_int*dz_int*0.01_rp;A_dz=dx_int*dy_int*0.01_rp
            do kc=1,10
                z_int=z-0.5_rp*dz_int+real(kc,kind=rp)/10._rp*dzf(k)! since we start from the cell center 
                do jc=1,10
                    y_int=y-0.5_rp*dy_int+real(jc,kind=rp)/10._rp*dy_int
                    do ic=1,10
                        x_int=x-0.5_rp*dx_int+real(ic,kind=rp)/10._rp*dx_int
                        if(isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,x_int,y_int,z_int,n,l,&
                                    hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then
                            mask_subcell(ic,jc,kc)=.true. ! this subcell is in solid
                            co=co+1
                        endif
                    end do 
                end do 
            end do
            do kc=1,10
                do jc=1,10
                    do ic=1,10
                            if(mask_subcell(ic,jc,kc))then
                                if(ic<10)then
                                    if(.not.mask_subcell(ic+1,jc,kc))then ! neigbour isnt in solid
                                        s_xp=s_xp+1! thus its area is what we look
                                    endif
                                endif
                                if(ic>1)then
                                    if(.not.mask_subcell(ic-1,jc,kc))then
                                        s_xm=s_xm+1
                                    endif
                                endif
                                if(jc<10)then
                                    if(.not.mask_subcell(ic,jc+1,kc))then
                                        s_yp=s_yp+1
                                    endif
                                endif
                                if(jc>1)then
                                    if(.not.mask_subcell(ic,jc-1,kc))then
                                        s_ym=s_ym+1
                                    endif
                                endif
                                if(kc<10)then
                                    if(.not.mask_subcell(ic,jc,kc+1))then
                                        s_zp=s_zp+1
                                    endif
                                endif
                                if(kc>1)then
                                    if(.not.mask_subcell(ic,jc,kc-1))then
                                        s_zm=s_zm+1
                                    endif
                                endif
                            endif                       
                    end do 
                end do 
            end do
            !print*,co
            if(co==0)then ! this means they are band but they sit perfectly
                if(i>lbound(mask_s,1))then
                    if(mask_s(i-1,j,k))then
                        normal=[1._rp,0._rp,0._rp]
                        A_xsi=dy_int*dz_int
                        A_tot=A_xsi
                    endif
                endif
                if(i<ubound(mask_s,1))then
                    if(mask_s(i+1,j,k))then
                        normal=[-1._rp,0._rp,0._rp]
                        A_xsi=dy_int*dz_int
                        A_tot=A_xsi
                    endif
                endif
                if(j>lbound(mask_s,2))then
                    if(mask_s(i,j-1,k))then
                        normal=[0._rp,1._rp,0._rp]
                        A_ysi=dx_int*dz_int
                        A_tot=A_ysi
                    endif
                endif
                if(j<ubound(mask_s,2))then
                    if(mask_s(i,j+1,k))then
                        normal=[0._rp,-1._rp,0._rp]
                        A_ysi=dx_int*dz_int
                        A_tot=A_ysi
                    endif
                endif
                if(k>lbound(mask_s,3))then
                    if(mask_s(i,j,k-1))then
                        normal=[0._rp,0._rp,1._rp]
                        A_zsi=dx_int*dy_int
                        A_tot=A_zsi
                    endif
                endif
                if(k<ubound(mask_s,3))then
                    if(mask_s(i,j,k+1))then
                        normal=[0._rp,0._rp,-1._rp]
                        A_zsi=dx_int*dy_int
                        A_tot=A_zsi
                    endif
                endif
                ! disable this for flat wall... 
                ! otherwise overestimates for curvy structure
            endif
            !lets get the inlet values
            x_in=count(.not.mask_subcell(1,:,:))*A_dx
            x_out=count(.not.mask_subcell(10,:,:))*A_dx
            y_in=count(.not.mask_subcell(:,1,:))*A_dy
            y_out=count(.not.mask_subcell(:,10,:))*A_dy
            z_in=count(.not.mask_subcell(:,:,1))*A_dz
            z_out=count(.not.mask_subcell(:,:,10))*A_dz
            ! surface area
            if(co>0)then
                A_x=(s_xp+s_xm)*A_dx;A_y=(s_yp+s_ym)*A_dy;A_z=(s_zp+s_zm)*A_dz;
                A_xsi=(s_xp-s_xm)*A_dx;A_ysi=(s_yp-s_ym)*A_dy;A_zsi=(s_zp-s_zm)*A_dz;
                A_tot=sqrt(A_xsi**2+A_ysi**2+A_zsi**2)
            endif
            if(co>0)call nomrlalize_vec([A_xsi,A_ysi,A_zsi],normal) ! we already fill the otherones before
            ! now lets calc the vol.
            cell_vol=dx_int*dy_int*dz_int 
            solid_vol=(real(co,rp)/1000._rp)*cell_vol
            fluid_vol=(1*cell_vol)-solid_vol
            if(inloop)then
                A_inlets(1,1)=x_in;A_inlets(2,1)=x_out;
                A_inlets(1,2)=y_in;A_inlets(2,2)=y_out;
                A_inlets(1,3)=z_in;A_inlets(2,3)=z_out;
            endif
            if(co==0)then
                A_tot     = 0._rp
                solid_vol = 0._rp
                fluid_vol = cell_vol
                normal    = 0._rp
            endif
            if(diag)then 
                print*,"cell volume is:",cell_vol
                print*,"cell fluid volume is:",fluid_vol
                print*,"cell solid volume is:",solid_vol
                print*,"total surface area in x:",A_x
                print*,"total surface area in y:",A_y
                print*,"total surface area in z:",A_z
                print*,"signed surface area in x:",A_xsi
                print*,"signed surface area in y:",A_ysi
                print*,"signed surface area in z:",A_zsi
                print*,"approx. surface area:",A_tot
                print*,"inlet area in x:",x_in
                print*,"outlet area in x:",x_out
                print*,"inlet area in y:",y_in
                print*,"outlet area in y:",y_out
                print*,"inlet area in z:",z_in
                print*,"outlet area in z:",z_out
                print*,"normal is x,y,z:",normal
            endif
            output(i,j,k,1)=cell_vol
            output(i,j,k,2)=fluid_vol
            output(i,j,k,3)=solid_vol
            output(i,j,k,4)=A_x
            output(i,j,k,5)=A_y
            output(i,j,k,6)=A_z
            output(i,j,k,7)=A_xsi
            output(i,j,k,8)=A_ysi
            output(i,j,k,9)=A_zsi
            output(i,j,k,10)=A_tot
            output(i,j,k,11)=x_in
            output(i,j,k,12)=x_out
            output(i,j,k,13)=y_in
            output(i,j,k,14)=y_out
            output(i,j,k,15)=z_in
            output(i,j,k,16)=z_out
            output(i,j,k,17)=normal(1)
            output(i,j,k,18)=normal(2)
            output(i,j,k,19)=normal(3)            
        endif
    end subroutine get_fric_cubes
    subroutine calc_fric_cubes(lo,u,v,w,p,dxi,dyi,dzci,visc,band_s,mom_saved,bforce,dti,cubeInput,output)
        real(rp),intent(in),dimension(0:,0:,0:)             :: u,v,w
        integer,intent(in),dimension(3)                     :: lo
        real(rp),intent(in),dimension(0:)                   :: dzci
        real(rp),intent(in)                                 :: dxi,dyi
        real(rp),intent(in)                                 :: visc,dti 
        logical, intent(in),dimension(0:,0:,0:)             :: band_s
        real(rp),intent(in),dimension(0:,0:,0:)             :: p 
        real(rp),intent(in)                                 :: bforce
        real(rp),intent(in),dimension(0:,0:,0:)             :: mom_saved
        real(rp)                                            :: A_tot,cell_vol,fluid_vol,solid_vol
        real(rp),dimension(3)                               :: normal
        real(rp)                                            :: conv_x_in,conv_x_out,p_x_in,p_x_out
        real(rp)                                            :: diff_x_in,diff_x_out
        real(rp)                                            :: F_wall,sum_fwall,sum_awall,real_fwall,sum_fvol
        real(rp),dimension(3)                               :: F_wall_c
        real(rp)                                            :: uuip,uuim,vujp,vujm,wukp,wukm
        real(rp)                                            :: dudxp,dudxm,dudyp,dudym,dudzp,dudzm
        integer                                             :: count
        real(rp)                                            :: mom_in,mom_out,shear_st,mom
        integer                                             :: i,j,k
        real(rp),intent(in),dimension(0:,0:,0:,:)        :: cubeInput
        real(rp)                                            :: xin,xout,yin,yout,zin,zout
        real(rp),intent(out),dimension(:,:)                 :: output
        logical :: diag
        integer :: ig,jg,kg
        diag=.false.
        sum_fwall = 0._rp;sum_awall = 0._rp;count = 0;sum_fvol=0._rp;shear_st=0._rp;
        do k=lbound(band_s,3)+1,ubound(band_s,3)-1
            do j=lbound(band_s,2)+1,ubound(band_s,2)-1
                do i=lbound(band_s,1)+1,ubound(band_s,1)-1
                    F_wall=0._rp;A_tot=0._rp;
                    if(band_s(i,j,k))then
                        ig=0;jg=0;kg=0;
                        conv_x_in=0._rp;conv_x_out=0._rp;p_x_in=0._rp;p_x_out=0._rp;
                        diff_x_in=0._rp;diff_x_out=0._rp;
                        ig=lo(1)+i-1;jg=lo(2)+j-1;kg=lo(3)+k-1;
                        fluid_vol=cubeInput(i,j,k,2);A_tot=cubeInput(i,j,k,10);
                        xin=cubeInput(i,j,k,11);xout=cubeInput(i,j,k,12);
                        yin=cubeInput(i,j,k,13);yout=cubeInput(i,j,k,14);
                        zin=cubeInput(i,j,k,15);zout=cubeInput(i,j,k,16);
                        if(A_tot<=tiny(1._rp))cycle
                        count=count+1
                        call get_momx_a(i,j,k,u,v,w,&
                                        uuip,uuim,vujp,vujm,wukp,wukm)
                        call get_momx_d(dxi,dyi,i,j,k,dzci,visc,u,&
                                        dudxp,dudxm,dudyp,dudym,dudzp,dudzm)  
                        p_x_in=p(i,j,k);p_x_out=p(i+1,j,k)
                        mom=u(i,j,k)-mom_saved(i,j,k)
                        F_wall= uuip*xout-uuim*xin&
                               +vujp*yout-vujm*yin&
                               +wukp*zout-wukm*zin&
                               +p_x_out*xout-p_x_in*xin&
                               -dudxp*xout+dudxm*xin&
                               -dudyp*yout+dudym*yin&
                               -dudzp*zout+dudzm*zin&
                               -bforce*fluid_vol&
                               +mom*fluid_vol*dti                                      
                        ! now we have a force lets get apply the normal to it
                        sum_fwall=sum_fwall+F_wall
                        sum_awall=sum_awall+A_tot
                        sum_fvol=sum_fvol+fluid_vol
                        shear_st=F_wall/A_tot
                        output(count,1)=real(ig,kind=rp)
                        output(count,2)=real(jg,kind=rp)
                        output(count,3)=real(kg,kind=rp)
                        output(count,4)=A_tot
                        output(count,5)=F_wall
                        output(count,6)=shear_st
                        if(diag)then
                            !conv
                            print*,"Convection fluxex in +x: ",uuip*xout
                            print*,"Convection fluxex in -x: ",uuim*xin
                            print*,"Convection fluxex in +y: ",vujp*yout
                            print*,"Convection fluxex in -y: ",vujm*yin
                            print*,"Convection fluxex in +z: ",wukp*zout
                            print*,"Convection fluxex in -z: ",wukm*zin
                            !diff
                            print*,"Diffusion fluxes in +x: ",dudxp*xout
                            print*,"Diffusion fluxes in -x: ",dudxm*xin
                            print*,"Diffusion fluxes in +y: ",dudyp*yout
                            print*,"Diffusion fluxes in -y: ",dudym*yin
                            print*,"Diffusion fluxes in +z: ",dudzp*zout
                            print*,"Diffusion fluxes in -z: ",dudzm*zin
                            !p
                            print*,"Pressre in fluxes +x: ",p_x_in*xout
                            print*,"Pressre in fluxes -x: ",p_x_out*xin
                            !bforce and mom cv
                            print*,"Body force flux",bforce*fluid_vol
                            print*,"Mometum fluxes: ",mom*fluid_vol*dti
                            !
                            print*,"Wall force in x dir: ",F_wall
                            print*,"shear stress for a cell:",shear_st
                        endif
                    endif
                end do 
            end do
        enddo
    end subroutine calc_fric_cubes
    !
    subroutine writeFVfric(fname,output)
        implicit none
        real(rp),dimension(:,:),intent(in)       :: output
        character(len=*)                         :: fname   
        integer                                  :: iunit,wr
        character(len=*), parameter              :: fmt_dp = '(*(es24.16e3,1x))', &
                                                    fmt_sp = '(*(es15.8e2,1x))'
#if !defined(_SINGLE_PRECISION)
        character(len=*), parameter              :: fmt_rp = fmt_dp
#else
        character(len=*), parameter              :: fmt_rp = fmt_sp
#endif
        open(newunit=iunit,file=fname)
            do wr=1,ubound(output,1)
                write(iunit,fmt_rp) output(wr,1),output(wr,2),output(wr,3),&
                                    output(wr,4),output(wr,5),output(wr,6)
            end do
        close(iunit)
    end subroutine writeFVfric
    subroutine save_mom(vel_id,band_id,mom_saved)
        implicit none
        real(rp),intent(in),dimension(0:,0:,0:)                 :: vel_id
        logical, intent(in),dimension(0:,0:,0:)                 :: band_id
        real(rp),intent(inout),dimension(0:,0:,0:)                :: mom_saved
        integer :: i,j,k
        !$acc parallel loop collapse(3) default(present) async(1)
        !$omp parallel do collapse(3) default(shared)
        do k=lbound(band_id,3)+1,ubound(band_id,3)-1
            do j=lbound(band_id,2)+1,ubound(band_id,2)-1
                do i=lbound(band_id,1)+1,ubound(band_id,1)-1
                    if(band_id(i,j,k))then
                        mom_saved(i,j,k)=vel_id(i,j,k)
                    endif
                end do 
            end do
        enddo 
    end subroutine save_mom
    !
    subroutine get_momx_a(i,j,k,u,v,w,uuip,uuim,vujp,vujm,wukp,wukm)
    implicit none
    real(rp), dimension(0:,0:,0:), intent(in   ) :: u,v,w
    integer :: i,j,k
    real(rp),intent(out) :: uuip,uuim,vujp,vujm,wukp,wukm

    uuip = 0.25*( u(i,j  ,k  )+u(i+1,j  ,k  ) )*( u(i,j,k)+u(i+1,j,k) )
    uuim = 0.25*( u(i,j  ,k  )+u(i-1,j  ,k  ) )*( u(i,j,k)+u(i-1,j,k) )
    vujp = 0.25*( v(i,j  ,k  )+v(i+1,j  ,k  ) )*( u(i,j,k)+u(i,j+1,k) )
    vujm = 0.25*( v(i,j-1,k  )+v(i+1,j-1,k  ) )*( u(i,j,k)+u(i,j-1,k) )
    wukp = 0.25*( w(i,j  ,k  )+w(i+1,j  ,k  ) )*( u(i,j,k)+u(i,j,k+1) )
    wukm = 0.25*( w(i,j  ,k-1)+w(i+1,j  ,k-1) )*( u(i,j,k)+u(i,j,k-1) )
  

    end subroutine get_momx_a
    !
    subroutine get_momx_d(dxi,dyi,i,j,k,dzci,visc,u,dudxp,dudxm,dudyp,dudym,dudzp,dudzm)
        implicit none
        real(rp), intent(in) :: dxi,dyi
        real(rp), intent(in), dimension(0:) :: dzci
        real(rp), intent(in) :: visc
        real(rp), dimension(0:,0:,0:), intent(in   ) :: u
        real(rp),intent(out) :: dudxp,dudxm,dudyp,dudym,dudzp,dudzm
        integer,intent(in) :: i,j,k
        !
        dudxp = (u(i+1,j,k)-u(i,j,k))*visc*dxi
        dudxm = (u(i,j,k)-u(i-1,j,k))*visc*dxi
        dudyp = (u(i,j+1,k)-u(i,j,k))*visc*dyi
        dudym = (u(i,j,k)-u(i,j-1,k))*visc*dyi
        dudzp = (u(i,j,k+1)-u(i,j,k))*visc*dzci(k)
        dudzm = (u(i,j,k)-u(i,j,k-1))*visc*dzci(k-1)

 
    end subroutine get_momx_d
    !
    subroutine write_data(fname,myid,grad_id)
        real(rp),dimension(:,:),intent(in)       :: grad_id
        integer,intent(in)                       :: myid
        character(len=*)                         :: fname   
        integer                                  :: iunit,wr
        character(len=*), parameter              :: fmt_dp = '(*(es24.16e3,1x))', &
                                                    fmt_sp = '(*(es15.8e2,1x))'
#if !defined(_SINGLE_PRECISION)
        character(len=*), parameter              :: fmt_rp = fmt_dp
#else
        character(len=*), parameter              :: fmt_rp = fmt_sp
#endif
        if(myid == 0) then
        open(newunit=iunit,file=fname)
            do wr=1,ubound(grad_id,1)
                write(iunit,fmt_rp) grad_id(wr,1),grad_id(wr,2),grad_id(wr,3),&
                                    grad_id(wr,4),grad_id(wr,5),grad_id(wr,6),&
                                    grad_id(wr,7),grad_id(wr,8),grad_id(wr,9),&
                                    grad_id(wr,10),grad_id(wr,11),grad_id(wr,12)
            end do
        close(iunit)
        end if
    end subroutine write_data
    subroutine apply_ibm_staircase(field,mask_id,dt)
        implicit none
        real(rp),intent(inout),dimension(0:,0:,0:)  :: field
        logical,intent(in),dimension(0:,0:,0:)      :: mask_id
        real(rp),intent(in)                         :: dt
        integer :: i,j,k
        !$acc parallel loop collapse(3) default(present) async(1)
        !$OMP parallel do   collapse(3) DEFAULT(shared)
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
    subroutine apply_ibm_staircase_scalar(scalar_field,scalar_mask_id,scalar_bc)
        implicit none
        real(rp),intent(inout),dimension(0:,0:,0:)  :: scalar_field
        logical,intent(in),dimension(0:,0:,0:)      :: scalar_mask_id
        real(rp),intent(in)                         :: scalar_bc
        integer :: i,j,k
        !$acc parallel loop collapse(3) default(present) async(1)
        !$OMP parallel do   collapse(3) DEFAULT(shared)
        do k = lbound(scalar_field,3)+1,ubound(scalar_field,3)-1
            do j = lbound(scalar_field,2)+1,ubound(scalar_field,2)-1
                do i = lbound(scalar_field,1)+1,ubound(scalar_field,1)-1
                    if (scalar_mask_id(i,j,k).eqv..true.)then
                        scalar_field(i,j,k) = scalar_bc 
                    endif
                end do 
            end do 
        end do
    end subroutine apply_ibm_staircase_scalar
end module mod_ibm