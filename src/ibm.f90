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
        xyz = [x,y,z]
        do side = 0,1
            do t = 1,3
                i=modulo(t,3)+1
                ii=modulo(t+1,3)+1
                if(ibm_direction(side,t))then
                    ! use wave wall  
                    if(.not. use_hmap)then
                        height(side,t)=amp_l(side,i)*0.5_rp*(1._rp+sin(2._rp*pi*&
                                        real(n_wave(side,i)*xyz(i)/l(i),rp)+phase_l(side,i)))+&
                                        amp_l(side,ii)*0.5_rp*(1._rp+sin(2._rp*pi*&
                                        real(n_wave(side,ii)*xyz(ii)/l(ii),rp)+phase_l(side,ii)))
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
        ibm_diagnostic=.true.
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
                    ! we create the real location of each velocity here
                    x = (real(ii,rp) -0.5d0+ real(dix,rp)*0.5d0)*dl(1)
                    y = (real(jj,rp) -0.5d0+ real(diy,rp)*0.5d0)*dl(2)
                    ! this part is due to grid stretching
                    ! and ofc staggered grid 
                    ! diz indicates if we are looking for w cells if we dont they are in the center
                    ! if we look into w cells they are in the face
                    if(diz/=1)then
                        z = zc(k)
                        zp = zc(k+1)
                        zm = zc(k-1)
                    else
                        z = zf(k)
                        zp = zf(k+1)
                        zm = zf(k-1)
                    endif
                    dzf_l=0._rp;dzc_l=0._rp
                    xp=x+dl(1);xm=x-dl(1);yp=y+dl(2);ym=y-dl(2);
                    do n_dir=1,6
                        select case(n_dir)
                            case(1)
                                ! xp
                                if(.not.mask_id(i,j,k).and.isInbody(ibm_direction,amp_l,n_wave,l_0,phase_l,xp,y,z,n,l,&
                                hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap))then
                                            call calc_lambda(x,y,z,xp,1,lambda,ibm_direction,amp_l,n_wave,l_0,&
                                                            phase_l,n,l,dl,hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap)
                                            laplacian_id(i,j,k)=laplacian_id(i,j,k)+lambda
                                            band_id(i,j,k) = .true. ! this means on band and fluid
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
                                band_id,visc,vel_id,dix,diy,diz,hmap,l1_hmap,l2_hmap,n1_hmap,n2_hmap,dzf,dzc)
        ! routine creates the tangent and normal plane for 3D IBM
        real(rp),dimension(:,:),intent(inout)                   :: grad_dist_id
        real(rp), dimension(0:,0:,0:), intent(in   )            :: vel_id
        real(rp), intent(in)                                    :: visc                        
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
                        if(wc<=3)then
                            ! if we have not enough points cycle
                            cycle
                        endif
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
                        grad_dist_id(count,4)=i!these are the locations of the velocities
                        grad_dist_id(count,5)=j!these are the locations of the velocities
                        grad_dist_id(count,6)=k!these are the locations of the velocities
                        grad_dist_id(count,7)=out_plane
                        grad_dist_id(count,8)=0._rp
                        grad_dist_id(count,9)=0._rp
                        ! here we only compute the geometrical distances! 
                        ! since we have a stat. wall we dont need to compute the distances
                        ! in each time step!
                    endif
                end do 
            end do 
        end do
    end subroutine calc_grad_dist
    subroutine calc_shear_st(fname,myid,grad_dist_id,vel_id,band_id)
        ! we assume the walls arent moving therefore we do the distance calc.
        ! time and compute the gradient here at each wanted time step
        logical, intent(in), dimension(0:,0:,0:)                :: band_id
        integer,intent(in)                                      :: myid
        character(len=*), intent(in)                            :: fname
        real(rp),dimension(0:,0:,0:),intent(in)                 :: vel_id
        real(rp),dimension(:,:),intent(inout)                   :: grad_dist_id
        integer                                                 :: i,j,k,count,kk
        real(rp)                                                :: grad,dist_grad,shear_st
        integer                                                 :: iunit
        character(len=*), parameter :: fmt_dp = '(*(es24.16e3,1x))', &
                                 fmt_sp = '(*(es15.8e2,1x))'
#if !defined(_SINGLE_PRECISION)
        character(len=*), parameter :: fmt_rp = fmt_dp
#else
        character(len=*), parameter :: fmt_rp = fmt_sp
#endif
        grad=0._rp;i=0;j=0;k=0;count=0
        do kk=lbound(grad_dist_id,1),ubound(grad_dist_id,1)
            count=count+1
            grad=0._rp
            ! now we go through the data
            i=grad_dist_id(kk,4)
            j=grad_dist_id(kk,5)
            k=grad_dist_id(kk,6)
            dist_grad=grad_dist_id(kk,7)
            !we have the location of the i,j,k for the velocity
            call comp_grad(dist_grad,vel_id(i,j,k),grad)
            shear_st=visc*grad
            grad_dist_id(kk,8)=grad
            grad_dist_id(kk,9)=shear_st
        end do 
        
        if(myid == 0) then
        open(newunit=iunit,file=fname)
            do k=1,count
                write(iunit,fmt_rp) grad_dist_id(k,1),grad_dist_id(k,2),grad_dist_id(k,3),&
                                grad_dist_id(k,4),grad_dist_id(k,5),grad_dist_id(k,6),&
                                grad_dist_id(k,7),grad_dist_id(k,8),grad_dist_id(k,9)
            end do
        close(iunit)
        end if
    end subroutine calc_shear_st

    subroutine comp_grad(d,vel,grad)
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
        real(rp),dimension(3)                               :: vec1,vec2,vec3,n1,n2,n3,n_av,n_ref
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
        n_av(:) = 0._rp;n_ref(:) = 0._rp;
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
        plane(:,1)=center(:)
        plane(:,2)=n_av(:)
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
end module mod_ibm