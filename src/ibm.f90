module mod_ibm
    use :: mod_types 
    use :: mod_param
    implicit none
    type :: ibm_type
        real(rp) ::  phase_x,phase_y
        real(rp) ::  amp_x,amp_y  
        real(rp) :: z0     
        real(rp)::solid = 1.0d30   
        real(rp),allocatable:: mask_u(:,:,:),mask_v(:,:,:),&
                                mask_w(:,:,:) 
        integer  :: n_wave_x, n_wave_y
    end type
    contains 
    ! we initalize the ibm coef. here
    ! at the main.f90
    subroutine init_ibm(ibm,l,n)
        implicit none
        type(ibm_type),intent(inout)                :: ibm
        real(rp), intent(in   ), dimension(3)   :: l
        integer , intent(in   ), dimension(3)       :: n
        ibm%amp_x = 0.15_rp*l(3)
        ibm%amp_y = 0.1_rp*l(3)
        ibm%z0 = 0.1_rp*l(3)
        ibm%phase_x = 0._rp
        ibm%phase_y = 0._rp
        ibm%n_wave_x = 2
        ibm%n_wave_y = 1
        allocate(ibm%mask_u(0:n(1)+1,0:n(2)+1,0:n(3)+1))
        allocate(ibm%mask_v(0:n(1)+1,0:n(2)+1,0:n(3)+1))
        allocate(ibm%mask_w(0:n(1)+1,0:n(2)+1,0:n(3)+1))
        ! isInbody ---> false 
        ! so they start w/ all fluid
        print*, "***IBM coefficients are initializing***"
        ibm%mask_u = 0._rp
        ibm%mask_v = 0._rp
        ibm%mask_w = 0._rp
    end subroutine init_ibm
    ! check if the real location(x,y,z) is in the given body shape 
    logical function isInbody(ibm,x,y,z,n,l)
        implicit none
        type(ibm_type),intent(inout) :: ibm
        integer , intent(in   ), dimension(3)   :: n
        real(rp), intent(in   ), dimension(3)   :: l
        real(rp),intent(in)                     :: x,y,z 
        real(rp),parameter                      :: pi = 3.141592653589793_rp
        real(rp)                                :: z_body
        z_body = ibm%z0+ibm%amp_x * 0.5d0 * &
                 (1.0d0 + sin(2.0d0*pi*real(ibm%n_wave_x,rp)*x/l(1) + ibm%phase_x))  + &
                 ibm%amp_y * 0.5d0 * &
                (1.0d0 + sin(2.0d0*pi*real(ibm%n_wave_y,rp)*y/l(2) + ibm%phase_y))
        isInBody = .false.
        if(z<z_body)then
            isInBody = .true.
        endif 
    end function isInbody
    ! 1st order IBM 
    ! we change the diL depending the velocity mask we are handling 
    ! e.g. we need to apply 1,0,0 for mask_u and we need to apply 0,1,0 for mask_v 
    ! 0,0,1 for mask_w
    subroutine set_ibm_staircase(lo,ibm,mask_id,dix,diy,diz,n,l,dl)
        implicit none
        type(ibm_type),intent(inout)                :: ibm
        real(rp), intent(in   ), dimension(3)       :: l
        real(rp), intent(in   ), dimension(3)       :: dl
        integer , intent(in   ), dimension(3)       :: n
        integer , intent(in   ), dimension(3)       :: lo
        real(rp),intent(inout),dimension(0:,0:,0:)   :: mask_id
        integer,intent(in)                          :: dix,diy,diz 
        integer                                     :: i,j,k
        integer                                     :: ii,jj,kk
        real(rp)                                    :: x,y,z
        print*, "***IBM coefficients are deploying***"
        do k = lbound(mask_id,3),ubound(mask_id,3)
            do j = lbound(mask_id,2),ubound(mask_id,2)
                do i = lbound(mask_id,1),ubound(mask_id,1)
                    ii = lo(1)+i-1
                    jj = lo(2)+j-1
                    kk = lo(3)+k-1
                    ! we create the real location of each velocity here
                    x = (real(ii,rp) - real(dix,rp)*0.5d0)*dl(1)
                    y = (real(jj,rp) - real(diy,rp)*0.5d0)*dl(2)
                    z = (real(kk,rp) - real(diz,rp)*0.5d0)*dl(3)
                    if(isInbody(ibm,x,y,z,n,l).eqv..true.)then
                        mask_id(i,j,k) = ibm%solid
                    endif 
                end do 
            end do 
        end do
    end subroutine set_ibm_staircase
    subroutine apply_ibm_staircase(ibm,field,mask_id,dt)
        implicit none
        type(ibm_type)                              :: ibm
        real(rp),intent(inout),dimension(0:,0:,0:)  :: field
        real(rp),intent(in),dimension(0:,0:,0:)     :: mask_id
        real(rp),intent(in)                         :: dt
        integer :: i,j,k
        do k = lbound(field,3),ubound(field,3)
            do j = lbound(field,2),ubound(field,2)
                do i = lbound(field,1),ubound(field,1)
                    if (mask_id(i,j,k)>=.5_rp*ibm%solid)then
                        field(i,j,k) = 0._rp 
                        !field(i,j,k)/(1._rp+dt*mask_id(i,j,k))
                    endif
                end do 
            end do 
        end do
    end subroutine apply_ibm_staircase
end module mod_ibm