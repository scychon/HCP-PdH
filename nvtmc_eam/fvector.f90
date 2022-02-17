!****************************************************************
! 
! this module contains utility functions for vector calculation
!
!****************************************************************

MODUlE fvector
    implicit none
 
contains

    ! Subroutines to expand the size of allocatable arrays
    subroutine expand1D ( array, nsize )
        implicit none
        integer, intent(in) :: nsize
        real*8, intent(inout), allocatable, dimension(:) :: array
        real*8, allocatable, dimension(:) :: temp
        integer :: isize
        
        isize = size(array(:))
        allocate(temp(isize+nsize))
        temp=0
        temp(:isize)=array
        deallocate(array)
        allocate(array(isize+nsize))
        array=temp
        deallocate(temp)
    end subroutine expand1D

    subroutine expand1Dint ( array, nsize )
        implicit none
        integer, intent(in) :: nsize
        integer, intent(inout), allocatable, dimension(:) :: array
        integer, allocatable, dimension(:) :: temp
        integer :: isize
        
        isize = size(array(:))
        allocate(temp(isize+nsize))
        temp=0
        temp(:isize)=array
        deallocate(array)
        allocate(array(isize+nsize))
        array=temp
        deallocate(temp)
    end subroutine expand1Dint

    subroutine expand2D ( array, nsize1, nsize2 )
        implicit none
        integer, intent(in) :: nsize1, nsize2
        real*8, intent(inout), allocatable, dimension(:,:) :: array
        real*8, allocatable, dimension(:,:) :: temp
        integer :: isize1, isize2
        
        isize1 = size(array(:,1))
        isize2 = size(array(1,:))
        allocate(temp(isize1+nsize1,isize2+nsize2))
        temp=0
        temp(:isize1,:isize2)=array
        deallocate(array)
        allocate(array(isize1+nsize1,isize2+nsize2))
        array=temp
        deallocate(temp)
    end subroutine expand2D

    subroutine expand2Dint ( array, nsize1, nsize2 )
        implicit none
        integer, intent(in) :: nsize1, nsize2
        integer, intent(inout), allocatable, dimension(:,:) :: array
        integer, allocatable, dimension(:,:) :: temp
        integer :: isize1, isize2
        
        isize1 = size(array(:,1))
        isize2 = size(array(1,:))
        allocate(temp(isize1+nsize1,isize2+nsize2))
        temp=0
        temp(:isize1,:isize2)=array
        deallocate(array)
        allocate(array(isize1+nsize1,isize2+nsize2))
        array=temp
        deallocate(temp)
    end subroutine expand2Dint

    subroutine expand3D ( array, nsize1, nsize2, nsize3 )
        implicit none
        integer, intent(in) :: nsize1, nsize2, nsize3
        real*8, intent(inout), allocatable, dimension(:,:,:) :: array
        real*8, allocatable, dimension(:,:,:) :: temp
        integer :: isize1, isize2, isize3
        
        isize1 = size(array(:,1,1))
        isize2 = size(array(1,:,1))
        isize3 = size(array(1,1,:))
        allocate(temp(isize1+nsize1,isize2+nsize2,isize3+nsize3))
        temp=0
        temp(:isize1,:isize2,:isize3)=array
        deallocate(array)
        allocate(array(isize1+nsize1,isize2+nsize2,isize3+nsize3))
        array=temp
        deallocate(temp)
    end subroutine expand3D

    subroutine expand4Dint ( array, nsize1, nsize2, nsize3, nsize4 )
        implicit none
        integer, intent(in) :: nsize1, nsize2, nsize3, nsize4
        integer, intent(inout), allocatable, dimension(:,:,:,:) :: array
        integer, allocatable, dimension(:,:,:,:) :: temp
        integer :: isize1, isize2, isize3, isize4
        
        isize1 = size(array(:,1,1,1))
        isize2 = size(array(1,:,1,1))
        isize3 = size(array(1,1,:,1))
        isize4 = size(array(1,1,1,:))
        allocate(temp(isize1+nsize1,isize2+nsize2,isize3+nsize3,isize4+nsize4))
        temp=-1
        temp(:isize1,:isize2,:isize3,:isize4)=array
        deallocate(array)
        allocate(array(isize1+nsize1,isize2+nsize2,isize3+nsize3,isize4+nsize4))
        array=temp
        deallocate(temp)
    end subroutine expand4Dint


    ! Subroutines to schrink the size of allocatable arrays to fit the dimension
    subroutine shrink1D ( array, nsize )
        implicit none
        integer, intent(in) :: nsize
        real*8, intent(inout), allocatable, dimension(:) :: array
        real*8, allocatable, dimension(:) :: temp
        
        allocate(temp(nsize))
        temp=0
        temp(:)=array(:nsize)
        deallocate(array)
        allocate(array(nsize))
        array=temp
        deallocate(temp)
    end subroutine shrink1D

    subroutine shrink1Dint ( array, nsize )
        implicit none
        integer, intent(in) :: nsize
        integer, intent(inout), allocatable, dimension(:) :: array
        integer, allocatable, dimension(:) :: temp
        
        allocate(temp(nsize))
        temp=0
        temp(:)=array(:nsize)
        deallocate(array)
        allocate(array(nsize))
        array=temp
        deallocate(temp)
    end subroutine shrink1Dint

    subroutine shrink2D ( array, nsize1, nsize2 )
        implicit none
        integer, intent(in) :: nsize1, nsize2
        real*8, intent(inout), allocatable, dimension(:,:) :: array
        real*8, allocatable, dimension(:,:) :: temp
        
        allocate(temp(nsize1,nsize2))
        temp=array(:nsize1,:nsize2)
        deallocate(array)
        allocate(array(nsize1,nsize2))
        array=temp
        deallocate(temp)
    end subroutine shrink2D

    subroutine shrink2Dint ( array, nsize1, nsize2 )
        implicit none
        integer, intent(in) :: nsize1, nsize2
        integer, intent(inout), allocatable, dimension(:,:) :: array
        integer, allocatable, dimension(:,:) :: temp
        
        allocate(temp(nsize1,nsize2))
        temp=array(:nsize1,:nsize2)
        deallocate(array)
        allocate(array(nsize1,nsize2))
        array=temp
        deallocate(temp)
    end subroutine shrink2Dint

    subroutine shrink3D ( array, nsize1, nsize2, nsize3 )
        implicit none
        integer, intent(in) :: nsize1, nsize2, nsize3
        real*8, intent(inout), allocatable, dimension(:,:,:) :: array
        real*8, allocatable, dimension(:,:,:) :: temp
        
        allocate(temp(nsize1,nsize2,nsize3))
        temp=array(:nsize1,:nsize2,:nsize3)
        deallocate(array)
        allocate(array(nsize1,nsize2,nsize3))
        array=temp
        deallocate(temp)
    end subroutine shrink3D

    subroutine shrink4Dint ( array, nsize1, nsize2, nsize3, nsize4 )
        implicit none
        integer, intent(in) :: nsize1, nsize2, nsize3, nsize4
        integer, intent(inout), allocatable, dimension(:,:,:,:) :: array
        integer, allocatable, dimension(:,:,:,:) :: temp
        
        allocate(temp(nsize1,nsize2,nsize3,nsize4))
        temp=array(:nsize1,:nsize2,:nsize3,:nsize4)
        deallocate(array)
        allocate(array(nsize1,nsize2,nsize3,nsize4))
        array=temp
        deallocate(temp)
    end subroutine shrink4Dint

    function vecAdd( vec1, vec2 )
        implicit none
        real*8, intent(in), dimension(3) :: vec1, vec2
        real*8, dimension(3) :: vecAdd
        vecAdd(:) = vec1(:) - vec2(:)
    end function vecAdd

    function vecSquare( vec )
        implicit none
        real*8, intent(in), dimension(3) :: vec
        real*8, dimension(3) :: vecSquare
        vecSquare(1) = vec(1)*vec(1)
        vecSquare(2) = vec(2)*vec(2)
        vecSquare(3) = vec(3)*vec(3)
    end function vecSquare

    function vecCross( vec1, vec2 )
        implicit none
        real*8, intent(in), dimension(3) :: vec1, vec2
        real*8, dimension(3) :: vecCross
        vecCross(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
        vecCross(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
        vecCross(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)
    end function vecCross

    function vecUnitNorm( vec1, vec2 )
        implicit none
        real*8, intent(in), dimension(3) :: vec1, vec2
        real*8, dimension(3) :: vecUnitNorm
        real*8, dimension(3) :: vecNorm
        vecNorm = vecCross(vec1,vec2)
        vecUnitNorm(:) = vecNorm(:)/norm(vecNorm)
    end function vecUnitNorm
    

    real*8 function getdr( vec_dr, vec_box )
        implicit none
        real*8, intent(in), dimension(3) :: vec_box, vec_dr
        real*8, dimension(3) :: tempvec
    
        tempvec(:) = vec_dr(:) - (nint(vec_dr(:)/vec_box(:)))*vec_box(:)
        getdr = norm(tempvec)
    end function getdr
    
    real*8 function norm( vec )
        implicit none
        real*8, intent(in), dimension(3) :: vec
    
        norm = sqrt(vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3))
    
    end function norm
end module fvector

MODUlE eam
    implicit none

    real*8,parameter :: fpdup0=0.2515463239137249, fhup0=-0.0296601759808718
    real*8,parameter :: rho0_h=7.98909, rho0_pd = 10.261
    real*8,parameter :: ah=9.99780, bh=60.0155, ch=0.000197047
    real*8,parameter :: dh=1.18860, eh=0.0540638
    real*8,parameter :: ch_rho=11.0025, delta_h=1.30927
    real*8,parameter :: D_hh=0.0661496, alpha_hh=3.67263, beta_hh=1.47797, r0_hh=2.51980
    real*8,parameter :: D_pdh=0.2494540,alpha_pdh=4.82613,beta_pdh=2.13158,r0_pdh=1.50964
    real*8,parameter :: PI=4.D0*DATAN(1.D0)
 
contains


!*******************************
! functions to shift rho and phi at rcut-rs
!*******************************
    real*8 function rcscale( r ) 
        use variables
        implicit none
        real*8, intent(in) :: r
        real*8             :: r5

        if (r .ge. rcut_all) then
            rcscale = 0
        elseif (r .le. rcut_all-rcut_rs) then
            rcscale = 1.d0
        else
            rcscale = (1.d0+dcos(PI*(r-rcut_all+rcut_rs)/rcut_rs))/2.d0
        endif
    end function rcscale

!*******************************
! functions to caclate rho
!*******************************
    real*8 function rhoa( r ) 
        implicit none
        real*8, intent(in) :: r
        real*8             :: r5

        r5 = r/5.d0
        rhoa = r5*(r5*(r5-2.7267629107325706) + 1.8716766113599643) &
                 *(r5*(r5-2.50290548851635) + 1.668549182690922) &
                 *(r5*(r5-2.0924467509943674) + 1.3150372774478005) &
                 *(r5*(r5-1.564328475106985) + 0.8987511149780485) &
                 *(r5*(r5-1.009780903403673) + 0.5124363774128722) &
                 *(r5*(r5-0.5304054524800665) + 0.2169886022464641) &
                 *(r5*(r5-0.1356566408715063) + 0.035852347523891395)

    end function rhoa

    real*8 function rhoa_pd( r ) 
        implicit none
        real*8, intent(in) :: r
        real*8             :: r5

        r5=r/5.d0
        rhoa_pd = -0.02972698211669922 &
                 +r5*(0.6676807403564453 + r5* &
                       (-255.8965835571289 + r5* &
                         (14673.409149169922-(2.597301181336601e7)*rhoa(r))))
        rhoa_pd = rhoa_pd * rcscale(r) 
    end function rhoa_pd

    real*8 function rhoa_h( r ) 
        implicit none
        real*8, intent(in) :: r

        rhoa_h = ch_rho*exp(-delta_h*r)*rcscale(r)

    end function rhoa_h

!*******************************
! functions to caclate F_alpha
!*******************************
    real*8 function fpdu( rho ) 
        implicit none
        real*8, intent(in) :: rho
        real*8             :: rhop

        rhop=rho/50.d0
        fpdu = 295878.9003038662 * (rhop-0.20581955357385892) &
              *(rhop-0.081228755904399)*rhop*(rhop+0.05298811034615951) &
              *(rhop*(rhop-2.4242616904962846) + 1.4791899886249564) &
              *(rhop*(rhop-2.1376274623740064) + 1.2169215689822592) &
              *(rhop*(rhop-1.6486007989726832) + 0.8159825255339774) &
              *(rhop*(rhop-1.0749204110338482) + 0.42007491336688396) &
              *(rhop*(rhop-0.5128056047933808) + 0.12468685331167456)

    end function fpdu

    real*8 function fpd( rho ) 
        implicit none
        real*8, intent(in) :: rho

        fpd = fpdu(rho) - fpdup0*rho
    end function fpd

    real*8 function fhu( rho )
        implicit none
        real*8, intent(in) :: rho
        real*8             :: rhoeh, rhoehdh

        rhoeh = rho+eh
        rhoehdh = rhoeh**dh
        fhu = -ch*rhoehdh*((rhoeh**2)/(2+dh) - (ah+bh)/(1+dh)*rhoeh + (ah*bh)/dh)

    end function fhu

    real*8 function fh( rho ) 
        implicit none
        real*8, intent(in) :: rho

        fh = fhu(rho) - fhup0*rho
    end function fh

!*******************************
! functions to caclate F_alpha
!*******************************
    real*8 function phi_pdu( r ) 
        implicit none
        real*8, intent(in) :: r
        real*8             :: r5

        r5=r/5.d0
        phi_pdu = -79415.24035137112 &
               *(r5-1.0699996145674568)*(r5-1.06015072612581) &
               *(r5-0.42433991011376526)*(r5+0.06169160085238687) &
               *(r5*(r5-2.0586473420376348) + 1.0683922574015199) &
               *(r5*(r5-1.6696359816422877) + 0.7337878627470482) &
               *(r5*(r5-1.1690370066230809) + 0.3909805777737639) &
               *(r5*(r5-0.2635598721249787) + 0.033551116514910245)
    end function phi_pdu

    real*8 function phi_pd( r ) 
        implicit none
        real*8, intent(in) :: r

        phi_pd = (phi_pdu(r)+2*fpdup0*rhoa_pd(r))*rcscale(r)
    end function phi_pd

    real*8 function phi_h( r ) 
        implicit none
        real*8, intent(in) :: r

        phi_h = (D_hh*(beta_hh*exp(-alpha_hh*(r-r0_hh))-alpha_hh*exp(-beta_hh*(r-r0_hh))))*rcscale(r)
    end function phi_h

    real*8 function phi_pdh( r ) 
        implicit none
        real*8, intent(in) :: r

        phi_pdh = (D_pdh*(beta_pdh*exp(-alpha_pdh*(r-r0_pdh))-alpha_pdh*exp(-beta_pdh*(r-r0_pdh))))*rcscale(r)
    end function phi_pdh

end module eam
