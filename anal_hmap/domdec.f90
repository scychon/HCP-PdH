!****************************************************************
! 
! this module contains functions to perform domain decomposition
!
!****************************************************************

MODUlE domdec
    implicit none
 
contains

    ! Subroutine to get command arguments and assign input and output file name
    subroutine init_domdec
        use variables
        implicit none
        integer :: iAtom,i,j,k,l,idxDomain,iBoxVol
        real*8,dimension(3) :: rmax
        integer,dimension(3) :: idxBox

        allocate(PdBox(3,nPdAtoms))

        rmax = maxval(PdPos,dim=2)
        nDomainVec(:) = int(rmax(:)/(rcut_pdpd*3))+1
        numDomainTot = nDomainVec(1)*nDomainVec(2)*nDomainVec(3)
        numBoxTot = nDomainVec(1)*nDomainVec(2)*nDomainVec(3)*27
        nBoxVec(:) = nDomainVec(:)*3
        write(*,*) 'nDomain : ', nDomainVec
        write(*,*) 'nBox : ', nBoxVec
        write(*,*) 'rmax : ', rmax
        write(*,*) 'rmax/rcut*3 : ', rmax/(rcut_pdpd*3)

        ! calculate number of Pd atoms in each domdec subdomain
        allocate(nPdBox(nBoxVec(1),nBoxVec(2),nBoxVec(3)))
        allocate(nPdDomain(numDomainTot))
        nPdBox=0
        nPdDomain=0
        nTotAtoms=nPdAtoms

        do iAtom=1,nPdAtoms
            idxBox(:) = PdPos(:,iAtom)/(rcut_pdpd)+1
            PdBox(:,iAtom) = idxBox(:)
            nPdBox(idxBox(1),idxBox(2),idxBox(3)) = nPdBox(idxBox(1),idxBox(2),idxBox(3))+1
            idxDomain = ((idxBox(1)-1)/3)*(nDomainVec(2)*nDomainVec(3)) + ((idxBox(2)-1)/3)*nDomainVec(3) + ((idxBox(3)-1)/3) + 1
            nPdDomain(idxDomain) = nPdDomain(idxDomain)+1
        enddo

        nPdBoxMax = maxval(nPdBox)
        nPdDomainMax = maxval(nPdDomain)
        write(*,*) 'max nPdBox : ',nPdBoxMax
        write(*,*) 'max nPdDomain : ',nPdDomainMax
        allocate(idxPdAtomInBox(nPdBoxMax,nBoxVec(1),nBoxVec(2),nBoxVec(3)))
        idxPdAtomInBox=-1
        nPdBox=0
        do iAtom=1,nPdAtoms
            idxBox(:) = PdBox(:,iAtom)
            nPdBox(idxBox(1),idxBox(2),idxBox(3)) = nPdBox(idxBox(1),idxBox(2),idxBox(3))+1
            idxPdAtomInBox(nPdBox(idxBox(1),idxBox(2),idxBox(3)),idxBox(1),idxBox(2),idxBox(3))=iAtom
        enddo

        call get_NL_Pd

        !get the volume of the boxes that contain Pd atoms
        iBoxVol = 0
        do i=1, nBoxVec(1)
          do j=1, nBoxVec(2)
            do k=1, nBoxVec(3)
                iBoxVol = iBoxVol + 1
            enddo
          enddo
        enddo
        volume=rcut_pdpd**3 * iBoxVol
        nqvol = nqunit* (1.008d0*temperature)**1.5d0 * volume

        !maximum nStepNLUpdate new H atoms can be added to each box
        !nHBoxMax = nStepNLUpdate*3
        nHBoxMax = 60
        allocate(idxHAtomInBox(nHBoxMax,nBoxVec(1),nBoxVec(2),nBoxVec(3)))
        idxHAtomInBox=-1

        !maximum nPdDomainMax new H atoms can be added to each domain (initially)
        nHMax = 2*nPdDomainMax*nDomainVec(1)*nDomainVec(2)*nDomainVec(3)
        write(*,*) 'nHMax : ',nHMax
        nHMax = nHBoxMax*nDomainVec(1)*nDomainVec(2)*nDomainVec(3)*27
        write(*,*) 'nHMax : ',nHMax
        allocate(HBox(3,nHMax)) 
        allocate(HPos(3,nHMax))
        allocate(nHDomain(numDomainTot))
        allocate(nHBox(nBoxVec(1),nBoxVec(2),nBoxVec(3)))
        allocate(nHNL(nBoxVec(1),nBoxVec(2),nBoxVec(3)))
        write(*,*) 'nhnl'
        HPos = -1
        nHDomain = 0
        nHBox=0
        HBox=-1     !unassigned H
        nHAtoms=0
        nHNL=0
        write(*,*) 'hbox, hpos nhdomain allocated'

        !assume maximum 300 H in NL
        nHNLMax = 300
        call get_NL_H

        ! allocate rho & f matrices
        allocate(rho_pdh(nPdAtoms))
        allocate(rho_hpd(nHMax))
        allocate(rho_hh(nHMax))
        allocate(phisum_pdh(nPdAtoms))
        allocate(phisum_hh(nHMax))
        allocate(f_h(nHMax))

        rho_pdh=0
        rho_hpd=0
        rho_hh=0
        phisum_pdh=0
        phisum_hh=0
        phitot=0
        phitot_pdh=0
        phitot_hh=0
        f_h=0

    end subroutine init_domdec

    ! subroutine to get NL for Pd atoms
    subroutine get_NL_Pd
        use variables
        use fvector
        implicit none
        integer :: i,j,k,iAtom
        integer,dimension(3) :: idxBox

        !assume maximum 100 pd in NL
        nPdNLMax = 100
        if (allocated(nPdNL)) deallocate(nPdNL)
        if (allocated(idxNLPd)) deallocate(idxNLPd)
        allocate(nPdNL(nBoxVec(1),nBoxVec(2),nBoxVec(3)))
        allocate(idxNLPd(nPdNLMax,nBoxVec(1),nBoxVec(2),nBoxVec(3)))

        nPdNL=0
        idxNLPd = -1
        do iAtom=1,nPdAtoms
          idxBox(:) = PdBox(:,iAtom)
          do i=max(idxBox(1)-1,1),min(idxBox(1)+1,nBoxVec(1))
            do j=max(idxBox(2)-1,1),min(idxBox(2)+1,nBoxVec(2))
              do k=max(idxBox(3)-1,1),min(idxBox(3)+1,nBoxVec(3))
                nPdNL(i,j,k)=nPdNL(i,j,k)+1
                if (nPdNL(i,j,k) .gt. nPdNLMax) then
                    call expand4Dint(idxNLPd,nPdNLMax,0,0,0)
                    nPdNLMax = nPdNLMax*2
                endif
                idxNLPd(nPdNL(i,j,k),i,j,k)=iAtom
              enddo
            enddo
          enddo
        enddo

        nPdNLMax = maxval(nPdNL)
        call shrink4Dint(idxNLPd,nPdNLMax,nBoxVec(1),nBoxVec(2),nBoxVec(3))

        write(*,*) 'max nPdNL : ',nPdNLMax

    end subroutine get_NL_Pd

    ! subroutine to get NL for H atoms
    subroutine get_NL_H
        use variables
        use fvector
        implicit none
        integer :: i,j,k,iAtom,tempmax
        integer,dimension(3) :: idxBox,idxTemp

        write(*,*) 'start get_NL_H'
        if (.not. allocated(nHNL)) allocate(nHNL(nBoxVec(1),nBoxVec(2),nBoxVec(3)))
        if (.not. allocated(idxNLH)) allocate(idxNLH(nHNLMax,nBoxVec(1),nBoxVec(2),nBoxVec(3)))
        nHNLMax = size(idxNLH(:,1,1,1))

        nHNL=0
        idxNLH = -1
        do iAtom=1,nHMax
          idxBox(:) = HBox(:,iAtom)
          if (idxBox(1) .lt. 0) then
            if (any(HPos(:,iAtom) .ge. 0)) write(*,*) 'ERROR : box is not assigned but hpos is assigned !'
            cycle
          endif
          idxTemp(:) = HPos(:,iAtom)/(rcut_pdpd)+1
          if (any(idxBox .ne. idxTemp)) then
            write(*,*) 'error box not match:', idxBox, idxTemp, iAtom,HPos(:,iAtom)
          endif
        
          do i=max(idxBox(1)-1,1),min(idxBox(1)+1,nBoxVec(1))
            do j=max(idxBox(2)-1,1),min(idxBox(2)+1,nBoxVec(2))
              do k=max(idxBox(3)-1,1),min(idxBox(3)+1,nBoxVec(3))
                !skip box that has no Pd atoms
                if (nPdBox(i,j,k)==0) cycle

                nHNL(i,j,k)=nHNL(i,j,k)+1
                if (nHNL(i,j,k) .gt. nHNLMax) then
                    write(*,*) 'tot nHNL : ', sum(nHNL)
                    call expand4Dint(idxNLH,nHNLMax,0,0,0)
                    nHNLMax = nHNLMax*2
                endif
                idxNLH(nHNL(i,j,k),i,j,k)=iAtom
              enddo
            enddo
          enddo
        enddo

        tempmax = maxval(nHNL)
        write(*,*) 'max nHNL : ',tempmax

        !after nStepNLUpdate mc moves, max number of H atoms within the NL (for corner box, 8 neighboring domains)
        if (tempmax .gt. nHNLMax/2) then
            call expand4Dint(idxNLH,nHNLMax,0,0,0)
            nHNLMax = nHNLMax*2
            write(*,*) 'max nHNL to be increased : ',nHNLMax
        endif

    end subroutine get_NL_H

    ! subroutine to calculate rho_pdh
    subroutine calc_pdh
        use variables
        use fvector
        use eam
        use omp_lib
        implicit none
        integer :: i,j
        real*8  :: dr
        integer,dimension(3) :: idxBox

        allocate(rho_pdh(nPdAtoms))
        rho_pdh=0.d0

        !$OMP PARALLEL &
        !$OMP   DEFAULT (FIRSTPRIVATE) &
        !$OMP   SHARED (rho_pdh)
        !$OMP DO
        do i=1,nPdAtoms
            do j=1,nHAtoms
                dr = norm(PdPos(:,i)-HPos(:,j))
                if (dr .le. rcut_pdh) then
                    rho_pdh(i) = rho_pdh(i)+rhoa_h(dr)
                endif
            enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        OPEN(unit=7,file='pdpos.dat')
        do i=1,nPdAtoms
            write(7,*) i, PdPos(:,i)
        enddo
        close(7)

        OPEN(unit=7,file='rho_pdh.dat')
        do i=1,nPdAtoms
            write(7,*) PdPos(:,i), rho_pdh(i)
        enddo
        close(7)

    end subroutine calc_pdh

    ! subroutine to calculate rho_pdpd & phisum_pdpd & f_pd
    subroutine calc_pdpd
        use variables
        use fvector
        use eam
        use omp_lib
        implicit none
        integer :: i,j
        real*8  :: dr
        integer,dimension(3) :: idxBox
        real*8, allocatable,dimension(:) :: phisum_pdpd

        allocate(f_pd(nPdAtoms))
        allocate(rho_pdpd(nPdAtoms))
        allocate(phisum_pdpd(nPdAtoms))
        f_pd = 0.d0
        rho_pdpd=0.d0
        phisum_pdpd=0.d0

        !$OMP PARALLEL &
        !$OMP   DEFAULT (FIRSTPRIVATE) &
        !$OMP   SHARED (rho_pdpd,phisum_pdpd,f_pd)
        !$OMP DO
        do i=1,nPdAtoms
            do j=1,nPdAtoms
              if ( i .ne. j ) then
                dr = norm(PdPos(:,i)-PdPos(:,j))
                if (dr .le. rcut_pdpd) then
                    rho_pdpd(i) = rho_pdpd(i)+rhoa_pd(dr)
                    phisum_pdpd(i) = phisum_pdpd(i)+phi_pd(dr)
                    !write(*,*) phi_pd(dr)
                endif
              endif
            enddo
            phisum_pdpd(i) = phisum_pdpd(i)*0.5d0
            f_pd(i) = fpd(rho_pdpd(i)+rho_pdh(i))
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        phitot_pdpd=0
        OPEN(unit=7,file='rho_phi_pdpd.dat')
        do i=1,nPdAtoms
            phitot_pdpd = phitot_pdpd + phisum_pdpd(i)
            write(7,*) i, rho_pdpd(i), phisum_pdpd(i), PdGroup(i)
        enddo
        close(7)
        write(*,*) 'phitot_pdpd : ', phitot_pdpd

        OPEN(unit=7,file='rhophi.dat')
        do i=1,5000
            dr = i*1.d-3
            write(7,*) dr, rhoa_pd(dr), rhoa_h(dr), phi_pd(dr), phi_pdh(dr), phi_h(dr)
        enddo
        close(7)

        OPEN(unit=7,file='fval.dat')
        do i=1,5000
            dr = i*1.d-2
            write(7,*) dr, fpd(dr), fh(dr)
        enddo
        close(7)
    end subroutine calc_pdpd

    ! subroutine to calculate rho_pdh
    subroutine map_slice
        use variables
        use fvector
        use eam
        use omp_lib
        implicit none
        integer :: i,j,k,l,idx,iAtom
        integer :: nx,ny,nz
        integer :: dn,npd,nh,nslice
        real*8  :: dr,sig,gauss,rho_curr,sum_gauss
        real*8  :: z_kn,z_k,z_kp,sum_pdrho
        integer,dimension(3) :: idxBox
        real*8, dimension(3) :: pos
        !real*8,dimension(:,:),allocatable :: sum_pdrho

        nx=230
        ny=230
        nz=size(SliceIdx)
        allocate(rho_hmap(nx,ny,nz))
        !allocate(sum_pdrho(nx,ny))
        rho_hmap=0.d0
        sig=5.d0
        dn=ceiling(10.d0/slice_res)
        z_k=(nx-SliceIdx(1))*slice_res
        z_kp=(nx-SliceIdx(2))*slice_res

        !$OMP PARALLEL &
        !$OMP   DEFAULT (FIRSTPRIVATE) &
        !$OMP   SHARED (rho_hmap)
        !$OMP DO
        do idx=2,nz-1
          k=SliceIdx(idx)
          z_kn=(nx-SliceIdx(idx-1))*slice_res
          z_k=(nx-SliceIdx(idx))*slice_res
          z_kp=(nx-SliceIdx(idx+1))*slice_res
          sum_pdrho=0.d0
          sum_gauss=0.d0
          nh = 0
          do iAtom=1,nHAtoms
            if ((HPos(3,iAtom).ge.(z_k+z_kn)/2.d0) .and. ((HPos(3,iAtom).le.(z_k+z_kp)/2.d0))) then
              nh=nh+1
            endif
          enddo
          npd = 0
          do iAtom=1,nPdAtoms
            if ((PdPos(3,iAtom).ge.(z_k+z_kn)/2.d0) .and. ((PdPos(3,iAtom).le.(z_k+z_kp)/2.d0))) then
              npd=npd+1
            endif
          enddo
          nslice=0
          do i=1,nx
            pos(1)=(nx-i+0.5d0)*slice_res
            do j=1,ny
              pos(2)=(ny-j+0.5d0)*slice_res
              rho_curr = 0.d0
              do iAtom=1,nHAtoms
                if ((HPos(3,iAtom).ge.(z_k+z_kn)/2.d0) .and. ((HPos(3,iAtom).le.(z_k+z_kp)/2.d0))) then
                  pos(3)=HPos(3,iAtom)
                  dr = norm(HPos(:,iAtom)-pos)
                  if (dr .le. 3*sig) then
                    gauss = exp(-dr**2/(2.d0*sig))
                    rho_curr = rho_curr+gauss
                  endif
                endif
              enddo
              if (rho_curr.gt.0.d0) then
                nslice=nslice+1
                sum_gauss=sum_gauss+rho_curr
                !$OMP CRITICAL
                rho_hmap(i,j,idx)=rho_curr
                !$OMP END CRITICAL
              endif
              !do iAtom=1,nPdAtoms
              !  if ((PdPos(3,iAtom).ge.(z_k+z_kn)/2.d0) .and. ((PdPos(3,iAtom).le.(z_k+z_kp)/2.d0))) then
              !    pos(3)=PdPos(3,iAtom)
              !    dr = norm(PdPos(:,iAtom)-pos)
              !    if (dr .le. 3*sig) then
              !      gauss = exp(-dr**2/(2.d0*sig))
              !      sum_pdrho=sum_pdrho+gauss
              !    endif
              !  endif
              !enddo
              !if (sum_pdrho .gt.0) then
                !rho_hmap(i,j,idx)=rho_hmap(i,j,idx)/sum_pdrho*nx*ny
              !endif
            enddo
          enddo
          write(*,*) idx,npd,nh, z_k, z_kn, z_kp,nslice,sum_gauss,sum_gauss/nslice,maxval(rho_hmap(:,:,idx)),minval(rho_hmap(:,:,idx)),sum(rho_hmap(:,:,idx))
          !!$OMP CRITICAL
          !rho_hmap(:,:,idx)=(rho_hmap(:,:,idx)*nslice*nh)/(sum_gauss*npd)
          !!$OMP END CRITICAL
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        OPEN(unit=7,file='rho_hmap.dat')
        do i=1,nz
            write(7,*) rho_hmap(:,:,i)
        enddo
        close(7)

    end subroutine map_slice

end module domdec
