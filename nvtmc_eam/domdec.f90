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
        integer :: i,j,k,l,iBoxVol
        integer :: iAtom
        integer :: nPdAtomsBak
        integer :: nOuterSurfAtoms, nInnerSurfAtoms
        real*8, DIMENSION(:,:),ALLOCATABLE      :: PdPosBak          ! array of Pd positions (xyz, nPdAtoms)

        !initial call to set the PdBox, NL, SurfPd
        nPdSurfAtoms=0
        call set_PdBox
        call get_NL_Pd
        call get_Surf_Pd

        !backup arrays and remove surface atoms
        nPdAtomsBak = nPdAtoms
        nOuterSurfAtoms = nPdSurfAtoms
        allocate(PdPosBak(3,nPdAtoms))
        PdPosBak = PdPos


        !remove surface atoms
        deallocate(PdPos)
        nPdAtoms = nPdAtoms-nPdSurfAtoms
        allocate(PdPos(3,nPdAtoms))
        do i=1,nPdAtoms
            PdPos(:,i) = PdPosBak(:,idxPdInAtom(i))
        enddo

        !call second surface Pd routine
        call set_PdBox
        call get_NL_Pd
        call get_Surf_Pd

        !reset PdPos and recall nearest neighbour Pd routine
        nPdAtoms = nPdAtomsBak
        deallocate(PdPos)
        allocate(PdPos(3,nPdAtoms))
        PdPos = PdPosBak
        call set_PdBox
        call get_NL_Pd
        call get_NN_Pd

        OPEN(unit=7,file='surfpd.xyz')
        write(7,*) nPdSurfAtoms
        write(7,*) 'Surface Pd XYZ'
        do i=1,nPdSurfAtoms
            if (i.le.nOuterSurfAtoms) then
                write(7,'(A3,3f13.6)') ' S1', PdPos(:,idxPdSurfAtom(i))
            else
                write(7,'(A3,3f13.6)') ' S2', PdPos(:,idxPdSurfAtom(i))
            endif
        enddo
        write(7,*) ' -1'
        close(7)

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

        ! Precalculate pd-pd interaction and rho_pdpd
        call calc_pdpd

        !initialize H matrices
        call init_Hmat

        !initialize rho,f,phi matrices
        call init_rhophi
    end subroutine init_domdec

    !initialize rho, f, phi matrices
    subroutine init_rhophi
        use variables
        use fvector
        use eam
        implicit none
        integer :: i,j
        real*8 :: dr

        ! allocate rho & f matrices
        allocate(rho_pdh(nPdAtoms))
        allocate(rho_hpd(nHMax))
        allocate(rho_hh(nHMax))
        allocate(phisum_pdh(nPdAtoms))
        allocate(phisum_hh(nHMax))
        allocate(f_h(nHMax))
        allocate(drho_pdh(nPdAtoms))
        allocate(drho_hpd(nHMax))
        allocate(drho_hh(nHMax))
        allocate(dphisum_pdh(nPdAtoms))
        allocate(dphisum_hh(nHMax))
        allocate(df_h(nHMax))

        rho_pdh=0
        rho_hpd=0
        rho_hh=0
        phisum_pdh=0
        phisum_hh=0
        phitot=0
        phitot_pdh=0
        phitot_hh=0
        f_h=0

        drho_pdh=0
        drho_hpd=0
        drho_hh=0
        dphisum_pdh=0
        dphisum_hh=0
        df_h=0

        !$OMP PARALLEL &
        !$OMP   DEFAULT (FIRSTPRIVATE) &
        !$OMP   SHARED (rho_pdh,phisum_pdh,f_pd)
        !$OMP DO
        do i=1,nPdAtoms
            do j=1,nHAtoms
                dr = norm(PdPos(:,i)-HPos(:,j))
                if (dr .le. rcut_pdh) then
                    rho_pdh(i) = rho_pdh(i)+rhoa_h(dr)
                    phisum_pdh(i) = phisum_pdh(i)+phi_pdh(dr)
                endif
            enddo
            f_pd(i) = fpd(rho_pdpd(i)+rho_pdh(i))
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        do i=1,nHAtoms
            do j=1,nPdAtoms
                dr = norm(HPos(:,i)-PdPos(:,j))
                if (dr .le. rcut_pdh) then
                    rho_hpd(i) = rho_hpd(i)+rhoa_pd(dr)
                endif
            enddo
            do j=1,nHAtoms
                if(i==j) cycle
                dr = norm(HPos(:,i)-HPos(:,j))
                if (dr .le. rcut_hh) then
                    rho_hh(i) = rho_hh(i)+rhoa_h(dr)
                    phisum_hh(i) = phisum_hh(i)+0.5*phi_h(dr)
                endif
            enddo
            f_h(i) = fh(rho_hpd(i)+rho_hh(i))
        enddo

        phitot_pdh=0
        OPEN(unit=7,file='rho_phi_pdh.dat')
        do i=1,nPdAtoms
            phitot_pdh = phitot_pdh + phisum_pdh(i)
            write(7,*) i, rho_pdh(i), phisum_pdh(i)
        enddo
        close(7)
        phitot_hh=0
        OPEN(unit=7,file='rho_phi_h.dat')
        do i=1,nHAtoms
            phitot_hh = phitot_hh + phisum_hh(i)
            write(7,*) i, rho_hh(i), phisum_hh(i)
        enddo
        close(7)
        write(*,*) 'phitot_pdpd, phitot_pdh, phitot_hh : ', phitot_pdpd, phitot_pdh, phitot_hh

    end subroutine init_rhophi


    !initialize matrices related with HPos
    subroutine init_Hmat
        use variables
        use fvector
        use gen_pos
        implicit none

        write(*,*) 'start Initializing HMat'
        !generate HPos
        if (rHPerPd .gt. 0) then
            nHAtoms = int(rHPerPd*nPdAtoms+.5)
            nHMax = nHAtoms
            call generate_HPos
            write(*,*) 'rHPerPd is set. rHPerPd : ',rHPerPd, ', nHAtoms : ',nHAtoms
            call set_HBox
        elseif (bHTrajAnal) then
            call set_HBox
        else
            nHAtoms = 0
            call allocHmat
        endif
    end subroutine init_Hmat

    subroutine set_HBox
        use variables
        use gen_pos
        implicit none
        integer :: idxDomain,iAtom
        integer,dimension(3) :: idxBox
        !maximum nStepNLUpdate new H atoms can be added to each box
        !nHBoxMax = nStepNLUpdate*3
        nHBoxMax = 60
        allocate(idxHAtomInBox(nHBoxMax,nBoxVec(1),nBoxVec(2),nBoxVec(3)))
        idxHAtomInBox=-1

        allocate(nHDomain(numDomainTot))
        allocate(nHBox(nBoxVec(1),nBoxVec(2),nBoxVec(3)))
        allocate(nHNL(nHAtoms))

        write(*,*) 'nHBox, nHNL, nHDomain allocated'
        nHDomain = 0
        nHBox=0
        nHNL=0

        do iAtom=1,nHAtoms
            idxBox(:) = HBox(:,iAtom)
            nHBox(idxBox(1),idxBox(2),idxBox(3)) = nHBox(idxBox(1),idxBox(2),idxBox(3))+1
            idxHAtomInBox(nHBox(idxBox(1),idxBox(2),idxBox(3)),idxBox(1),idxBox(2),idxBox(3))=iAtom
            idxDomain = ((idxBox(1)-1)/3)*(nDomainVec(2)*nDomainVec(3)) + ((idxBox(2)-1)/3)*nDomainVec(3) + ((idxBox(3)-1)/3) + 1
            nHDomain(idxDomain) = nHDomain(idxDomain)+1
        enddo
        !assume maximum 300 H in NL
        nHNLMax = 300
        call get_NL_H
        call get_NN_H
    end subroutine set_HBox

    subroutine allocHmat
        use variables
        implicit none
        !maximum nPdDomainMax new H atoms can be added to each domain (initially)
        nHMax = 2*nPdDomainMax*nDomainVec(1)*nDomainVec(2)*nDomainVec(3)
        write(*,*) 'nHMax : ',nHMax
        nHMax = nHBoxMax*nDomainVec(1)*nDomainVec(2)*nDomainVec(3)*27
        write(*,*) 'nHMax : ',nHMax
        allocate(HBox(3,nHMax)) 
        allocate(HPos(3,nHMax))
        write(*,*) 'nhnl'

        HPos = -1
        HBox=-1     !unassigned H
        nHAtoms=0
        write(*,*) 'hbox, hpos allocated'

        call set_HBox        
    end subroutine allocHmat


    ! Subroutine to get command arguments and assign input and output file name
    subroutine set_PdBox
        use variables
        implicit none
        integer :: iAtom,i,j,k,l,idxDomain
        real*8,dimension(3) :: rmax
        integer,dimension(3) :: idxBox

        if (.not. allocated(PdBox))   allocate(PdBox(3,nPdAtoms))

        rmax = maxval(PdPos,dim=2)+3
        nDomainVec(:) = int(rmax(:)/(rcut_pdpd*3))+1
        numDomainTot = nDomainVec(1)*nDomainVec(2)*nDomainVec(3)
        numBoxTot = nDomainVec(1)*nDomainVec(2)*nDomainVec(3)*27
        nBoxVec(:) = nDomainVec(:)*3
        write(*,*) 'nDomain : ', nDomainVec
        write(*,*) 'nBox : ', nBoxVec
        write(*,*) 'rmax : ', rmax
        write(*,*) 'rmax/rcut*3 : ', rmax/(rcut_pdpd*3)

        ! calculate number of Pd atoms in each domdec subdomain
        if (.not. allocated(nPdBox))    allocate(nPdBox(nBoxVec(1),nBoxVec(2),nBoxVec(3)))
        if (.not. allocated(nPdDomain)) allocate(nPdDomain(numDomainTot))
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
        if (allocated(idxPdAtomInBox))   deallocate(idxPdAtomInBox)
        allocate(idxPdAtomInBox(nPdBoxMax,nBoxVec(1),nBoxVec(2),nBoxVec(3)))
        idxPdAtomInBox=-1
        nPdBox=0
        do iAtom=1,nPdAtoms
            idxBox(:) = PdBox(:,iAtom)
            nPdBox(idxBox(1),idxBox(2),idxBox(3)) = nPdBox(idxBox(1),idxBox(2),idxBox(3))+1
            idxPdAtomInBox(nPdBox(idxBox(1),idxBox(2),idxBox(3)),idxBox(1),idxBox(2),idxBox(3))=iAtom
        enddo
    end subroutine set_PdBox

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
          do k=max(idxBox(3)-1,1),min(idxBox(3)+1,nBoxVec(3))
            do j=max(idxBox(2)-1,1),min(idxBox(2)+1,nBoxVec(2))
              do i=max(idxBox(1)-1,1),min(idxBox(1)+1,nBoxVec(1))
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

    ! subroutine to get nearest neighbour Pd atom for H atoms
    subroutine get_NN_H
        use variables
        use fvector
        implicit none
        integer :: iAtom,jAtom,i
        integer,dimension(3) :: idxBox
        real*8 :: dr

        if (.not. allocated(idxHNNAtom)) allocate(idxHNNAtom(nHAtoms+1))
        if (.not. allocated(drHNNAtom)) allocate(drHNNAtom(nHAtoms+1))

        !construct nearest neighbour list for Pd atoms
        idxHNNAtom = -1
        drHNNAtom = rcut_new
        !$OMP PARALLEL &
        !$OMP   DEFAULT (FIRSTPRIVATE) &
        !$OMP   SHARED (idxHNNAtom,drHNNAtom,HBox,nPdNL,idxNLPd)
        !$OMP DO
        do iAtom=1,nHAtoms
            idxBox(:) = HBox(:,iAtom)
            do i=1,nPdNL(idxBox(1),idxBox(2),idxBox(3))
                jAtom = idxNLPd(i,idxBox(1),idxBox(2),idxBox(3))
                dr = norm(PdPos(:,jAtom)-HPos(:,iAtom))
                if (dr .lt. drHNNAtom(iAtom)) then
                    idxHNNAtom(iAtom) = jAtom
                    drHNNAtom(iAtom) = dr
                endif
            enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        write(*,*) 'Pd nearest neighbour of H atoms are set !'
    end subroutine get_NN_H

    ! subroutine to get nearest neighbour Pd atoms
    subroutine get_NN_Pd
        use variables
        use fvector
        implicit none
        integer :: iAtom,jAtom,i
        integer,dimension(3) :: idxBox

        if (allocated(idxPdNNAtoms)) deallocate(idxPdNNAtoms)
        if (allocated(drPdNNAtoms)) deallocate(drPdNNAtoms)
        allocate(idxPdNNAtoms(nPdNNAtoms,nPdAtoms))
        allocate(drPdNNAtoms(nPdNNAtoms,nPdAtoms))

        !construct nearest neighbour list for Pd atoms
        idxPdNNAtoms = -1
        drPdNNAtoms = rcut_pdpd
        !$OMP PARALLEL &
        !$OMP   DEFAULT (FIRSTPRIVATE) &
        !$OMP   SHARED (idxPdNNAtoms,drPdNNAtoms,PdBox,nPdNL,idxNLPd)
        !$OMP DO
        do iAtom=1,nPdAtoms
            idxBox(:) = PdBox(:,iAtom)
            do i=1,nPdNL(idxBox(1),idxBox(2),idxBox(3))
                jAtom = idxNLPd(i,idxBox(1),idxBox(2),idxBox(3))
                if (iAtom==jAtom) cycle
                call add_Pd_NN(iAtom,jAtom)
            enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        write(*,*) 'Pd nearest neighbours are set !'
    end subroutine get_NN_Pd

    ! subroutine to get Surface Pd atoms
    subroutine get_Surf_Pd
        use variables
        use fvector
        implicit none
        integer :: iAtom,jAtom
        integer :: i,j,k,l,idxDomain
        integer :: nPdSurf,nPdNN
        logical :: bExistsOuter
        real*8  :: scaleDrNorm
        real*8  :: drNNAvg
        real*8,dimension(3) :: drCOM,drNN
        !integer,dimension(3) :: idxCOMBox
        logical,dimension(nPdAtoms) :: bSurfAtom
        integer,dimension(:),allocatable :: idxPdInAtomCurr

        !idxCOMBox(:) = rCOM(:)/rcut_pdpd+1
        !assume maximum 100 pd on surface
        if (.not. allocated(nPdNL)) call get_NL_Pd

        call get_NN_Pd

        if (nPdSurfAtoms.eq.0) then
            scaleDrNorm = 0.5
        else
            scaleDrNorm = 0.2
        endif
        bSurfAtom = .false.
        nPdSurf = 0
        !calculate dr vector to COM
        !$OMP PARALLEL &
        !$OMP   DEFAULT (FIRSTPRIVATE) &
        !$OMP   SHARED (idxPdNNAtoms,drPdNNAtoms,bSurfAtom,nPdSurf,PdPos)
        !$OMP DO
        do iAtom=1,nPdAtoms
            drCOM = rCOM(:)-PdPos(:,iAtom)
            drNNAvg = norm(drCOM)
            !drCOM = 0.d0
            !nPdNN = 0
            !do i=1,nPdNNAtoms
            !    jAtom = idxPdNNAtoms(i,iAtom)
            !    if (jAtom .gt.0) then
            !        drCOM(:) = drCOM(:) + PdPos(:,jAtom)
            !        nPdNN = nPdNN +1
            !    endif
            !enddo
            !drCOM(:) = drCOM(:)/nPdNN - PdPos(:,iAtom)
            !drNNAvg = norm(drCOM)
            !drCOM(:) = rCOM(:)-PdPos(:,iAtom)
            bExistsOuter = .false.
            do i=1,nPdNNAtoms
                jAtom = idxPdNNAtoms(i,iAtom)
                if (jAtom .lt. 0) then
                    exit
                else
                    drNN = PdPos(:,jAtom)-PdPos(:,iAtom)
                    if (dot_product(drCOM,drNN)/norm(drNN) .lt. -scaleDrNorm*drNNAvg) then
                    !if (dot_product(drCOM,drNN)/norm(drNN) .lt. -0.25*drNNAvg) then
                        bExistsOuter = .true.
                    endif
                endif
            enddo
            if (.not. bExistsOuter) then
                !$OMP CRITICAL
                bSurfAtom(iAtom)= .true.
                nPdSurf = nPdSurf + 1
                !$OMP END CRITICAL
            endif
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        write(*,*) nPdSurf, ' Pd Surface Atoms are found !'
        if (allocated(idxPdSurfAtom)) then
            call expand1Dint(idxPdSurfAtom,nPdSurf)
        else
            allocate(idxPdSurfAtom(nPdSurf))
            allocate(bPdSurfAtom(nPdAtoms))
            idxPdSurfAtom=0
            bPdSurfAtom = .false.
        endif

        jAtom=0
        allocate(idxPdInAtomCurr(nPdAtoms-nPdSurf))
        do iAtom=1,nPdAtoms
            if (bSurfAtom(iAtom)) then
                nPdSurfAtoms = nPdSurfAtoms+1
                idxPdSurfAtom(nPdSurfAtoms) = idxPdInAtom(iAtom)
                bPdSurfAtom(idxPdInAtom(iAtom)) = .true.
            else
                jAtom = jAtom+1
                idxPdInAtomCurr(jAtom) = idxPdInAtom(iAtom)
            endif
        enddo
        idxPdInAtom(:jAtom) = idxPdInAtomCurr(:)
        call shrink1Dint(idxPdInAtom,jAtom)
        if (jAtom .ne. nPdAtoms-nPdSurf) then
            write(*,*) 'ERROR IN idxPdInAtom !',jAtom, nPdAtoms,nPdSurf, nPdAtoms-nPdSurf
        endif
    end subroutine get_Surf_Pd

    ! subroutine to check and add Pd atom to NN
    subroutine add_Pd_NN (iAtom, jAtom)
        use variables
        use fvector
        implicit none
        integer, intent(in) :: iAtom,jAtom
        logical :: bSmaller
        integer :: i,j
        real*8  :: dr
        integer, dimension(nPdNNAtoms) :: idxPdNN
        real*8, dimension(nPdNNAtoms)  :: drPdNN

        idxPdNN(:) = idxPdNNAtoms(:,iAtom)
        drPdNN(:) = drPdNNAtoms(:,iAtom)
        dr = norm(PdPos(:,iAtom)-PdPos(:,jAtom))
        if (dr .gt. rcut_pdpd)  return

        bSmaller = .false.
        do i= nPdNNAtoms,1,-1
            if ((idxPdNN(i) .lt. 0) .or. (drPdNN(i) .gt. dr)) then
                if (i .lt. nPdNNAtoms) then
                    idxPdNN(i+1) = idxPdNN(i)
                    drPdNN(i+1) = drPdNN(i)
                endif
                idxPdNN(i) = jAtom
                drPdNN(i) = dr
                if (i.gt.1) then
                    cycle
                endif
            elseif (i==nPdNNAtoms) then
                exit
            endif

            !$OMP CRITICAL
            idxPdNNAtoms(:,iAtom) = idxPdNN(:)
            drPdNNAtoms(:,iAtom) = drPdNN(:)
            !$OMP END CRITICAL
        enddo
    end subroutine add_Pd_NN

    ! subroutine to calculate rho_pdpd & phisum_pdpd & f_pd
    subroutine calc_pdpd
        use variables
        use fvector
        use eam
        use omp_lib
        implicit none
        integer :: i,j,idx
        real*8  :: dr
        integer,dimension(3) :: idxBox
        real*8,dimension(nRbins) :: rho_pdpd_rbin, rho_pdh_rbin
        real*8,dimension(nPdNNAtoms) :: drlist
        real*8,dimension(:),allocatable :: rho_pdpd_new
        integer,dimension(nRbins) :: nPdRbins

        allocate(f_pd(nPdAtoms))
        allocate(rho_pdpd(nPdAtoms))
        allocate(rho_pdpd_new(nPdAtoms))
        allocate(phisum_pdpd(nPdAtoms))
        allocate(df_pd(nPdAtoms))

        f_pd = 0.d0
        rho_pdpd=0.d0
        rho_pdpd_new=0.d0
        phisum_pdpd=0.d0

        !$OMP PARALLEL &
        !$OMP   DEFAULT (FIRSTPRIVATE) &
        !$OMP   SHARED (rho_pdpd,phisum_pdpd,f_pd,rho_pdpd_new)
        !$OMP DO
        do i=1,nPdAtoms
            do j=1,nPdAtoms
              if ( i .ne. j ) then
                dr = norm(PdPos(:,i)-PdPos(:,j))
                if (dr .le. rcut_pdpd) then
                    rho_pdpd(i) = rho_pdpd(i)+rhoa_pd(dr)
                    rho_pdpd_new(i) = rho_pdpd_new(i)+rhoa_h(dr)
                    phisum_pdpd(i) = phisum_pdpd(i)+phi_pd(dr)
                    !write(*,*) phi_pd(dr)
                endif
              endif
            enddo
            phisum_pdpd(i) = phisum_pdpd(i)*0.5d0
            f_pd(i) = fpd(rho_pdpd(i))
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        rho_pdpd_rbin = 0
        nPdRbins = 0
        do i=1,nPdAtoms
            idx=idxPdRbin(i)
            nPdRbins(idx) = nPdRbins(idx) + 1
            rho_pdpd_rbin(idx) = rho_pdpd_rbin(idx) + rho_pdpd_new(i)
        enddo


        OPEN(unit=7,file='rho_pdpd_rhoah_rbin.dat')
        write(7,*) '#ibin,nPdRbin, rho_pdpd_rhoah'
        do i=1,nRbins
            write(7,*) i, nPdRbins(i), rho_pdpd_rbin(i)
        enddo
        close(7)

        OPEN(unit=7,file='dr_pdnn.dat')
        write(7,*) '# iAtom, pdgroup, drPdNNAtoms(1:12)'
        do i=1,nPdAtoms
            drlist = drPdNNAtoms(:,i)
            write(7,*) idxPdRbin(i), PdGroup(i), drlist(:)
        enddo
        close(7)

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

        df_pd = 0
    end subroutine calc_pdpd

    ! subroutine to calculate fpd & fh after H initialization
    subroutine calc_h
        use variables
        use fvector
        use eam
        use omp_lib
        implicit none
        integer :: i,j
        real*8  :: dr
        integer,dimension(3) :: idxBox

        allocate(f_pd(nPdAtoms))
        allocate(rho_pdpd(nPdAtoms))
        allocate(phisum_pdpd(nPdAtoms))
        allocate(df_pd(nPdAtoms))

        f_pd = 0.d0
        
    end subroutine calc_h
end module domdec
