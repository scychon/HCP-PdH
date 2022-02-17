!****************************************************************
! 
! this module contains functions to perform H Trajectory analysos
!
!****************************************************************

MODUlE func_anal
    implicit none
 
contains

    ! subroutine to initialize trajectory analysis
    subroutine init_hanal (funit)
        use variables
        implicit none
        integer,intent(in) :: funit

        call calc_pdconc

        !read first.xyz file and initialize to read traj.xyz file
        call init_readH(funit)

        !initialize f, phi, rho, hconc matrices
        call init_traj_mat

    end subroutine init_hanal

    ! subroutine to cacluate Pd concentration map from RCOM
    subroutine calc_pdconc
        use variables
        use fvector
        implicit none
        real*8  :: dr
        integer :: i,idx

        if (.not. allocated(idxPdRbin)) allocate(idxPdRbin(nPdAtoms))
        if (.not. allocated(nPdConcRbin)) allocate(nPdConcRbin(nRbins))
        if (.not. allocated(pPdConcRbin)) allocate(pPdConcRbin(nRbins))

        idxPdRbin = -1
        nPdConcRbin = 0
        pPdConcRbin = 0
        !$OMP PARALLEL &
        !!$OMP   DEFAULT (FIRSTPRIVATE) &
        !$OMP   PRIVATE (i,dr,idx)
        !$OMP DO
        do i=1,nPdAtoms
            dr = norm(PdPos(:,i)-rCOM)
            idx = floor(dr+1)
            idxPdRbin(i) = idx
            !$OMP CRITICAL
            nPdConcRbin(idx) = nPdConcRbin(idx) +1
            pPdConcRbin(idx) = pPdConcRbin(idx) +1.d0/(dr**2)
            !$OMP END CRITICAL
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        OPEN(unit=7,file='pdconc.dat')
        do i=1,nRbins
            !pPdConcRbin(i) =  (1.d0*nPdConcRbin(i))/(1.d0*nPdAtoms* i**2)
            pPdConcRbin(i) =  pPdConcRbin(i)/(1.d0*nPdAtoms)
            write(7,*) i, nPdConcRbin(i), pPdConcRbin(i)
        enddo
        close(7)
        write(*,*) 'pd concentration analysis done'
    end subroutine calc_pdconc

    ! subroutine to initialize matrices related to H positions & coordinates
    subroutine init_readH (funit)
        use variables
        implicit none
        integer, intent(in) :: funit
        integer :: i,iAtom,nskipline
        real*8, dimension(3) :: posnew
        character(len=4) :: atname

        OPEN(unit=funit,file='first.xyz',status='old')
        read(funit, *) nHAtoms
        nHMax = nHAtoms
        if (.not. allocated(HPos)) then
            allocate(HPos(3,nHAtoms))
            allocate(HBox(3,nHAtoms))
            HPos = -1
            HBox = -1
            write(*,*) 'HPos, HBox are allocated', nHAtoms

            allocate(nHCoordPdAtoms(nHAtoms))
            allocate(idxHCoordPdAtoms(nMaxHCoord,nHAtoms))
            allocate(drHCoordPdAtoms(nMaxHCoord,nHAtoms))
            allocate(pHConcRbin(nRbins))
            allocate(nHConcRbin(nRbins))
            allocate(idxHRbin(nHAtoms))
            nHCoordPdAtoms=0
            idxHCoordPdAtoms=-1
            drHCoordPdAtoms=rcut_pdh
            pHConcRbin=0
            nHConcRbin=0
            idxHRbin=-1
        endif
        read(funit,*) !text line
        do iAtom=1,nHAtoms
            read(funit,*) atname, posnew(:)
            HPos(:,iAtom) = posnew(:)
            HBox(:,iAtom) = posnew(:)/(rcut_pdpd) + 1
        enddo
        close(funit)
        write(*,*) 'Read H Pos from first frame first.xyz'

        OPEN(unit=funit,file=strTrajFile,status='old')
        write(*,*) 'Initiallizing to read ',strTrajFile,' trajectory file'
        iFrame = 0
        if (nSkipFrame .gt. 0) then
            nskipline=(nHAtoms+2)*nSkipFrame
            write(*,*) 'Skipping first ',nSkipFrame,' frames.'
            do i=1,nskipline
                read(funit,*)
            enddo
        endif
        write(*,*) 'Skiped first ',nSkipFrame,' frames'

    end subroutine init_readH

    !initialize rho, f, phi, hconc matrices
    subroutine init_traj_mat
        use variables
        implicit none

        ! allocate rho & f matrices
        allocate(nHConcRbinTraj(nRbins,nFrames))
        allocate(nHRbinSum(nRbins))
        allocate(rho_pdh_traj(nPdAtoms,nFrames))
        allocate(rho_hpd_rbin_traj(nRbins, nFrames))
        allocate(rho_hh_rbin_traj(nRbins, nFrames))
        allocate(rho_pdh_avg(nPdAtoms))
        allocate(rho_hpd_rbin_avg(nRbins))
        allocate(rho_hh_rbin_avg(nRbins))
        allocate(f_pd_traj(nPdAtoms,nFrames))
        allocate(f_h_rbin_traj(nRbins, nFrames))
        allocate(f_pd_avg(nPdAtoms))
        allocate(f_h_rbin_avg(nRbins))
        allocate(phisum_pdh_traj(nPdAtoms,nFrames))
        allocate(phisum_pdh_avg(nPdAtoms))
        allocate(phisum_hpd_rbin_traj(nRbins, nFrames))
        allocate(phisum_hh_rbin_traj(nRbins, nFrames))
        allocate(phisum_hpd_rbin_avg(nRbins))
        allocate(phisum_hh_rbin_avg(nRbins))

        allocate(typeHSite(nHAtoms))
        allocate(pdGroupHSite(nHAtoms))
        allocate(phisum_hpd(nHAtoms))

        nHConcRbinTraj       = 0 
        nHRbinSum            = 0 
        rho_pdh_traj         = 0 
        rho_hpd_rbin_traj    = 0 
        rho_hh_rbin_traj     = 0 
        rho_pdh_avg          = 0 
        rho_hpd_rbin_avg     = 0 
        rho_hh_rbin_avg      = 0 
        f_pd_traj            = 0 
        f_h_rbin_traj        = 0 
        f_pd_avg             = 0 
        f_h_rbin_avg         = 0 
        phisum_pdh_traj      = 0 
        phisum_pdh_avg       = 0 
        phisum_hpd_rbin_traj = 0 
        phisum_hh_rbin_traj  = 0 
        phisum_hpd_rbin_avg  = 0 
        phisum_hh_rbin_avg   = 0 

        typeHSite=-1
        pdGroupHSite=-1
        phisum_hpd=0

    end subroutine init_traj_mat

    ! subroutine to read one frame from the H trajectory file
    subroutine readHFrame ( funit, istat )
        use variables
        implicit none
        integer,intent(in)  :: funit
        integer,intent(out) :: istat
        integer :: iAtom,nAtom
        character(len=4) :: atname

        read(funit, *,IOSTAT=istat) nAtom
        if (istat .eq. 0) then
            if (nAtom .ne. nHAtoms) then
                write(*,*) 'ERROR number of atoms ', nAtom, ' does not match with nHAtoms ',nHAtoms
                stop
            endif
            read(funit, *) ! header text
            do iAtom=1,nHAtoms
                read(funit,*) atname, HPos(:,iAtom)
            enddo
        endif

    end subroutine readHFrame
        
    ! subroutine to analize current frame
    subroutine analHFrame
        use variables
        use domdec
        use gen_pos
        implicit none
        integer :: iAtom,f,nAtom
        character(len=256)         :: str, line

        iFrame = iFrame + 1
        call get_HCoord

        call calc_HH

        call update_traj_mat

    end subroutine analHFrame
        
    ! subroutine to analize current frame
    subroutine sumHMap
        use variables
        use domdec
        use gen_pos
        implicit none
        integer :: idx,i,j,k
        integer :: ix,iy
        integer :: dnmax = 20
        real*8  :: drsq, drmaxsq, twosig

        drmaxsq = 20.1**2
        twosig = 10./(2.*rScale)

        !$OMP PARALLEL &
        !$OMP   DEFAULT (FIRSTPRIVATE) &
        !$OMP   SHARED (rhoh_sum)
        !$OMP DO
        do idx=1,ngridvec*ngridvec
            ix=(idx-1)/ngridvec +1
            iy=idx-(ix-1)*ngridvec
            do i=max(1,ix-dnmax),min(ngridvec,ix+dnmax)
              do j=max(1,iy-dnmax),min(ngridvec,iy+dnmax)
                do k=1,nslice
                  if (num_h_grid(i,j,k) > 0) then
                    drsq = (i-ix)**2+(j-iy)**2
                    if (drsq < drmaxsq) then
                        !$OMP CRITICAL
                        rhoh_sum(ix,iy,k) = rhoh_sum(ix,iy,k)+exp(-drsq/(twosig))*num_h_grid(i,j,k)
                        !$OMP END CRITICAL
                        write(*,*) 'test write'
                    endif
                  endif
                enddo
              enddo
            enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
    end subroutine sumHMap
        
    ! subroutine to analize current frame
    subroutine analHMap
        use variables
        use domdec
        use gen_pos
        implicit none
        integer :: iAtom,i
        integer :: ix,iy
        real*8  :: iSlice

        iFrame = iFrame + 1
        do iAtom=1,nHAtoms
            ix = nint(HPos(1,iAtom)/rScale)
            iy = nint(HPos(2,iAtom)/rScale)
            iSlice = HPos(3,iAtom)/rScale
            if (iSlice >= sliceIdxs(1)) then
                num_h_grid(ix,iy,1) = num_h_grid(ix,iy,1)+1
            elseif (iSlice < sliceIdxs(nslice-1)) then
                num_h_grid(ix,iy,nslice) = num_h_grid(ix,iy,nslice)+1
            else
                do i=2,nslice-1
                    if (iSlice < sliceIdxs(i)) then
                        cycle
                    else
                        num_h_grid(ix,iy,i) = num_h_grid(ix,iy,i)+1
                        exit
                    endif
                    write(*,*) 'ERROR : unknown sliceidx ',iSlice
                enddo
            endif
        enddo
    end subroutine analHMap

    ! subroutine to get coordinating Pd atoms for each H atoms
    ! also calculates rho_hpd and phi_hpd
    subroutine get_HCoord
        use variables
        use fvector
        use eam
        implicit none
        integer :: iAtom,jAtom,i,idx
        integer :: iType, iGroup,ibin,iTypeGroup
        integer :: nCoordAtoms,idxNN
        integer,dimension(3) :: idxBox
        logical :: bGroupCoord
        real*8 :: dr,drNN

        if (.not. allocated(nHCoordPdAtoms)) allocate(nHCoordPdAtoms(nHAtoms))

        !construct nearest neighbour list for Pd atoms
        nHCoordPdAtoms = 0
        idxHCoordPdAtoms = -1
        idxHNNAtom = -1
        drHNNAtom = rcut_new
        typeHSite=-1
        pdGroupHSite=-1
        rho_hpd = 0
        phisum_hpd = 0
        !$OMP PARALLEL &
        !$OMP   PRIVATE (i,iAtom,jAtom,dr,idxBox,nCoordAtoms,idxNN,drNN)
        !$OMP DO
        do iAtom=1,nHAtoms
            HBox(:,iAtom) = HPos(:,iAtom)/(rcut_pdpd) + 1
            idxBox(:) = HBox(:,iAtom)
            nCoordAtoms = 0
            drNN = rcut_pdh
            do i=1,nPdNL(idxBox(1),idxBox(2),idxBox(3))
                jAtom = idxNLPd(i,idxBox(1),idxBox(2),idxBox(3))
                dr = norm(PdPos(:,jAtom)-HPos(:,iAtom))
                if (dr .le. rcut_pdh) then
                    rho_hpd(iAtom) = rho_hpd(iAtom)+rhoa_pd(dr)
                    phisum_hpd(iAtom) = phisum_hpd(iAtom)+phi_pdh(dr)
                    if (dr .lt. drNN) then
                        idxNN = jAtom
                        drNN = dr
                    endif
                    if (dr .lt. rcut_hcoord) then
                        nCoordAtoms = nCoordAtoms+1
                        if (nCoordAtoms .gt. nMaxHCoord) then
                            write(*,*) 'ERROR ! Too many coordinating Pd Atoms in ',(iFrame-1)*nStride+nSkipFrame+1,' th frame for ',iAtom,' th H atom'
                            cycle
                        endif
                        idxHCoordPdAtoms(nCoordAtoms,iAtom) = jAtom
                    endif
                endif
            enddo
            nHCoordPdAtoms(iAtom) = nCoordAtoms
            idxHNNAtom(iAtom) = idxNN
            drHNNAtom(iAtom) = drNN
            phisum_hpd(iAtom) = phisum_hpd(iAtom)*0.5d0     ! put only half of the phi_hpd since we will separately calculate phi_pdh

            if (nCoordAtoms .lt. 4) then
                if (bPdSurfAtom(idxNN)) then
                    typeHSite(iAtom) = 4        ! surface H atom
                else
                    typeHSite(iAtom) = 6        ! etc
                endif
            elseif (nCoordAtoms .eq. 4) then
                typeHSite(iAtom) = 2            ! Tetrahedral site
            elseif (nCoordAtoms .eq. 5) then
                typeHSite(iAtom) = 3            ! In between of Tetrahedral and Octahedral site
            elseif (nCoordAtoms .eq. 6) then
                typeHSite(iAtom) = 1            ! Octahedral site
            elseif (nCoordAtoms .gt. 6) then
                typeHSite(iAtom) = 5            ! Crowded site (7 or 8 coordinating Pd atoms)
            endif

            if (bUsePdGroup) then
                pdGroupHSite(iAtom) = PdGroup(idxNN)
            endif
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

!    ! subroutine to assign Htype matrices
!    subroutine stat_htype
        nHType = 0
        nHPdGroupNN = 0
        nHPdGroupCoord = 0
        do iAtom=1,nHAtoms
            iType = typeHSite(iAtom)
            nHType(iType) = nHType(iType)+1
            if (bUsePdGroup) then
                iGroup = pdGroupHSite(iAtom)
                if ((iGroup ==1) .or.(iGroup==2)) then
                    nHPdGroupNN(iGroup) = nHPdGroupNN(iGroup)+1
                    iTypeGroup = 3+(iGroup-1)*6+iType
                    nHPdGroupNN(iTypeGroup) = nHPdGroupNN(iTypeGroup)+1
                    bGroupCoord = .TRUE.
                    do i=1,min(nHCoordPdAtoms(iAtom),8)
                        if( PdGroup(idxHCoordPdAtoms(i,iAtom)).ne.iGroup) bGroupCoord = .FALSE.
                    enddo
                    if (bGroupCoord) then
                        nHPdGroupCoord(iGroup) = nHPdGroupCoord(iGroup)+1
                        nHPdGroupCoord(iTypeGroup) = nHPdGroupCoord(iTypeGroup)+1
                    else
                        nHPdGroupCoord(3) = nHPdGroupCoord(3)+1
                    endif
                else
                    nHPdGroupNN(3) = nHPdGroupNN(3)+1
                    nHPdGroupCoord(3) = nHPdGroupCoord(3)+1
                endif
            endif
        enddo
!    end subroutine stat_htype

    end subroutine get_HCoord

    ! subroutine to calculate HH interactions
    ! first construct HH neighbor list, then loop to get rho_hh and phi_hh
    subroutine calc_HH
        use variables
        use fvector
        use eam
        implicit none
        integer :: i,j
        real*8  :: dr
        real*8  :: rhoval,phival

        !$OMP PARALLEL &
        !$OMP   PRIVATE (i,j,dr,rhoval,phival)
        !$OMP DO
        do i=1,nPdAtoms
            rhoval=0
            phival=0
            do j=1,nHAtoms
                dr = norm(PdPos(:,i)-HPos(:,j))
                if (dr .le. rcut_pdh) then
                    rhoval = rhoval + rhoa_h(dr)
                    phival = phival + phi_pdh(dr)
                endif
            enddo
            rho_pdh(i) = rhoval
            phisum_pdh(i) = 0.5d0*phival
            f_pd(i) = fpd(rho_pdpd(i)+rho_pdh(i))
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        !$OMP PARALLEL &
        !$OMP   PRIVATE (i,j,dr,rhoval,phival)
        !$OMP DO
        do i=1,nHAtoms
            rhoval=0
            phival=0
            do j=1,nHAtoms
                if(i==j) cycle
                dr = norm(HPos(:,i)-HPos(:,j))
                if (dr .le. rcut_hh) then
                    rhoval = rhoval + rhoa_h(dr)
                    phival = phival + phi_h(dr)
                endif
            enddo
            rho_hh(i) = rhoval
            phisum_hh(i) = 0.5d0*phival
            f_h(i) = fh(rho_hpd(i)+rho_hh(i))
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

    end subroutine calc_HH

    ! subroutine to update trajectory matrices
    ! also perform concentration mapping
    subroutine update_traj_mat
        use variables
        use fvector
        use eam
        implicit none
        integer :: i,j,idx
        real*8  :: dr
        real*8  :: rhoval,phival
        real*8,dimension(nRbins) :: rho_pdpd_rbin, rho_pdh_rbin, rho_hpd_rbin, rho_hh_rbin, f_h_rbin, f_pd_rbin
        real*8,dimension(nRbins) :: phisum_pdh_rbin, phisum_pdpd_rbin, phisum_hpd_rbin, phisum_hh_rbin

        if (iFrame .gt. nFrames) then
            call expand2Dint(nHConcRbinTraj,0,nFrames)
            call expand2D(rho_pdh_traj,0,nFrames)
            call expand2D(rho_hpd_rbin_traj,0,nFrames)
            call expand2D(rho_hh_rbin_traj,0,nFrames)
            call expand2D(f_pd_traj,0,nFrames)
            call expand2D(f_h_rbin_traj,0,nFrames)
            call expand2D(phisum_pdh_traj,0,nFrames)
            call expand2D(phisum_hpd_rbin_traj,0,nFrames)
            call expand2D(phisum_hh_rbin_traj,0,nFrames)
            nFrames = 2*nFrames
        endif

        idxHRbin = -1
        nHConcRbin = 0
        !$OMP PARALLEL &
        !$OMP   PRIVATE (i,dr,idx)
        !$OMP DO
        do i=1,nHAtoms
            dr = norm(HPos(:,i)-rCOM)
            idx = floor(dr+1)
            idxHRbin(i) = idx
            !$OMP CRITICAL
            nHConcRbin(idx) = nHConcRbin(idx) +1
            !$OMP END CRITICAL
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        nHConcRbinTraj(:,iFrame)    = nHConcRbin(:)
        nHRbinSum(:)                = nHRbinSum(:) + nHConcRbin(:)

        rho_pdh_traj(:,iFrame)      = rho_pdh(:)
        rho_pdh_avg(:)              = rho_pdh_avg(:)+rho_pdh(:)
        f_pd_traj(:,iFrame)         = f_pd(:)
        f_pd_avg(:)                 = f_pd_avg(:) + f_pd(:)
        phisum_pdh_traj(:,iFrame)   = phisum_pdh(:)
        phisum_pdh_avg(:)           = phisum_pdh_avg(:) + phisum_pdh(:)

        rho_hpd_rbin = 0
        rho_hh_rbin = 0
        f_h_rbin = 0
        phisum_hpd_rbin = 0
        phisum_hh_rbin = 0
        do i=1,nHAtoms
            idx=idxHRbin(i)
            rho_hpd_rbin(idx) = rho_hpd_rbin(idx) + rho_hpd(i)
            rho_hh_rbin(idx) = rho_hh_rbin(idx) + rho_hh(i)
            f_h_rbin(idx) = f_h_rbin(idx) + f_h(i)
            phisum_hpd_rbin(idx) = phisum_hpd_rbin(idx) + phisum_hpd(i)
            phisum_hh_rbin(idx) = phisum_hh_rbin(idx) + phisum_hh(i)
        enddo

        rho_pdh_rbin = 0
        rho_pdpd_rbin = 0
        f_pd_rbin = 0
        phisum_pdpd_rbin = 0
        phisum_pdh_rbin = 0
        do i=1,nPdAtoms
            idx=idxPdRbin(i)
            rho_pdh_rbin(idx) = rho_pdh_rbin(idx) + rho_pdh(i)
            rho_pdpd_rbin(idx) = rho_pdpd_rbin(idx) + rho_pdpd(i)
            f_pd_rbin(idx) = f_pd_rbin(idx) + f_pd(i)
            phisum_pdh_rbin(idx) = phisum_pdh_rbin(idx) + phisum_pdh(i)
            phisum_pdpd_rbin(idx) = phisum_pdpd_rbin(idx) + phisum_pdpd(i)
        enddo

        write(20,*) ''
        do i=1,nRbins
            ftot_pd = f_pd_rbin(i)
            ftot_h = f_h_rbin(i)
            ftot = ftot_h + ftot_pd
            Etot_pd = ftot_pd + phisum_pdpd_rbin(i) + phisum_pdh_rbin(i)
            Etot_h = ftot_h + phisum_hpd_rbin(i) + phisum_hh_rbin(i)
            write(20,*) i,nPdConcRbin(i),nHConcRbin(i), Etot_pd, Etot_h, (Etot_pd+Etot_h), &
                    rho_pdpd_rbin(i)+rho_pdh_rbin(i), rho_hpd_rbin(i)+rho_hh_rbin(i), rho_pdpd_rbin(i), rho_pdh_rbin(i), rho_hpd_rbin(i), rho_hh_rbin(i), &
                    f_pd_rbin(i), f_h_rbin(i), phisum_pdpd_rbin(i)+phisum_pdh_rbin(i), phisum_hpd_rbin(i)+phisum_hh_rbin(i), &
                    phisum_pdpd_rbin(i), phisum_pdh_rbin(i), phisum_hpd_rbin(i), phisum_hh_rbin(i)
        enddo
        flush(20)

        rho_hpd_rbin_traj(:,iFrame)   = rho_hpd_rbin(:)
        rho_hpd_rbin_avg(:)           = rho_hpd_rbin_avg(:) + rho_hpd_rbin(:)
        rho_hh_rbin_traj(:,iFrame)   = rho_hh_rbin(:)
        rho_hh_rbin_avg(:)           = rho_hh_rbin_avg(:) + rho_hh_rbin(:)
        f_h_rbin_traj(:,iFrame)   = f_h_rbin(:)
        f_h_rbin_avg(:)           = f_h_rbin_avg(:) + f_h_rbin(:)
        phisum_hpd_rbin_traj(:,iFrame)   = phisum_hpd_rbin(:)
        phisum_hpd_rbin_avg(:)           = phisum_hpd_rbin_avg(:) + phisum_hpd_rbin(:)
        phisum_hh_rbin_traj(:,iFrame)   = phisum_hh_rbin(:)
        phisum_hh_rbin_avg(:)           = phisum_hh_rbin_avg(:) + phisum_hh_rbin(:)

        ! sum f, phi terms
        ftot_pd = sum(f_pd)
        ftot_h = sum(f_h)
        ftot = ftot_pd+ftot_h
        phitot_pdh = sum(phisum_pdh)
        phitot_hpd = sum(phisum_hpd)
        phitot_hh = sum(phisum_hh)
        if ((phitot_pdh-phitot_hpd)**2 .gt. .01) then
            write(*,*) 'ERROR phitot_pdh and phitot_hpd does not match ! ',phitot_pdh, phitot_hpd
        endif
        phitot = phitot_pdpd + phitot_pdh + phitot_hpd + phitot_hh
        Etot_pd = ftot_pd + phitot_pdpd + phitot_pdh
        Etot_h = ftot_h + phitot_hpd + phitot_hh
        E_pd = Etot_pd/(1.d0*nPdAtoms)
        E_h = Etot_h/(1.d0*nHAtoms)
        Etot = ftot+phitot

    end subroutine update_traj_mat

    ! subroutine to wrap up analysis and print final results
    subroutine finalize_anal
        use variables
        use fvector
        use eam
        implicit none
        integer :: i,j,idx
        real*8  :: dr
        real*8,dimension(nRbins) :: rho_pdh_rbin_avg, f_pd_rbin_avg, rho_pd_rbin_avg

        rho_pdh_avg(:)              = rho_pdh_avg(:)/(1.d0*iFrame)
        f_pd_avg(:)                 = f_pd_avg(:) /(1.d0*iFrame)
        phisum_pdh_avg(:)           = phisum_pdh_avg(:) /(1.d0*iFrame)
        do i=1,nPdAtoms
            idx = idxPdRbin(i)
            
        enddo

        rho_hpd_rbin_avg(:)           = rho_hpd_rbin_avg(:) / (1.d0*nHRbinSum(:))
        rho_hh_rbin_avg(:)           = rho_hh_rbin_avg(:) / (1.d0*nHRbinSum(:))
        f_h_rbin_avg(:)           = f_h_rbin_avg(:) / (1.d0*nHRbinSum(:))
        phisum_hpd_rbin_avg(:)           = phisum_hpd_rbin_avg(:) / (1.d0*nHRbinSum(:))
        phisum_hh_rbin_avg(:)           = phisum_hh_rbin_avg(:) / (1.d0*nHRbinSum(:))


        !OPEN(unit=7,file='rbin_all.dat')
        !write(7,*) '# ibin, nPdConc, pPdConc, rho_pd_avg, rho_pdh_avg, rho_h_avg, rho_hpd_avg, rho_hh_avg, f_h_avg, f_pd_avg, '
        !do i=1,nRbins
        !    !pPdConcRbin(i) =  (1.d0*nPdConcRbin(i))/(1.d0*nPdAtoms* i**2)
        !    pPdConcRbin(i) =  pPdConcRbin(i)/(1.d0*nPdAtoms)
        !    write(7,*) i, nPdConcRbin(i), pPdConcRbin(i), rho_pdh_avg_rbin_avg
        !enddo
        !close(7)

    end subroutine finalize_anal

end module func_anal
