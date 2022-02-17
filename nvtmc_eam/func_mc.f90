!****************************************************************
! 
! this module contains functions to perform MC moves
!
!****************************************************************

MODUlE func_mc
    implicit none
 
contains

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

    ! Subroutine to get command arguments and assign input and output file name
    subroutine metropolis (Ediff, bAccept, iAddPtl)
        use variables
        implicit none
        real*8, intent(in)   :: Ediff
        logical, intent(out) :: bAccept
        integer, intent(in)  :: iAddPtl

        real*8 :: rd,compval

        if (iAddPtl .eq. 1) then
            compval = betamu + log(nqvol/(nHAtoms+1))
        elseif (iAddPtl .eq. -1) then
            compval = -betamu + log(nHAtoms/nqvol)
        elseif (iAddPtl .eq. 0) then
            compval = 0
        else
            write(*,*) 'ERROR: unknown MC attempt, iAddPtl : ',iAddPtl
            stop
        endif
        compval = -Ediff/kT + compval

        if (compval .gt. 0) then
            bAccept=.true.
        else
            call random_number(rd)
            if (rd .lt. exp(compval)) then
                bAccept=.true.
            else
                bAccept=.false.
            endif
        endif
    end subroutine


    ! subroutine to do exchange move
    subroutine mcmove (bAccept)
        use variables
        use fvector
        use gen_pos
        use eam
        implicit none
        logical, intent(out) :: bAccept     !whether the mcexch is accepted or not

        logical :: bOverlap

                                            !neighbor list (NL) runs over 27(3*3*3) cells with the chosen box at center
                                            !each box has dimension of rcut_pdpd**3
        integer :: i,j,k,l,idx              !counters
        integer :: idxHMove                 !index of H atom to be moved (idx in hpos)
        integer :: idxHMoveInBox            !index of H atom to be moved (idx in box)
        integer :: idxH                     !index to be looped over all hydrogen neighbor list (HNL)
        integer :: idxPd                    !index to be looped over all Pd neighbor list (PdNL)
        integer :: numPdNL                  !number of Pd atoms within NL
        integer :: numHNL                   !number of H atoms within NL
        integer :: numH                     !number of H atoms within chosen box


        real*8  :: drold                    !temp variable to save old distance
        real*8  :: drnew                    !temp variable to save new distance
        real*8  :: drmin                    !minimum distance to nearest pd atom
        real*8  :: rdH                      !temp variable to save random number for h index
        real*8  :: rhoahval                 !temp variable to save value of rho^a_h
        real*8  :: dphival                  !temp variable to save value of phi_a
        real*8  :: rhoself_hpd              !value of rho_hpd for the selected H atom (contributions from Pd in numPdNL)
        real*8  :: rhoself_hh               !value of rho_hh for the selected H atom (contributions from H in numHNL)
        real*8  :: phiself_hh               !value of phi_hh for the selected H atom (contributions from H in numHNL)
        real*8  :: Ediff                    !energy difference after mc move

        real*8  :: dphitotbox_pdh, dphitotbox_hh, dphitotbox

        real*8  :: df_sum
        integer,dimension(3) :: idxBox
        integer,allocatable,dimension(:) :: idxNLPdbox              !list of indices of Pd atoms in NL (numPdNL)
        integer,allocatable,dimension(:) :: idxNLHbox               !list of indices of H atoms in NL (numHNL)

        real*8, dimension(3) :: rdvec                   !temp vector to save random numbers (3)
        real*8, dimension(3) :: posold                  !old position of selected ptl to be moved
        real*8, dimension(3) :: posnew                  !new position of selected ptl to be moved
        real*8, dimension(3) :: pos0                    !position of origin of chosen box

        real*8, allocatable,dimension(:) :: rhobox_pdh              !array of rho_pdh in NL (numPdNL)
        real*8, allocatable,dimension(:) :: rhobox_hh               !array of rho_hh in NL (numHNL); last index for newly inserted H atom

        real*8, allocatable,dimension(:) :: fbox_pd,fbox_h
        real*8, allocatable,dimension(:) :: phibox_pdh,phibox_hh

        integer :: chknum, chknum2, iadd
        real*8 :: phitotbak, phitotbak_pdh, phitotbak_hh

        bAccept = .false.
        call random_number(rdH)
        idxHMove = int(rdH*nHAtoms)+1
        posold = HPos(:,idxHMove)

        ! get new position
        do i=1,10
            call get_new_HPos(posnew)
            call checkOverlapAll(posnew,bOverlap,idxHMove)
            if (.not. bOverlap) exit
        enddo
        if (bOverlap) return                ! return immediately if new Position overlaps with other atom

        phitotbak_pdh= sum(phisum_pdh)
        phitotbak_hh = sum(phisum_hh)
        phitotbak = phitotbak_pdh+phitotbak_hh
        drho_hh = 0
        drho_hpd = 0
        drho_pdh = 0
        dphisum_pdh=0
        dphitot_pdh=0

        !call get_NL_H
        ! calculate rhophi due to previous H atom position
        call calc_rhophipdh(posold,idxHMove,.false.)
        !dphitotbox_pdh = dphitot_pdh

        ! calculate rhophi due to new H atom position
        call calc_rhophipdh(posnew,idxHMove,.true.)
        !dphitot_pdh = dphitot_pdh + dphitotbox_pdh

        call calc_rhophihh(posnew,idxHMove)
        dphitot = dphitot_pdh + dphitot_hh

        ! calculate changes in f_pd and f_h
        call calc_df(df_sum)

        !calculate sum of dF and dPhi 
        Ediff = df_sum + dphitot
        !Do Metropolis check on new energy; iAddPtl = 0 for particle move
        call metropolis(Ediff,bAccept, 0)

        !write(*,*) Ediff,df_sum,dphitot,dphitot_pdh, dphitot_hh, bAccept
        if (.not. bAccept) return

        !if accepted, update all corresponding matrices
        if (bAccept) then
            call update_newH(posold,posnew,idxHMove)
        endif
        phitot_pdh = sum(phisum_pdh)
        phitot_hh = sum(phisum_hh)
        phitot = phitot_pdh+phitot_hh
        if (abs(phitot-phitotbak-dphitot).gt.0.001) then
            write(*,*) 'pitot', phitot, phitot_pdh, phitot_hh,phitotbak,phitotbak_pdh,phitotbak_hh,phitotbak+dphitot,phitotbak_pdh+dphitot_pdh,phitotbak_hh+dphitot_hh
        endif
    end subroutine mcmove

    !calc changes in rho,phi given the idxHAtom
    subroutine calc_rhophipdh (pos, idxHMove, bAdd)
        use variables
        use fvector
        use gen_pos
        use eam
        implicit none
        real*8,dimension(3), intent(in) :: pos
        integer,intent(in)   :: idxHMove
        logical,intent(in)   :: bAdd
        integer,dimension(3) :: boxidx

        integer :: numPdNL, numHNL  !number of Pd-neighbour, H-neighbour, H-in-box
        integer :: iscale           !1 for addition -1 for deletion
        integer :: i,j,k            !counters
        integer :: idxPd, idxNNPd,idxH
        real*8  :: dr
        real*8  :: drmin                    !minimum distance to nearest pd atom
        real*8  :: dphival                  !temp variable to save value of phi_a
        real*8  :: rhoahval                 !temp variable to save value of rho^a_h
        real*8  :: dphitotcheck

        if (bAdd) then
            iscale = 1
        else
            iscale = -1
        endif

        boxidx(:) = pos(:)/(rcut_pdpd) + 1
        numPdNL = nPdNL(boxidx(1),boxidx(2),boxidx(3))
        numHNL = nHNL(idxHMove)
        if (bAdd) numHNL = nHAtoms

        !dphitot_pdh=0
        dphitotcheck=0

        drmin = rcut_new
        idxNNPd = -1
        !Calculate rho_pd and f_pd for numPdNL neighbors
        do i=1,numPdNL
          idxPd = idxNLPd(i,boxidx(1),boxidx(2),boxidx(3))
          rhoahval = 0.d0
          dphival = 0.d0
        
          !subtract rho & phi values for the old position
          dr=norm(PdPos(:,idxPd)-pos)
          if ( dr .le. rcut_pdh ) then
            if (dr .le. rcut_min) write(*,*) 'ERROR IN DR: ',dr,' is smaller than rcut'
            drho_pdh(idxPd) = drho_pdh(idxPd) + iscale*rhoa_h(dr)
            dphival = iscale*phi_pdh(dr)
            dphisum_pdh(idxPd) = dphisum_pdh(idxPd) + dphival
            if (bAdd) then
                drho_hpd(idxHMove) = drho_hpd(idxHMove) + iscale*rhoa_pd(dr)
                if (dr .lt. drmin) then
                    drmin = dr
                    idxNNPd = idxPd
                endif
            endif
          endif

          if ((drho_pdh(idxPd)+rho_pdh(idxPd)).lt.-0.001) then
            write(*,*) 'ERROR DRHOPDH PDH PDH PDH ',i,idxPd,idxHMove,drho_pdh(idxPd),rho_pdh(idxPd),idxNLPd(i,boxidx(1),boxidx(2),boxidx(3)),'test'
          endif
          !update dphitot
          dphitot_pdh = dphitot_pdh + dphival
        enddo
        if (bAdd) then
            idxHNNAtom(nHAtoms+1) = idxPd        ! last index for trial insertion
            drHNNAtom(nHAtoms+1) = dr
        else
            drho_hpd(idxHMove) = iscale*rho_hpd(idxHMove)
        endif

        if (abs(sum(dphisum_pdh)-dphitot_pdh).gt.0.001) then
            write(*,*) 'ERROR DPHISUM_PDH',sum(dphisum_pdh),dphitot_pdh
        endif

    end subroutine calc_rhophipdh


    subroutine calc_rhophihh (posnew, idxHMove)
        use variables
        use fvector
        use gen_pos
        use eam
        implicit none
        real*8,dimension(3), intent(in) :: posnew
        integer,intent(in)   :: idxHMove
        integer,dimension(3) :: boxidx

        integer :: numPdNL, numHNL  !number of Pd-neighbour, H-neighbour, H-in-box
        integer :: iscale           !1 for addition -1 for deletion
        integer :: i,j,k            !counters
        integer :: idxPd, idxNNPd,idxH
        real*8  :: drnew,drold
        real*8  :: drmin                    !minimum distance to nearest pd atom

        drho_hh = 0
        dphisum_hh = 0
        dphitot_hh = 0
        
        !Calculate rho_h for numHNL neighbors
        do idxH=1,nHAtoms
          if (idxHMove==idxH) cycle        ! skip self interaction

          !subtract rho & phi values for the old position
          drnew=norm(HPos(:,idxH)-posnew)
          if ( drnew .le. rcut_hh ) then
              drho_hh(idxH) = rhoa_h(drnew)
              dphisum_hh(idxH) = phi_h(drnew)
              drho_hh(idxHMove) = drho_hh(idxHMove) + rhoa_h(drnew)           !change of rho_hh for the selected ptl
          endif
          drold=norm(HPos(:,idxH)-HPos(:,idxHMove))
          if ( drold .le. rcut_hh ) then
              drho_hh(idxH) = drho_hh(idxH)-rhoa_h(drold)
              dphisum_hh(idxH) = dphisum_hh(idxH) - phi_h(drold)
              drho_hh(idxHMove) = drho_hh(idxHMove) - rhoa_h(drold)           !change of rho_hh for the selected ptl
          endif

          if ((drho_hh(idxH)+rho_hh(idxH).lt.-0.0001)) then
            write(*,*) 'ERROR DRHOHH HH HH HH ',i,idxH,idxHMove,drho_hh(idxH),rho_hh(idxH),idxNLH(:,idxHMove),'test'
          endif
          !update dphitot
        enddo
        dphitot_hh = sum(dphisum_hh)

    end subroutine calc_rhophihh

    !calculate changes in f after calc_rhophibox
    subroutine calc_df (dftot)
        use variables
        use eam
        implicit none
        integer :: i            !counters
        real*8,intent(out) :: dftot

        dftot = 0
        do i=1,nPdAtoms
            if (drho_pdh(i) .ne. 0) then
                df_pd(i) = fpd(rho_pdh(i)+drho_pdh(i)+rho_pdpd(i))-f_pd(i)
                dftot = dftot + df_pd(i)
            endif
        enddo
        !write(*,*) 'DFTot : ',dftot
        do i=1,nHAtoms
            if ((drho_hpd(i) .ne. 0) .or. (drho_hh(i) .ne. 0)) then
                if (rho_hpd(i)+drho_hpd(i)+rho_hh(i)+drho_hh(i) .le. 0) then
                    write(*,*) 'ERROR IN RHO_H : ',i,rho_hpd(i),drho_hpd(i),rho_hh(i),drho_hh(i)
                    stop
                endif
                df_h(i) = fh(rho_hpd(i)+drho_hpd(i)+rho_hh(i)+drho_hh(i))-f_h(i)
                dftot = dftot + df_h(i)
            endif
        enddo
        !write(*,*) 'DFTotH : ',dftot
    end subroutine calc_df

    !update H position to newly selected one
    subroutine update_newH (posold, posnew, idxHMove)
        use variables
        use gen_pos
        implicit none
        real*8,dimension(3), intent(in) :: posold,posnew
        integer,intent(in)   :: idxHMove
        integer,dimension(3) :: boxidx
        integer :: i
        integer :: numH,idxDomain
        real*8  :: rhpd,rhh

        do i=1,nPdAtoms
            if (drho_pdh(i) .ne. 0) then
                rho_pdh(i) = rho_pdh(i) + drho_pdh(i)
                drho_pdh(i) = 0
                f_pd(i) = f_pd(i) + df_pd(i)
                df_pd(i) = 0
                phisum_pdh(i) = phisum_pdh(i) + dphisum_pdh(i)
                dphisum_pdh(i) = 0
            endif
        enddo
          
        do i=1,nHAtoms 
            if ((drho_hpd(i) .ne. 0) .or. (drho_hh(i) .ne. 0)) then
                rhpd = rho_hpd(i)
                rhh = rho_hh(i)
                rho_hpd(i) = rho_hpd(i) + drho_hpd(i)
                rho_hh(i) = rho_hh(i) + drho_hh(i)
                if ((rho_hpd(i) .lt. -0.001) .or. (rho_hh(i) .lt. -0.001)) then
                    write(*,*) 'ERROR HPD HH : ',i, rho_hpd(i), rho_hh(i), drho_hpd(i), drho_hh(i),rhpd,rhh
                    stop
                endif
                drho_hpd(i) = 0
                drho_hh(i) = 0
                f_h(i) = f_h(i) + df_h(i)
                df_h(i) = 0
                phisum_hh(i) = phisum_hh(i) + dphisum_hh(i)
                dphisum_hh(i) = 0
            endif
        enddo   

        !delete neighbour list of old atom
        boxidx(:) = posold(:)/rcut_pdpd+1
        numH = nHBox(boxidx(1),boxidx(2),boxidx(3))
        do i=1,numH
            if (idxHMove==idxHAtomInBox(i,boxidx(1),boxidx(2),boxidx(3))) then
                idxHAtomInBox(i,boxidx(1),boxidx(2),boxidx(3))=idxHAtomInBox(numH,boxidx(1),boxidx(2),boxidx(3))
                idxHAtomInBox(numH,boxidx(1),boxidx(2),boxidx(3))=-1
                exit
            elseif (i==numH) then
                write(*,*) 'ERROR idxHAtomInBox DOES NOT HAVE idxHMOVE',idxHMove,idxHAtomInBox(1:numH,boxidx(1),boxidx(2),boxidx(3))
            endif
        enddo
        numH = numH-1
        idxDomain = ((boxidx(1)-1)/3)*(nDomainVec(2)*nDomainVec(3)) + ((boxidx(2)-1)/3)*nDomainVec(3) + ((boxidx(3)-1)/3) + 1
        nHBox(boxidx(1),boxidx(2),boxidx(3))=numH
        nHDomain(idxDomain) = nHDomain(idxDomain)-1
        !call updateNL(boxidx,idxHMove,1)

        !add neighbour list of new atom
        boxidx(:) = posnew(:)/rcut_pdpd+1
        numH = nHBox(boxidx(1),boxidx(2),boxidx(3))
        idxHAtomInBox(numH+1,boxidx(1),boxidx(2),boxidx(3))=idxHMove
        nHBox(boxidx(1),boxidx(2),boxidx(3))=numH+1
        idxDomain = ((boxidx(1)-1)/3)*(nDomainVec(2)*nDomainVec(3)) + ((boxidx(2)-1)/3)*nDomainVec(3) + ((boxidx(3)-1)/3) + 1
        nHDomain(idxDomain) = nHDomain(idxDomain)+1
        call update_NL_H(idxHMove)
        !call updateNL(boxidx,idxHMove,-1)

        !update HBox, HPos, idxHNNAtom
        HPos(:,idxHMove) = posnew(:)
        HBox(:,idxHMove) = boxidx(:)
        idxHNNAtom(idxHMove) = idxHNNAtom(nHAtoms+1)
        drHNNAtom(idxHMove) = drHNNAtom(nHAtoms+1)
        idxHNNAtom(nHAtoms+1) = -1
        drHNNAtom(nHAtoms+1) = rcut_new
        !call get_NL_H

    end subroutine update_newH
end module func_mc
