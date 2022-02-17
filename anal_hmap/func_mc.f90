!****************************************************************
! 
! this module contains functions to perform MC moves
!
!****************************************************************

MODUlE func_mc
    implicit none
 
contains

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
    subroutine mcmove (boxidx, bAccept)
        use variables
        use fvector
        use eam
        implicit none
        integer,dimension(3), intent(in) :: boxidx
        logical, intent(out) :: bAccept     !whether the mcexch is accepted or not

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

        real*8  :: df_sum, dphitot
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

        chknum = sum(nHNL)
        iadd = 0

        bAccept = .false.
        pos0(:)=(boxidx(:)-1)*rcut_pdpd

        numPdNL = nPdNL(boxidx(1),boxidx(2),boxidx(3))
        numHNL = nHNL(boxidx(1),boxidx(2),boxidx(3))
        numH = nHBox(boxidx(1),boxidx(2),boxidx(3))

        allocate(rhobox_pdh(numPdNL))
        allocate(rhobox_hh(numHNL))
        allocate(phibox_pdh(numPdNL))
        allocate(phibox_hh(numHNL))
        allocate(idxNLPdbox(numPdNL))
        allocate(idxNLHbox(numHNL))

        allocate(fbox_pd(numPdNL))
        allocate(fbox_h(numHNL))

        idxNLPdbox(:)=idxNLPd(1:numPdNL,boxidx(1),boxidx(2),boxidx(3))
        idxNLHbox(:)=idxNLH(1:numHNL,boxidx(1),boxidx(2),boxidx(3))

        do i=1,numHNL
          if (idxNLH(i,boxidx(1),boxidx(2),boxidx(3)).le.0) then
            write(*,*) 'ERROR idxNLH NOT ASSIGNED : ',numHNL, i, idxNLH(i,boxidx(1),boxidx(2),boxidx(3)), idxNLH(:,boxidx(1),boxidx(2),boxidx(3))
            write(*,*) 'numHNL, i, idxNLH(i), idxNLH'
            write(*,*) boxidx
            stop
          endif
        enddo
        do i=numHNL+1,nHNLMax
          if (idxNLH(i,boxidx(1),boxidx(2),boxidx(3)).ne.-1) then
            write(*,*) 'ERROR idxNLH assigned after numHNL : ',numHNL, i, idxNLH(i,boxidx(1),boxidx(2),boxidx(3)), idxNLH(:,boxidx(1),boxidx(2),boxidx(3))
            write(*,*) 'numHNL, i, idxNLH(i), idxNLH'
            write(*,*) boxidx
            stop
          endif
        enddo

        rhobox_pdh=0
        rhobox_hh=0
        rhoself_hpd=0
        rhoself_hh=0
        phibox_pdh=0
        phibox_hh=0
        phiself_hh=0
        fbox_pd=0
        fbox_h=0

        df_sum=0
        dphitot=0


        ! choose h atom to move
        call random_number(rdH)
        idxHMoveInBox = int(rdH*numH)+1
        idxHMove = idxHAtomInBox(idxHMoveInBox,boxidx(1),boxidx(2),boxidx(3))
        if (idxHMove .lt. 0) then
            write(*,*) 'error idxHMove : ', boxidx,idxHMoveInBox, idxHMove, numH
            write(*,*) 'boxidx, idxHMoveInBox, idxHMove, numH'
        endif
        posold(:)=HPos(:,idxHMove)

        ! get new position
        !posnew(:)=pos0 + rdvec(:)*rcut_pdpd
        do i=1,10
            call random_number(rdvec)
            posnew(:)=posold(:) + rdvec(:)*rcut_mcmove
            if (all(boxidx(:).eq.(int(posnew(:)/(rcut_pdpd)+1)))) exit
        enddo
        !write(*,*) boxidx, posnew, int(posnew(:)/(rcut_pdpd)+1)
        if (any(boxidx(:).ne.(int(posnew(:)/(rcut_pdpd)+1)))) return        ! return immediately if new postion is outside of box

        !write(*,*) 'box set is successful'
        

        drmin = rcut_pdh+1
        !Calculate rho_pd and f_pd for numPdNL neighbors
        do i=1,numPdNL
          idxPd = idxNLPdbox(i)
          rhoahval = 0.d0
          dphival = 0.d0
        
          !subtract rho & phi values for the old position
          drold=norm(PdPos(:,idxPd)-posold)
          if ( drold .le. rcut_pdh ) then
            rhoahval = -rhoa_h(drold)
            dphival = -phi_pdh(drold)
          endif

          !add rho & phi values for the new position
          drnew=norm(PdPos(:,idxPd)-posnew)
          if ( drnew .le. rcut_pdh ) then
            if ( drnew .lt. rcut_min ) return                    ! return immediately if too close to existing atoms
            rhoahval = rhoahval + rhoa_h(drnew)
            dphival  = dphival  + phi_pdh(drnew)
            rhoself_hpd = rhoself_hpd + rhoa_pd(drnew)
            if ( drnew .lt. drmin) drmin = drnew
          endif

          !update rhobox_pdh & phibox_pdh
          rhobox_pdh(i)=rho_pdh(idxPd)+rhoahval
          phibox_pdh(i)=phisum_pdh(idxPd)+dphival
          dphitot = dphitot + dphival

          fbox_pd(i)=fpd(rhobox_pdh(i)+rho_pdpd(idxPd))
          df_sum = df_sum + fbox_pd(i)-f_pd(idxPd)
        enddo

        if (drmin .gt. rmin_pdh) return                       ! return immediately if no pd atom is within the cutoff

        !Calculate rho_h for numHNL neighbors
        do i=1,numHNL
          idxH = idxNLHbox(i)
          if (idxHMove==idxH) cycle        ! skip self interaction

          rhoahval = 0.d0
          dphival = 0.d0

          !subtract rho & phi values for the old position
          drold=norm(HPos(:,idxH)-posold)
          if ( drold .le. rcut_hh ) then
              rhoahval = -rhoa_h(drold)
              dphival = -phi_h(drold)
          endif

          !add rho & phi values for the new position
          drnew=norm(HPos(:,idxH)-posnew)
          if ( drnew .le. rcut_hh ) then
            if ( drnew .lt. rcut_min ) return                    ! return immediately if too close to existing atoms
            rhoahval = rhoahval + rhoa_h(drnew)
            dphival  = dphival  + phi_h(drnew)
            rhoself_hh = rhoself_hh + rhoa_h(drnew)           !change of rho_hh for the selected ptl
            phiself_hh = phiself_hh + phi_h(drnew)           !change of phi_hh for the selected ptl
          endif

          !update rhobox_hh & phibox_h
          rhobox_hh(i)=rho_hh(idxH)+rhoahval
          phibox_hh(i)=phisum_hh(idxH)+dphival
          dphitot = dphitot + dphival

          fbox_h(i)=fh(rhobox_hh(i)+rho_hpd(idxH))
          df_sum = df_sum + fbox_h(i)-f_h(idxH)
        enddo

        !save the rho_hh & f_h value for the selected H
        rhobox_hh(idxHMoveInBox)=rhoself_hh
        phibox_hh(idxHMoveInBox)=phiself_hh
        fbox_h(idxHMoveInBox) = fh(rhoself_hpd+rhoself_hh)
        df_sum = df_sum + fbox_h(idxHMoveInBox)-f_h(idxHMove)
        !dphitot does not need to be updated since phi interaction only needs to be counted once

        !calculate sum of dF and dPhi
        Ediff = df_sum + dphitot

        !Do Metropolis check on new energy; iAddPtl = 0 for particle move
        call metropolis(Ediff,bAccept, 0)

        if (.not. bAccept) return

        !if accepted, update all corresponding matrices
        if (bAccept) then
        !update rho_i & f_i & phi_ij
          do i=1,numPdNL
            idxPd = idxNLPdbox(i)
            rho_pdh(idxPd) = rhobox_pdh(i)
            f_pd(idxPd) = fbox_pd(i)
            phisum_pdh(idxPd) = phibox_pdh(i)
          enddo
          do i=1,numHNL
            idxH = idxNLHbox(i)
            rho_hh(idxH) = rhobox_hh(i)
            f_h(idxH) = fbox_h(i)
            phisum_hh(idxH) = phibox_hh(i)
          enddo
          rho_hpd(idxHMove) = rhoself_hpd
          HPos(:,idxHMove) = posnew
        endif

    end subroutine mcmove

    ! subroutine to do exchange move
    subroutine mcexch (boxidx,bAccept)
        use variables
        use fvector
        use eam
        implicit none
        integer,dimension(3), intent(in) :: boxidx
        logical, intent(out) :: bAccept     !whether the mcexch is accepted or not

        integer :: icnt  ! integer to count number of H atoms
        integer :: nHAtomsOld, nHAtomsNew

                                            !neighbor list (NL) runs over 27(3*3*3) cells with the chosen box at center
                                            !each box has dimension of rcut_pdpd**3
        integer :: i,j,k,l,idx              !counters
        integer :: idxDomain                !index of the domain (calculated from boxidx), between 1 and numDomainTot
        integer :: idxHDel                  !index of H atom to be deleted
        integer :: idxHDelInBox             !index of H atom to be deleted
        integer :: idxH                     !index to be looped over all hydrogen neighbor list (HNL)
        integer :: idxPd                    !index to be looped over all Pd neighbor list (PdNL)
        integer :: numPdNL                  !number of Pd atoms within NL
        integer :: numHNL                   !number of H atoms within NL
        integer :: numH                     !number of H atoms within chosen box
        integer :: numHDomain               !number of H atoms within chosen domain
        integer :: idxBoxVal                !index of the chosen box


        real*8  :: dr                       !temp variable to save distance
        real*8  :: drmin                    !minimum distance to nearest pd atom
        real*8  :: rd                       !temp variable to save random number for delete/insert
        real*8  :: rdcut                    !temp variable to save random number for delete/insert
        real*8  :: rdH                      !temp variable to save random number for h index
        real*8  :: rhoahval                 !temp variable to save value of rho^a_h
        real*8  :: dphival                  !temp variable to save value of phi_a
        real*8  :: rhobox_hpd               !value of rho_hpd for the newly inserted H atom (contributions from Pd in numPdNL)
        real*8  :: Ediff                    !energy difference after insert/delete mc move

        real*8  :: df_sum, dphitot
        integer,dimension(3) :: idxBox
        integer,allocatable,dimension(:) :: idxNLPdbox              !list of indices of Pd atoms in NL (numPdNL)
        integer,allocatable,dimension(:) :: idxNLHbox               !list of indices of H atoms in NL (numHNL)

        real*8, dimension(3) :: rdvec                   !temp vector to save random numbers (3)
        real*8, dimension(3) :: posvec                  !temp vector to save position of selected ptl to be inserted or deleted
        real*8, dimension(3) :: pos0                    !position of origin of chosen box

        real*8, allocatable,dimension(:) :: rhobox_pdh              !array of rho_pdh in NL (numPdNL)
        real*8, allocatable,dimension(:) :: rhobox_hh               !array of rho_hh in NL (numHNL+1); last index for newly inserted H atom

        real*8, allocatable,dimension(:) :: fbox_pd,fbox_h
        real*8, allocatable,dimension(:) :: phibox_pdh,phibox_hh

        integer :: chknum, chknum2, iAddPtl

        chknum = sum(nHNL)
        iAddPtl = 0

        bAccept = .false.
        pos0(:)=(boxidx(:)-1)*rcut_pdpd
        idxDomain = ((boxidx(1)-1)/3)*(nDomainVec(2)*nDomainVec(3)) + ((boxidx(2)-1)/3)*nDomainVec(3) + ((boxidx(3)-1)/3) + 1
        idxBoxVal = (boxidx(1)-1)*(nDomainVec(2)*nDomainVec(3)*9) + (boxidx(2)-1)*nDomainVec(3)*3 + boxidx(3)

        numPdNL = nPdNL(boxidx(1),boxidx(2),boxidx(3))
        numHNL = nHNL(boxidx(1),boxidx(2),boxidx(3))
        numH = nHBox(boxidx(1),boxidx(2),boxidx(3))
        numHDomain = nHDomain(idxDomain)

        allocate(rhobox_pdh(numPdNL))
        allocate(rhobox_hh(numHNL+1))
        allocate(phibox_pdh(numPdNL))
        allocate(phibox_hh(numHNL+1))
        allocate(idxNLPdbox(numPdNL))
        allocate(idxNLHbox(numHNL))

        allocate(fbox_pd(numPdNL))
        allocate(fbox_h(numHNL+1))

        idxNLPdbox(:)=idxNLPd(1:numPdNL,boxidx(1),boxidx(2),boxidx(3))
        idxNLHbox(:)=idxNLH(1:numHNL,boxidx(1),boxidx(2),boxidx(3))

        do i=1,numHNL
          if (idxNLH(i,boxidx(1),boxidx(2),boxidx(3)).le.0) then
            write(*,*) 'ERROR idxNLH NOT ASSIGNED : ',numHNL, i, idxNLH(i,boxidx(1),boxidx(2),boxidx(3)), idxNLH(:,boxidx(1),boxidx(2),boxidx(3))
            write(*,*) 'numHNL, i, idxNLH(i), idxNLH'
            write(*,*) boxidx
            stop
          endif
        enddo
        do i=numHNL+1,nHNLMax
          if (idxNLH(i,boxidx(1),boxidx(2),boxidx(3)).ne.-1) then
            write(*,*) 'ERROR idxNLH assigned after numHNL : ',numHNL, i, idxNLH(i,boxidx(1),boxidx(2),boxidx(3)), idxNLH(:,boxidx(1),boxidx(2),boxidx(3))
            write(*,*) 'numHNL, i, idxNLH(i), idxNLH'
            write(*,*) boxidx
            stop
          endif
        enddo

        rhobox_pdh=0
        rhobox_hh=0
        rhobox_hpd=0
        phibox_pdh=0
        phibox_hh=0
        fbox_pd=0
        fbox_h=0

        df_sum=0
        dphitot=0

        rdcut = .5
        ! choose delete vs. insert; insert only when numH=0
        if (numH .eq. 0) then
            rd=1.0
        else
            call random_number(rd)
        endif


        ! delete when rd < rdcut
        if (rd .lt. rdcut) then
            call random_number(rdH)
            idxHDelInBox = int(rdH*numH)+1
            idxHDel = idxHAtomInBox(idxHDelInBox,boxidx(1),boxidx(2),boxidx(3))
            if (idxHDel .lt. 0) then
                write(*,*) 'error idxHDel : ', boxidx,idxHDelInBox, idxHDel, numH
                write(*,*) 'boxidx, idxHDelInBox, idxHDel, numH'
            endif
            posvec(:)=HPos(:,idxHDel)
            iAddPtl = -1
        ! insert when rd >= .5
        else
            call random_number(rdvec)
            posvec(:)=pos0 + rdvec(:)*rcut_pdpd
            iAddPtl = 1
        endif

        drmin = rcut_pdh+1
        !Calculate rho_pd and f_pd for numPdNL neighbors
        do i=1,numPdNL
          idxPd = idxNLPdbox(i)
          dr=norm(PdPos(:,idxPd)-posvec)
          if ( dr .le. rcut_pdh ) then
            if (rd .lt. rdcut) then
              rhoahval = -rhoa_h(dr)
              dphival = -phi_pdh(dr)
            else
              if ( dr .lt. rcut_min ) return                    ! return immediately if too close to existing atoms
              rhoahval = rhoa_h(dr)
              rhobox_hpd = rhobox_hpd + rhoa_pd(dr)
              dphival = phi_pdh(dr)
            endif
            if ( dr .lt. drmin) drmin = dr
          else
            rhoahval = 0.d0
            dphival = 0.d0
          endif
          rhobox_pdh(i)=rho_pdh(idxPd)+rhoahval
          phibox_pdh(i)=phisum_pdh(idxPd)+dphival
          dphitot = dphitot + dphival

          fbox_pd(i)=fpd(rhobox_pdh(i)+rho_pdpd(idxPd))
          df_sum = df_sum + fbox_pd(i)-f_pd(idxPd)
        enddo

        if (drmin .gt. rmin_pdh) return                       ! return immediately if no pd atom is within the cutoff

        !Calculate rho_h for numHNL neighbors
        do i=1,numHNL
          idxH = idxNLHbox(i)
          if ((rd .lt. rdcut) .and. (idxHDel==idxH)) cycle        ! skip self interaction for delete move

          dr=norm(HPos(:,idxH)-posvec)
          if ( dr .le. rcut_hh ) then
            if (rd .lt. rdcut) then
              rhoahval = -rhoa_h(dr)
              dphival = -phi_h(dr)
            else
              if ( dr .lt. rcut_min ) return                    ! return immediately if too close to existing atoms
              rhoahval = rhoa_h(dr)
              rhobox_hh(numHNL+1)=rhobox_hh(numHNL+1)+rhoahval
              dphival = phi_h(dr)
              phibox_hh(numHNL+1)=phibox_hh(numHNL+1)+dphival
            endif
          else
            rhoahval = 0.d0
            dphival = 0.d0
          endif
          rhobox_hh(i)=rho_hh(idxH)+rhoahval
          phibox_hh(i)=phisum_hh(idxH)+dphival
          dphitot = dphitot + dphival

          fbox_h(i)=fh(rhobox_hh(i)+rho_hpd(idxH))
          df_sum = df_sum + fbox_h(i)-f_h(idxH)
        enddo

        !save the rho_hpd value for the newly inserted H
        if (rd .ge. rdcut) then
            fbox_h(numHNL+1)=fh(rhobox_hpd+rhobox_hh(numHNL+1))
            df_sum = df_sum + fbox_h(numHNL+1)
            !dphitot does not need to be updated since phi interaction only needs to be counted once
        endif

        

        !calculate sum of dF and dPhi

        Ediff = df_sum + dphitot
        !Do Metropolis check on new energy; iAddPtl = +1 for insertion and -1 for deletion
        call metropolis(Ediff,bAccept, iAddPtl)

        if (.not. bAccept) return


        !check for the hpos & hnl list that does not update properly
        icnt = 0
        do i=1,nHMax
            if (HPos(1,i) .gt. 0) icnt=icnt+1
        enddo
        nHAtomsOld = icnt
        if (icnt .ne. sum(nHDomain)) then
            write(*,*) 'ERROR NHATOMS : ', icnt, sum(nHDomain)
            write(*,*) 'missing atoms in xth domain'
        endif

        !if accepted, update all corresponding matrices
        if (bAccept) then
        !update rho_i & f_i & phi_ij
          do i=1,numPdNL
            idxPd = idxNLPdbox(i)
            rho_pdh(idxPd) = rhobox_pdh(i)
            f_pd(idxPd) = fbox_pd(i)
            phisum_pdh(idxPd) = phibox_pdh(i)
          enddo
          do i=1,numHNL
            idxH = idxNLHbox(i)
            rho_hh(idxH) = rhobox_hh(i)
            f_h(idxH) = fbox_h(i)
            phisum_hh(idxH) = phibox_hh(i)
          enddo
          !update deleted or added h atom indices in NL
          if (rd .lt. rdcut) then
            idxH = idxHAtomInBox(numH,boxidx(1),boxidx(2),boxidx(3))    !idx of last h atom in the box
            !move last atom (idxH) to idxHDel
            if (idxHDel .ne. idxH) then
                !write(*,*) 'deleted posold : ', HPos(:,idxHDel), HPos(:,idxH)
                f_h(idxHDel) = f_h(idxH)
                rho_hh(idxHDel) = rho_hh(idxH)
                rho_hpd(idxHDel) = rho_hpd(idxH)
                phisum_hh(idxHDel) = phisum_hh(idxH)
                HPos(:,idxHDel) = HPos(:,idxH)
                if (any(HBox(:,idxHDel).ne.HBox(:,idxH))) then
                  write(*,*) 'error boxidx mistmatch', hBox(:,idxHDel),hBox(:,idxH),idxHDel,idxH
                endif
                !idxHAtomInBox(idxHDelInBox,boxidx(1),boxidx(2),boxidx(3))=idxHAtomInBox(numH,boxidx(1),boxidx(2),boxidx(3))
            endif
            !clear last atom in the box (idxH)
            f_h(idxH)=0             !clear f_h
            rho_hh(idxH)=0          !clear rho_hh
            rho_hpd(idxH)=0         !clear rho_hpd
            phisum_hh(idxH)=0       !clear phi_hh
            if (any(HBox(:,idxH).le.0)) then
                write(*,*) 'ERROR DELETING ATOM HBox is not assigned', HBox(:,idxH),idxH
            endif
            if (any(HPos(:,idxH).le.0)) then
                write(*,*) 'ERROR DELETING ATOM HPos is not assigned', HPos(:,idxH),idxH
            endif
            HBox(:,idxH)=-1         !unassign last atom
            HPos(:,idxH)=-1         !unassign last atom
            idxHAtomInBox(numH,boxidx(1),boxidx(2),boxidx(3))=-1

            !reduce number of h atoms
            numH = numH-1
            numHDomain = numHDomain-1
            nHBox(boxidx(1),boxidx(2),boxidx(3))=numH
            nHDomain(idxDomain) = numHDomain
            !write(*,*) 'deleted an atom', idxDomain, nHDomain(idxDomain),numH, sum(nHDomain), sum(nHBox),idxH,idxHDel
            !write(*,*) 'deleted pos : ', HPos(:,idxHDel), HPos(:,idxH)
            tempIdx0=idxH
            tempIdx1=idxHDel
            !update neighborlist
            call updateNL(boxidx,idxH,idxHDel)
          else
            idxH = numDomainTot*numHDomain + idxDomain            !index of newly added atom
            idxH = numBoxTot*numH + idxBoxVal            !index of newly added atom
            HBox(:,idxH) = boxidx(:)
            HPos(:,idxH) = posvec(:)
            rho_hh(idxH) = rhobox_hh(numHNL+1)
            rho_hpd(idxH) = rhobox_hpd
            f_h(idxH) = fbox_h(numHNL+1)
            phisum_hh(idxH) = phibox_hh(numHNL+1)

            idxHAtomInBox(numH+1,boxidx(1),boxidx(2),boxidx(3))=idxH
            nHBox(boxidx(1),boxidx(2),boxidx(3))=numH+1
            nHDomain(idxDomain) = nHDomain(idxDomain) + 1
            !if (idxDomain==47) write(*,*) 'added an atom',nHDomain(idxDomain),numH+1,sum(nHDomain), sum(nHBox), idxH
            !update neighborlist
            call updateNL(boxidx,idxH,-1)
          endif
        endif
        !update


        !check for the hpos & hnl list that does not update properly
        icnt = 0
        do i=1,nHMax
            if (HPos(1,i) .gt. 0)   icnt=icnt+1
        enddo
        if (icnt .ne. sum(nHDomain)) then
            write(*,*) 'ERROR NHATOMS AFTER MCMOVE : ', icnt, sum(nHDomain),nHAtomsOld,idxH, idxDomain, numH, nHDomain(idxDomain), posvec(:), HPos(:,idxH), nHMax
            write(*,*) 'missing atoms in xth domain'
            write(*,*) 'icnt, sumH, nHAtomsOld, idxH, idxdomain, numh, nhdomain(idxdomain),posvec, hpos'
            write(*,*) 'tempidx0, tempidx1, pos0, pos1', tempIdx0, tempIdx1, HPos(:,tempIdx0), HPos(:,tempIdx1)
            do i=1,nHMax
                idxDomain = mod(i,numDomainTot)
                !if (HPos(1,i) .gt. 0)   write(*,*) i, idxDomain, nHDomain(idxDomain),HPos(:,i)
            enddo
            
            do i=1,numDomainTot
                do j=0,nHDomain(i)-1
                    tempIdx0 = numDomainTot*j + i
                    if (HPos(1,tempIdx0) .le. 0) then
                        write(*,*) 'FOUND error atom in ', i,j+1,tempIdx0, HPos(:,tempIdx0), HBox(:,tempIdx0), numDomainTot
                        write(*,*) 'domainIdx, ithatom in domain, idxh, hpos,hbox'
                    endif
                enddo
            enddo
            stop
        endif

        !Return accepted or rejected;

    end subroutine mcexch

    ! subroutine to update NL
    subroutine updateNL (boxidx, idxH, idxHDel)
        use variables
        implicit none
        integer, intent(in) :: idxH, idxHDel
        integer,dimension(3), intent(in) :: boxidx

        integer     :: i,j,k,idx        !counters
        integer     :: chknum, chknum2, iadd

        chknum = sum(nHNL)
        iadd = 0
        do i=max(boxidx(1)-1,1),min(boxidx(1)+1,nBoxVec(1))
          do j=max(boxidx(2)-1,1),min(boxidx(2)+1,nBoxVec(2))
            do k=max(boxidx(3)-1,1),min(boxidx(3)+1,nBoxVec(3))
              !skip box that has no Pd atoms
              if (nPdBox(i,j,k)==0) cycle

              iadd = iadd +1
              !delete idx of removed H atom from NL
              if (idxHDel .gt. 0) then
                do idx=1,nHNL(i,j,k)
                  if (idxNLH(idx,i,j,k)==idxH) then
                    idxNLH(idx,i,j,k) = idxNLH(nHNL(i,j,k),i,j,k)
                    idxNLH(nHNL(i,j,k),i,j,k) = -1
                    nHNL(i,j,k)=nHNL(i,j,k)-1
                    exit
                  endif
                  if (idx .eq. nHNL(i,j,k)) then
                        write(*,*) 'error idx', boxidx, idxH, idxHDel, idx,i,j,k,nHNL(i,j,k)
                        write(*,*) 'boxidx, idxh, idxhdel, idx, i,j,k, nHNL'
                        !write(*,*) 'idxNLH', idxNLH(:,i,j,k)
                  endif
                enddo
              !append idx of newly added H atom to NL
              else
                nHNL(i,j,k) = nHNL(i,j,k)+1
                idxNLH(nHNL(i,j,k),i,j,k)=idxH
              endif
            enddo
          enddo
        enddo

        chknum2 = sum(nHNL)
        if (idxHDel .lt. 0) then
            if (chknum+iadd .ne. chknum2) write(*,*) 'addition NL update failed'
        else
            !write(*,*) 'number of nHNL deleted : ',iadd
            if (chknum-iadd .ne. chknum2) write(*,*) 'deletion NL update failed'
        endif
    end subroutine updateNL

end module func_mc
