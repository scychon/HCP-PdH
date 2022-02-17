!****************************************************************
! 
! this module contains functions to generate positions & NL lists
!
!****************************************************************

MODUlE gen_pos
    implicit none
 
contains

    ! subroutine to generate H positions for NVT MC
    subroutine generate_HPos
        use variables
        implicit none
        integer :: iiter,iAtom,iFail
        logical :: bOverlap
        real*8, dimension(3) :: posnew

        if (.not. allocated(HPos)) then
            allocate(HPos(3,nHAtoms))
            allocate(HBox(3,nHAtoms))
            HPos = -1
            HBox = -1
            write(*,*) 'HPos, HBox are allocated'
        endif

        iFail = 0
        do iAtom=1,nHAtoms
            ! iterate up to 10 times
            do iiter=1,20
                call get_new_HPos(posnew)
                call checkOverlapAll(posnew,bOverlap,-1)
                if (.not. bOverlap) exit
            enddo
            if (bOverlap) then
                write(*,*) 'Failed to add non-overlapping H atom : ',iAtom
                iFail = iFail + 1
            endif
            HPos(:,iAtom) = posnew(:)
            HBox(:,iAtom) = posnew(:)/(rcut_pdpd) + 1
        enddo
        write(*,*) 'Inserted ',nHAtoms,' H atom'
        write(*,*) 'Failed to insert ',iFail,' H atom'
    end subroutine generate_HPos

    ! get random position for HAtom near a random Pd atom
    subroutine get_new_HPos (posnew)
        use variables
        implicit none
        real*8,dimension(3),intent(out) :: posnew
        integer :: idxPdAtom
        integer :: nPdCurrAtoms
        real*8  :: rdIdx
        real*8,dimension(3) :: rdvec

        nPdCurrAtoms = nPdAtoms
        if (bExcludeSurfPd) nPdCurrAtoms = nPdCurrAtoms-nPdSurfAtoms
        ! choose Pd atom to associate with the new H atom
        call random_number(rdIdx)
        idxPdAtom = int(rdIdx*nPdCurrAtoms)+1
        if (bExcludeSurfPd) idxPdAtom = idxPdInAtom(idxPdAtom)

        ! get new position
        call random_number(rdvec)
        posnew(:)=PdPos(:,idxPdAtom) + rcut_minpdh + rdvec(:)*(rcut_new-rcut_minpdh)
    end subroutine get_new_HPos

    !check whether new atom overlaps with any Pd atom (within rcut_min)
    subroutine checkOverlapPd (posnew, bOverlap)
        use variables
        use fvector
        implicit none
        real*8,dimension(3),intent(in) :: posnew
        logical,intent(out) :: bOverlap
        integer :: i
        integer :: numPdNL
        real*8 :: dr
        integer,dimension(3) :: boxidx
        
        boxidx(:) = posnew(:)/(rcut_pdpd)+1
        numPdNL = nPdNL(boxidx(1),boxidx(2),boxidx(3))
        bOverlap = .false.
        !check if it overlaps with any Pd atoms in the box
        do i=1,numPdNL
            dr=norm(PdPos(:,idxNLPd(i,boxidx(1),boxidx(2),boxidx(3)))-posnew)
            if ( dr .le. rcut_minpdh ) then
                bOverlap = .true.
                return
            endif
        enddo
    end subroutine checkOverlapPd

    !check whether new atom overlaps with any H atom (within rcut_min)
    subroutine checkOverlapH (posnew, bOverlap,idxHDel)
        use variables
        use fvector
        implicit none
        real*8,dimension(3),intent(in) :: posnew
        logical,intent(out) :: bOverlap
        integer,intent(in)  :: idxHDel
        integer :: i
        real*8 :: dr
        
        bOverlap = .false.
        !check if it overlaps with any H atoms in the box
        do i=1,nHAtoms
            if ((i==idxHDel) .or. (HPos(1,i).le.0)) cycle
            dr=norm(HPos(:,i)-posnew)
            if ( dr .le. rcut_minhh ) then
                bOverlap = .true.
                return
            endif
        enddo
    end subroutine checkOverlapH

    !check whether new atom overlaps with any atom (within rcut_min)
    subroutine checkOverlapAll (posnew, bOverlap,idxHDel)
        use variables
        use fvector
        implicit none
        real*8,dimension(3),intent(in) :: posnew
        logical,intent(out) :: bOverlap
        integer,intent(in)  :: idxHDel
        
        call checkOverlapPd(posnew,bOverlap)
        if (.not. bOverlap) call checkOverlapH(posnew,bOverlap,idxHDel)
    end subroutine checkOverlapAll

    ! subroutine to get NL for H atoms
    subroutine update_NL_H (idxH)
        use variables
        use fvector
        implicit none
        integer, intent(in) :: idxH
        integer :: i,j,k,iAtom,jAtom,tempmax
        integer :: idxTemp
        integer,dimension(3) :: idxBox
        real*8 :: dr

        !write(*,*) 'start get_NL_H'
        nHNLMax = size(idxNLH(:,1))

        if (nHNL(idxH).gt.0) then
          do iAtom=1,nHNL(idxH)
            idxTemp = idxNLH(iAtom,idxH)
            if (idxTemp .lt. 0) write(*,*) 'ERROR idxNLH is not assigned but nHNL is not zero',idxH,iAtom,idxTemp,nHNL(idxH),idxNLH(:,idxH)
            do jAtom=1, nHNL(idxTemp)
                if (idxNLH(jAtom,idxTemp).eq.idxH) then
                    idxNLH(jAtom,idxTemp)=idxNLH(nHNL(idxTemp),idxTemp)
                    idxNLH(nHNL(idxTemp),idxTemp) = -1
                    nHNL(idxTemp)=nHNL(idxTemp)-1
                    exit
                elseif (jAtom==nHNL(idxTemp)) then
                    write(*,*) 'ERROR nHNL reached END IDX WITHOUT FINDING NN',idxTemp,idxH
                endif
            enddo
          enddo
        endif
          nHNL(idxH)=0
          idxNLH(:,idxH) = -1
          do iAtom=1,nHAtoms
            if (iAtom==idxH) cycle
            dr=norm(HPos(:,iAtom)-HPos(:,idxH))
            if (dr .lt. rcut_hh) then
                nHNL(idxH) = nHNL(idxH)+1
                nHNL(iAtom) = nHNL(iAtom)+1
                if ((nHNL(iAtom) .gt. nHNLMax).or.(nHNL(idxH) .gt. nHNLMax)) then
                    write(*,*) 'tot nHNL : ', sum(nHNL)
                    call expand2Dint(idxNLH,nHNLMax,0)
                    nHNLMax = nHNLMax*2
                endif
                idxNLH(nHNL(idxH),idxH) = iAtom
                idxNLH(nHNL(iAtom),iAtom) = idxH
            endif
          enddo
        !tempmax = maxval(nHNL)
        !write(*,*) 'max nHNL : ',tempmax

    end subroutine update_NL_H

    ! subroutine to get NL for H atoms
    subroutine get_NL_H
        use variables
        use fvector
        implicit none
        integer :: i,j,k,iAtom,jAtom,tempmax
        integer,dimension(3) :: idxBox,idxTemp
        real*8 :: dr

        !write(*,*) 'start get_NL_H'
        if (.not. allocated(nHNL)) allocate(nHNL(nHAtoms))
        if (.not. allocated(idxNLH)) allocate(idxNLH(nHNLMax,nHAtoms))
        nHNLMax = size(idxNLH(:,1))

        nHNL=0
        idxNLH = -1
        do iAtom=1,nHAtoms-1
          do jAtom=iAtom+1,nHAtoms
            dr=norm(HPos(:,iAtom)-HPos(:,jAtom))
            if (dr .lt. rcut_hh) then
                nHNL(iAtom) = nHNL(iAtom)+1
                nHNL(jAtom) = nHNL(jAtom)+1
                if ((nHNL(iAtom) .gt. nHNLMax).or.(nHNL(jAtom) .gt. nHNLMax)) then
                    write(*,*) 'tot nHNL : ', sum(nHNL)
                    call expand2Dint(idxNLH,nHNLMax,0)
                    nHNLMax = nHNLMax*2
                endif
                idxNLH(nHNL(iAtom),iAtom) = jAtom
                idxNLH(nHNL(jAtom),jAtom) = iAtom
            endif
          enddo
        enddo
        tempmax = maxval(nHNL)
        !write(*,*) 'max nHNL : ',tempmax

    end subroutine get_NL_H

end module gen_pos
