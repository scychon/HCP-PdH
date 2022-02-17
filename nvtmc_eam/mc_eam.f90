program mc_eam
    use variables
    use filehandle
    use fvector
    use eam
    use func_mc
    use domdec
    use gen_pos
    use omp_lib
    use dcd, only: openmmdcdfile

    implicit none
    type(openmmdcdfile) :: openmmdcd

    integer, parameter  :: niter = 10                !max # of MCtrials in each step
    integer :: nMove                                 !# of H ptl move MC trial (==numH)
    integer, parameter  :: nExch = 20                !max # of MCtrials in each step
    integer             :: i,j,k,ii,jj,kk,l,icnt,idx,ibox     !counters
    integer             :: idxpd,idxh,idxmol         !atom indices
    integer             :: iiter,istep, rdint
    integer             :: narg, cptArg              !#of arg & counter of arg
    integer             :: nsize                     !size of Hatom array
    integer             :: numthread,threadid
    integer,dimension(3) :: domainGroupIdx, boxidx         ! index of the current box (domainGroupIdx*3 + ii/jj/kk)

    logical :: bAccept                          !whether the mc move is accepted or not
    integer :: nAccept, nReject, nTry           !# of mc moves that is accepted, rejected, tried in each step
    integer :: nAccTot, nRejTot, nTryTot        !total # of mc moves that is accepted, rejected, tried
    integer :: nAccMVTot, nRejMVTot, nTryMVTot        !total # of mc moves that is accepted, rejected, tried
    real*8  :: rAccept,rAccTot                  !acceptance ratio

    real*8              :: dr    !distance
    real*8              :: rd    !random number
    real*8, dimension(3) :: rd_val
    integer,allocatable,dimension(:) :: pdlist,hlist

    character(len=256)         :: str, line
    character(len=256),allocatable,dimension(:) :: aStrArgs
! ----------------------------------------------------------------

    write(*,*) 'get omp num threads'
    !$OMP PARALLEL
    numthread = omp_get_num_threads()
    !$OMP END PARALLEL
    write(6,*) 'Using ',numthread,' number of threads'

    !Check if any arguments are found
    narg=command_argument_count()
    !Loop over the arguments
    if(narg>0)then
      !loop across options
      allocate(aStrArgs(narg))
      do cptArg=1,narg
        call get_command_argument(cptArg,str)
        aStrArgs(cptArg) = str
      enddo
      !assign in and out files
      call getargs(aStrArgs)
    else
      write(*,*) 'no arguments are given'
      write(*,*) 'usage : calc_cond -f param.dat -o outfile.dat'
      stop
    endif

    ! Read param file to locate strXtcfiles
    call readparam

    ! Read xyz file to save Pd coordinates
    call readxyz

    ! Initialize rbin matrices
    call calc_pdconc

    ! Initialize the domain decomposition
    call init_domdec

    ! initialize dcd trajectory file
    call openmmdcd % init(strDCDFile,nHMax)

!    do i=1,15
!      do j=1,15
!        do k=1,15
!            write(*,*) 'nPdDomain : ',i,j,k,nPdDomain(i,j,k)
!        enddo
!      enddo
!    enddo
    write(*,*) 'total : ', sum(nPdBox)

    OPEN(unit=17,file='traj.xyz')
    OPEN(unit=18,file=strEnerFile)
    write(18,*) '#istep,   nHAtoms,     nAccpt,      nReject,      rAccpt,  ftot, phitot, ftot+phitot, ftot_pd, ftot_h, phitot_pdh, phitot_hh'

    rAccTot = 0.d0
    nAccTot = 0
    nRejTot = 0
    nTryTot = 0
    nAccMVTot = 0
    nRejMVTot = 0
    nTryMVTot = 0

    iFrame = 0
    openmmdcd % pos = HPos
    call openmmdcd % write
    nHAtoms = sum(nHDomain)
    OPEN(unit=16,file='first.xyz')
    write(16,'(I0)') nHAtoms
    write(16,*) iFrame,'th Frame'
    do i=1, nHMax
        if (HPos(1,i) .gt.0) write(16,'(A3,3f13.6)') '  H', HPos(:,i)
    enddo
    close(16)

    rAccept = 0.d0
    nAccept = 0
    nReject = 0
    nTry = 0
  do istep=1,nSimStep
    call random_number(rd)
    ibox = floor(27*rd)
    ii=ibox/9
    jj=(ibox-9*ii)/3
    kk=mod(ibox,3)


    !if (mod(istep,1) .eq. 0) then
    if (mod(istep,nStepNLUpdate) .eq. 0) then
        nHAtoms = sum(nHDomain)
        write(*,*) 'step : ', istep,' nHAtoms : ',nHAtoms, ', sum(nHDomain) : ',sum(nHDomain), ', sum(nHNL) : ', sum(nHNL)
        call get_NL_H
        write(*,*) 'step : ', istep,' nHAtoms : ',nHAtoms, ', sum(nHDomain) : ',sum(nHDomain), ', sum(nHNL) : ', sum(nHNL)
        write(*,*) 'step : ', istep, ' nHAtoms : ', nHAtoms, ' total ', nAccTot, ' accepted, ', nRejTot, ' rejected, ', nTryTot, ' tried, ', rAccTot, ' ratio' 
    endif

        nMove = 20
        !call get_neighbor(boxidx,pdlist,hlist)
        do iiter=1,niter
            call random_number(rd)
            rdint = int(rd*(nMove))+1
            !rdint = int(rd*(nMove+nExch))+1
            if (rdint .le. nMove) then
                call mcmove(bAccept)
                if (bAccept) nAccMVTot = nAccMVTot+1
                nTryMVTot = nTryMVTot+1
            else
                write(*,*) 'MCEXCH move is done WHY?',rdint,nMove
                !call mcexch(boxidx,bAccept)
            endif
            !stop iteration if accepted
            if (bAccept) then
                nAccept = nAccept + 1
                nReject = nReject + iiter-1
                nTry = nTry + iiter
                exit
            endif
        enddo
        if (bAccept) then
            bAccept = .false.
        else
            nReject = nReject + niter
            nTry = nTry + niter
            !write(*,*) 'domain ',idx, 'got rejected'
        endif

    rAccept = float(nAccept)/float(nTry)
    nAccTot = nAccTot + nAccept
    nRejTot = nRejTot + nReject
    nTryTot = nTryTot + nTry
    rAccTot = float(nAccTot)/float(nTryTot)
    if (mod(istep,nSaveEner)==0) then
        nHAtoms = sum(nHDomain)
        write(*,*) nAccept, ' accepted, ', nReject, ' rejected, ', nTry, ' tried, ', rAccept, ' ratio' 
        write(*,*) 'step : ', istep, ' nHAtoms : ', nHAtoms, ' total ', nAccTot, ' accepted, ', nRejTot, ' rejected, ', nTryTot, ' tried, ', rAccTot, ' ratio' 
        ftot_pd = sum(f_pd)
        ftot_h = sum(f_h)
        ftot = ftot_pd+ftot_h
        phitot_pdh = sum(phisum_pdh)
        phitot_hh = sum(phisum_hh)
        phitot = phitot_pdpd+phitot_pdh+phitot_hh
        write(18,*) istep, nHAtoms, nAccTot, nRejTot, rAccTot, ftot, phitot, ftot+phitot, ftot_pd, ftot_h, phitot_pdh, phitot_hh,phitot_pdpd
        flush(18)
        rAccTot = float(nAccTot-nAccMVTot)/float(nTryTot-nTryMVTot)
        write(*,*) 'stepEx : ', istep, ' nHAtoms : ', nHAtoms, ' total ', nAccTot-nAccMVTot, ' accepted, ', nRejTot-(nTryMVTot-nAccMVTot), ' rejected, ', nTryTot-nTryMVTot, ' tried, ', rAccTot, ' ratio' 
        rAccTot = float(nAccMVTot)/float(nTryMVTot)
        write(*,*) 'stepMV : ', istep, ' nHAtoms : ', nHAtoms, ' total ', nAccMVTot, ' accepted, ', nTryMVTot-nAccMVTot, ' rejected, ', nTryMVTot, ' tried, ', rAccTot, ' ratio' 
        rAccept = 0.d0
        nAccept = 0
        nReject = 0
        nTry = 0
    endif

    if (mod(istep,nSaveDcd)==0) then
        openmmdcd % pos = HPos
        call openmmdcd % write
        nHAtoms = sum(nHDomain)
        iFrame = iFrame + 1
        write(17,'(I0)') nHAtoms
        write(17,*) iFrame,'th Frame'
        do i=1, nHMax
            if (HPos(1,i) .gt.0) write(17,'(A3,3f13.6)') '  H', HPos(:,i)
        enddo
        flush(17)
    endif

  enddo
    call openmmdcd % close

    close(18)
    close(17)
end program mc_eam
