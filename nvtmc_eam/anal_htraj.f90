program anal_htraj
    use variables
    use filehandle
    use fvector
    use eam
    use func_anal
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

    integer :: istat                !file handle flag

    real*8              :: dr    !distance
    real*8              :: rd    !random number
    real*8              :: rho_pd_avgall, rho_h_avgall, rho_pdh_avgall, rho_hpd_avgall
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
      call getargsAnal(aStrArgs)
    else
      write(*,*) 'no arguments are given'
      write(*,*) 'usage : calc_cond -f param.dat -o outfile.dat'
      stop
    endif

    ! Read param file to locate strXtcfiles
    call readparam

    nFrames = nSimStep/nSaveDcd/nStride
    ! Read xyz file to save Pd coordinates
    call readxyz

    ! Initialize nHAtoms and skip first nSkipFrame
    call init_hanal(17)

    ! Initialize the domain decomposition
    call init_domdec

    write(*,*) 'total : ', sum(nPdBox)

    iFrame = 0
    if (bUseHMap) then
        open(unit=21,file='hmap.dat')
        open(unit=22,file='rhohsum.dat')
        num_h_grid =0
        rhoh_sum = 0
    else
        call backupfile(strEnerFile)
        OPEN(unit=18,file=strEnerFile)
        OPEN(unit=19,file='htype.dat')
        OPEN(unit=20,file='rbin_traj.dat')
        write(18,*) '#istep, nPdAtoms,  nHAtoms, Etot, E_pd, E_h, ftot, ftot_pd, ftot_h, rho_pd_avg, rho_h_avg, rho_pdh_avg, rho_hpd_avg, phitot, phitot_pdpd, phitot_pdh, phitot_hpd, phitot_hh'
        write(19,*) '#istep, nPdAtoms, nHAtoms, nHType(1:6), nHPdGroupNN(1:3), nHPdGroupCoord(1:3)'
        write(20,*) '# ibin, nPdConc, nHConc, E_pd, E_h, E_tot, rho_pd, rho_h, rho_pdpd, rho_pdh, rho_hpd, rho_hh, f_pd, f_h, phi_pd, phi_h, phi_pdpd, phi_pdh, phi_hpd, phi_hh'
    endif

    write(*,*) 'nFrames:',nFrames, nSimStep, nSaveDcd, nStride
  do istep=1,nSimStep/nSaveDcd
    call readHFrame(17,istat)
    if (istat .gt. 0) then
        write(*,*) 'ERROR in trajectory file ! The program will terminate'
        stop
    elseif (istat .lt. 0) then
        write(*,*) 'The trajectory is shorter than the given number of frames in parameter file'
        exit
    endif

    if (mod(istep-1,nStride) .eq. 0) then
        if (bUseHMap) then
            call analHMap
            write(*,*) istep, nPdAtoms, maxval(num_h_grid)
        else
            call analHFrame
            rho_pdh_avgall = (sum(rho_pdh))/(1.d0*nPdAtoms)
            rho_hpd_avgall = (sum(rho_hpd))/(1.d0*nHAtoms)
            rho_pd_avgall = (sum(rho_pdpd))/(1.d0*nPdAtoms) + rho_pdh_avgall
            rho_h_avgall = (sum(rho_hh))/(1.d0*nHAtoms) + rho_hpd_avgall
            write(18,*) istep, nPdAtoms,  nHAtoms, Etot, E_pd, E_h, ftot, ftot_pd, ftot_h, rho_pd_avgall, rho_h_avgall, rho_pdh_avgall, rho_hpd_avgall, phitot, phitot_pdpd, phitot_pdh, phitot_hpd, phitot_hh
            flush(18)
            write(*,*) istep, nPdAtoms,  nHAtoms, Etot, E_pd, E_h, ftot, ftot_pd, ftot_h, rho_pd_avgall, rho_h_avgall, rho_pdh_avgall, rho_hpd_avgall, phitot, phitot_pdpd, phitot_pdh, phitot_hpd, phitot_hh
            write(19,*) istep, nPdAtoms, nHAtoms, nHType(1:6), nHPdGroupNN(1:3), nHPdGroupCoord(1:3), nHPdGroupNN(4:15), nHPdGroupCoord(4:15)
            flush(19)
        endif
    endif
  enddo
    close(17)
    if (bUseHMap) then
        call sumHMap
        do i=1,ngridvec
          do j=1,ngridvec
            write(21,*) i,j,num_h_grid(i,j,:)
            write(22,*) i,j,rhoh_sum(i,j,:)
          enddo
          write(21,*) ''
          write(22,*) ''
        enddo
        close(21)
        close(22)
    else
        close(18)
        close(19)
        close(20)
        call finalize_anal
    endif


end program anal_htraj
