program anal_hmap
    use variables
    use filehandle
    use fvector
    use eam
    use func_mc
    use domdec
    use omp_lib

    implicit none
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

    ! Read xyz file to load Pd coordinates
    call readxyz

    ! Read xyz file to load H coordinates
    call readlast

    ! Initialize the domain decomposition
    !call init_domdec

    ! Precalculate pd-h interaction and rho_h
    call calc_pdh

    ! Precalculate pd-pd interaction and rho_pdpd
    call calc_pdpd

    allocate(SliceIdx(31))
    SliceIdx = (/ 212,206,199,193,186,179,173,166,160,153,147,140,133,127,120,114,107,101,94,88,81,75,68,61,55,48,42,36,30,23,16 /)
    write(*,*) SliceIdx

    call map_slice


end program anal_hmap
