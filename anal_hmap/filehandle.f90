!****************************************************************
! 
! this module contains utility functions for vector calculation
!
!****************************************************************

MODUlE filehandle
    use variables
    implicit none
 
contains

    ! Subroutine to get command arguments and assign input and output file name
    subroutine getargs ( aStrArgs )
        use variables
        implicit none

        character(len=256),intent(in),allocatable,dimension(:) :: aStrArgs
        integer :: narg,iarg
        logical::bLookForInp=.FALSE.
        logical::bLookForOut=.FALSE.
        logical::bFileExist
        character(len=256)         :: str, line

        ! Initialization
        strInFile = ""
        strOutFile = ""
        narg=size(aStrArgs)
        !loop across options
        do iarg=1,narg
          str = aStrArgs(iarg)
          select case(adjustl(str))
            case("--help","-h")
               write(*,*)"This is program TestArg : Version 0.1"
          
            !First known args
            case("-f")
               bLookForInp=.TRUE. !change logical value
            case("-o")
               bLookForOut=.TRUE.
          
            case default
            !Treat the second arg of a serie
              if(bLookForInp)then
                strInFile=adjustl(str) !assign a value to pedfile
                bLookForInp=.FALSE. !put the logical variable to its initial value
              elseif(bLookForOut)then
                strOutFile=adjustl(str)
                inquire(file=strOutFile,exist=bFileExist)
                if(bFileExist)then
                 write(*,*)'file ',strOutFile,' exist'
                endif
                bLookForOut=.FALSE.
              else
                write(*,*)"Option ",adjustl(str),"unknown"
              endif
          end select
        end do

        if(strInFile .eq. "") then
          write(*,*) 'input file has not set'
          write(*,*) 'usage : conductivity -f param.dat -o outfile.dat'
          stop
        endif
        
        inquire(file=strInFile,exist=bFileExist)!check if it exist
        if(.not.bFileExist)then
          write(*,*)'file ',strInFile,' not found'
          write(*,*) 'usage : conductivity -f param.dat -o outfile.dat'
          stop
        endif
        
        if(strOutFile .eq. "") then
          write(*,*) 'usage : conductivity -f param.dat -o outfile.dat'
          write(*,*) 'output file has not set'
          write(*,*) 'will use default outfile name conductivity_out.dat'
          strOutFile = 'conductivity_out.dat'
        endif

    end subroutine getargs

    ! subroutine to read the parameter file
    subroutine readparam
        use variables
        implicit none
        character(len=256)         :: str, line

        OPEN(unit=7,file=strInFile,status='old')
        read(7, '(A)') line
        write(*,*) trim(line)
        read(7,'(A)') strPDXYZFile
        write(*,*) trim(strPDXYZFile)
        read(7, '(A)') line
        write(*,*) trim(line)
        read(7,'(A)') strDCDFile
        write(*,*) trim(strDCDFile)
        read(7, '(A)') line
        write(*,*) trim(line)
        read(7,'(A)') strEnerFile
        write(*,*) trim(strEnerFile)
        read(7, '(A)') line
        read(7, *) temperature
        kT= temperature*8.617281e-5
        read(7, '(A)') line
        read(7, *) betamu
        read(7, '(A)') line
        read(7, *) nSimStep
        read(7, '(A)') line
        read(7, *) nSaveDcd
        read(7, '(A)') line
        read(7, *) nSaveEner
        read(7, '(A)') line
        read(7, *) nStepNLUpdate
        read(7, '(A)') line
        read(7, *) rcut_mcmove
        close(7)

    end subroutine readparam

    ! subroutine to read the h coordinate file (last frame)
    subroutine readlast
        use variables
        implicit none
        integer :: iAtom
        character(len=256)         :: str, line

        OPEN(unit=7,file=strDCDFile,status='old')
        read(7, *) nHAtoms
        write(*,*) 'number of H atoms : ',nHAtoms
        allocate(HPos(3,nHAtoms))
        do iAtom=1,nHAtoms
            !read(7, '(A3,3F13.6)') str,PdPos(:,iAtom)
            read(7,*) str,HPos(:,iAtom)
            !read(7, '(A3,3F13.6,I13)') str,PdPos(:,iAtom),PdGroup(iAtom)
            !write(*,*) PdPos(:,iAtom)
        enddo
        !read(7, '(A)') line
        close(7)
        
    end subroutine readlast

    ! subroutine to read the pd coordinate file
    subroutine readxyz
        use variables
        implicit none
        integer :: iAtom
        character(len=256)         :: str, line

        OPEN(unit=7,file=strPDXYZFile,status='old')
        read(7, *) nPdAtoms
        write(*,*) 'number of Pd atoms : ',nPdAtoms
        read(7, '(A)') line
        allocate(PdPos(3,nPdAtoms))
        allocate(PdGroup(nPdAtoms))
        do iAtom=1,nPdAtoms
            !read(7, '(A3,3F13.6)') str,PdPos(:,iAtom)
            read(7,*) str,PdPos(:,iAtom)
            !read(7, '(A3,3F13.6,I13)') str,PdPos(:,iAtom),PdGroup(iAtom)
            !write(*,*) PdPos(:,iAtom)
        enddo
        !read(7, '(A)') line
        close(7)
        
    end subroutine readxyz

    subroutine backupfile ( strFile )
      implicit none
      character(len=256), intent(in)   :: strFile
      character(len=256)               :: strFileIn, strFileBak
      LOGICAL :: file_exists
      integer :: i
      INQUIRE(FILE=strFile, EXIST=file_exists)
      if(file_exists) then
        strFileBak = trim(strFile) // '.bak'
        strFileIn = strFileBak
        i = 0
        INQUIRE(FILE=strFileBak, EXIST=file_exists)
        do while(file_exists)
          write(strFileBak,'(A,I0)') trim(strFileIn) // '.',i
          i = i+1
          INQUIRE(FILE=strFileBak, EXIST=file_exists)
        enddo
        call system('cp ' // trim(strFile) // ' ' // trim(strFileBak))
        write(*,*) 'file ' // trim(strFile) // ' exists and is backed up to file ' // trim(strFileBak)
      endif
    end subroutine backupfile

end module filehandle

