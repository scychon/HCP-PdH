!  DCD Fortran Interface DCD Example Program with Wrapper
!  2017 (c) Chang Yun Son <cson@chem.wisc.edu>
!
 
module dcd
  implicit none
  private
  public :: openmmdcdfile

  integer :: fp = 10
  real    :: prec = 1e-4

  type openmmdcdfile
    integer :: STAT, STEP, NATOMS, NFRAMES, INTERVAL, FIRSTSTEP
    real :: time, prec, dt, dtframe
    logical           :: bBoxFlag
    real*8,dimension(3,2) :: box
    real,dimension(:,:),allocatable :: pos
    
    contains
      procedure :: init => dcd_init
      procedure :: write => dcd_write
      procedure :: close => dcd_close
  end type openmmdcdfile
contains
  subroutine dcd_init(this, fname, natom)
    class(openmmdcdfile), intent(inout) :: this
    CHARACTER(len=256), intent(in) :: fname
    integer, intent(in) :: natom
    integer           :: istat,iBoxFlag
    integer           :: dummyi,i,testi
    integer           :: dummyiar(8)
    real              :: dummyr
    character*4       :: dummyc
    character*80      :: dummystr,dummystr2

    dummyc='CORD'
    dummystr='Created by dcdmod'
    dummystr2='Created by dcdmod'
    this % FIRSTSTEP = 0
    this % INTERVAL = 1
    this % dt = 1.0/0.04888821
    iBoxFlag = 0
    dummyi = 0
    testi = 0
    dummyiar = 0
    write(*,*) 'start reading input file ',trim(fname)
    open(fp,file=trim(fname),form='unformatted')
    write(fp) dummyc, this % NFRAMES, this % FIRSTSTEP, this % INTERVAL, (dummyi,i=1,6), this % dt, &
            iBoxFlag, (dummyiar(i),i=1,8),testi

    this % STEP = this % FIRSTSTEP
    this % dtframe = this % INTERVAL * this % dt
    this % time = this % FIRSTSTEP * this % dt

    write(*,*) 'dummyc NFRAMES FIRSTSTEP STEP INTERVAL dt',dummyc, this % NFRAMES, &
                this % FIRSTSTEP, this % STEP, this % INTERVAL, this % dt
    !read(fp) iBoxFlag, (dummyi,i=1,8),testi
    !read(fp) iBoxFlag, (dummyiar(i),i=1,8),testi
    write(*,*) 'iBoxFlag',iBoxFlag, (dummyiar(i),i=1,8),testi
    write(fp) dummyi,dummystr,dummystr2
    write(*,*) 'dummystr',dummyi,dummystr,dummystr2
    this % NATOMS = natom
    write(fp) this % NATOMS
    write(*,*) 'NATOMS', this % NATOMS
    
    if(iBoxFlag .ne. 0) this % bBoxFlag = .True.
    this % box(:,1) = 400
    this % box(:,2) = 0
    if(.not. allocated(this % pos)) then
        allocate(this % pos(3,this % NATOMS))
        write(*,*) 'position array is allocated '
    endif
  end subroutine dcd_init

  subroutine dcd_write(this)
    class(openmmdcdfile), intent(inout) :: this
    integer           :: istat
    integer           :: dummyi, ilength
    real              :: dummyr
    integer           :: natom,i,j

    natom = this % NATOMS
    this % STEP = this % STEP + this % INTERVAL
    this % time = this % time +  this % dtframe

    write(fp) (this % box(1,j), j=1,2),(this % box(2,j),j=1,2), & 
            (this % box(3,j),j=2,1,-1)
    !write(*,*) 'Box length',(this % box(i,1), i=1,3),(this % box(i,2),i=1,3)
    write(fp) (this % pos(1,j),j=1,natom)
    write(fp) (this % pos(2,j),j=1,natom)
    write(fp) (this % pos(3,j),j=1,natom)
    flush(fp)
  end subroutine dcd_write

  subroutine dcd_close(this)
    class(openmmdcdfile), intent(inout) :: this
    close(fp)
  end subroutine dcd_close

end module dcd
