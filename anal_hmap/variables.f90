!**************************************************
! these are global variables used by many subroutines throughout
! the program
!**************************************************

MODULE variables

  integer :: nSaveDcd, nSaveLog, nSaveEner         !output file frequency
  integer :: nSimStep               !# of simulation steps
  integer :: nPdAtoms               !# of Pd Atoms in system; fixed value, read from xyz file
  integer :: nHMax                  !max # of H Atoms in system (3*nPdAtoms)
  integer :: nHAtoms                !# of H Atoms in system; updated after each mc_exch (add/del) move
  integer :: nTotAtoms              !# of total atoms in system; nPdAtoms + nHAtoms

  integer :: nPdBoxMax                  !max # of Pd Atoms in chosen box; calculated during initialization from xyz file
  integer :: nPdDomainMax               !max # of Pd Atoms in the domain; calculated during initialization from xyz file
  integer :: nHBoxMax                   !max # of H Atoms in chosen box (3*nPdBoxMax)
  integer :: nHDomainMax                !max # of H Atoms in the domain
  integer :: nStepNLUpdate              !number of mc steps between NL update
  integer :: nPdNLMax                   !max # of Pd Atoms in neighbor list;calculated during initialization from xyz file
  integer :: nHNLMax                    !max # of Pd Atoms in neighbor list (3*nPdNLMax)
                                        !NL is composed of 27 (3*3*3) boxes around the chosen box (domain)
  integer :: numDomainTot               !# of domain groups; calculated during initialization from xyz file
  integer :: numBoxTot                  !# of boxes; calculated during initialization from xyz file
  integer,dimension(3) :: nDomainVec    !# of domain groups in xyz (each domain group has 27 boxes)
  integer,dimension(3) :: nBoxVec       !# of boxes in xyz (nDomainVec*3)

  real*8  :: temperature            !temperature of the system
  real*8  :: kT,beta                !kT in eV unit
  real*8  :: betamu                 !mu/kT (dimensionless)
  real*8  :: volume                 !volume of the box (in angstrom**3 unit)
  real*8  :: nqvol                  !nq*volume = volume * (m*T)**1.5*nqunit (dimensionless)
  real*8  :: rcut_mcmove            !maximum length to move ptl for mcmove
  real*8, parameter    :: rcut_pdpd=5.35, rcut_hh=5.35, rcut_pdh=4.95     !cutoff distances for each pair interactions phi_ij
  real*8, parameter    :: rcut_min=0.5                                    !cutoff distance to immediately reject insertion
  !real*8, parameter    :: rcut_mcmove=2                                   !maximum length to move ptl for mcmove
  real*8, parameter    :: rmin_pdh=2.2                                    !maximum distance to the closest pd atom for new h atom
  real*8, parameter    :: nqunit=0.00018793204250346798                   !unit of nq=(2pi*m*kB*T/h)**1.5 in (amu*kelvin)**1.5/angstrom**3 unit
                                    !obtained from ((2*np.pi*0.001*kilogram/AVOGADRO_CONSTANT_NA/mole*BOLTZMANN_CONSTANT_kB*kelvin)**1.5/(6.62607015e-34*joules*second)**3)*angstrom**3

  real*8  :: ftot                   !sum of f energy; ftot = ftot_pd+ftot_h
  real*8  :: ftot_pd                !sum of f_pd
  real*8  :: ftot_h                 !sum of f_h
  real*8  :: phitot                 !sum of pair interaction energy phi_ij; sum of phitot_pdpd + phitot_pdh + phitot_hh
  real*8  :: phitot_pdpd            !contribution from pd-pd interaction phi_pdpd; calculated during initialization from xyz file
  real*8  :: phitot_pdh             !contribution from pd-h interaction; can be obtained by sum(phisum_pdh) == sum(phisum_hpd)
  real*8  :: phitot_hh              !contribution from h-h interaction; sum(phisum_hh)
  real*8, allocatable,dimension(:) :: phisum_pdh       !sum of phi_pd for pd ptl i due to h in system (nPdAtoms); update after each MC move
  real*8, allocatable,dimension(:) :: phisum_hh        !sum of phi_h for h ptl i due to h in system (nHMax); update after each MC move
    
  real*8, parameter     :: slice_res=0.3479          !size of each slice
  integer,DIMENSION(:),ALLOCATABLE    :: SliceIdx   !Idxs of slices in z direction to draw H concentration map
  real*8,DIMENSION(:,:,:),ALLOCATABLE :: rho_hmap   !density map at the ith slice in z direction to draw H concentration map

  integer,DIMENSION(:),ALLOCATABLE    :: PdGroup   !Pd group assigned from tomography data (nPdAtoms); -1~5, read from xyz file
  integer,DIMENSION(:),ALLOCATABLE    :: nHDomain  ! array of number of H ptls in domain (idxDomain, between 1 and numDomainTot)
  integer,DIMENSION(:),ALLOCATABLE    :: nPdDomain ! array of number of Pd ptls in domain (idxDomain, between 1 and numDomainTot)
  integer,DIMENSION(:,:),ALLOCATABLE  :: PdBox, HBox  ! array of box idxs of idxAtom th ptl  (3, idxAtom)
  integer,DIMENSION(:,:,:),ALLOCATABLE:: nPdBox, nHBox,nPdNL,nHNL  ! array of number of Pd/H ptls in box (i,j,k)
  integer,DIMENSION(:,:,:,:),ALLOCATABLE:: idxPdAtomInBox, idxHAtomInBox  ! array of atom idxs in box (idxAtom,i,j,k)
  integer,DIMENSION(:,:),ALLOCATABLE  :: idxNLPdH  ! array of neighbor idxs (idx, idxAtom)
  integer,DIMENSION(:,:,:,:),ALLOCATABLE  :: idxNLPd, idxNLH  ! array of neighbor idxs (idxNLAtom,i,j,k)

  real*8,DIMENSION(:,:),ALLOCATABLE   :: PdPos,HPos    ! array of Pd positions (xyz, nPdAtoms)

  real*8, allocatable,dimension(:) :: f_pd,f_h
  real*8, allocatable,dimension(:) :: rho_pdpd         !array of rho_pdpd in system (nPdAtoms); calculate only once during initialization
  real*8, allocatable,dimension(:) :: rho_pdh          !array of rho_pd due to h in system (nPdAtoms); update after each MC move
  real*8, allocatable,dimension(:) :: rho_hpd          !array of rho_h due to pd in system (nHMax); update after each MC move
  real*8, allocatable,dimension(:) :: rho_hh           !array of rho_h due to h in system (nHMax); update after each MC move

  character(len=256)         :: strInFile, strOutFile, strTopFile
  character(len=256)         :: strPDXYZFile, strDCDFile, strEnerFile

  integer :: tempIdx0, tempIdx1

END MODULE variables
