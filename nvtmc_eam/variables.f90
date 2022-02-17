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
  integer :: nPdSurfAtoms           !# of Pd Atoms on the surface; fixed value, calculated from xyz file
  integer :: nHSurfAtoms            !# of H Atoms on surface; updated if old/new position is on the surface
  integer :: nHSurfAtomsMax         !max # of H Atoms on surface; (5*nPdSurfAtoms)
  integer,parameter :: nPdNNAtoms=12    !# of nearest neighbour atoms for each Pd
  logical :: bExcludeSurfPd         !true if exclude surface Pd atoms when generating new H position nearby

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
  real*8, dimension(3) :: rCOM      !position vector of COM for all Pd atoms
  real*8, parameter    :: rcut_pdpd=5.35, rcut_hh=5.35, rcut_pdh=5.35     !cutoff distances for each pair interactions phi_ij
  real*8, parameter    :: rcut_all=5.35, rcut_rs=0.3                      !distance from cutoff distance to start shifting value
  real*8, parameter    :: rcut_min=0.5,rcut_minhh=1.5,rcut_minpdh=1.0     !cutoff distance to immediately reject insertion
  real*8, parameter    :: rcut_new=2.5                                    !maximum distance of PdH bond for new atom
  real*8, parameter    :: rcut_hcoord=2.7                                 !maximum distance of PdH bond for coordination
  !real*8, parameter    :: rcut_mcmove=2                                   !maximum length to move ptl for mcmove
  real*8, parameter    :: rmin_pdh=2.2                                    !maximum distance to the closest pd atom for new h atom
  real*8, parameter    :: nqunit=0.00018793204250346798                   !unit of nq=(2pi*m*kB*T/h)**1.5 in (amu*kelvin)**1.5/angstrom**3 unit
                                    !obtained from ((2*np.pi*0.001*kilogram/AVOGADRO_CONSTANT_NA/mole*BOLTZMANN_CONSTANT_kB*kelvin)**1.5/(6.62607015e-34*joules*second)**3)*angstrom**3

  real*8  :: rHPerPd                !ratio of H atoms per Pd atoms
  real*8  :: ftot                   !sum of f energy; ftot = ftot_pd+ftot_h
  real*8  :: ftot_pd                !sum of f_pd
  real*8  :: ftot_h                 !sum of f_h
  real*8  :: phitot                 !sum of pair interaction energy phi_ij; sum of phitot_pdpd + phitot_pdh + phitot_hh
  real*8  :: phitot_pdpd            !contribution from pd-pd interaction phi_pdpd; calculated during initialization from xyz file
  real*8  :: phitot_pdh             !contribution from pd-h interaction; can be obtained by sum(phisum_pdh) == sum(phisum_hpd)
  real*8  :: phitot_hpd             !contribution from pd-h interaction; can be obtained by sum(phisum_pdh) == sum(phisum_hpd)
  real*8  :: phitot_hh              !contribution from h-h interaction; sum(phisum_hh)
  real*8  :: dphitot                !changes in sum of pair interaction energy phi_ij; sum of phitot_pdpd + phitot_pdh + phitot_hh
  real*8  :: dphitot_pdh            !changes in contribution from pd-h interaction; can be obtained by sum(phisum_pdh) == sum(phisum_hpd)
  real*8  :: dphitot_hh             !changes in contribution from h-h interaction; sum(phisum_hh)
  real*8, allocatable,dimension(:) :: phisum_pdpd      !sum of phi_pd for pd ptl i due to pd in system (nPdAtoms)
  real*8, allocatable,dimension(:) :: phisum_pdh       !sum of phi_pd for pd ptl i due to h in system (nPdAtoms); update after each MC move
  real*8, allocatable,dimension(:) :: phisum_hh        !sum of phi_h for h ptl i due to h in system (nHMax); update after each MC move
  real*8, allocatable,dimension(:) :: dphisum_pdh    !backup sum of phi_pd for pd ptl i due to h in system (nPdAtoms); update after each MC move
  real*8, allocatable,dimension(:) :: dphisum_hh     !backup sum of phi_h for h ptl i due to h in system (nHMax); update after each MC move
    

  integer,DIMENSION(:),ALLOCATABLE    :: PdGroup   !Pd group assigned from tomography data (nPdAtoms); -1~5, read from xyz file
  integer,DIMENSION(:),ALLOCATABLE    :: nHDomain  ! array of number of H ptls in domain (idxDomain, between 1 and numDomainTot)
  integer,DIMENSION(:),ALLOCATABLE    :: nHNL      ! array of number of NLs per H atom (nHAtoms)
  integer,DIMENSION(:),ALLOCATABLE    :: nPdDomain ! array of number of Pd ptls in domain (idxDomain, between 1 and numDomainTot)
  integer,DIMENSION(:,:),ALLOCATABLE  :: PdBox, HBox  ! array of box idxs of idxAtom th ptl  (3, idxAtom)
  integer,DIMENSION(:,:,:),ALLOCATABLE:: nPdBox, nHBox,nPdNL  ! array of number of Pd/H ptls in box (i,j,k)
  integer,DIMENSION(:,:,:,:),ALLOCATABLE:: idxPdAtomInBox, idxHAtomInBox  ! array of atom idxs in box (idxAtom,i,j,k)
  integer,DIMENSION(:,:),ALLOCATABLE  :: idxNLPdH  ! array of neighbor idxs (idx, idxAtom)
  integer,DIMENSION(:,:,:,:),ALLOCATABLE  :: idxNLPd  ! array of neighbor idxs (idxNLAtom,i,j,k)
  integer,DIMENSION(:,:),ALLOCATABLE  :: idxNLH  ! array of neighbor idxs (idxNLAtom,nHAtoms)

  logical,DIMENSION(:),ALLOCATABLE    :: bPdSurfAtom   !true if the Pd atom is on surface (nPdAtoms)
  integer,DIMENSION(:),ALLOCATABLE    :: idxPdSurfAtom !idx of surface Pd atoms (nPdSurfAtoms)
  integer,DIMENSION(:),ALLOCATABLE    :: idxPdInAtom   !idx of inner Pd atoms (nPdAtoms - nPdSurfAtoms)
  integer,DIMENSION(:),ALLOCATABLE    :: idxHSurfAtom  !idx of surface H atoms (nHSurfAtoms)
  integer,DIMENSION(:),ALLOCATABLE    :: idxHInAtom    !idx of inner H atoms (nHAtoms - nHSurfAtoms)
  integer,DIMENSION(:,:),ALLOCATABLE  :: idxPdNNAtoms  !idxs of nearest Pd neighbours for each Pd atoms (nPdNNAtoms,nPdAtoms);  idx of 1-12the neighbouring atoms, idx of Pd atom
  real*8, DIMENSION(:,:),ALLOCATABLE  :: drPdNNAtoms   !distance to the ith nearest Pd neighbour for each Pd atoms (nPdNNAtoms,nPdAtoms)
  integer,DIMENSION(:),ALLOCATABLE    :: idxHNNAtom    !idx of nearest Pd neighbour for each H atom (nHAtoms+1), last index is for trial insert
  real*8, DIMENSION(:),ALLOCATABLE    :: drHNNAtom     !distance to the nearest Pd neighbour for each H atom (nHAtoms+1), last index for trial insert

  real*8,DIMENSION(:,:),ALLOCATABLE   :: PdPos,HPos    ! array of Pd positions (xyz, nPdAtoms)

  real*8, allocatable,dimension(:) :: f_pd,f_h
  real*8, allocatable,dimension(:) :: df_pd,df_h
  real*8, allocatable,dimension(:) :: rho_pdpd         !array of rho_pdpd in system (nPdAtoms); calculate only once during initialization
  real*8, allocatable,dimension(:) :: rho_pdh          !array of rho_pd due to h in system (nPdAtoms); update after each MC move
  real*8, allocatable,dimension(:) :: rho_hpd          !array of rho_h due to pd in system (nHMax); update after each MC move
  real*8, allocatable,dimension(:) :: rho_hh           !array of rho_h due to h in system (nHMax); update after each MC move
  real*8, allocatable,dimension(:) :: drho_pdh        !backup array of rho_pd due to h in system (nPdAtoms); update after each MC move
  real*8, allocatable,dimension(:) :: drho_hpd        !backup array of rho_h due to pd in system (nHMax); update after each MC move
  real*8, allocatable,dimension(:) :: drho_hh         !backup array of rho_h due to h in system (nHMax); update after each MC move

  character(len=256)         :: strInFile, strOutFile, strTopFile
  character(len=256)         :: strPDXYZFile, strDCDFile, strEnerFile, strTrajFile

  integer :: tempIdx0, tempIdx1

  integer :: nFrames                !# of frames in trajectory analysis
  integer :: nSkipFrame             !# of frames to skip from the start of the trajectory
  integer :: iFrame                 !idx of frames of H trajectory
  integer :: nStride                !read every # of frames from the H trajectory
  logical :: bHTrajAnal             !true for anal_hpos program
  logical :: bUsePdGroup            !true for nanoptl pd group assignment

  integer :: nRbins=50
  integer :: nMaxHCoord=8
  real*8  :: Etot, Etot_pd, Etot_h, E_pd, E_h
  integer,DIMENSION(:),ALLOCATABLE      :: nHCoordPdAtoms      ! array of number of coordinating Pd atoms per H atom (nHAtoms)
  integer,DIMENSION(:,:),ALLOCATABLE    :: idxHCoordPdAtoms    !idx of nearest Pd neighbour for each H atom (nMaxHCoord,nHAtoms)
  real*8, DIMENSION(:,:),ALLOCATABLE    :: drHCoordPdAtoms     !distance to the nearest Pd neighbour for each H atom (nMaxHCoord,nHAtoms)
  real*8, DIMENSION(:),ALLOCATABLE      :: pHConcRbin          ! probability density of H concentration on Rbin (nRbins)
  integer,DIMENSION(:),ALLOCATABLE      :: nHConcRbin          ! # of H atoms in the Rbin distance from COM (nRbins)
  real*8, DIMENSION(:),ALLOCATABLE      :: pPdConcRbin         ! probability density of Pd concentration on Rbin (nRbins)
  integer,DIMENSION(:),ALLOCATABLE      :: nPdConcRbin         ! # of Pd atoms in the Rbin distance from COM (nRbins)
  integer,DIMENSION(:),ALLOCATABLE      :: idxPdRbin           ! rbin idx of Pd atoms = drCOM/nRbins (nPdAtoms)
  integer,DIMENSION(:),ALLOCATABLE      :: idxHRbin            ! rbin idx of H atoms = drCOM/nRbins (nHAtoms)
  real*8, allocatable,dimension(:)      :: phisum_hpd          !sum of phi_pd for h ptl i due to pd in system (nHAtoms)

  integer,DIMENSION(:,:),ALLOCATABLE    :: nHConcRbinTraj       ! # of H atoms in the Rbin distance from COM (nRbins,nFrames)
  integer,DIMENSION(:),ALLOCATABLE      :: nHRbinSum            ! # of H atoms in the Rbin distance from COM (nRbins)
  real*8, allocatable,dimension(:,:)    :: rho_pdh_traj         ! neighboring H concentration on each Pd atom (nPdAtoms,nFrames)
  real*8, allocatable,dimension(:,:)    :: rho_hpd_rbin_traj    ! rho_hpd trajectory (nRbins, nFrames)
  real*8, allocatable,dimension(:,:)    :: rho_hh_rbin_traj     ! rho_hh trajectory (nRbins, nFrames)
  real*8, allocatable,dimension(:)      :: rho_pdh_avg          ! average neighboring H concentration on each Pd atom (nPdAtoms)
  real*8, allocatable,dimension(:)      :: rho_hpd_rbin_avg     ! rho_hpd trajectory (nRbins)
  real*8, allocatable,dimension(:)      :: rho_hh_rbin_avg      ! rho_hh trajectory (nRbins)
  real*8, allocatable,dimension(:,:)    :: f_pd_traj            ! f_pd trajectory (nPdAtoms,nFrames)
  real*8, allocatable,dimension(:,:)    :: f_h_rbin_traj        ! f_h trajectory (nRbins, nFrames)
  real*8, allocatable,dimension(:)      :: f_pd_avg             ! average f_pd (nPdAtoms)
  real*8, allocatable,dimension(:)      :: f_h_rbin_avg         ! average f_h (nRbins)
  real*8, allocatable,dimension(:,:)    :: phisum_pdh_traj      ! sum of phi_pd for pd ptl i due to h in system (nPdAtoms,nFrames)
  real*8, allocatable,dimension(:)      :: phisum_pdh_avg       ! average sum of phi_pd for pd ptl i due to h in system (nPdAtoms)
  real*8, allocatable,dimension(:,:)    :: phisum_hpd_rbin_traj ! phisum_hpd trajectory (nRbins, nFrames)
  real*8, allocatable,dimension(:,:)    :: phisum_hh_rbin_traj  ! phisum_hh trajectory (nRbins, nFrames)
  real*8, allocatable,dimension(:)      :: phisum_hpd_rbin_avg  ! phisum_hpd trajectory (nRbins)
  real*8, allocatable,dimension(:)      :: phisum_hh_rbin_avg   ! phisum_hh trajectory (nRbins)

  integer,allocatable,dimension(:)      :: typeHSite            ! type of Hsite (nHAtoms) (1:Octahedral, 2:Tetrahedral, 3:InBetween, 4:Surface, 5:Crowded, 6:etc)
  integer,allocatable,dimension(:)      :: pdGroupHSite         ! Pd group of NN Pd per H atom (nHAtoms) (1,2, or else)
  integer,dimension(6)      :: nHType                    ! # of H atoms in each type (1:Octahedral, 2:Tetrahedral, 3:InBetween, 4:Surface, 5:Crowded, 6:etc)
  integer,dimension(15)      :: nHPdGroupNN               ! # of H atoms NN to pd atoms in each pd group (1,2, or else, grp1:type1-6, grp2:type1-6)
  integer,dimension(15)      :: nHPdGroupCoord            ! # of H atoms coordinated with pd atoms in each pd group (1,2, or else, grp1:type1-6, grp2:type1-6) 1 or 2 only if all coordinating Pd atoms are in group 1,2

  real*8 :: rScale = 0.3479
  integer :: ngridvec = 230, nslice = 31
  real*8, parameter :: sliceIdxs(30) = (/209. , 202.5, 196. , 189.5, 182.5, 176. , 169.5, 163. , 156.5, &
       150. , 143.5, 136.5, 130. , 123.5, 117. , 110.5, 104. ,  97.5, &
        91. ,  84.5,  78. ,  71.5,  64.5,  58. ,  51.5,  45. ,  39. , &
        33. ,  26.5,  19.5/)
  integer,dimension(230,230,31) :: num_h_grid ! number of H atoms on this grid point
  real*8, dimension(230,230,31) :: rhoh_sum   ! output htrajectory
  logical :: bUseHMap

END MODULE variables
