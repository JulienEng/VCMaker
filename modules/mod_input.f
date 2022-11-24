!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module mod_input
!!
!! Contains the variables that required for 
!! any VCMaker runs
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   !! General Module
   !!!   Loop integers
      integer::rest
      integer::soft
   

   !!! Job controle variables --
      logical:: log_grad                                     !inp Grad calc.
      logical:: log_disp                                     !inp Disp. calc.
      logical:: log_gap                                      !inp GFI Anal.
      logical:: log_lambda                                   !inp vibronic coupling calc.
      logical:: log_vsa                                      !inp VSA calc.
      logical:: log_dncscan                                  !inp scan along DNC
      logical:: log_dncdiag                                  !inp scan along DNC
      logical:: log_dncgrid                                  !inp scan along DNC
      logical:: log_genxyz                                   !inp Generate XYZ structure distorted along N DNC.
      logical:: log_quantics                                 !inp Generate Quantics output
      logical:: log_min                                      !inp Generate "theoretical minimum

   !!! Molecular Geometry --
      integer::n_at                                          !#of Atoms
      character(len=666)::file_gs                            ! File containing the ground state geometry
      double precision,dimension(:,:),allocatable::coord_gs  ! Ground state structure coordinates
      character(len=2),dimension(:),allocatable::label       ! Label of the Atoms
      double precision,dimension(3)::com                     ! Centre of Mass
      double precision,dimension(3,3)::inertia_axis          ! Axis of Inertia

   !!! Molecular Hessian --
      character(len=666)::hessian                            ! File containing the hessian
      double precision,dimension(:,:),allocatable::hess      ! Hessian matrix
      double precision,dimension(:,:),allocatable::mw_hess   ! Mass-weighted Hessian matrix
      double precision,dimension(:,:),allocatable::hess_tmp
      double precision,dimension(:,:),allocatable::mwh_evec  ! Eigenvectors of the mass-weighted Hessian matrix
      double precision,dimension(:,:),allocatable::nm        ! Mass-weighted normal modes
      double precision,dimension(:,:),allocatable::dnc       ! Dimensionless normal coordinates
      double precision,dimension(:,:),allocatable::tdnc      ! Transposed matrix of hess_dnc
   !!! Masses ---
      double precision,dimension(:),allocatable::mass        ! Atomic mass
      double precision,dimension(:),allocatable::mass3n      ! Atomic mass -- Format for x,y,z/atom
      double precision,dimension(:),allocatable::red_mass    ! Reduced mass of the normal modes.
   !!! Frequencies ---
      double precision,dimension(:),allocatable::freq        ! Frequencies of the normal modes



      !Logical
      logical::vibronic

      end module mod_input

