!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   module mod_disp                                        !!!
!!                                                        !!! 
!! Variables requires if log_disp=.TRUE.                  !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer::n_geo                                              !#of geometries
character(len=666),dimension(:),allocatable::file_geo       ! Files containing the excited states geometries
double precision,dimension(:,:),allocatable::coord_tmp      ! TEMP matrix to read coordinates
character(len=2),dimension(:),allocatable::labels_tmp       ! TEMP matrix to read labels
double precision,dimension(:),allocatable::kdisp_tmp        ! TEMP matrix to store Kappas
double precision,dimension(:),allocatable::edisp_tmp        ! TEMP matrix to store energy shifts
double precision,dimension(:,:,:),allocatable::coord_es     ! Excited States structures geometries
double precision,dimension(:,:),allocatable::kappa_dist     ! Intrastate coupling
double precision,dimension(:,:),allocatable::kappa_disps    ! Displacements associated with Kappa
double precision,dimension(:,:),allocatable::eshift_disp    ! Energy shift associated with Kappa
double precision,dimension(:,:),allocatable::kappa_star     ! ???
double precision,dimension(:,:),allocatable::omega_star     ! ???

end module mod_disp
