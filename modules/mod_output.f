!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  module mod_output                                         !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer::n_mode,n_vmd
integer,dimension(:),allocatable::list_mode,list_vmd
!!! Output Control ---
logical::log_out_avo                                          !out Avogadro file 
logical::log_out_quant                                        !out quantics op file
logical::log_out_vmd                                          !out VMD files
logical::log_out_avo_vsa                                      !out Avogadro VSA files
integer::out_prtlvl                                           ! Not yet implemented
logical::log_out_scannm                                       !out NM distort. geo
logical::log_out_scanmwnm                                     !out MWNM distort.  geo
logical::log_out_scandnc                                      !out DNC distort. geo
!!! Output integer ---
!!integer::n_scannm    !#of nm distorted geo to print
!!integer::n_scanmwnm  !#of mw-nm distorted geo to print
!!integer::n_scandnc   !#of dnc distorted geo to print
double precision,dimension(:,:),allocatable::list_scannm   ! List for NM distort.
double precision,dimension(:,:),allocatable::list_scanmwnm ! List for MWNM distort.
double precision,dimension(:,:),allocatable::list_scandnc  ! List for DNC distort.

!Writing 
integer,dimension(:),allocatable::write_nm
integer,dimension(:),allocatable::write_dnc
integer,dimension(:),allocatable::write_vsa_nm
integer,dimension(:),allocatable::write_vsa_dnc

end module mod_output
