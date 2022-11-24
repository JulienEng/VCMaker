SUBROUTINE alignxyz(coords_init,coords_fin,labels_init,n_at)
      IMPLICIT NONE
      !Cart Geometry
      integer,intent(in)::n_at
      double precision,dimension(3,n_at),intent(in)::coords_init
      double precision,dimension(3,n_at),intent(inout)::coords_fin
      character(len=2),dimension(n_at),intent(in)::labels_init
!      character(len=2),dimension(n_at)::labels_fin
!      character(len=666)::file_init
!      character(len=666)::file_fin

      !Internal Coordinates
      integer::n_bl  ! #Bond lengths
      integer::n_ba  ! #Bond angles
      integer::n_di  ! #Dihedral
      integer::n_int ! #Internal coordinates
      integer,dimension(:,:),allocatable::internal_def
      double precision,dimension(:),allocatable::internal_init
      double precision,dimension(:),allocatable::internal_fin
      character(len=16),dimension(:),allocatable::internal_labels

      !Internal Internal to Cartesian Back Transformation
      double precision,dimension(:,:),allocatable::bmat
      double precision,dimension(:,:),allocatable::amat
      double precision::rmsd
      double precision,dimension(:),allocatable::dQ
      double precision,dimension(:),allocatable::internal_tmp
      double precision,dimension(:),allocatable::delta_init
      double precision,dimension(:),allocatable::delta_tmp
      double precision,dimension(:),allocatable::xyz_tmp
      double precision,dimension(:),allocatable::xyz_tmpa
      double precision,dimension(:,:),allocatable::coords_tmp
      double precision,dimension(:,:),allocatable::coords_rmsd
      integer::counter

      integer::i  !,j
      double precision,parameter::pi=4.0d0*atan(1.0d0)

      call get_rmsd(coords_init,coords_fin,n_at,rmsd)
!
!      ! Output Cartesian Coordinates
!        write(*,*) '+------------------------+'
!        write(*,*) '| Cartesian Coordinates: |'
!        write(*,*) '+------------------------+'
!        write(*,*)
!        write(*,'(2(A12,A20,A3))') 'Geometry 1: ',adjustl(file_init),' | ', 'Geometry 2: ', adjustl(file_fin),'   '
!        write(*,'(A70)') '----------------------------------------------------------------------'
!        do i=1,n_at
!        write(*,'(2(A2,3F10.5,A3))') labels_init(i), (coords_init(j,i),j=1,3),' | ', labels_fin(i), (coords_fin(j,i),j=1,3),'   '
!        end do
!        write(*,*)
!        write(*,'(A6,E11.4)') 'RMSD: ',rmsd
!        write(*,*)

      ! Determine Internal coordinates
      ! What if a molecule is linear?

      n_bl=n_at-1
      n_ba=n_at-2
      n_di=n_at-3
     
      n_int=n_bl+n_ba+n_di ! n_int=3*n_at-6 !

      allocate(internal_def(4,n_int))
      allocate(internal_init(n_int))
      allocate(internal_fin(n_int))
      allocate(dQ(n_int))
      allocate(internal_labels(n_int))

      allocate(bmat(n_int,3*n_at))
      allocate(amat(3*n_at,n_int))

      ! Get the automatic definition of the internal coordinates
      call get_defint(n_at,coords_init,internal_def,internal_labels,n_bl,n_ba,n_di)
      
      ! Get the Internal coordinates for the initial geometry
      call get_internals(n_at,coords_init,internal_def,internal_init,n_bl,n_ba,n_di)
      ! Get the Internal coordinates for the final geometry
      call get_internals(n_at,coords_fin,internal_def,internal_fin,n_bl,n_ba,n_di)

!      ! Print Initial State
!        write(*,*) '+------------------------+'
!        write(*,*) '| Internal  Coordinates: |'
!        write(*,*) '+------------------------+'
!        write(*,*)
!        write(*,*) 'This set of coordinates is automatically built in the form of '
!        write(*,*) 'a non-redundant Z-Matrice.'
!        write(*,*) 'Atoms may be swapped to avoid having three aligned atoms.                              '
!        write(*,*)
!        write(*,'(A16,A3,A10,A10,A13,A8,A10)') 'Internal Coords.',' | ','Geo 1:     ','     | ',&
!                &'Geo 2:       ','      | ', 'Diff:     '
!        write(*,'(A91)') '-------------------------------------------------------------------------------------------'
!
!      do i=1,n_int
!      if (i.le.n_bl) then
!      write(*,'(A16,A3,2X,F10.3,A7,4X,F10.3,A7,F10.3,A4,A3,F7.2,A2)') internal_labels(i),' | ', internal_init(i),' Ang. | ',&
!                                            & internal_fin(i),' Ang. | ',&
!                                            & internal_fin(i)-internal_init(i),' Ang.','  (',&
!                                            & ((internal_init(i)-internal_fin(i))/internal_init(i))*100, '%)'
!      else 
!      write(*,'(A16,A3,2X,F10.3,A7,4X,F10.3,A7,F10.3,A4,A3,F7.2,A2)') internal_labels(i),' | ',&
!                                                & internal_init(i)*180.0d0/pi,' Deg. | ',&
!                                                & internal_fin(i)*180.0d0 /pi,' Deg. | ',&
!                                                & (internal_fin(i)-internal_init(i))*180.0d0 /pi,' Deg.',&
!                                                &'  (',((internal_init(i)-internal_fin(i))/internal_init(i))*100, '%)'
!      end if
!      end do
!
!!        write(*,*)
!        write(*,*)
!
!        write(*,*) '+-------------------------------------------+'
!        write(*,*) '| Internal to Cartesian Back transformation |'
!        write(*,*) '+-------------------------------------------+'
!        write(*,*)
!        write(*,*)
        allocate(delta_init(n_int))
        allocate(delta_tmp(n_int))
        allocate(internal_tmp(n_int))
        allocate(xyz_tmp(3*n_at))
        allocate(xyz_tmpa(3*n_at))
        allocate(coords_tmp(3,n_at))
        allocate(coords_rmsd(3,n_at))
        
        
        call xyz_nto3n(coords_init,xyz_tmp,n_at,1) ! Transform a (3,n_at) xyz matrix to a (3*n_at) one
        counter=0
        coords_tmp=coords_init

        do while (rmsd.gt.1E-6.AND.counter.le.20) !!!! Starting the optimisation loop
        counter=counter+1
        !Get the internal coordinates of the new structure
        call get_internals(n_at,coords_tmp,internal_def,internal_tmp,n_bl,n_ba,n_di) 

        if (counter.eq.1) then
        do i=1,n_int  ! First step
           delta_init(i)= internal_fin(i)-internal_init(i) ! Sq in Trygve paper
           delta_tmp=0.0d0
        end do
        else
        do i=1,n_int  ! all the remaining steps
           delta_tmp(i)= internal_tmp(i)-internal_init(i)
        end do
        end if

    
        dQ=delta_init-delta_tmp
        do i=1,n_int  ! all the remaining steps
            if (i.gt.n_bl+n_ba) then ! Correction of going over the -pi:pi range
                    if (ABS(dQ(i)).gt.&
                     &(2*pi-ABS(dQ(i)))) then
                dQ(i)=-sign(2*pi-ABS(dQ(i)),dQ(i))
             end if
            end if
        end do

        !! Generate a new B-Matrix and A-Matrix from the new geometry
        call mk_bmat(bmat,n_at,n_bl,n_ba,n_di,xyz_tmp,internal_def)   !! xyz_tmp needs to be set
        call mk_amat(amat,bmat,n_at,labels_init,n_int)

        xyz_tmp=xyz_tmp+MATMUL(amat,dQ)


        call xyz_nto3n(coords_tmp,xyz_tmp,n_at,-1) ! Transform a (3*n_at) xyz matrix to a (3,n_at) one
        call get_rmsd(coords_tmp,coords_rmsd,n_at,rmsd)

        coords_rmsd=coords_tmp !! We save the new geometry for the RMSD calculation of next step

!!         write(*,'(A5,I3,A9,E18.9)') 'Step ',counter,'  RMSD:  ',rmsd
        end do

        call get_internals(n_at,coords_tmp,internal_def,internal_tmp,n_bl,n_ba,n_di)
        coords_fin=coords_tmp

!!         write(*,*)
!!         write(*,*) 'Convergence reached!'
!! 
!!         write(*,*) 'Last geometry: Internal coordinates'
!!         write(*,*)
!!         write(*,'(A16,A3,A10,A12,A16)') 'Internal Coords.',' | ','Geo 2:     ','        |   ','Aligned Geo:    '
!!         write(*,'(A69)') '---------------------------------------------------------------------'
!!         do i=1,n_int
!!         if (i.le.n_bl) then
!!         write(*,'(A16,A3,2X,F10.3,A7,4X,F10.3,A4,A3,F7.2,A2)') internal_labels(i),' | ', internal_fin(i),' Ang. | ',&
!!                                                              & internal_tmp(i),' Ang.',&
!!                                               &'  (',((internal_fin(i)-internal_tmp(i))/internal_fin(i))*100, '%)'
!!         else
!!         write(*,'(A16,A3,2X,F10.3,A7,4X,F10.3,A4,A3,F7.2,A2)') internal_labels(i),' | ', internal_fin(i)*180.0d0/pi,' Deg. | ',&
!!                                                                             & internal_tmp(i)*180.0d0 /pi,' Deg.',&
!!                                                 &'  (',((internal_fin(i)-internal_tmp(i))/internal_tmp(i))*100, '%)'
!!         end if
!!         end do
!! 
!!         write(*,*)
!!         write(*,*)
!!         write(*,*) 'Last geometry: Cartesian coordinates'
!!         write(*,*)
!!         write(*,*)
!!         write(*,'(A12,A20,A3,A13)') 'Geometry 2: ',adjustl(file_init),' | ', 'Aligned Geo: '
!!         write(*,'(A70)') '----------------------------------------------------------------------'
!!         do i=1,n_at
!!         write(*,'(2(A2,3F10.5,A3))') labels_init(i), (coords_fin(j,i),j=1,3),' | ', labels_fin(i), (coords_tmp(j,i),j=1,3),'   '
!!         end do

      deallocate(internal_def)
      deallocate(internal_init)
      deallocate(internal_fin)
      deallocate(dQ)
      deallocate(internal_labels)

      deallocate(bmat)
      deallocate(amat)

      deallocate(delta_init)
      deallocate(delta_tmp)
      deallocate(internal_tmp)
      deallocate(xyz_tmp)
      deallocate(xyz_tmpa)
      deallocate(coords_tmp)
      deallocate(coords_rmsd)

END SUBROUTINE
