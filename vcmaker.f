!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                          !!
!!    VCMaker -- Vibronic Coupling Maker                    !!
!!                                                          !!
!!    Extraction of parameters for the construction of      !! 
!!    vibronic coupling Hamiltonians.                       !!
!!                                                          !!
!!         Julien Eng                                       !!
!!                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        PROGRAM VCMAKER
                                                               
        use mod_input
        use mod_constants
        use mod_grad
        use mod_disp
        use mod_gap
        use mod_dncscan
        use mod_lambda
        use mod_min
        use mod_err
        IMPLICIT NONE
        integer::i,j
        character(len=666)::input

!!!!!!!!!!!!
!! HEADER !!
!!!!!!!!!!!!

        CALL HEADER()  ! Printing the fancy header!

!!!!!!!!!!!!!!!!!!!!!!!
!! READING THE INPUT !!
!!!!!!!!!!!!!!!!!!!!!!!

        call getarg(1,input)   ! Reading the input file name
        if (input.eq."") then
                write(err_msg,*) 'No input file specified. Exiting now.'
                call mk_error()
        end if
        call get_input(input)  ! Reading the job specifications

        call get_atn(trim(adjustl(file_gs)),n_at)                ! Reading the No of atoms of the molecule
        allocate(coord_gs(3,n_at)); coord_gs(:,:)=0.0d0
        allocate(label(n_at)); label(:)="AA"
        call get_xyz(coord_gs,label,n_at,trim(adjustl(file_gs))) ! Reading the coordinates of the Ground State geometry 
        allocate(mass(n_at)); mass(:)=0.0d0
        call get_mass(label,mass,n_at)                           ! Reading the mass for each atom
        call get_com(coord_gs,mass,n_at,com)                     ! Computing the centre of mass
        call get_inertia(coord_gs,mass,n_at,com,inertia_axis)    ! computing the axes of Inertia
        allocate(hess(3*n_at,3*n_at)); hess(:,:)=0.0d0
        allocate(mw_hess(3*n_at,3*n_at)); mw_hess(:,:)=0.0d0
        call get_hessian(trim(adjustl(hessian)),hess)            ! Reading the Cartesian Hessian

        ! Printing time
        call out_input 

!! The following is always performed
! Hessian mass-weighting
        allocate(freq(3*n_at)); freq(:)=0.0d0
        allocate(mwh_evec(3*n_at,3*n_at)); mwh_evec(:,:)=0.0d0
        allocate(nm(3*n_at,3*n_at)); nm(:,:)=0.0d0
        allocate(red_mass(3*n_at)); red_mass(:)=0.0d0
        call mk_hess(hess,mass,mw_hess,coord_gs,com,freq,mwh_evec,n_at)  ! Transform the Cartesian Hessian to the Mass-weighted Cartesian Hessian
        call mk_nm(mwh_evec,nm,red_mass,mass,n_at)                       ! Generate the unitary displacements of mass-weighted normal modes
        write(*,*)
        write(*,*) '================================='
        Write(*,*) '-- Mass-weighted Normal Modes'
        write(*,*) '================================='
        write(*,*)
        call out_sixcol(nm,n_at)                                        ! Printing the normal modes
        allocate(dnc(3*n_at,3*n_at)); dnc(:,:)=0.0d0
        allocate(tdnc(3*n_at,3*n_at)); tdnc(:,:)=0.0d0
        call mk_dnc(dnc,tdnc,freq,mwh_evec,mass,n_at)                   ! Transform m.w.-normal modes to dimensionless normal coordinates (DNC)
        write(*,*)
        write(*,*) '================================='
        Write(*,*) '-- Dimensionless Normal Coordinates (DNC)'
        write(*,*) '================================='
        write(*,*) ' Cartesian displacements in Bohr associated to each DNC'
        write(*,*)
        call out_sixcol(tdnc,n_at)                                      ! Printing the DNC


        !! Specific jobs
        if (log_grad) then                                              ! Estimation of Kappa from dE/dx
           ! Allocating the arrays to be used
           allocate(grad(3*n_at,n_grad))
           allocate(kappa_grad(3*n_at,n_grad))
           allocate(eshift_grad(n_grad))
           allocate(grad_tmp(3*n_at))
           allocate(kgrad_tmp(3*n_at))
        
           do i=1,n_grad                                                ! Loop over the number of gradients provided
               grad_tmp(:)=0.0d0
               call get_grad(trim(adjustl(file_grad(i))),grad_tmp,n_at) ! Read the gradient from formatted file
               grad(:,i)=grad_tmp(:)
               kgrad_tmp(:)=0.0d0
               call mk_kappa_grad(freq,grad_tmp,kgrad_tmp,tdnc,degrad,n_at) ! Project the Cartesian gradient on the DNC
               kappa_grad(:,i)=kgrad_tmp(:)
               
               if (i.eq.1) then                                         ! Printing Loop. One table/electronic state
                  write(*,*)
                  write(*,*) '================================='
                  Write(*,*) '-- LVC Analysis from Gradients'
                  write(*,*) '================================='
                  write(*,*)
               end if
               write(*,*) 'File:  ', TRIM(ADJUSTL(file_grad(i)))
               write(*,*)
               call out_kgrad(freq,kgrad_tmp,degrad,red_mass,n_at)      ! Printing subroutine
               write(*,*)
               if (log_min) then                                        ! If theoretical minimum is requested
                if (min_n.eq.-1) then
                        min_n=3*n_at-6
                        allocate(min_mode(min_n))
                        do j=1,min_n ; min_mode(j)=j ; end do
                end if
                call mk_min(coord_gs,label,tdnc,kgrad_tmp,freq,n_at,'MIN-'//TRIM(ADJUSTL(file_grad(i))))
               end if
          end do
        end if

        if (log_disp) then                                              ! Estimation of Kappa from dx
        ! Allocating the arrays to be used
           allocate(coord_es(3,n_at,n_geo))
           allocate(kappa_dist(3*n_at,n_geo))
           allocate(kappa_disps(3*n_at,n_geo))
           allocate(eshift_disp(3*n_at,n_geo))
           allocate(coord_tmp(3,n_at))
           allocate(kdisp_tmp(3*n_at))
           allocate(edisp_tmp(3*n_at))
           allocate(labels_tmp(n_at))

                do i=1,n_geo                                            ! Loop over the number of geometries
                        coord_tmp(:,:)=0.0d00                           ! Reinitialising the temp mat
                        kdisp_tmp(:)=0.0d0
                        edisp_tmp(:)=0.0d0
                        call get_xyz(coord_tmp,labels_tmp,n_at,trim(adjustl(file_geo(i)))) ! Read the geometries 
        
                        do j=1,n_at                                     ! Checking for consistency
                                if (labels_tmp(j).ne.label(j)) then
                                        write(err_msg,*) 'Atoms read in file: ',TRIM(ADJUSTL(file_geo(i))),&
                                        &' are not in the same order!'
                                        call mk_error()
                                end if
                        end do
                        CALL alignxyz(coord_gs,coord_tmp,label,n_at)    ! Align the molecules using the Pulay approach
                        coord_es(:,:,i)=coord_tmp(:,:)

                        CALL mk_kappa_disp(freq,coord_tmp,coord_gs,kdisp_tmp,dnc,edisp_tmp,n_at) ! Projecting the dx to dQ
                        !! Copy back kappas to big matrice
                        kappa_dist(:,i)=kdisp_tmp(:)
                        eshift_disp(:,i)=edisp_tmp(:)
                        if (i.eq.1) then                                ! Printing loop. One table/Geometry
                                write(*,*)
                                write(*,*) '=========================================='
                                Write(*,*) '-- LVC Analysis from Displaced Geometries '
                                write(*,*) '=========================================='
                                write(*,*)
                        end if
                        write(*,*) 'File:  ', TRIM(ADJUSTL(file_geo(i)))
                        write(*,*)
                call out_kdisp(freq,kdisp_tmp,edisp_tmp,red_mass,n_at) ! Printing subroutine
                        write(*,*)

                        coord_tmp(:,:)=0.0d00                           ! Reinitialising the temp mat
                end do 
        end if

        if (log_gap) then
        !!! Allocation section
        allocate(gap(n_gap))
        allocate(rsma(3*n_at,n_gap))
        !!!!
                do i=1,n_gap

                        if (i.eq.1) then                                ! Printing loop. One table/Geometry
                                write(*,*)
                                write(*,*) '=========================================='
                                Write(*,*) '-- Anharmonicity Analysis  '
                                write(*,*) '=========================================='
                                write(*,*)
                        end if
                        write(*,*) 'Gradiant file : ', TRIM(ADJUSTL(file_grad(gap_grad(i))))
                        write(*,*) 'Geometry file : ', TRIM(ADJUSTL(file_geo(gap_disp(i))))
                        write(*,*)
                        
                        call mk_gap(kappa_grad(:,gap_grad(i)),&
                                   &kappa_dist(:,gap_disp(i)),&
                                   &gap(i),&
                                   &rsma(:,i),&
                                   &n_at)
                        write(*,'(A,F8.5)') 'Global Anharmonicity Parameter:   GAP =',gap(i)
                        write(*,*)
                        call out_rsma(freq,rsma(:,i),red_mass,n_at) 
                end do
        end if

        if (log_lambda) then                                            ! Estimation of lambdas from FC Hessians
        allocate(freq_es(3*n_at)) ; freq_es(:)=0.0d0
        allocate(hess_es_tmp(3*n_at,3*n_at)) ; hess_es_tmp(:,:)=0.0d0
        allocate(scratch_a(3*n_at,3*n_at)) ; scratch_a(:,:)=0.0d0
        allocate(scratch_b(3*n_at,3*n_at)) ; scratch_b(:,:)=0.0d0
        allocate(hess_es(3*n_at,3*n_at,n_hess)) ; hess_es_tmp(:,:)=0.0d0
        allocate(lambda(3*n_at,n_lambda))
        allocate(lambda_up(3*n_at,n_lambda))
        allocate(lambda_down(3*n_at,n_lambda))
                do i=1,n_hess
                        call get_hessian(trim(adjustl(file_eshess(i))),hess_es_tmp) ! Reading the Cartesian Hessian
! correction of the excited-state hessian
                call mk_hess(hess_es_tmp,mass,scratch_a,coord_gs,com,freq_es,scratch_b,n_at) ! Transforming the Cartesian Hessian to the GS DNC.
                        hess_es(:,:,i)=hess_es_tmp(:,:) 
                end do
        deallocate(scratch_a)
        deallocate(scratch_b)
               do i=1,n_lambda                                          ! Loop over the lambdas
        call mk_lambda(hess_es(:,:,lambda_pairs(i,1))&
                     &,hess_es(:,:,lambda_pairs(i,2))&
                     &,hess& 
                     &,vib_nrj(lambda_pairs(i,1))/hartoeV&              ! We convert back the energies to Eh...
                     &,vib_nrj(lambda_pairs(i,2))/hartoeV&              ! Otherwise the total unit is wrong!
                     &,tdnc,lambda(:,i)&
                     &,lambda_up(:,i)&
                     &,lambda_down(:,i)&
                     &,n_at)

                       if (i.eq.1) then                                 ! Printing moment. One table/coupling
                          write(*,*)
                          write(*,*) '================================='
                          Write(*,*) '-- Interstate Vibronic Coupling'
                          write(*,*) '================================='
                          write(*,*)
                        end if
                        write(*,*) 'File:  ', TRIM(ADJUSTL(file_eshess(lambda_pairs(i,1))))
                        write(*,*) 'File:  ', TRIM(ADJUSTL(file_eshess(lambda_pairs(i,2))))
                        write(*,*)
                call out_lambda(freq,lambda(:,i),lambda_up(:,i),lambda_down(:,i),red_mass,n_at)
                        write(*,*)


                end do
        end if


!!! Structure generation
        if (log_dncscan) then
                do i=1,n_dnc
! Scan along a given  DNC
                   call mk_dncscan(coord_gs&             ! Initial Geometry
                        &,label&                         ! Atom labels
                        &,tdnc&                          ! DNC displacements
                        &,dncscan_ptleft(i)&             ! #Points -DNC
                        &,dncscan_stepleft(i)&           ! Increment
                        &,dncscan_ptright(i)&            ! #Points +DNC
                        &,dncscan_stepright(i)&          ! Increment
                        &,dncscan_mode(i)&               ! DNC label
                        &,n_at)                          ! #Atom
                        if (log_dncdiag.AND.i.lt.n_dnc) then
                                do j=i+1,n_dnc
                        call mk_dncscan_diag(coord_gs&   ! Initial Geometry
                        &,label&                         ! Atom labels
                        &,tdnc&                          ! DNC displacements
                        &,dncscan_ptleft(i)&             ! #Points -DNC
                        &,dncscan_stepleft(i)&           ! Increment
                        &,dncscan_ptright(i)&            ! #Points +DNC
                        &,dncscan_stepright(i)&          ! Increment
                        &,dncscan_mode(i)&               ! 1st DNC label
                        &,dncscan_ptleft(j)&             ! #Points -DNC
                        &,dncscan_stepleft(j)&           ! Increment
                        &,dncscan_ptright(j)&            ! #Points +DNC
                        &,dncscan_stepright(j)&          ! Increment
                        &,dncscan_mode(j)&               ! 2nd DNC label 
                        &,n_at)                          ! #Atom
                                end do
                        end if
                end do
        end if


        if (log_dncgrid) then
                call mk_dncscan_grid(coord_gs&          ! Initial Geometry
                        &,label&                        ! Atom labels
                        &,tdnc&                         ! DNC displacements
                        &,grd_mode&                     ! 1st DNC label
                        &,grd_ptleft&                   ! #Points -DNC
                        &,grd_stepleft&                 ! Increment
                        &,grd_ptright&                  ! #Points +DNC
                        &,grd_stepright&                ! Increment
                        &,n_dim&
                        &,n_at) 
        end if

        if (log_genxyz) then
                call mk_genxyz(coord_gs&                ! Initial Geometry
                 &,label&                               ! Atom labels
                 &,tdnc&                                ! DNC displacements
                 &,gen_mode&                            ! 1st DNC label
                 &,gen_disto&                           ! #Points -DNC
                 &,n_gen&
                 &,n_at&
                 &,'genxyz_out.xyz') 
        end if

        if (log_quantics) then
                if (log_grad) then
                call out_quantics(n_grad)
                else if (log_disp) then
                call out_quantics(n_geo)
                end if
        end if

        call endoftimes()                               ! All done !

END
