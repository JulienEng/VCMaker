Subroutine mk_dncscan_diag(mkscan_coord&      ! Initial Geometry
                    &,mkscan_label&           ! Atom labels
                    &,mkscan_tdnc&            ! DNC displacements
                    &,mkscan_ptleft&          ! #Points -DNC
                    &,mkscan_stleft&          ! Increment
                    &,mkscan_ptright&         ! #Points +DNC
                    &,mkscan_stright&         ! Increment
                    &,mkscan_modea&           ! DNC label
                    &,mkscan_ptleftb&         ! #Points -DNC for 2nd mode
                    &,mkscan_stleftb&         ! Increment for 2nd mode
                    &,mkscan_ptrightb&        ! #Points +DNC for 2nd mode
                    &,mkscan_strightb&        ! Increment for 2nd mode
                    &,mkscan_modeb&           ! DNC label for 2nd mode
                    &,mkscan_n)               ! #Atom
use mod_constants
            implicit none
            integer::k,l,m                                                              !! Loop integers
            integer,intent(in)::mkscan_n                                                !! Number of modes to scan
            character(len=2),dimension(mkscan_n),intent(in)::mkscan_label               !! Atom labels
            double precision,dimension(3,mkscan_n),intent(in)::mkscan_coord             !! Molecular Cartesian coordinates
            double precision,dimension(3*mkscan_n,3*mkscan_n),intent(in)::mkscan_tdnc   !! Transpose DNC matrice
            integer,intent(in)::mkscan_ptleft                                           !! Number of points in negative distortion for mode A
            integer,intent(in)::mkscan_ptright                                          !! Number of points in positive distortion for mode A
            double precision,intent(in)::mkscan_stleft                                  !! Increment in negative distortion for mode A
            double precision,intent(in)::mkscan_stright                                 !! Increment in positive distortion for mode A
            integer,intent(in)::mkscan_modea                                            !! Index of mode A
            ! Second mode
            integer,intent(in)::mkscan_ptleftb                                          !! Number of points in negative distortion for mode B
            integer,intent(in)::mkscan_ptrightb                                         !! Number of points in positive distortion for mode B
            double precision,intent(in)::mkscan_stleftb                                 !! Increment in negative distortion for mode B
            double precision,intent(in)::mkscan_strightb                                !! Increment in positive distortion for mode B
            integer,intent(in)::mkscan_modeb                                            !! Index of mode B

            double precision,dimension(3,mkscan_n)::mkscan_outcoord                     !! Cartesian coordinates of the distorted geometry
            character(len=666)::mkscan_out                                              !! Output path 
            integer::mkscan_counter                                                     !! Loop Counter 
            character(len=100)::scratch                                                 !! Just a string.
            character(len=100)::DIR                                                     !! Output directory
            character(len=100)::XYZFILE                                                 !! Output file
            integer::mkscan_ptl                                                         !! Harmonised number of point in the negative distortion
            integer::mkscan_ptr                                                         !! Harmonised number of point in the positive distortion

            !! The number of points in the negative distortion must 
            !! be the same for both modes
            if (mkscan_ptleft.ne.mkscan_ptleftb) then
            !! If not, we select the larger number of points.
                mkscan_ptl=MAX(mkscan_ptleft,mkscan_ptleftb)
            else
                mkscan_ptl=mkscan_ptleft
            end if
            !! The number of points in the positive distortion must 
            !! be the same for both modes
            if (mkscan_ptright.ne.mkscan_ptrightb) then
                !! If not, we select the larger number of points.
                mkscan_ptr=MAX(mkscan_ptright,mkscan_ptrightb)
            else
                mkscan_ptr=mkscan_ptright
            end if

            !! Initialisation of the counter to give the output its name.
            mkscan_counter=0 

            !! Creation of the Directory where the geometries will be saved.
            write(DIR,'(A9,i0,A1,i0)') 'DNC_SCAN_',mkscan_modea,'_',mkscan_modeb
            write(scratch,*) 'mkdir -p '//TRIM(ADJUSTL(DIR))
            call system(trim(adjustl(scratch))) 

            !! Loop over the number of points on the negative distortion
            do k=1,mkscan_ptl
            !! Reinitialisation of the output coordinates.
            mkscan_outcoord(:,:)=0.0d0
            !! Incrementation of the counter.
            mkscan_counter=mkscan_counter+1
                        
            !! Generating the filename
            write(XYZFILE,'(A3,i0,A4,i0,A1,i0,A4)') 'DNC',mkscan_modea,'_DNC',mkscan_modeb,'_',mkscan_counter,'.xyz'
            WRITE(mkscan_out,'(A)') TRIM(ADJUSTL(DIR))//'/'//TRIM(ADJUSTL(XYZFILE)) 
            !! Loop over the atoms
            do l=1,mkscan_n 
            !! Loop over X,Y and Z
                  do m=1,3     
            !! Generation of the distorted geometry. 
            !! By convention, mode number n in the VCMaker input corresponds to the (n+6)th column on the matrice
                       mkscan_outcoord(m,l)=mkscan_coord(m,l)+&
                       &(-mkscan_ptl+k-1)*(&
                       &mkscan_stleft*mkscan_tdnc(mkscan_modea+6,((l-1)*3+m))+&
                       &mkscan_stleftb*mkscan_tdnc(mkscan_modeb+6,((l-1)*3+m)))&
                       &/angtoa0
                  end do
               end do
            !! Writing the geometry!
               call out_xyz(mkscan_outcoord,mkscan_label,mkscan_n,trim(adjustl(mkscan_out)))
            end do

            !! Incrementation of the counter.
            mkscan_counter=mkscan_counter+1

            !! Generating the filename
            write(XYZFILE,'(A3,i0,A4,i0,A1,i0,A4)') 'DNC',mkscan_modea,'_DNC',mkscan_modeb,'_',mkscan_counter,'.xyz'
            WRITE(mkscan_out,'(A)') TRIM(ADJUSTL(DIR))//'/'//TRIM(ADJUSTL(XYZFILE)) 
            !! Writing the geometry!
            call out_xyz(mkscan_coord,mkscan_label,mkscan_n,trim(adjustl(mkscan_out)))

            !! Loop over the number of points on the positive distortion
            do k=1,mkscan_ptr
             !! Reinitialisation of the output coordinates.
            mkscan_outcoord(:,:)=0.0d0
            !! Incrementation of the counter.
            mkscan_counter=mkscan_counter+1

            !! Generating the filename
            write(XYZFILE,'(A3,i0,A4,i0,A1,i0,A4)') 'DNC',mkscan_modea,'_DNC',mkscan_modeb,'_',mkscan_counter,'.xyz'
            WRITE(mkscan_out,'(A)') TRIM(ADJUSTL(DIR))//'/'//TRIM(ADJUSTL(XYZFILE)) 
            !! Loop over the atoms
               do l=1,mkscan_n 
            !! Loop over X,Y and Z
                  do m=1,3     
            !! Generation of the distorted geometry. 
            !! By convention, mode number n in the VCMaker input corresponds to the (n+6)th column on the matrice
                       mkscan_outcoord(m,l)=mkscan_coord(m,l)+k*(&
                       &mkscan_stright*mkscan_tdnc(mkscan_modea+6,((l-1)*3+m))+&
                       &mkscan_strightb*mkscan_tdnc(mkscan_modeb+6,((l-1)*3+m)))&
                       &/angtoa0
                  end do
               end do
               !! Writing the geometry!
               call out_xyz(mkscan_outcoord,mkscan_label,mkscan_n,trim(adjustl(mkscan_out)))
            end do


end subroutine



