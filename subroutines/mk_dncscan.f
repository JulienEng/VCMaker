      Subroutine mk_dncscan(mkscan_coord&      ! Initial Geometry
                          &,mkscan_label&      ! Atom labels
                          &,mkscan_tdnc&       ! DNC displacements
                          &,mkscan_ptleft&     ! #Points -DNC
                          &,mkscan_stleft&     ! Increment
                          &,mkscan_ptright&    ! #Points +DNC
                          &,mkscan_stright&    ! Increment
                          &,mkscan_mode&       ! DNC label
                          &,mkscan_n)          ! #Atom
            use mod_constants
            implicit none
            integer::k,l,m
            integer,intent(in)::mkscan_n
            character(len=2),dimension(mkscan_n),intent(in)::mkscan_label
            double precision,dimension(3,mkscan_n),intent(in)::mkscan_coord
            double precision,dimension(3*mkscan_n,3*mkscan_n),intent(in)::mkscan_tdnc
            integer,intent(in)::mkscan_ptleft
            integer,intent(in)::mkscan_ptright
            double precision,intent(in)::mkscan_stleft
            double precision,intent(in)::mkscan_stright
            integer,intent(in)::mkscan_mode

            double precision,dimension(3,mkscan_n)::mkscan_outcoord
            character(len=666)::mkscan_out
            integer::mkscan_counter
            character(len=100)::scratch
            character(len=100)::DIR,XYZFILE

            !! Initialisation of the counter (for filename writing)
            mkscan_counter=0 
            !! Creation of the Directory where the geometries will be saved.
            write(DIR,'(A9,i0)') 'DNC_SCAN_',mkscan_mode
            write(scratch,*) 'mkdir -p '//TRIM(ADJUSTL(DIR))
            call system(trim(adjustl(scratch))) ! Make directory


            !! Loop over the number of point towards the negative direction.
            do k=1,mkscan_ptleft
            mkscan_outcoord(:,:)=0.0d0
            !! Incrementing the counter
            mkscan_counter=mkscan_counter+1
            
            !! Generating the filename
            write(XYZFILE,'(A3,i0,A1,i0,A4)') 'DNC',mkscan_mode,'_',mkscan_counter,'.xyz'
            WRITE(mkscan_out,'(A)') TRIM(ADJUSTL(DIR))//'/'//TRIM(ADJUSTL(XYZFILE)) 
            
            !! Loop over the number of atom
              do l=1,mkscan_n 
               !! Loop over X,Y and Z
                  do m=1,3 
                       mkscan_outcoord(m,l)=mkscan_coord(m,l)+&
                       &(-mkscan_ptleft+k-1)*mkscan_stleft*mkscan_tdnc(mkscan_mode+6,((l-1)*3+m))/angtoa0
                  end do
               end do
               !! Writing the Geometry
               call out_xyz(mkscan_outcoord,mkscan_label,mkscan_n,trim(adjustl(mkscan_out)))
            end do

            !! Write the non-distorted geometry.
            mkscan_counter=mkscan_counter+1
            write(XYZFILE,'(A3,i0,A1,i0,A4)') 'DNC',mkscan_mode,'_',mkscan_counter,'.xyz'
            WRITE(mkscan_out,'(A)') TRIM(ADJUSTL(DIR))//'/'//TRIM(ADJUSTL(XYZFILE)) 
            call out_xyz(mkscan_coord,mkscan_label,mkscan_n,trim(adjustl(mkscan_out)))

            !! Loop over the number of point towards the positive direction.
            do k=1,mkscan_ptright
            mkscan_outcoord(:,:)=0.0d0
            !! Incrementing the counter
            mkscan_counter=mkscan_counter+1
            !! Generating the filename
            write(XYZFILE,'(A3,i0,A1,i0,A4)') 'DNC',mkscan_mode,'_',mkscan_counter,'.xyz'
            WRITE(mkscan_out,'(A)') TRIM(ADJUSTL(DIR))//'/'//TRIM(ADJUSTL(XYZFILE)) 
               !! Loop over the number of atom
            do l=1,mkscan_n 
               !! Loop over X,Y and Z
                  do m=1,3     
                       mkscan_outcoord(m,l)=mkscan_coord(m,l)+&
                       &k*mkscan_stright*mkscan_tdnc(mkscan_mode+6,((l-1)*3+m))/angtoa0
                  end do
               end do
               !! Write the non-distorted geometry.
               call out_xyz(mkscan_outcoord,mkscan_label,mkscan_n,trim(adjustl(mkscan_out)))
            end do


end subroutine



