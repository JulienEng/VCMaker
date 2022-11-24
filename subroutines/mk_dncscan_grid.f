Subroutine mk_dncscan_grid(grd_coord&      ! Initial Geometry
                    &,grd_label&           ! Atom labels
                    &,grd_tdnc&            ! DNC displacements
                    &,grd_mode&            ! DNC label for 2nd mode
                    &,grd_ptleft&          ! #Points -DNC
                    &,grd_stleft&          ! Increment
                    &,grd_ptright&         ! #Points +DNC
                    &,grd_stright&         ! Increment
                    &,grd_dim&             ! #Atom
                    &,grd_n)               ! #Atom
            use mod_constants
            implicit none
            integer::i,j
            integer::k,m
            integer,intent(in)::grd_n
            integer,intent(in)::grd_dim
            character(len=2),dimension(grd_n),intent(in)::grd_label
            double precision,dimension(3,grd_n),intent(in)::grd_coord
            double precision,dimension(3*grd_n,3*grd_n),intent(in)::grd_tdnc
            integer,dimension(grd_dim),intent(in)::grd_ptleft
            integer,dimension(grd_dim),intent(in)::grd_ptright
            double precision,dimension(grd_dim),intent(in)::grd_stleft
            double precision,dimension(grd_dim),intent(in)::grd_stright
            integer,dimension(grd_dim),intent(in)::grd_mode


            integer::n_pt
            integer::n_tmp,n_tmp2,n_tmp3
            double precision,dimension(:,:),allocatable::grd_disto
            double precision,dimension(grd_dim)::grd_indivdisto
           
            character(len=666)::grd_out
            character(len=100)::FMT_B,scratch
            character(len=100)::DIR

            !! The first step of the generation of a grd_dim-dimension grid of distorted geometries
            !! is to generate an array that contains on each line the distortion along each DNC
            !! for a given geometry
            !! Example: 3D grid, with distortion -1 and 1 along each DNC:
            !!
            !!    Geo:     DNC1   DNC2   DNC3
            !!      1      -1      -1      -1
            !!      2      -1      -1       0
            !!      3      -1      -1       1
            !!      4      -1       0      -1
            !!      5      -1       0       0
            !!      6      -1       0       1
            !!      .......
            !!      .......
            !!

            n_pt=1
            !! Computing the total number of points to generate
            do i=1,grd_dim
                n_pt=n_pt*(grd_ptleft(i)+grd_ptright(i)+1)
            end do
            !! Allocation of the distortion matrice
            !!     Number of modes v     v Total number of points
            allocate(grd_disto(grd_dim,n_pt))
        
            !! Initialisation of the grid
            grd_disto=0.0d0
           
                n_tmp=n_pt
                n_tmp2=1
            do i=1,grd_dim
            !! Number of points between to incrementation of the displacement along mode i
                n_tmp=n_tmp/(grd_ptleft(i)+grd_ptright(i)+1)

            !! If Mode 1
                if (i.eq.1) then
            !! Loop over the number of point along the first DNC
                       do j=1,(grd_ptleft(i)+grd_ptright(i)+1)
              !! If j< number of points in the negative direction: Negative distortion
                            if (j.le.grd_ptleft(i)) then
                           grd_disto(i,(j-1)*n_tmp+1)=-(grd_ptleft(i)-j+1)*grd_stleft(i)
              !! If j= number of points in the negative direction+1: No distortion
                             else if (j.eq.grd_ptleft(i)+1) then
                           grd_disto(i,(j-1)*n_tmp+1)=0.0d0
              !! If j> number of points in the negative direction +1: Positive distortion            
                             else if (j.gt.grd_ptleft(i)+1) then
                           grd_disto(i,(j-1)*n_tmp+1)=(j-(grd_ptleft(i)+1))*grd_stright(i)
                             end if
              !! We copy the distortion along mode i for all the points between two incrementations.
                           do m=1,n_tmp
                           grd_disto(i,(j-1)*n_tmp+m)=grd_disto(i,(j-1)*n_tmp+1)
                           end do
                        end do
                 else
            !! Number of "Cycle" we have to do on dimension i (== number of incrementation
            !! of mode i-1)
                         n_tmp2=n_tmp2*(grd_ptleft(i-1)+grd_ptright(i-1)+1)
            !! Spacing between two incrementations of mode i
                         n_tmp3=n_tmp*(grd_ptleft(i)+grd_ptright(i)+1)
            !! Loop over each incrementations of mode i-1
                         do k=1,n_tmp2
            !! Loop over the number of point along the DNC i
                                do j=1,(grd_ptleft(i)+grd_ptright(i)+1)
              !! If j< number of points in the negative direction: Negative distortion
                                if (j.le.grd_ptleft(i)) then
                                grd_disto(i,(k-1)*n_tmp3+((j-1)*n_tmp+1))=-(grd_ptleft(i)-j+1)*grd_stleft(i)
              !! If j= number of points in the negative direction+1: No distortion
                                else if (j.eq.grd_ptleft(i)+1) then
                                grd_disto(i,(k-1)*n_tmp3+((j-1)*n_tmp+1))=0.0d0
              !! If j> number of points in the negative direction +1: Positive distortion            
                                else if (j.gt.grd_ptleft(i)+1) then
                                grd_disto(i,(k-1)*n_tmp3+((j-1)*n_tmp+1))=(j-(grd_ptleft(i)+1))*grd_stright(i)

                                end if
                                if (i.ne.grd_dim) then
              !! We copy the distortion along mode i for all the points between two incrementations, but only
              !! if we have NOT reached the last mode of the grid.
                                do m=1,n_tmp
                                grd_disto(i,(k-1)*n_tmp3+(j-1)*n_tmp+m)=grd_disto(i,(k-1)*n_tmp3+(j-1)*n_tmp+1)
                                end do
                                end if
                                end do
                        end do
                 end if
            end do


            !! Making the directory where the Cartesian coordinates files will be stored
            write(DIR,*) 'DNC_GRID'
            write(scratch,*) 'mkdir -p '//TRIM(ADJUSTL(DIR))
            call system(trim(adjustl(scratch)))
            !! Generating the file where the coordinate grid will be printed.
                write(grd_out,*) TRIM(ADJUSTL(DIR))//'/DNC_GRID_OUTPUT.vcm'
                OPEN(667,file=trim(adjustl(grd_out)))
            !! Determining the write format for the coordinate grid.
                write(FMT_B,'(A4,i0,A6)') '(I4,',grd_dim,'F12.6)'
            !! Loop over the number of geometeries to build.
               do j=1,n_pt
            !! Copying the jth line of the grid matrice to a single vector
                grd_indivdisto(:)=grd_disto(:,j)   
            !! Generating the filename
                write(scratch,*) j             
                write(grd_out,*) TRIM(ADJUSTL(DIR))//'/'//TRIM(ADJUSTL(scratch))//'.xyz'
            !! Building the distorted geometry.
                call mk_genxyz(grd_coord& 
                 &,grd_label&                
                 &,grd_tdnc&                 
                 &,grd_mode&             
                 &,grd_indivdisto&            
                 &,grd_dim&
                 &,grd_n&
                 &,trim(adjustl(grd_out)))
            !! Writing the Coordinate grid.
                write(667,FMT_B) j, (grd_indivdisto(k),k=1,grd_dim)
               end do
                close(667)



end subroutine



