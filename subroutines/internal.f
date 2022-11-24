        SUBROUTINE get_defint(gi_n,gi_xyz,gi_def,gi_labels,gi_nbl,gi_nba,gi_ndi)
        !input
                implicit none
                integer,INTENT(IN)::gi_nbl
                integer,INTENT(IN)::gi_nba
                integer,INTENT(IN)::gi_ndi
                integer,INTENT(IN)::gi_n
                double precision,dimension(3,gi_n),INTENT(IN)::gi_xyz
                integer,dimension(4,gi_nbl+gi_nba+gi_ndi),INTENT(OUT)::gi_def
                character(len=16),dimension(gi_nbl+gi_nba+gi_ndi),INTENT(OUT)::gi_labels
        !Looping Integer
                integer::i,j,k
        !Internal variables
                double precision,dimension(3)::vec_a
                double precision,dimension(3)::vec_b
                double precision,dimension(3)::vec_c
                double precision,dimension(3)::n_a
                double precision,dimension(3)::n_b
                integer::c_bl ! Just Itsy bitsy counters  
                integer::c_ba ! Just Itsy bitsy counters
                integer::c_di ! Just Itsy bitsy counters
                integer,dimension(gi_n)::indices
                integer::tampon
                double precision::normvec
        !Just the number pi. Don't mind it!
                double precision,parameter::pi=4.0d0*atan(1.0d0)
        
        
                ! Initialisation: 
                do i=1,gi_n
                   indices(i)=i
                end do
                c_bl=1
                c_ba=1
                c_di=1
                gi_def=0
        
                !! Start from i=1 or i=n_at ?
          if (gi_ndi.gt.0) then
                do i=1,gi_n-3
                   vec_a(:)=gi_xyz(:,indices(i+1))-gi_xyz(:,indices(i))
                   vec_b(:)=gi_xyz(:,indices(i+2))-gi_xyz(:,indices(i+1))
        
                   !Check if vec_a and vec_b are colinear
                   !If yes, switch atom i and i+4
                   call vecto(vec_a,vec_b,n_a)
                   if (normvec(n_a).lt.1E-5) then
                        tampon=indices(i)
                        indices(i)=indices(i+4)
                        indices(i+4)=tampon
                   end if
        
                   vec_a(:)=gi_xyz(:,indices(i+1))-gi_xyz(:,indices(i))
                   vec_b(:)=gi_xyz(:,indices(i+2))-gi_xyz(:,indices(i+1))
                   vec_c(:)=gi_xyz(:,indices(i+3))-gi_xyz(:,indices(i+2))
        
                   !Check if vec_b and vec_c are colinear
                   !If yes, switch atom i+3 and i+4
                   call vecto(vec_b,vec_c,n_b)
                   if (normvec(n_b).lt.1E-5) then
                        tampon=indices(i+3)
                        indices(i+3)=indices(i+4)
                        indices(i+4)=tampon
                   end if
        
                   do j=1,4
                     if (j.le.2) then
                        gi_def(j,c_bl)=indices(i+(j-1))
                     end if
                     if (j.le.3) then
                        gi_def(j,gi_nbl+c_ba)=indices(i+(j-1))
                     end if
                     if (j.le.4) then
                        gi_def(j,gi_nbl+gi_nba+c_di)=indices(i+(j-1))
                     end if
                   end do
        
                        write(gi_labels(c_bl),'(A1,I3,A1,I3,8X)')            'B',indices(i+(1-1)),&
                                                                              &'-',indices(i+(2-1))
                        write(gi_labels(gi_nbl+c_ba),'(3(A1,I3),4X)')        'A',indices(i+(1-1)),&
                                                                              &'-',indices(i+(2-1)),&
                                                                              &'-',indices(i+(3-1))
                        write(gi_labels(gi_nbl+gi_nba+c_di),'(4(A1,I3))') 'D',indices(i+(1-1)),&
                                                                              &'-',indices(i+(2-1)),&
                                                                              &'-',indices(i+(3-1)),&
                                                                              &'-',indices(i+(4-1))
        
                        !Initialisation of the next step:
                        if (i.lt.gi_n-3) then
                           c_bl=c_bl+1
                           c_ba=c_ba+1
                           c_di=c_di+1
                        !If it is the last step, we need to add the two remaining Bond Length
                        else if (i.eq.gi_n-3) then
                           do k=1,2  ! Loop over the two remaining angles
                              do j=1,2  ! Loop over the two atoms involved in the bond
                                 gi_def(j,c_bl+k)=indices(i+k+(j-1))
                              end do
                              write(gi_labels(c_bl+k),'(A1,I3,A1,I3,8X)')            'B',indices(i+k+(1-1)),&
                                                                                    &'-',indices(i+k+(2-1))
                            end do
                           ! And the remaining valence angle.
                           do j=1,3 ! Loop over the three atoms involved in the angle
                              gi_def(j,gi_nbl+c_ba+1)=indices(i+1+(j-1))
                           end do
                           write(gi_labels(gi_nbl+c_ba+1),'(3(A1,I3),4X)')        'A',indices(i+1+(1-1)),&
                                                                              &'-',indices(i+1+(2-1)),&
                                                                              &'-',indices(i+1+(3-1))
                        end if
        
        
        
                 end do
         else if (gi_ndi.eq.0) then
                                 gi_def(1,1)=indices(1)
                                 gi_def(2,1)=indices(2)
                                 gi_def(1,2)=indices(2)
                                 gi_def(2,2)=indices(3)
                              write(gi_labels(1),'(A1,I3,A1,I3,8X)')            'B',1,'-',2
                              write(gi_labels(2),'(A1,I3,A1,I3,8X)')            'B',2,'-',3
                                 gi_def(1,3)=indices(1)
                                 gi_def(2,3)=indices(2)
                                 gi_def(3,3)=indices(3)
                           write(gi_labels(3),'(3(A1,I3),4X)')        'A',1,'-',2,'-',3
         end if
        
        END SUBROUTINE
        
        
        SUBROUTINE get_internals(gi_n,gi_xyz,gi_def,gi_int,gi_nbl,gi_nba,gi_ndi)
        IMPLICIT NONE
        !Input
        integer,INTENT(IN)::gi_n
        integer,INTENT(IN)::gi_nbl
        integer,INTENT(IN)::gi_nba
        integer,INTENT(IN)::gi_ndi
        double precision,dimension(3,gi_n),INTENT(IN)::gi_xyz
        integer,dimension(4,gi_nbl+gi_nba+gi_ndi),INTENT(IN)::gi_def
        !output
        double precision,dimension(gi_nbl+gi_nba+gi_ndi),INTENT(OUT)::gi_int
        !Local variables
        integer::i
        double precision,dimension(3)::vec_a
        double precision::alpha
        double precision::normvec
        double precision,parameter::pi=4.0d0*atan(1.0d0)
        
        do i=1,gi_nbl  !! Bond lengths
        vec_a(:)=gi_xyz(:,gi_def(2,i))-gi_xyz(:,gi_def(1,i))
        gi_int(i)=normvec(vec_a)
        end do
        
        do i=1,gi_nba
        call valence_angle(gi_xyz(:,gi_def(1,i+gi_nbl)),&
                          &gi_xyz(:,gi_def(2,i+gi_nbl)),&
                          &gi_xyz(:,gi_def(3,i+gi_nbl)),&
                          &alpha)
        gi_int(i+gi_nbl)=alpha
        end do
        
        do i=1,gi_ndi
        call sub_died(gi_xyz(:,gi_def(1,i+gi_nbl+gi_nba)),&
                     &gi_xyz(:,gi_def(2,i+gi_nbl+gi_nba)),&
                     &gi_xyz(:,gi_def(3,i+gi_nbl+gi_nba)),&
                     &gi_xyz(:,gi_def(4,i+gi_nbl+gi_nba)),&
                     &alpha)
        gi_int(i+gi_nbl+gi_nba)=alpha
        end do
        END SUBROUTINE


               SUBROUTINE sub_died(ca,cb,cc,cd,died)  !! As in GeomClass
               implicit none
               DOUBLE PRECISION::cos_died
               DOUBLE PRECISION,DIMENSION(3),INTENT(IN)::ca,cb,cc,cd
               DOUBLE PRECISION,DIMENSION(3)::va,vb,vc
               DOUBLE PRECISION,DIMENSION(3)::na,nb,nc
               DOUBLE PRECISION::normvec,phi,tau,scalar
               DOUBLE PRECISION::check
               DOUBLE PRECISION,INTENT(OUT)::died
               DOUBLE PRECISION,PARAMETER::pi=4.0d0*atan(1.0d0)
               INTEGER::signe
               va(:)=cb(:)-ca(:)
               vb(:)=cc(:)-cb(:)
               vc(:)=cd(:)-cc(:)
               va=va/normvec(va)
               vb=vb/normvec(vb)
               vc=vc/normvec(vc)

               phi=acos(scalar(-va,vb))
               tau=acos(scalar(-vb,vc))

               call vecto(va,vb,na)
               call vecto(vb,vc,nb)
               cos_died=scalar(na,nb)/ (sin(phi)*sin(tau))
               call vecto(na,nb,nc)

               ! Checking the sign of the diedral
!               check=scalar( nc,va)
               check=scalar( nc,vb)
               if ( check .gt. 0.0d0 ) then
                        signe=1                       
               else
                        signe=-1                       
               end if


               ! Taking care of numerical problems
               if (cos_died .gt. 1) then
                       died=0.0d0
               else if (cos_died.lt.-1) then
                       died=DBLE(signe)*pi
               else
               died=DBLE(signe)*acos(cos_died)
               end if


             END SUBROUTINE

             SUBROUTINE valence_angle(ca,cb,cc,angle) ! Angle from Cart coords
              IMPLICIT NONE
              DOUBLE PRECISION,DIMENSION(3),INTENT(IN)::ca,cb,cc
              DOUBLE PRECISION,DIMENSION(3)::va,vb
              DOUBLE PRECISION,INTENT(OUT)::angle
              DOUBLE PRECISION::normvec,scalar

        va(:)=(ca(:)-cb(:))
        vb(:)=(cc(:)-cb(:))
        angle=acos(scalar(va,vb)/(normvec(va)*normvec(vb)))

       END SUBROUTINE


        !!!!!!!! SUBROUTINES FOR THE !!!!!!!!
        !!!!! GENERATION OF THE WILSON !!!!!!
        !!!!!!!! B MATRIX ELEMENTS !!!!!!!!!!
        SUBROUTINE mk_bmat(mat,n,n_bl,n_ba,n_di,&
                          &coords,int_def) !!! Check what to put here.
                IMPLICIT NONE
                integer,INTENT(IN)::n_bl,n_ba,n_di,n
                integer,dimension(4,n_bl+n_ba+n_di),INTENT(IN)::int_def
                double precision,dimension(3,n),INTENT(IN)::coords

                double precision,dimension(n_bl+n_ba+n_di,3*n),INTENT(OUT)::mat

                integer::i,j !! Fruity Loops
                double precision,dimension(:),allocatable::tmp_bmat

          do i=1,n_bl              !! Loop over bonds
          allocate(tmp_bmat(2*3))  !! XYZ for each atom of the bond
          tmp_bmat=0.0d0
          call bmat_bond(tmp_bmat,coords(:,int_def(1,i))&
                                &,coords(:,int_def(2,i)))            !returns a tmp matrice, allocate to the entire Bmat here

          do j=1,3
          mat(i,(int_def(1,i)-1)*3+j)=tmp_bmat(j)
          mat(i,(int_def(2,i)-1)*3+j)=tmp_bmat(j+3)
          end do

          deallocate(tmp_bmat)
          end do

          do i=1,n_ba !! Loop over angles
          allocate(tmp_bmat(3*3))
          tmp_bmat=0.0d0
          call bmat_angle(tmp_bmat,coords(:,int_def(1,i+n_bl))&
                                 &,coords(:,int_def(2,i+n_bl))&
                                 &,coords(:,int_def(3,i+n_bl)))
                         !returns a tmp matrice, allocate to the entire Bmat here
          do j=1,3
          mat(i+n_bl,(int_def(1,i+n_bl)-1)*3+j)=tmp_bmat(j)
          mat(i+n_bl,(int_def(2,i+n_bl)-1)*3+j)=tmp_bmat(j+3)
          mat(i+n_bl,(int_def(3,i+n_bl)-1)*3+j)=tmp_bmat(j+6)
          end do
          deallocate(tmp_bmat)
          end do

          do i=1,n_di !! Loop over dihedral
          allocate(tmp_bmat(4*3))
          !call bmat_dihedral_alt(tmp_bmat&
          call bmat_dihedral(tmp_bmat&
                            &,coords(:,int_def(1,i+n_bl+n_ba))&
                            &,coords(:,int_def(2,i+n_bl+n_ba))&
                            &,coords(:,int_def(3,i+n_bl+n_ba))&
                            &,coords(:,int_def(4,i+n_bl+n_ba)))
          do j=1,3
          mat(i+n_bl+n_ba,(int_def(1,i+n_bl+n_ba)-1)*3+j)=tmp_bmat(j)
          mat(i+n_bl+n_ba,(int_def(2,i+n_bl+n_ba)-1)*3+j)=tmp_bmat(j+3)
          mat(i+n_bl+n_ba,(int_def(3,i+n_bl+n_ba)-1)*3+j)=tmp_bmat(j+6)
          mat(i+n_bl+n_ba,(int_def(4,i+n_bl+n_ba)-1)*3+j)=tmp_bmat(j+9)
          end do
          deallocate(tmp_bmat)
          end do

        END SUBROUTINE

        SUBROUTINE mk_Amat(amat,bmat,n,labels,ni)
        IMPLICIT NONE
        INTEGER::i,j
        INTEGER,intent(in)::n
        INTEGER,intent(in)::ni             ! n_int
        double precision,dimension(:,:),allocatable::m
        double precision,dimension(3*n,ni),intent(out)::amat
        double precision,dimension(ni,3*n),intent(in)::bmat
        double precision,dimension(:,:),allocatable::bmbt
        character(len=2),dimension(n),intent(in)::labels
        double precision,dimension(:),allocatable::mass
        !! DGETRF
        integer,dimension(:),allocatable::dgetrf_ipiv
        integer::dgetrf_info
        integer::dgetri_lwork
        double precision,dimension(:),allocatable::dgetri_work

        ! We compute A such as A= m.B^T.(B.m.B^T)^(-1)
        ! where m is any non singular 3Nx3N matrice --> We can take mii=1/mass
        allocate(m(3*n,3*n))
        allocate(mass(n))
        call get_mass(labels,mass,n)
        m=0.0d0
        do i=1,n
          do j=1,3
        m((i-1)*3+j,(i-1)*3+j)=1.0d0/mass(i)
          end do
        end do

        !We comput the second factor: B.m.B^T
        allocate(bmbt(ni,ni))
        amat=0.0d0
        bmbt=0.0d0
        bmbt=MATMUL(bmat,MATMUL(m,TRANSPOSE(bmat)))
        !Inversion.
        allocate(dgetrf_ipiv(ni))
        call dgetrf(ni,ni,bmbt,ni,dgetrf_ipiv,dgetrf_info)
        if (dgetrf_info.ne.0) then
           write(*,*) 'Error in DGETRF routine (inversion of B.M.Bt)'
           write(*,*) 'Info= ',dgetrf_info
           stop
        end if
        dgetri_lwork=-1
        allocate(dgetri_work(1))
        call dgetri(ni,bmbt,ni,dgetrf_ipiv,dgetri_work&
                   &,dgetri_lwork,dgetrf_info)
        dgetri_lwork=INT(dgetri_work(1))
        deallocate(dgetri_work)
        allocate(dgetri_work(dgetri_lwork))
        call dgetri(ni,bmbt,ni,dgetrf_ipiv,dgetri_work&
                   &,dgetri_lwork,dgetrf_info)
        deallocate(dgetri_work)
        if (dgetrf_info.ne.0) then
            write(*,*) 'Error in DGETRI routine (inversion of B.M.Bt)'
            write(*,*) 'Info= ',dgetrf_info
            stop
        end if
        amat=MATMUL( MATMUL(m,TRANSPOSE(bmat)),bmbt)
        END SUBROUTINE


        SUBROUTINE bmat_bond(bmat,ca,cb) !! ca=coords atom a !! cb=coords atom b
        IMPLICIT NONE
        double precision,dimension(3*2)::bmat
        double precision,dimension(3)::ca,cb
        double precision,dimension(3)::vec
        double precision::normvec
        integer::j

        ! A-----B

        vec(:)=cb(:)-ca(:)

        do j=1,3
        ! Atom A
          bmat(j)=-vec(j)/normvec(vec)
        ! Atom B
        bmat(j+3)=vec(j)/normvec(vec)
        end do
        END SUBROUTINE



        SUBROUTINE bmat_angle(bmat,ca,cb,cc) !! ca=coords atom a !! cb=coords atom b !! cc=coords atomc !! B is the APEX atom.
        IMPLICIT NONE
        double precision,dimension(3*3)::bmat
        double precision,dimension(3)::ca,cb,cc
        double precision,dimension(3)::veca,vecb
        double precision,dimension(3)::ea,eb ! Unitary vectors
        double precision::normvec,scalar,alpha,sina,cosa
        integer::j

        !   A-----B
        !       |  \
        !       a-- \
        !            C

        !Definition of the vectors and the associated
        !Unitary displacements
        veca(:)=ca(:)-cb(:)
        vecb(:)=cc(:)-cb(:)
        ea(:)=veca(:)/normvec(veca)
        eb(:)=vecb(:)/normvec(vecb)

        !Description of alpha=angle(veca,vecb)
        cosa=scalar(veca,vecb)/(normvec(veca)*normvec(vecb))
        sina=dsqrt(1-cosa**2)
        alpha=acos(cosa)

        do j=1,3
        ! Atom A
        bmat(j)=(cosa*ea(j)-eb(j))/(normvec(veca)*sina)
        ! Atom C
        bmat(j+6)=(cosa*eb(j)-ea(j))/(normvec(vecb)*sina)
        ! Atom B (APEX)
        bmat(j+3)=-(bmat(j)+bmat(j+6))
        end do
        END SUBROUTINE


        SUBROUTINE bmat_dihedral(bmat,ca,cb,cc,cd) !! ca=coords atom a !! cb=coords atom b !! cc=coords atomc !! B is the APEX atom.
        IMPLICIT NONE
        double precision,dimension(3*4),intent(OUT)::bmat
        double precision,dimension(3),intent(IN)::ca,cb,cc,cd
        double precision,dimension(3)::veca,vecb,vecc
        double precision,dimension(3)::e12,e23,e34 ! Unitary vectors
        double precision,dimension(3)::h2,h3
        double precision,dimension(3)::f2,f3
        double precision::normvec,scalar
        double precision::alpha,sina,cosa
        double precision::beta,sinb,cosb
        integer::j

        !   A-----B
        !       |  \---b
        !       a-- \  \
        !            C-----D

        !Definition of the vectors and the associated
        !Unitary displacements
        veca(:)=cb(:)-ca(:)
        vecb(:)=cc(:)-cb(:)
        vecc(:)=cd(:)-cc(:)
        e12(:)=veca(:)/normvec(veca)
        e23(:)=vecb(:)/normvec(vecb)
        e34(:)=vecc(:)/normvec(vecc)

        !Description of alpha=angle(-veca,vecb)
        cosa=scalar(-veca,vecb)/(normvec(veca)*normvec(vecb))
        alpha=acos(cosa)
        sina=sin(alpha)
        !sina=dsqrt(1-cosa**2)
        !Description of beta=angle(-vecb,vecc)
        cosb=scalar(-vecb,vecc)/(normvec(vecb)*normvec(vecc))
        beta=acos(cosb)
        sinb=sin(beta)
        !sinb=dsqrt(1-cosb**2)

        call vecto(e12,e23,h2)
        call vecto(e23,e34,h3)
        ! we use:
        ! (-u)^v=u^(-v)=-(u^v)
        ! (-u)^(-v)=u^v

        do j=1,3
        h2(j)=h2(j)/sina
        h3(j)=h3(j)/sinb
        end do

        call vecto(h2,e23,f2)
        call vecto(h3,e23,f3)

        do j=1,3
        bmat(0+j)=-h2(j)/(normvec(veca)*sina) ! St1
        bmat(9+j)=+h3(j)/(normvec(vecc)*sinb) ! St4
        end do

        do j=1,3
        bmat(3+j)=((normvec(veca)*cosa-normvec(vecb))/normvec(vecb))*bmat(0+j)&
                 &-((normvec(vecc)*cosb)/(normvec(vecb)))*bmat(9+j)
        bmat(6+j)=((normvec(vecc)*cosb-normvec(vecb))/normvec(vecb))*bmat(9+j)&
                 &-((normvec(veca)*cosa)/(normvec(vecb)))*bmat(0+j)
        end do

        END SUBROUTINE

                SUBROUTINE get_rmsd(ca,cb,n,rmsd)
                IMPLICIT NONE
                double precision,dimension(3,n),intent(IN)::ca,cb
                double precision,intent(OUT)::rmsd
                integer,intent(IN)::n
                integer::i,j

                rmsd=0.0d0
                do i=1,n
                    do j=1,3
                rmsd=rmsd+(ca(j,i)-cb(j,i))**2
                    end do
                end do
                rmsd=dsqrt((1.0d0/DBLE(n))*rmsd)
        END SUBROUTINE

