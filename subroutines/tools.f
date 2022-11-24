!!!! This file contains all the general tools needed for VCMAKER:
      ! upper subroutine: Convert strings in all capital letters
      ! diag subroutine: general subroutine to diagonalise real
      ! symmetric matrices
      ! centerofmass: centrer a set of coordinates on its center of mass
      ! align: compute the angle/axe and generate the rotation matrice
      ! to align a vector on another one
      ! scalar function: Returns the scalar vector
      ! vecto subroutine: Returns the vectorial product of two vectors
      ! normvec function: returns the norm of a vector
      ! Rotation subroutine: generate the rotation matrix
!!! Maybe you should bloody order it alphabetically. FFS


      SUBROUTINE upper(string)
        implicit none
        integer::i,ic
        character(*)::string
    
        character(len=26), Parameter :: &
        & cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        character(len=26), Parameter :: &
        & low = 'abcdefghijklmnopqrstuvwxyz'

        do i = 1, LEN_TRIM(string)
            ic = INDEX(low, string(i:i))
            if (ic > 0) string(i:i) = cap(ic:ic)
        end do


      END SUBROUTINE


      SUBROUTINE diag(init_mat,d_n,eigenvec,eigenval)
      implicit none
      integer,intent(IN)::d_n
      double precision,dimension(3*d_n,3*d_n),intent(in)::init_mat
      double precision,dimension(3*d_n,3*d_n),intent(out)::eigenvec
      double precision,dimension(3*d_n),intent(out)::eigenval
      ! Diagonalization routine
        integer::diag_lwork
        integer::diag_info
        double precision,dimension(:),allocatable::diag_work
      
      eigenvec=init_mat
      
      diag_lwork=-1
      allocate(diag_work(9*d_n))
       
      call DSYEV('V','U',3*d_n,eigenvec,3*d_n,eigenval,&
      &diag_work,diag_lwork,diag_info)
      
      diag_lwork=INT(diag_work(1))
      deallocate(diag_work)
      allocate(diag_work(diag_lwork))
      
      call DSYEV('V','U',3*d_n,eigenvec,3*d_n,eigenval,&
      &diag_work,diag_lwork,diag_info)
      deallocate(diag_work)
      
      END SUBROUTINE

      SUBROUTINE centermass(com_coords,re_com)
        use mod_input
        use mod_constants
        implicit none
        integer::i,j
        double precision,dimension(3,n_at)::com_coords
        double precision,dimension(3)::re_com
        double precision::mtot

        re_com=0.0d0
        mtot=0.0d0
        do i=1,n_at
          do j=1,3
            com(j)=com(j)+mass(i)*com_coords(j,i)
          end do
          mtot=mtot+mass(i)
        end do
        com=com/mtot

        do i=1,n_at
          do j=1,3
            com_coords(j,i)=com_coords(j,i)-com(j)
          end do
        end do

       END SUBROUTINE


       function normvec(vec)
               implicit none
               double precision::normvec
               double precision,dimension(3)::vec
       
       normvec=sqrt(vec(1)**2+vec(2)**2+vec(3)**2)
       return
       end function

       subroutine vecto(veca,vecb,vectorial)
               implicit none
               double precision,dimension(3)::veca,vecb,vectorial
       vectorial(1)=(veca(2)*vecb(3)-veca(3)*vecb(2))
       vectorial(2)=(veca(3)*vecb(1)-veca(1)*vecb(3))
       vectorial(3)=(veca(1)*vecb(2)-veca(2)*vecb(1))
       return
       end subroutine

       function scalar(veca,vecb)
               implicit none
               double precision::scalar
               double precision,dimension(3)::veca,vecb
       
               scalar=veca(1)*vecb(1)+veca(2)*vecb(2)+veca(3)*vecb(3)
               return
        end function

        subroutine align(veca,vecb,rot,angle,axe,orient)
        implicit none
        integer::j
        integer::orient
        double precision,dimension(3)::veca,vecb,axe
        double precision::scalar
        double precision::angle
        double precision::normvec
        double precision::tmp
        double precision,dimension(3,3)::rot
        double precision,parameter::pi=4.0d0*atan(1.0d0)
        double precision,parameter::angle_thresh=1E-5

        tmp=scalar(veca,vecb)/(normvec(veca)*normvec(vecb))
        if (tmp.lt.-1.AND.tmp.gt.-1-angle_thresh) THEN
                tmp=-1
        elseif (tmp.gt.1.AND.tmp.lt.1+angle_thresh) THEN
                tmp=1
        end if
        angle=acos(tmp)
        if (angle.eq.0.OR.angle.eq.pi/2.0d0) then
        axe=0.0d0
        axe(orient)=1.0d0
        else
        call vecto(veca,vecb,axe)
        end if

        if ((orient.eq.1).AND.(axe(1).lt.0.d0)) then
                do j=1,3
                axe(j)=-1*axe(j)
                angle=-angle
                end do
        else if ((orient.eq.2).AND.(axe(2).lt.0.d0)) then
                do j=1,3
                axe(j)=-1*axe(j)
                angle=-angle
                end do
        else if ((orient.eq.3).AND.(axe(3).lt.0.d0)) then
                do j=1,3
                axe(j)=-1*axe(j)
                angle=-angle
                end do
        end if
        call rotation(rot,axe,angle)
        end subroutine

        subroutine rotation(rot_mat,axe_tau,tau)
                implicit none
                integer::i
                double precision,dimension(3,3)::rot_mat
                double precision,dimension(3)::axe_tau
                double precision::tau
                double precision::normvec,norm
                norm=normvec(axe_tau)
                do i=1,3
                  axe_tau(i)=axe_tau(i)/norm
                end do
        
        rot_mat(1,1)= cos(tau) + axe_tau(1)**2*(1-cos(tau)) !!
        rot_mat(1,2)= axe_tau(1) * axe_tau(2) * (1-cos(tau)) - axe_tau(3)*sin(tau)
        rot_mat(1,3)= axe_tau(1) * axe_tau(3) * (1-cos(tau)) + axe_tau(2)*sin(tau)
        
        rot_mat(2,1)= axe_tau(2) * axe_tau(1) * (1-cos(tau)) + axe_tau(3)*sin(tau)
        rot_mat(2,2)= cos(tau) + axe_tau(2)**2*(1-cos(tau))
        rot_mat(2,3)= axe_tau(2) * axe_tau(3) * (1-cos(tau)) - axe_tau(1)*sin(tau)
        
        rot_mat(3,1)= axe_tau(3) * axe_tau(1) * (1-cos(tau)) - axe_tau(2)*sin(tau)
        rot_mat(3,2)= axe_tau(3) * axe_tau(2) * (1-cos(tau)) + axe_tau(1)*sin(tau)
        rot_mat(3,3)= cos(tau) + axe_tau(3)**2*(1-cos(tau)) !!
        
        end subroutine

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine xyz_nto3n(ca,cb,n,job)       !!
        ! Convert xyz(3,n) -> xyx(3*n) if job=1 !!
        ! Convert xyz(3*n) -> xyx(3,n) if job=-1!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        implicit none
        integer::n
        integer::i,j
        integer::job
        double precision,dimension(3,n)::ca
        double precision,dimension(3*n)::cb

        if (job.eq.1) then ! n to 3n
                do i=1,n
                        do j=1,3
                                cb((i-1)*3+j)=ca(j,i)
                        end do
                end do
        else if (job.eq.-1) then
                do i=1,n
                        do j=1,3
                                ca(j,i)=cb((i-1)*3+j)
                        end do
                end do
        end if

        END SUBROUTINE


