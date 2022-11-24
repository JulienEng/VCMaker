subroutine mk_hess(hess_hess,hess_mass,hess_mwhess,hess_coord,hess_com,hess_freq,hess_evec,hess_n)
implicit none
integer::i,j
integer,intent(in)::hess_n
double precision,dimension(hess_n),intent(in)::hess_mass
double precision,dimension(3,hess_n),intent(in)::hess_coord
double precision,dimension(3),intent(in)::hess_com

double precision,dimension(3*hess_n,3*hess_n),intent(inout)::hess_hess
double precision,dimension(3*hess_n),intent(inout)::hess_freq
double precision,dimension(3*hess_n,3*hess_n),intent(inout)::hess_mwhess
double precision,dimension(3*hess_n,3*hess_n),intent(inout)::hess_evec
double precision,dimension(:,:),allocatable::mass_mat
allocate(mass_mat(3*hess_n,3*hess_n))

do i=1,3*hess_n ! Make the hessian totally sym
    do j=1,i
        hess_hess(i,j) = (hess_hess(i,j)+hess_hess(j,i))/2.0d0
        hess_hess(j,i) =  hess_hess(i,j)
    end do
end do

!generate the mass matrix
mass_mat=0.0d0
do i=1,hess_n
   do j=1,3
mass_mat((i-1)*3+j,(i-1)*3+j)=1/dsqrt(hess_mass(i))
   end do
end do

hess_mwhess=MATMUL(MATMUL(mass_mat,hess_hess),mass_mat)
call transfo(hess_coord,hess_mass,hess_com,hess_mwhess,hess_n)
call diag(hess_mwhess,hess_n,hess_evec,hess_freq)

hess_evec=TRANSPOSE(hess_evec)

do i=1,6
hess_freq(i)=0.0d0
end do

deallocate(mass_mat)
end subroutine

SUBROUTINE transfo(tr_coord,tr_mass,tr_com,tr_mwhess,tr_n)
use mod_input
use mod_constants
implicit none
integer::i,j,m

integer,intent(IN)::tr_n
double precision,dimension(3,tr_n),intent(IN)::tr_coord
double precision,dimension(tr_n),intent(IN)::tr_mass
double precision,dimension(3),intent(IN)::tr_com

double precision,dimension(3*tr_n,3*tr_n),intent(INOUT)::tr_mwhess


double precision,dimension(6,3*n_at)::TROT
double precision,dimension(3*n_at,3*n_at)::proj

double precision::norm

TROT=0.d0
do i=1,3
   do j=1,tr_n
TROT(i,(j-1)*3+i)=dsqrt(tr_mass(j))
   end do
end do
do j=1,tr_n
TROT(4,(j-1)*3+3)= sqrt(tr_mass(j))*(tr_coord(2,j)-tr_com(2))
TROT(4,(j-1)*3+2)=-sqrt(tr_mass(j))*(tr_coord(3,j)-tr_com(3))

TROT(5,(j-1)*3+1)= sqrt(tr_mass(j))*(tr_coord(3,j)-tr_com(3))
TROT(5,(j-1)*3+3)=-sqrt(tr_mass(j))*(tr_coord(1,j)-tr_com(1))

TROT(6,(j-1)*3+2)= sqrt(tr_mass(j))*(tr_coord(1,j)-tr_com(1))
TROT(6,(j-1)*3+1)=-sqrt(tr_mass(j))*(tr_coord(2,j)-tr_com(2))
end do
norm=0.d0

  do i=1,6
norm=0.d0
        do j=1,3*tr_n
norm=norm+TROT(i,j)**2
        end do
        do j=1,3*tr_n
TROT(i,j)=TROT(i,j)/dsqrt(norm)
        end do
norm=0.d0
        do j=1,3*tr_n
norm=norm+TROT(i,j)**2
        end do
 end do

proj=0.d0
do i=1,3*tr_n
proj(i,i)=1.0d0
    do m=1,6
proj(i,i)=proj(i,i)-TROT(m,i)**2
    end do
end do

do i=1,3*tr_n
   do j=i+1,3*tr_n
    do m=1,6
proj(i,j)=proj(i,j)-(TROT(m,i)*TROT(m,j))
    end do
proj(j,i)=proj(i,j)
   end do
end do

tr_mwhess=matmul(proj,MATMUL(tr_mwhess,TRANSPOSE(proj)))

end subroutine

