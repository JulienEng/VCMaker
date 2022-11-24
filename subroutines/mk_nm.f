subroutine mk_nm(nm_evec,nm_nmod,nm_redmass,nm_mass,nm_n)
implicit none
integer::i,j

integer,intent(IN)::nm_n
double precision,dimension(nm_n),intent(IN)::nm_mass
double precision,dimension(3*nm_n,3*nm_n),intent(INOUT)::nm_evec
double precision,dimension(3*nm_n,3*nm_n),intent(INOUT)::nm_nmod
double precision,dimension(3*nm_n),intent(INOUT)::nm_redmass

double precision,dimension(3*nm_n)::nm_mass3n

do i=1,nm_n
   do j=1,3
nm_mass3n((i-1)*3+j)=nm_mass(i)
end do
end do

do i=1,3*nm_n
    do j=1,3*nm_n
    nm_nmod(i,j)=nm_evec(i,j)/dsqrt(nm_mass( FLOOR((DBLE(j)+2)/3.0d0)   ) )
    end do
end do

nm_redmass=0.d00
do i=1,3*nm_n
    do j=1,3*nm_n
nm_redmass(i)=nm_redmass(i)+nm_nmod(i,j)**2
    end do
    do j=1,3*nm_n
nm_nmod(i,j)=nm_nmod(i,j)/dsqrt(nm_redmass(i))
    end do
end do
end subroutine

