subroutine mk_gap(gap_kgrad,gap_kdisp,gap_gap,gap_rsma,gap_n)
    implicit none
    integer::i
    integer,intent(in)::gap_n
    double precision,dimension(3*gap_n),intent(in)::gap_kgrad
    double precision,dimension(3*gap_n),intent(in)::gap_kdisp

    double precision::gap_K
    double precision::gap_Ks
    double precision::gap_normRsmf
    double precision,intent(out)::gap_gap
    double precision,dimension(3*gap_n),intent(out)::gap_rsma
    

    gap_K=0.0d0
    gap_Ks=0.0d0
    gap_normRsmf=0.0d0
do i=7,3*gap_n  ! First 6 modes have NAN gradients
     gap_K=gap_K + ABS(gap_kgrad(i))
    gap_ks=gap_ks+ ABS(gap_kdisp(i))
    gap_rsma(i)=(gap_kgrad(i)-gap_kdisp(i))*1000
    gap_normRsmf=gap_normRsmf+gap_rsma(i)**2
end do

gap_rsma=gap_rsma/dsqrt(gap_normRsmf)
gap_gap=(MAX(gap_K,gap_Ks)-MIN(gap_K,gap_Ks))/MAX(gap_K,gap_Ks)


end subroutine
