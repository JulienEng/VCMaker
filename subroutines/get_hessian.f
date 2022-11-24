      subroutine get_hessian(filen,hessmat)
      use mod_input
      use mod_constants
      use mod_err
      implicit none
      integer::i,j
      character(*)::filen
      integer::iost
      double precision,dimension(3*n_at,3*n_at)::hessmat

      !! Opening the VCM formatted file containing the Cartesian Hessian
      open(78,file=filen)
      !! Loop over the number of atoms*3
      do j=1,3*n_at
        !! Loop over the number of atoms*3
        read(78,*,iostat=iost)  (hessmat(j,i),i=1,3*n_at)
        !! If fewer lines than anticipated or not enough terms/line
        if (iost.ne.0) then
          !! Write the error
          write(err_msg,*) 'Error while reading Hessian file: ',TRIM(ADJUSTL(filen))
          call mk_error()
        end if
      end do
      !! We close everything before leaving.     
      close(78)
end subroutine
