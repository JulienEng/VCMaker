      subroutine get_grad(grad_file,grad_grad,grad_n)
      use mod_err
      IMPLICIT NONE
      integer::l
      integer,intent(IN)::grad_n
      integer::iost
      double precision,dimension(3*grad_n),intent(INOUT)::grad_grad
      character(len=*),intent(IN)::grad_file
      
      !! Opening the VCM formatted file containing the Cartesian Gradient
      OPEN(78,file=trim(adjustl(grad_file)))
      !! Loop over the number of atoms*3
      do l=1,3*grad_n
            read(78,*,iostat=iost) grad_grad(l)
            !! If fewer lines than anticipated
            if (iost.ne.0) then
                  !! Write the error
                  write(err_msg,*) 'Error while reading Gradient file: ',TRIM(ADJUSTL(grad_file))
                  call mk_error()
            end if
      end do
      !! We close everything before leaving.
      CLOSE(78)
end subroutine
