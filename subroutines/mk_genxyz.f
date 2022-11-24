        SUBROUTINE mk_genxyz(gen_coord&
              &,gen_label&
              &,gen_tdnc&
              &,gen_mode&
              &,gen_disto&
              &,n_gen&
              &,n_at&
              &,gen_out)
        use mod_constants
        implicit none
        integer,intent(in)::n_at
        integer,intent(in)::n_gen
        integer::i,j,k
        character(len=2),dimension(n_at),intent(in)::gen_label
        double precision,dimension(3,n_at),intent(in)::gen_coord
        double precision,dimension(3*n_at,3*n_at),intent(in)::gen_tdnc
        integer,dimension(n_gen),intent(in)::gen_mode
        double precision,dimension(n_gen),intent(in)::gen_disto
        double precision,dimension(3,n_at)::gen_coordin
        double precision,dimension(3,n_at)::gen_coordout
        character(len=*)::gen_out

        gen_coordin(:,:)=gen_coord(:,:)
        do i=1,n_gen
            do k=1,n_at
               do j=1,3 !x,y,z
                  gen_coordout(j,k)=gen_coordin(j,k)+gen_disto(i)*gen_tdnc(gen_mode(i)+6,((k-1)*3+j))/angtoa0 
               end do
            end do
            gen_coordin(:,:)=gen_coordout(:,:)
        end do

        call out_xyz(gen_coordout,gen_label,n_at,gen_out)
        
        
        END SUBROUTINE
