!! SUBROUTINES
!! ----------------------------------------------------------
!! get_atn(atn_file,atn_n)   get the atom number from a xyz file
!!
!! atn_file : character | xyz file
!! atn_n    : integer   | number of atom 
!! ----------------------------------------------------------
!! get_xyz(xyz_coords,xyz_labels,xyz_n,xyz_file)    get the labels/coordinates from a xyz file
!!
!! xyz_coords(3,xyz_n) : double precision | contains the coordinates
!! xyz_labels(xyz_n)   : character(len=2) | contains the label of each atom
!! xyz_n               : integer          | number of atoms
!! xyz_file            : character        | xyz file
!! ----------------------------------------------------------
!! get_mass(mass_labels,mass_mass,mass_n)
!!
!! mass_labels(mass_n) : character(len=2) | contains the labels of the atoms
!! mass_mass(mass_n)   : double precision | returns the mass of the atoms
!! mass_n              : integer          | number of atoms
!! ----------------------------------------------------------

SUBROUTINE get_atn(atn_file,atn_n)
        use mod_err
        IMPLICIT NONE
        integer,intent(OUT)::atn_n
        integer::iost
        character(len=*),intent(IN)::atn_file

        OPEN(67,file=atn_file)
        read(67,*,iostat=iost) atn_n
              if (iost.ne.0) then
                write(err_msg,*) 'Error while retrieving the number of atom from the XYZ file: ',TRIM(ADJUSTL(atn_file))
                call mk_error()
              end if

        CLOSE(67)
END SUBROUTINE

SUBROUTINE get_xyz(xyz_coords,xyz_labels,xyz_n,xyz_file)
        use mod_err
        IMPLICIT NONE
        integer::k,l  ! Loop Integers
        integer,intent(IN)::xyz_n  ! #Atoms
        integer::iost
        double precision,dimension(3,xyz_n),intent(OUT)::xyz_coords ! coordinates(xyz,atom)
        character(len=2),dimension(xyz_n),intent(OUT)::xyz_labels 
        character(len=*),intent(IN)::xyz_file

        open(67,file=xyz_file)
        read(67,*)
        read(67,*)
        do k=1,xyz_n
        read(67,*,iostat=iost) xyz_labels(k), (xyz_coords(l,k),l=1,3)
                if (iost.ne.0) then
                write(err_msg,*) 'Error while reading XYZ file: ',TRIM(ADJUSTL(xyz_file))
                call mk_error()
                end if
        end do
        
        close(67)
END SUBROUTINE

SUBROUTINE get_com(com_coords,com_mass,com_n,re_com)
        implicit none
!! Inputs
        integer,intent(IN)::com_n
        double precision,dimension(3,com_n),intent(IN)::com_coords
        double precision,dimension(com_n),intent(IN)::com_mass
!! Output
        double precision,dimension(3),intent(INOUT)::re_com
!! Internal variables
        integer::k,l
        double precision::mtot

        re_com=0.0d0
        mtot=0.0d0
        do k=1,com_n
          do l=1,3
            re_com(l)=re_com(l)+com_mass(k)*com_coords(l,k)
          end do
          mtot=mtot+com_mass(k)
        end do
        re_com=re_com/mtot

END SUBROUTINE

SUBROUTINE get_inertia(in_coords,in_mass,in_n,in_com,ivec)
implicit none
integer::l,m,k
integer,intent(IN)::in_n
double precision,dimension(3,in_n),intent(IN)::in_coords
double precision,dimension(3,in_n)::in_coords_tmp
double precision,dimension(in_n),INTENT(IN)::in_mass
double precision,dimension(3),INTENT(IN)::in_com

double precision,dimension(3)::in_vec
double precision,dimension(3)::ival
double precision,dimension(3,3)::in_inertia
double precision,dimension(3,3),INTENT(INOUT)::ivec
double precision::normvec

      do l=1,in_n
         do m=1,3
in_coords_tmp(m,l)=in_coords(m,l)-in_com(m)
         end do
      end do
      in_inertia=0.0d0
      do l=1,3
         do m=1,3
          if (l.eq.m) THEN
              do k=1,in_n
in_vec(:)=in_coords_tmp(:,k)
in_inertia(l,m)=in_inertia(l,m)+in_mass(k)*(normvec(in_vec)**2-in_coords_tmp(l,k)*in_coords_tmp(m,k))
              end do
          ELSE
              do k=1,in_n
in_inertia(l,m)=in_inertia(l,m)-in_mass(k)*in_coords_tmp(l,k)*in_coords_tmp(m,k)
              end do
          END IF
         end do
      end do
      call diag(in_inertia,1,ivec,ival)
END SUBROUTINE

        SUBROUTINE get_mass(mass_labels,mass_mass,mass_n)
                use mod_err
                implicit none
                integer::k,l
                integer::mass_n
                character(len=2),dimension(mass_n),intent(IN)::mass_labels
                double precision,dimension(mass_n),intent(INOUT)::mass_mass

                character(len=2),dimension(:),allocatable::dtbase_labels
                double precision,dimension(:),allocatable::dtbase_mass
                integer::dtbase_len
        
                ! To update if database is updated.
                dtbase_len=118
                allocate(dtbase_labels(dtbase_len))
                allocate(dtbase_mass(dtbase_len))

                mass_mass=0.0d0

                dtbase_labels(1)='H'
                dtbase_labels(2)='He'
                dtbase_labels(3)='Li'
                dtbase_labels(4)='Be'
                dtbase_labels(5)='B'
                dtbase_labels(6)='C'
                dtbase_labels(7)='N'
                dtbase_labels(8)='O'
                dtbase_labels(9)='F'
                dtbase_labels(10)='Ne'
                dtbase_labels(11)='Na'
                dtbase_labels(12)='Mg'
                dtbase_labels(13)='Al'
                dtbase_labels(14)='Si'
                dtbase_labels(15)='P'
                dtbase_labels(16)='S'
                dtbase_labels(17)='Cl'
                dtbase_labels(18)='Ar'
                dtbase_labels(19)='K'
                dtbase_labels(20)='Ca'
                dtbase_labels(21)='Sc'
                dtbase_labels(22)='Ti'
                dtbase_labels(23)='V'
                dtbase_labels(24)='Cr'
                dtbase_labels(25)='Mn'
                dtbase_labels(26)='Fe'
                dtbase_labels(27)='Co'
                dtbase_labels(28)='Ni'
                dtbase_labels(29)='Cu'
                dtbase_labels(30)='Zn'
                dtbase_labels(31)='Ga'
                dtbase_labels(32)='Ge'
                dtbase_labels(33)='As'
                dtbase_labels(34)='Se'
                dtbase_labels(35)='Br'
                dtbase_labels(36)='Kr'
                dtbase_labels(37)='Rb'
                dtbase_labels(38)='Sr'
                dtbase_labels(39)='Y'
                dtbase_labels(40)='Zr'
                dtbase_labels(41)='Nb'
                dtbase_labels(42)='Mo'
                dtbase_labels(43)='Tc'
                dtbase_labels(44)='Ru'
                dtbase_labels(45)='Rh'
                dtbase_labels(46)='Pd'
                dtbase_labels(47)='Ag'
                dtbase_labels(48)='Cd'
                dtbase_labels(49)='In'
                dtbase_labels(50)='Sn'
                dtbase_labels(51)='Sb'
                dtbase_labels(52)='Te'
                dtbase_labels(53)='I'
                dtbase_labels(54)='Xe'
                dtbase_labels(55)='Cs'
                dtbase_labels(56)='Ba'
                dtbase_labels(57)='La'
                dtbase_labels(58)='Ce'
                dtbase_labels(59)='Pr'
                dtbase_labels(60)='Nd'
                dtbase_labels(61)='Pm'
                dtbase_labels(62)='Sm'
                dtbase_labels(63)='Eu'
                dtbase_labels(64)='Gd'
                dtbase_labels(65)='Tb'
                dtbase_labels(66)='Dy'
                dtbase_labels(67)='Ho'
                dtbase_labels(68)='Er'
                dtbase_labels(69)='Tm'
                dtbase_labels(70)='Yb'
                dtbase_labels(71)='Lu'
                dtbase_labels(72)='Hf'
                dtbase_labels(73)='Ta'
                dtbase_labels(74)='W'
                dtbase_labels(75)='Re'
                dtbase_labels(76)='Os'
                dtbase_labels(77)='Ir'
                dtbase_labels(78)='Pt'
                dtbase_labels(79)='Au'
                dtbase_labels(80)='Hg'
                dtbase_labels(81)='Tl'
                dtbase_labels(82)='Pb'
                dtbase_labels(83)='Bi'
                dtbase_labels(84)='Po'
                dtbase_labels(85)='At'
                dtbase_labels(86)='Rn'
                dtbase_labels(87)='Fr'
                dtbase_labels(88)='Ra'
                dtbase_labels(89)='Ac'
                dtbase_labels(90)='Th'
                dtbase_labels(91)='Pa'
                dtbase_labels(92)='U'
                dtbase_labels(93)='Np'
                dtbase_labels(94)='Pu'
                dtbase_labels(95)='Am' 
                dtbase_labels(96)='Cm'
                dtbase_labels(97)='Bk'
                dtbase_labels(98)='Cf'
                dtbase_labels(99)='Es'
                dtbase_labels(100)='Fm'
                dtbase_labels(101)='Md'
                dtbase_labels(102)='No'
                dtbase_labels(103)='Lr'
                dtbase_labels(104)='Rf'
                dtbase_labels(105)='Db'
                dtbase_labels(106)='Sg'
                dtbase_labels(107)='Bh'
                dtbase_labels(108)='Hs'
                dtbase_labels(109)='Mt'
                dtbase_labels(110)='Ds'                                                                                                                      
                dtbase_labels(111)='Rg'
                dtbase_labels(112)='Cn'
                dtbase_labels(113)='Nh'
                dtbase_labels(114)='Fl'
                dtbase_labels(115)='Mc'
                dtbase_labels(116)='Lv'
                dtbase_labels(117)='Ts'
                dtbase_labels(118)='Og'         

                dtbase_mass(1)=1.0079
                dtbase_mass(2)=4.0026
                dtbase_mass(3)=6.941
                dtbase_mass(4)=9.0122
                dtbase_mass(5)=10.811
                dtbase_mass(6)=12.0107
                dtbase_mass(7)=14.0067
                dtbase_mass(8)=15.9994
                dtbase_mass(9)=18.9984
                dtbase_mass(10)=20.1797
                dtbase_mass(11)=22.9897
                dtbase_mass(12)=24.305
                dtbase_mass(13)=26.9815
                dtbase_mass(14)=28.0855
                dtbase_mass(15)=30.9738
                dtbase_mass(16)=32.065
                dtbase_mass(17)=35.453
                dtbase_mass(18)=39.0983
                dtbase_mass(19)=39.948
                dtbase_mass(20)=40.078
                dtbase_mass(21)=44.9559
                dtbase_mass(22)=47.867
                dtbase_mass(23)=50.9415
                dtbase_mass(24)=51.9961
                dtbase_mass(25)=54.938
                dtbase_mass(26)=55.845
                dtbase_mass(27)=58.6934
                dtbase_mass(28)=58.9332
                dtbase_mass(29)=63.546
                dtbase_mass(30)=65.39
                dtbase_mass(31)=69.723
                dtbase_mass(32)=72.64
                dtbase_mass(33)=74.9216
                dtbase_mass(34)=78.96
                dtbase_mass(35)=79.904
                dtbase_mass(36)=83.8
                dtbase_mass(37)=85.4678
                dtbase_mass(38)=87.62
                dtbase_mass(39)=88.9059
                dtbase_mass(40)=91.224
                dtbase_mass(41)=92.9064
                dtbase_mass(42)=95.94
                dtbase_mass(43)=98
                dtbase_mass(44)=101.07
                dtbase_mass(45)=102.9055
                dtbase_mass(46)=106.42
                dtbase_mass(47)=107.8682
                dtbase_mass(48)=112.411
                dtbase_mass(49)=114.818
                dtbase_mass(50)=118.71
                dtbase_mass(51)=121.76
                dtbase_mass(52)=126.9045
                dtbase_mass(53)=127.6
                dtbase_mass(54)=131.293
                dtbase_mass(55)=132.9055
                dtbase_mass(56)=137.327
                dtbase_mass(57)=138.9055
                dtbase_mass(58)=140.12
                dtbase_mass(59)=140.91
                dtbase_mass(60)=144.24
                dtbase_mass(61)=145.00
                dtbase_mass(62)=150.36
                dtbase_mass(63)=151.96
                dtbase_mass(64)=157.25
                dtbase_mass(65)=158.93
                dtbase_mass(66)=162.50
                dtbase_mass(67)=164.93
                dtbase_mass(68)=167.26
                dtbase_mass(69)=168.93
                dtbase_mass(70)=173.05
                dtbase_mass(71)=174.97
                dtbase_mass(72)=178.49
                dtbase_mass(73)=180.95
                dtbase_mass(74)=183.84
                dtbase_mass(75)=186.027
                dtbase_mass(76)=190.23
                dtbase_mass(77)=192.22
                dtbase_mass(78)=195.08
                dtbase_mass(79)=196.97
                dtbase_mass(80)=200.59
                dtbase_mass(81)=204.38
                dtbase_mass(82)=207.20
                dtbase_mass(83)=208.98
                dtbase_mass(84)=209
                dtbase_mass(85)=210
                dtbase_mass(86)=222
                dtbase_mass(87)=223
                dtbase_mass(88)=226
                dtbase_mass(89)=227
                dtbase_mass(90)=232.04
                dtbase_mass(91)=231.04
                dtbase_mass(92)=238.03
                dtbase_mass(93)=237
                dtbase_mass(94)=244
                dtbase_mass(95)=243
                dtbase_mass(96)=247
                dtbase_mass(97)=247
                dtbase_mass(98)=251
                dtbase_mass(99)=252
                dtbase_mass(100)=257
                dtbase_mass(101)=258
                dtbase_mass(102)=259
                dtbase_mass(103)=262
                dtbase_mass(104)=267
                dtbase_mass(105)=268
                dtbase_mass(106)=269
                dtbase_mass(107)=270
                dtbase_mass(108)=269
                dtbase_mass(109)=277
                dtbase_mass(110)=281
                dtbase_mass(111)=282
                dtbase_mass(112)=285
                dtbase_mass(113)=286
                dtbase_mass(114)=290
                dtbase_mass(115)=290
                dtbase_mass(116)=293
                dtbase_mass(117)=294
                dtbase_mass(118)=294

                do k=1,mass_n
                   do l=1,dtbase_len
                   if (mass_labels(k).eq.dtbase_labels(l)) then
                        mass_mass(k)=dtbase_mass(l)
                   end if
                   end do
                   if (mass_mass(k).eq.0.0d0) then
                             write(err_msg,*) 'Atom ',mass_labels(k),' is not recognised or not supported!'
                             call mk_error()
                     end if
                end do

                deallocate(dtbase_mass)
                deallocate(dtbase_labels)

end subroutine

