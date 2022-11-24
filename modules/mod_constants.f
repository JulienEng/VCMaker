module mod_constants

  double precision,parameter     :: pi      = 4.0*atan(1.0d0)
  double precision,parameter     :: hbar    = 1.05457162825177E-34   ! Reduced Planck Constant [J.s]
  double precision,parameter     :: hjs     = 6.62607015E-34         ! Planck Constant [J.s]
  double precision,parameter     :: c       = 299792458d2            ! Speed of light in [cm/s]
  double precision,parameter     :: bohrtom = 5.2917724900001E-11    ! Conversion factor: a0 to m
  double precision,parameter     :: amu     = 1.660538921E-27        ! Conversion factor: amu to kg
  double precision,parameter     :: me      = 9.10938356E-31         ! Conversion factor: me to kg
  double precision,parameter     :: hartoJ  = 4.3597482E-18          ! Conversion factor: Eh to J
  double precision,parameter     :: hartoeV = 27.21139664131         ! Conversion factor: Eh to eV
  double precision,parameter     :: angtoa0 = 1.8897259886           ! Conversion factor: Ang to a0
  double precision,parameter     :: tautos  = 2.42E-17               ! Conversion factor: u.a of time to s

end module mod_constants

