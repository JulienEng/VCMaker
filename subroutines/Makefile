FC          = gfortran
FLAGS_DEBUG =  -g -fcheck=all -Wall
FLAGS_OPT   =
FLAGS_LIB   = -llapack
FLAGS_FORM  = -ffixed-line-length-none -ffree-form

FLAGS       =  $(FLAGS_MPI) $(FLAGS_OPT) $(FLAGS_LIB) $(FLAGS_FORM) $(FLAGS_DEBUG)

.f.o:
	$(FC) -I../modules/ -c $(FLAGS) $(FLAGS_DEBUG) $< 

OBJS= tools.o \
      get_grad.o \
      get_hessian.o \
      get_input.o \
      get_xyz.o \
      header.o \
      mk_dnc.o \
      mk_dncscan.o \
      mk_dncscan_diag.o \
      mk_dncscan_grid.o \
      mk_hess.o \
      mk_kappa_grad.o \
      mk_kappa_disp.o \
      mk_gap.o \
      mk_lambda_work.o \
      mk_nm.o \
      mk_genxyz.o \
      mk_min.o \
      mk_convert.o \
      out_6col.o \
      out_input.o \
      out_kgrad.o \
      out_kdisp.o \
      out_rsma.o \
      out_lambda.o \
      get_alignxyz.o\
      out_quantics.o \
      out_xyz.o \
      mk_error.o \
      internal.o 

ai: $(OBJS)
        

