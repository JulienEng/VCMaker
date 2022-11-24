FC          = gfortran
FLAGS_DEBUG =  -g -fcheck=all -Wall
FLAGS_OPT   =
FLAGS_LIB   = -llapack
FLAGS_FORM  = -ffixed-line-length-none -ffree-form

FLAGS       =  $(FLAGS_MPI) $(FLAGS_OPT) $(FLAGS_LIB) $(FLAGS_FORM) $(FLAGS_DEBUG)

.f.o:
	$(FC) -c $(FLAGS) $(FLAGS_DEBUG) $< 

OBJS = mod_input.o \
       mod_constants.o \
       mod_disp.o \
       mod_dncscan.o \
       mod_grad.o \
       mod_gap.o \
       mod_lambda.o \
       mod_quantics.o \
       mod_output.o \
       mod_min.o \
       mod_err.o

ai: $(OBJS)

