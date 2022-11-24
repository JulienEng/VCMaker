# Makefile for VCMaker by J.Eng, 2021
#  With the help of T.Pope
FC          = gfortran
FLAGS_DEBUG = -g -fcheck=all -Wall -fbounds-check -fsanitize=address 
FLAGS_OPT   = -O2
FLAGS_MPI   =
FLAGS_FORM  = -ffixed-line-length-none -ffree-form
LIBRARY     = -llapack
MODULES     = ./modules
SUBDIRS     = ./subroutines
FLAGS       = $(FLAGS_MPI) $(FLAGS_DEBUG) $(FLAGS_OPT) $(FLAGS_FORM)
PROGRAM     = vcmaker.x
PRECOMP     = vcmaker.o
SOURCE      = vcmaker.f

SUBOBJS= $(SUBDIRS)/tools.o \
         $(SUBDIRS)/get_grad.o \
         $(SUBDIRS)/get_hessian.o \
         $(SUBDIRS)/get_input.o \
         $(SUBDIRS)/get_xyz.o \
         $(SUBDIRS)/header.o \
         $(SUBDIRS)/mk_dnc.o \
         $(SUBDIRS)/mk_dncscan.o \
         $(SUBDIRS)/mk_dncscan_diag.o \
         $(SUBDIRS)/mk_dncscan_grid.o \
         $(SUBDIRS)/mk_hess.o \
         $(SUBDIRS)/mk_kappa_grad.o \
         $(SUBDIRS)/mk_kappa_disp.o \
         $(SUBDIRS)/mk_gap.o \
         $(SUBDIRS)/mk_lambda_work.o \
         $(SUBDIRS)/mk_nm.o \
         $(SUBDIRS)/mk_genxyz.o \
         $(SUBDIRS)/mk_min.o \
         $(SUBDIRS)/mk_convert.o \
         $(SUBDIRS)/out_6col.o \
         $(SUBDIRS)/out_input.o \
         $(SUBDIRS)/out_kgrad.o \
         $(SUBDIRS)/out_kdisp.o \
         $(SUBDIRS)/out_rsma.o \
         $(SUBDIRS)/out_lambda.o \
         $(SUBDIRS)/internal.o \
         $(SUBDIRS)/mk_error.o \
         $(SUBDIRS)/get_alignxyz.o \
         $(SUBDIRS)/out_quantics.o \
         $(SUBDIRS)/out_xyz.o  

MODOBJS= $(MODULES)/mod_input.o \
         $(MODULES)/mod_constants.o \
         $(MODULES)/mod_disp.o \
         $(MODULES)/mod_dncscan.o \
         $(MODULES)/mod_grad.o \
         $(MODULES)/mod_gap.o \
         $(MODULES)/mod_lambda.o \
         $(MODULES)/mod_output.o \
         $(MODULES)/mod_quantics.o \
         $(MODULES)/mod_min.o \
         $(MODULES)/mod_err.o

.PHONY: all $(MODULES) $(SUBDIRS) $(PROGRAM)

$(MODULES):
	$(MAKE) --directory=$@
        
$(SUBDIRS):
	$(MAKE) --directory=$@

all: $(MODULES) $(SUBDIRS)
	$(FC)  -I$(SUBDIRS) -I$(MODULES) $(FLAGS) $(LIBRARY) -c $(SOURCE) -o $(PRECOMP)
	$(FC)  $(SUBOBJS) $(MODOBJS) $(PRECOMP) $(FLAGS) $(LIBRARY) -o $(PROGRAM)
clean:
	rm $(MODULES)/*.o $(MODULES)/*.mod $(SUBDIRS)/*.o $(PRECOMP)
                                                                                   
