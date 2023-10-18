#
#to use intel debugger, use -g flag and run the program by typing "idb -gui makeda" at the prompt.
#Ifort flags
#FFLAGS = -132 -r8 -save -w -noarg_temp_created -zero -fpconstant -traceback -g #-warn nousage -C #-shared-intel
#gfortran flags
FFLAGS = -ffree-line-length-132 -fdefault-real-8 -fdefault-double-8 -fno-automatic -w -finit-local-zero -fbacktrace -g

MESA_DIR = /home/agnes/mesa
include $(MESA_DIR)/utils/makefile_header

FC = gfortran
#FC = ifort

%.o :	%.f90
	$(FC) -c $(FFLAGS) -o $@ $<

% :	%.o
	$(FC) $(FFLAGS) -o $@ $<

OBJ1 = wd_eos_mod.o
OBJ2 = wd_opalz_mod.o
OBJ3v15 = block_data.o commonblocks_v15.o
OBJ3v20 = block_data.o commonblocks_v20.o
OBJ4 = utils_subroutines.o ocon_functions.o phase_functions.o nuax_subroutines.o
OBJ5dxv15 = chemprofiles_subroutines_v15.o
OBJ5dxv20 = chemprofiles_subroutines_v20.o
OBJ5he = chemprofiles_subroutines.heliumcore.o
OBJ6 = eprep_subroutines.o istat_subroutines.o
OBJ7v15 = evol_subroutines_v15.o pulse_subroutines.o
OBJ7v20 = evol_subroutines_v20.o pulse_subroutines.o
OBJ8v15 = evol_mainroutine.o pulse_mainroutine.o calcp_mainroutine_v15.o
OBJ8v20 = evol_mainroutine.o pulse_mainroutine.o calcp_mainroutine_v20.o
OBJ9v15 = getpar_grid_v15.o
OBJ9v20 = getpar_grid_v20.o

#OBJ1
#commonblocks.o: commonblocks.f90 
#	$(FC)  $(FFLAGS) -c $<
#block_data.o: block_data.f90
#	$(FC)  $(FFLAGS) -c $<

#OBJ2
#utils_subroutines.o: utils_subroutines.f90
#	$(FC)  $(FFLAGS) -c $<
#ocon_functions.o: ocon_functions.f90
#	$(FC)  $(FFLAGS) -c $<
#phase_functions.o: phase_functions.f90
#	$(FC)  $(FFLAGS) -c $<
#nuax_subroutines.o: nuax_subroutines.f90
#	$(FC)  $(FFLAGS) -c $<

#OBJ3
#eprep_subroutines.o: eprep_subroutines.f90
#	$(FC)  $(FFLAGS) -c $<
#istat_subroutines.o: istat_subroutines.f90
#	$(FC)  $(FFLAGS) -c $<

#OBJ4
#evol_subroutines.o: evol_subroutines.f90
#	$(FC)  $(FFLAGS) -c $<
#pulse_subroutines.o: pulse_subroutines.f90
#	$(FC)  $(FFLAGS) -c $<

#OBJ5
#chemprofiles_subroutines.o: chemprofiles_subroutines.f90
#	$(FC)  $(FFLAGS) -c $<

#OBJ5 
#evol_mainroutine.o: evol_mainroutine.f90
#	$(FC)  $(FFLAGS) -c $<
#pulse_mainroutine.o: pulse_mainroutine.f90
#	$(FC)  $(FFLAGS) -c $<
#calcp_mainroutine.o: calcp_mainroutine.f90
#	$(FC)  $(FFLAGS) -c $<

#OBJ6 
#getpar_grid.o: 	getpar_grid.f90 
#	$(FC) $(FFLAGS) -c $<

#For fixed oxygen profiles
makedx_v15: 	$(OBJ1) $(OBJ2) $(OBJ3v15) $(OBJ4) $(OBJ5dxv15) $(OBJ6) $(OBJ7v15) $(OBJ8v15) $(OBJ9v15)
	$(LOADER) $(FCopenmp) -o ../$@ $^ -L$(MESA_LIB_DIR) $(LOAD_MESA_MICRO) 

#For asteroseismic fitting of oxygen profiles (and other structure of the model)
makedx_v20: 	$(OBJ1) $(OBJ2) $(OBJ3v20) $(OBJ4) $(OBJ5dxv20) $(OBJ6) $(OBJ7v20) $(OBJ8v20) $(OBJ9v20)
	$(LOADER) $(FCopenmp) -o ../$@ $^ -L$(MESA_LIB_DIR) $(LOAD_MESA_MICRO) 

makeda_he: 	$(OBJ1) $(OBJ2) $(OBJ3v20) $(OBJ4) $(OBJ5he) $(OBJ6) $(OBJ7v20) $(OBJ8v20) $(OBJ9v20)
	$(LOADER) $(FCopenmp) -o ../$@ $^ -L$(MESA_LIB_DIR) $(LOAD_MESA_MICRO)

wd_eos_mod.o:	wd_eos_mod.f90
	$(FC) $(FCfree) $(MY_FC_FLAGS) -I$(MESA_INCLUDE_DIR) -c $<

wd_opalz_mod.o:	wd_opalz_mod.f90	
	$(FC) $(FCfree) $(MY_FC_FLAGS) -I$(MESA_INCLUDE_DIR) -c $<

wd_test.o:	wd_test.f90	
	$(FC) $(FCfree) $(MY_FC_FLAGS) -I$(MESA_INCLUDE_DIR) -c $<

wd_test:	wd_eos_mod.o wd_opalz_mod.o wd_test.o	
	$(LOADER) $(FCopenmp) $^ -L$(MESA_LIB_DIR) $(LOAD_MESA_MICRO) -o $@

clean:		
		rm *.o *.mod *~

eos_wd.o: 	eos_wd.f90
	$(FC) $(FCfree) $(MY_FC_FLAGS) -I$(MESA_INCLUDE_DIR) -c $<

