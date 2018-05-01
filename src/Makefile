#
#to use intel debugger, use -g flag and run the program by typing "idb -gui makeda" at the prompt.
#Ifort flags
#FFLAGS = -132 -r8 -save -w -noarg_temp_created -zero -fpconstant -traceback -g #-warn nousage -C #-shared-intel
#gfortran flags
FFLAGS = -ffree-line-length-132 -fdefault-real-8 -fdefault-double-8 -fno-automatic -w -finit-local-zero -fbacktrace -g

MESA_DIR = /home/axk55/mesa-r8118
include $(MESA_DIR)/utils/makefile_header

FC = gfortran
#FC = ifort

%.o :	%.f90
	$(FC) -c $(FFLAGS) -o $@ $<

% :	%.o
	$(FC) $(FFLAGS) -o $@ $<

OBJ1 = wd_eos_mod.o
OBJ2 = wd_opalz_mod.o
OBJ3 = block_data.o commonblocks.o
OBJ4 = utils_subroutines.o ocon_functions.o phase_functions.o nuax_subroutines.o
OBJ5dx = chemprofiles_subroutines.o
OBJ5he = chemprofiles_subroutines.heliumcore.o
OBJ6 = eprep_subroutines.o istat_subroutines.o
OBJ7 = evol_subroutines.o pulse_subroutines.o
OBJ8 = evol_mainroutine.o pulse_mainroutine.o calcp_mainroutine.o
OBJ9 = getpar_grid.o

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

makeda: 	$(OBJ3) $(OBJ4) $(OBJ5dx) $(OBJ6) $(OBJ7) $(OBJ8) $(OBJ9)
	$(FC)	$(FFLAGS) -o ../$@ $^ 

makedx_mesa: 	$(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(OBJ5dx) $(OBJ6) $(OBJ7) $(OBJ8) $(OBJ9)
	$(LOADER) $(FCopenmp) -o ../$@ $^ -L$(MESA_LIB_DIR) $(LOAD_MESA_MICRO) 

makeda_he: 	$(OBJ1) $(OBJ2) $(OBJ3) $(OBJ4) $(OBJ5he) $(OBJ6) $(OBJ7) $(OBJ8) $(OBJ9)
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

