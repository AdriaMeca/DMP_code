#Author: Adria Meca Montserrat.
#Last modified date: 29/05/22.

#Function that filters the strings that contain the substring 'v' from a given
#list of strings.
my_filter = $(foreach v,$(2),$(if $(findstring $(1),$(v)),,$(v)))

#Folder containing the modules whose procedures are used in the main simulations.
modules := ../modules/
#Folder containing the mod files that result from the compilation of the modules.
mod_dir := ../mod_dir/
#Folder containing the 'o' files that result from the compilation of the modules.
obj_dir := ../obj_dir/

#Compilation flags.
options := -I $(mod_dir) -J $(mod_dir) -g -Wall -Wextra -pedantic -std=f2003 -O3

#Modules of interest.
sources := $(call my_filter,properties,$(wildcard $(modules)*.f90))
#Resulting output files.
objects := $(patsubst $(modules)%.f90,$(obj_dir)%.o,$(sources))

#List of commonly used 'o' files.
in_list := $(filter $(obj_dir)dmp% $(obj_dir)mc%,$(objects))
out_list := $(filter-out %procedures.o %generator.o,$(objects))

#Main program and its byproducts.
src := ./patient_zero.f90
obj := $(patsubst %.f90,%.o,$(src))
exe := $(patsubst %.o,%.exe,$(obj))

#Compilation instructions.
$(exe): $(obj)
	gfortran $(objects) $< -o $@
$(obj): $(src) $(objects)
	gfortran $(options) -c $< -o $@
$(objects): $(obj_dir)%.o: $(modules)%.f90
	gfortran $(options) -c $< -o $@

#Module interdependencies.
$(out_list): $(obj_dir)array_procedures.o
$(filter %problem.o,$(out_list)): $(in_list)
$(call my_filter,dmp,$(out_list)): $(obj_dir)random_number_generator.o
$(call my_filter,generation,$(out_list)): $(obj_dir)network_generation.o

#Command that cleans up the byproducts of the main program's compilation.
clean:
	rm -f $(exe) $(obj)

#Command that removes all the files that result from the compilation process.
hard-clean:
	rm -f $(mod_dir)*.mod $(obj_dir)*.o $(exe) $(obj)

#Command that runs the executable of the main program.
run:
	@$(exe)
