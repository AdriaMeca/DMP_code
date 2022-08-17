# Author: Adria Meca Montserrat.
# Last modified date: 17/08/22.

# Function that filters out strings containing 'v' from a given list of strings.
my_filter = $(foreach v,$(2),$(if $(findstring $(1),$(v)),,$(v)))

modules := ../modules/
mod_dir := ../mod_dir/
obj_dir := ../obj_dir/

# Compilation flags.
options := -I $(mod_dir) -J $(mod_dir) -g -Wall -Wextra -pedantic -std=f2003 -O3

# Modules of interest and their resulting object files.
sources := $(wildcard $(modules)*.f90)
objects := $(patsubst $(modules)%.f90,$(obj_dir)%.o,$(sources))

# Lists of commonly used object files.
list1 := $(filter-out %types.o %generator.o,$(objects))
list2 := $(filter-out %procedures.o %algorihtms.o %simulations.o %properties.o,$(list1))

# Main program and its byproducts.
src := ./pz_simulation.f90
obj := $(patsubst %.f90,%.o,$(src))
exe := $(patsubst %.o,%.exe,$(obj))

# Compilation instructions.
$(exe): $(obj)
	gfortran $(objects) $< -o $@
$(obj): $(src) $(objects)
	gfortran $(options) -c $< -o $@
$(objects): $(obj_dir)%.o: $(modules)%.f90
	gfortran $(options) -c $< -o $@

# Module interdependencies.
$(list1): $(obj_dir)derived_types.o
$(list2) $(obj_dir)dmp_algorithms.o: $(obj_dir)array_procedures.o
$(list2) $(obj_dir)mc_simulations.o: $(obj_dir)random_number_generator.o

$(obj_dir)pz_simulation.o: $(filter %algorithms.o %simulations.o,$(list1))
$(obj_dir)network_generation.o $(obj_dir)rewiring_algorithms: $(obj_dir)network_properties.o

# Commands.
clean:
	rm -f $(exe) $(obj)

hard-clean:
	rm -f $(mod_dir)*.mod $(obj_dir)*.o $(exe) $(obj)

run:
	@$(exe)
