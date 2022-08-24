# Author: Adria Meca Montserrat.
# Last modified date: 24/08/22.

# Command line arguments.
ifeq ($(name),)
	name := probabilities
endif

ifeq ($(type),)
	type := probabilities
endif

# Commonly used folders.
modules := ../modules/
mod_dir := ../mod_dir/
obj_dir := ../obj_dir/

# Dependencies and attachments of the main program.
dependencies := $(obj_dir)derived_types.o
dependencies += $(obj_dir)network_generation.o
dependencies += $(obj_dir)random_number_generator.o
dependencies += $(obj_dir)rewiring_algorithms.o

ifeq ($(type), probabilities)
	dependencies += $(obj_dir)dmp_algorithms.o
	dependencies += $(obj_dir)mc_simulations.o

	attachments := $(dependencies)
	attachments += $(obj_dir)array_procedures.o
	attachments += $(obj_dir)network_properties.o
else
	ifeq ($(type), pz_simulation)
		dependencies += $(obj_dir)array_procedures.o
		dependencies += $(obj_dir)patient_zero_problem.o

		attachments := $(dependencies)
		attachments += $(obj_dir)dmp_algorithms.o
		attachments += $(obj_dir)mc_simulations.o
		attachments += $(obj_dir)network_properties.o
	endif
endif

# Compilation flags.
options := -I $(mod_dir) -J $(mod_dir) -g -Wall -Wextra -pedantic -std=f2003 -O3

# Modules of interest and their resulting object files.
sources := $(wildcard $(modules)*.f90)
objects := $(patsubst $(modules)%.f90,$(obj_dir)%.o,$(sources))

# Main program and its byproducts.
src := ./$(name).f90
obj := $(patsubst %.f90,%.o,$(src))
exe := $(patsubst %.o,%.exe,$(obj))

# Compilation instructions.
$(exe): $(obj)
	gfortran $(attachments) $< -o $@

$(obj): $(src) $(dependencies)
	gfortran $(options) -c $< -o $@

$(objects): $(obj_dir)%.o: $(modules)%.f90
	gfortran $(options) -c $< -o $@

# Module interdependencies.
$(obj_dir)array_procedures.o: $(obj_dir)derived_types.o

$(obj_dir)dmp_algorithms.o: $(obj_dir)array_procedures.o
$(obj_dir)dmp_algorithms.o: $(obj_dir)derived_types.o

$(obj_dir)mc_simulations.o: $(obj_dir)derived_types.o
$(obj_dir)mc_simulations.o: $(obj_dir)random_number_generator.o

$(obj_dir)network_generation.o: $(obj_dir)array_procedures.o
$(obj_dir)network_generation.o: $(obj_dir)derived_types.o
$(obj_dir)network_generation.o: $(obj_dir)network_properties.o
$(obj_dir)network_generation.o: $(obj_dir)random_number_generator.o

$(obj_dir)network_properties.o: $(obj_dir)derived_types.o

$(obj_dir)patient_zero_problem.o: $(obj_dir)array_procedures.o
$(obj_dir)patient_zero_problem.o: $(obj_dir)derived_types.o
$(obj_dir)patient_zero_problem.o: $(obj_dir)dmp_algorithms.o
$(obj_dir)patient_zero_problem.o: $(obj_dir)mc_simulations.o
$(obj_dir)patient_zero_problem.o: $(obj_dir)random_number_generator.o

$(obj_dir)rewiring_algorithms.o: $(obj_dir)array_procedures.o
$(obj_dir)rewiring_algorithms.o: $(obj_dir)derived_types.o
$(obj_dir)rewiring_algorithms.o: $(obj_dir)network_properties.o
$(obj_dir)rewiring_algorithms.o: $(obj_dir)random_number_generator.o

# Commands.
clean:
	rm -f $(exe) $(obj)

hard-clean:
	rm -f $(mod_dir)*.mod $(obj_dir)*.o

run:
	@$(exe)
