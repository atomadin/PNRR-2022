FC=ifort -module ./modules -diag-disable=10448
TROUT=2>&1 | head -n 10 | fold -s -w 72 | ccze -A
RM=rm -f
LIBS= 

# Default compilation flags.
FFLAGS=-fpp -warn all -heap-arrays -xHost -O3 -qopenmp
.PHONY: print-date

all: build print-date

build:
	$(FC) -o boltzmann_graphene.x $(FFLAGS) boltzmann_graphene.f90 $(TROUT)

print-date:
	@printf "\033[36m$$(date)\033[0m\n"


# Safe flags.  Run target with "make safe".
safe: FFLAGS=-fpp -O0 -warn all -check all -check bounds -debug extended -debug-parameters all -g -traceback
safe: all

# Optimization
# ------------
# -inline speed -unroll -ipo -parallel -qopenmp -no-prec-div
# -fno-inline-functions
# -unroll0 -fno-inline-function
# -qopenmp
# -qopt-report=1
# -qopt-report-phase=par
# -qopt-report-routine=calculate_pbe_gammas
#
# Profiling
# ---------
# With -p
# function analyzer creates file gmon.out
# then do the following to access report
# gprof boltzmann_graphene.x > gmon.txt
#
# Floating-point exceptions
# -------------------------
# -fpp: Fortran pre-processor
# -fpe[integer value]: floating-point exception
#      0: Abort on overflow and divide-by-zero, flush denormalized to zero.
#      1: All floating-point exceptions disabled. Some underflow set to zero.
#      3: All floating-point exceptions disabled.  Underflow gradual. (Default.)
# -ftz: flush denormal to zero
# -ipo: inter-procedural optimization
# -fast: -ipo -O3 -no-prec-div -static -xHost

.PHONY: clean
clean:
	@$(RM) ./modules/*.mod
	@$(RM) boltzmann_graphene.x









