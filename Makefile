#===============================================================================
# Makefile for HNCR-nCLJ: Hypernetted Chain Solver for Charged LJ Mixtures
#===============================================================================

# Compiler selection (gfortran, ifx, ifort, nvfortran)
ifeq ($(origin FC),default)
  FC = gfortran
endif
FC ?= gfortran

# Target executable name
TARGET = hncr_nclj

# Directories
SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin

# Source files (ordered for compilation / dependencies)
SRCS = $(SRC_DIR)/types_mod.f90 \
       $(SRC_DIR)/params_mod.f90 \
       $(SRC_DIR)/grid_mod.f90 \
       $(SRC_DIR)/fft_mod.f90 \
       $(SRC_DIR)/arrays_mod.f90 \
       $(SRC_DIR)/linalg_mod.f90 \
       $(SRC_DIR)/potentials_mod.f90 \
       $(SRC_DIR)/initial_guess_mod.f90 \
       $(SRC_DIR)/thermo_mod.f90 \
       $(SRC_DIR)/hnc_solver_mod.f90 \
       $(SRC_DIR)/main.f90

OBJS = $(patsubst $(SRC_DIR)/%.f90, $(OBJ_DIR)/%.o, $(SRCS))

# Linear Algebra Libraries (LAPACK and BLAS)
# Overridable via: make LIBS="-lflexiblas" or make LIBS="-qmkl"
ifeq ($(findstring ifx,$(FC)),ifx)
  LIBS ?= -qmkl
else ifeq ($(findstring ifort,$(FC)),ifort)
  LIBS ?= -mkl
else
  # Default for gfortran: OpenBLAS (includes LAPACK) or FlexiBLAS
  ifneq ($(EBROOTOPENBLAS),)
    LIBS ?= -L$(EBROOTOPENBLAS)/lib -Wl,-rpath,$(EBROOTOPENBLAS)/lib -lopenblas
  else ifneq ($(wildcard /usr/local/modules/x64_v4/software/OpenBLAS/0.3.30-GCC-14.3.0/lib),)
    LIBS ?= -L/usr/local/modules/x64_v4/software/OpenBLAS/0.3.30-GCC-14.3.0/lib -Wl,-rpath,/usr/local/modules/x64_v4/software/OpenBLAS/0.3.30-GCC-14.3.0/lib -lopenblas
  else
    LIBS ?= -lopenblas
  endif
endif

# Compiler flags
ifeq ($(FC),gfortran)
  FFLAGS_OPT   = -O3 -march=native -funroll-loops -ffast-math -finline-functions
  FFLAGS_DEBUG = -g -O0 -Wall -Wextra -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow
  MODULE_FLAG  = -J$(OBJ_DIR)
else ifeq ($(findstring if,$(FC)),if)
  FFLAGS_OPT   = -O3 -xHost -ipo -fp-model fast=2
  FFLAGS_DEBUG = -g -O0 -check all -traceback -fpe0
  MODULE_FLAG  = -module $(OBJ_DIR)
else
  FFLAGS_OPT   = -O3
  FFLAGS_DEBUG = -g -O0
  MODULE_FLAG  = -I$(OBJ_DIR)
endif

# Build mode: release (default) or debug
MODE ?= release
ifeq ($(MODE),debug)
  FFLAGS = $(FFLAGS_DEBUG) $(MODULE_FLAG)
else
  FFLAGS = $(FFLAGS_OPT) $(MODULE_FLAG)
endif

.PHONY: all release debug clean distclean test help

all: $(BIN_DIR)/$(TARGET)

release:
	$(MAKE) MODE=release all

debug:
	$(MAKE) MODE=debug all

$(BIN_DIR)/$(TARGET): $(OBJS) | $(BIN_DIR)
	$(FC) $(FFLAGS) $(OBJS) $(LIBS) -o $@
	@echo "Build successful: $@"

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -o $@

# Explicit module dependencies
$(OBJ_DIR)/params_mod.o: $(OBJ_DIR)/types_mod.o
$(OBJ_DIR)/grid_mod.o: $(OBJ_DIR)/types_mod.o $(OBJ_DIR)/params_mod.o
$(OBJ_DIR)/fft_mod.o: $(OBJ_DIR)/types_mod.o $(OBJ_DIR)/params_mod.o
$(OBJ_DIR)/arrays_mod.o: $(OBJ_DIR)/types_mod.o $(OBJ_DIR)/params_mod.o $(OBJ_DIR)/grid_mod.o $(OBJ_DIR)/fft_mod.o
$(OBJ_DIR)/linalg_mod.o: $(OBJ_DIR)/types_mod.o
$(OBJ_DIR)/potentials_mod.o: $(OBJ_DIR)/types_mod.o $(OBJ_DIR)/params_mod.o $(OBJ_DIR)/grid_mod.o $(OBJ_DIR)/fft_mod.o $(OBJ_DIR)/arrays_mod.o
$(OBJ_DIR)/initial_guess_mod.o: $(OBJ_DIR)/types_mod.o $(OBJ_DIR)/params_mod.o $(OBJ_DIR)/grid_mod.o $(OBJ_DIR)/fft_mod.o
$(OBJ_DIR)/thermo_mod.o: $(OBJ_DIR)/types_mod.o $(OBJ_DIR)/params_mod.o $(OBJ_DIR)/grid_mod.o $(OBJ_DIR)/fft_mod.o $(OBJ_DIR)/arrays_mod.o $(OBJ_DIR)/linalg_mod.o
$(OBJ_DIR)/hnc_solver_mod.o: $(OBJ_DIR)/types_mod.o $(OBJ_DIR)/params_mod.o $(OBJ_DIR)/grid_mod.o $(OBJ_DIR)/fft_mod.o $(OBJ_DIR)/arrays_mod.o $(OBJ_DIR)/linalg_mod.o
$(OBJ_DIR)/main.o: $(OBJS:$(OBJ_DIR)/main.o=)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

test: all
	@echo "Running test calculation with datosLJ.hncr..."
	@$(BIN_DIR)/$(TARGET) datosLJ.hncr
	@echo "Test completed successfully."

clean:
	rm -rf $(OBJ_DIR) *.mod

distclean: clean
	rm -rf $(BIN_DIR) Salida.hncr Thermohncr.dat RDFhncr.dat skhnc.dat lsol.hncr fort.* *.mod

help:
	@echo "HNCR-nCLJ Build System"
	@echo "Usage:"
	@echo "  make              - Build optimized release executable (bin/hncr_nclj)"
	@echo "  make MODE=debug   - Build with debugging flags and runtime checks"
	@echo "  make test         - Run test calculation with sample input"
	@echo "  make clean        - Remove object and module files"
	@echo "  make distclean    - Remove build artifacts, binaries, and calculation outputs"
