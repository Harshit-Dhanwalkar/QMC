CC     = gcc # or cc
CFLAGS = -Wall -Wextra -O2 -fPIC -fopenmp
# CFLAGS += -Wno-error=implicit-function-declaration
CFLAGS += -fsanitize=address -g
CFLAGS += -MMD -MP

OUTPUT_DIR  = output
CFLAGS     += -DQMC_OUTPUT_DIR=\"$(OUTPUT_DIR)\"

# Plot backend selection: GR (default) -> GNUPLOT -> MATPLOTLIB -> NONE
PLOT_BACKEND ?= GR
GR_PREFIX    ?= $(PWD)/third_party/gr/install

PLOT_SRC     =
PLOT_CFLAGS  =
PLOT_LDFLAGS =

ifeq ($(PLOT_BACKEND),GR)
    GR_INC   := $(GR_PREFIX)/include
    GR_LIB   := $(GR_PREFIX)/lib
    GR_AVAIL := $(shell test -f $(GR_LIB)/libGR.so && echo yes || echo no)
    ifeq ($(GR_AVAIL),no)
        $(warning GR library not found at $(GR_LIB)/libGR.so)
        override PLOT_BACKEND := GNUPLOT
    endif
endif

ifeq ($(PLOT_BACKEND),GNUPLOT)
    GNUPLOT_AVAIL := $(shell command -v gnuplot >/dev/null 2>&1 && echo yes || echo no)
    ifeq ($(GNUPLOT_AVAIL),no)
        $(warning gnuplot not found in PATH)
        override PLOT_BACKEND := MATPLOTLIB
    endif
endif

ifeq ($(PLOT_BACKEND),MATPLOTLIB)
    MATPLOTLIB_AVAIL := $(shell command -v python3 >/dev/null 2>&1 && python3 -c "import matplotlib" >/dev/null 2>&1 && echo yes || echo no)
    ifeq ($(MATPLOTLIB_AVAIL),no)
        $(warning matplotlib not found for python3)
        override PLOT_BACKEND := NONE
    endif
endif

# TODO: define macro -DHAVE_LATEX to conditionally compile LaTeX features
LATEX_AVAIL := $(shell command -v pdflatex >/dev/null 2>&1 && command -v pdftoppm >/dev/null 2>&1 && echo yes || echo no)
ifeq ($(LATEX_AVAIL),no)
    $(warning pdflatex or pdftoppm not found in PATH – LaTeX rendering will fail at runtime)
endif

# LaTeX backend selection : pdflatex and lualatex 
LATEX_COMPILER ?= pdflatex


ifeq ($(PLOT_BACKEND),GR)
    PLOT_SRC     = export/plot_gr.c \
                   export/gr/gr_plot.c \
                   export/gr/formats/png.c \
                   export/gr/formats/jpeg.c \
                   export/gr/formats/svg.c \
                   export/gr/formats/pdf.c \
                   export/gr/formats/interactive.c
    PLOT_LIBS    = -lGR -lGR3 -lqhull_r
    PLOT_CFLAGS  = -DUSE_GR -I$(GR_INC) -Iexport
    PLOT_LDFLAGS = -L$(GR_LIB) -Wl,-rpath,$(GR_LIB) -L/usr/lib/x86_64-linux-gnu
    $(info Plot backend: GR ($(GR_LIB)))
else ifeq ($(PLOT_BACKEND),GNUPLOT)
    PLOT_SRC     = export/plot_gnuplot.c \
                   export/gnuplot/gnuplot_pipe.c
    PLOT_CFLAGS  = -DUSE_GNUPLOT -Iexport
    PLOT_LDFLAGS =
    $(info Plot backend: GNUPLOT)
else ifeq ($(PLOT_BACKEND),MATPLOTLIB)
    PLOT_SRC     = export/plot_matplotlib.c \
                   export/matplotlib/matplotlib_pipe.c
    PLOT_CFLAGS  = -DUSE_MATPLOTLIB -Iexport
    PLOT_LDFLAGS =
    $(info Plot backend: MATPLOTLIB)
else
    PLOT_BACKEND := NONE
    PLOT_SRC     = export/plot_none.c
    PLOT_CFLAGS  = -Iexport
    PLOT_LDFLAGS =
    $(warning No usable plotting backend found - building NO-OP backend)
endif

CFLAGS += -DQMC_PLOT_BACKEND_NAME=\"$(PLOT_BACKEND)\"
CFLAGS += -DQMC_LATEX_COMPILER=\"$(LATEX_COMPILER)\"
BACKEND_SENTINEL := $(BUILD_DIR)/.plot_backend
CFLAGS  += $(PLOT_CFLAGS)
LDFLAGS  = -lm $(PLOT_LDFLAGS) $(PLOT_LIBS)

# Additive-only hooks for CI/ad-hoc flag injection
EXTRA_CFLAGS  ?=
EXTRA_LDFLAGS ?=
CFLAGS  += $(EXTRA_CFLAGS)
LDFLAGS += $(EXTRA_LDFLAGS)

# Directories
CORE_DIR     = core
PHYSICS_DIR  = physics
EXPORT_DIR   = export
LATEX_DIR    = latex
EXAMPLES_DIR = examples
TESTS_DIR    = tests
BUILD_DIR    = build

# Source files
CORE_SRCS    = $(CORE_DIR)/vector.c \
               $(CORE_DIR)/matrix.c \
               $(CORE_DIR)/utils.c \
               $(CORE_DIR)/random.c \
               $(CORE_DIR)/sparse.c \
               $(CORE_DIR)/fft/fft.c \
               $(CORE_DIR)/fft/fft2d.c \
               $(CORE_DIR)/fft/fft3d.c \
               $(CORE_DIR)/fft/fft_wrapper.c \
               $(CORE_DIR)/linalg/eigen_generic.c \
               $(CORE_DIR)/linalg/tridiag_eigen.c \
               $(CORE_DIR)/linalg/tridiag_eigh.c \
               $(CORE_DIR)/linalg/qr.c \
               $(CORE_DIR)/linalg/lu.c \
               $(CORE_DIR)/linalg/svd.c \
               $(CORE_DIR)/linalg/complex_eigh.c \
               $(CORE_DIR)/ode/numerov.c \
               $(CORE_DIR)/ode/rk4.c \
               $(CORE_DIR)/ode/crank_nicolson.c \
               $(CORE_DIR)/special/hermite.c \
               $(CORE_DIR)/special/laguerre.c \
               $(CORE_DIR)/special/legendre.c \
               $(CORE_DIR)/special/bessel.c \
               $(CORE_DIR)/special/polynomials.c \
               $(CORE_DIR)/special/spherical_harmonics.c

PHYSICS_SRCS = $(PHYSICS_DIR)/potentials.c \
               $(PHYSICS_DIR)/wavefn.c \
               $(PHYSICS_DIR)/schrodinger.c \
               $(PHYSICS_DIR)/uncertainty.c \
               $(PHYSICS_DIR)/angular.c \
               $(PHYSICS_DIR)/central_potential.c \
               $(PHYSICS_DIR)/hydrogen.c \
               $(PHYSICS_DIR)/helium.c \
               $(PHYSICS_DIR)/perturbation.c \
               $(PHYSICS_DIR)/variational.c \
               $(PHYSICS_DIR)/wkb.c \
               $(PHYSICS_DIR)/scattering.c \
               $(PHYSICS_DIR)/rabi.c \
               $(PHYSICS_DIR)/identical.c \
               $(PHYSICS_DIR)/relativistic.c \
               $(PHYSICS_DIR)/fine_structure.c \
               $(PHYSICS_DIR)/qubits.c \
               $(PHYSICS_DIR)/lindblad.c \
               $(PHYSICS_DIR)/hartree_fock.c \
               $(PHYSICS_DIR)/driven.c \
               $(PHYSICS_DIR)/soft.c \
               $(PHYSICS_DIR)/boson_sampling.c \
               $(PHYSICS_DIR)/zeeman.c \
               $(PHYSICS_DIR)/vmc.c \
               $(PHYSICS_DIR)/dmc.c \
               $(PHYSICS_DIR)/pimc.c \
               $(PHYSICS_DIR)/vqe.c \
               $(PHYSICS_DIR)/mp2.c \
               $(PHYSICS_DIR)/lattice.c \
               $(PHYSICS_DIR)/quantum_info.c \
               $(PHYSICS_DIR)/qec.c \
               $(PHYSICS_DIR)/quantum_algorithms.c \
               $(PHYSICS_DIR)/second_quant.c \
               $(PHYSICS_DIR)/molecular_integrals.c

LATEX_SRCS   = $(LATEX_DIR)/latex_gen.c

PLOT_SRCS    = $(PLOT_SRC)

# Object files
ALL_SRCS    = $(CORE_SRCS) $(PHYSICS_SRCS) $(LATEX_SRCS) $(PLOT_SRC)
ALL_OBJS    = $(patsubst %.c,$(BUILD_DIR)/%.o,$(ALL_SRCS))

-include $(ALL_OBJS:.o=.d)

# Targets
EXAMPLES    = $(BUILD_DIR)/eg_01_particle_box \
              $(BUILD_DIR)/eg_02_harmonic \
              $(BUILD_DIR)/eg_03_finite_well \
              $(BUILD_DIR)/eg_04_infinite_well \
              $(BUILD_DIR)/eg_05_uncertainty \
              $(BUILD_DIR)/eg_06_hydrogen \
              $(BUILD_DIR)/eg_07_central_potential \
              $(BUILD_DIR)/eg_08_helium \
              $(BUILD_DIR)/eg_09_identical_particles \
              $(BUILD_DIR)/eg_10_perturbation \
              $(BUILD_DIR)/eg_11_wkb \
              $(BUILD_DIR)/eg_12_tunnelling \
              $(BUILD_DIR)/eg_13_scattering \
              $(BUILD_DIR)/eg_14_rabi \
              $(BUILD_DIR)/eg_15_angular_coupling \
              $(BUILD_DIR)/eg_16_finestructure \
              $(BUILD_DIR)/eg_17_dirac \
              $(BUILD_DIR)/eg_18_qubits \
              $(BUILD_DIR)/eg_19_lindblad \
              $(BUILD_DIR)/eg_20_hartree_fock \
              $(BUILD_DIR)/eg_21_driven \
              $(BUILD_DIR)/eg_22_soft \
              $(BUILD_DIR)/eg_23_fermi_golden_rule \
              $(BUILD_DIR)/eg_24_central_potential_3d \
              $(BUILD_DIR)/eg_25_boson_sampling \
              $(BUILD_DIR)/eg_26_zeeman \
              $(BUILD_DIR)/eg_27_cap_tdse \
              $(BUILD_DIR)/eg_28_vmc_helium \
              $(BUILD_DIR)/eg_29_dmc_helium \
              $(BUILD_DIR)/eg_30_pimc_helium \
              $(BUILD_DIR)/eg_31_vqe \
              $(BUILD_DIR)/eg_32_mp2 \
              $(BUILD_DIR)/eg_33_lattice \
              $(BUILD_DIR)/eg_34_quantum_info \
              $(BUILD_DIR)/eg_35_qec \
              $(BUILD_DIR)/eg_36_quantum_algorithms \
              $(BUILD_DIR)/eg_37_second_quant \
              $(BUILD_DIR)/eg_38_landau_levels

TESTS       = $(BUILD_DIR)/test_complex \
              $(BUILD_DIR)/test_matrix \
              $(BUILD_DIR)/test_numerov \
              $(BUILD_DIR)/test_rk4 \
              $(BUILD_DIR)/test_fft \
              $(BUILD_DIR)/test_hydrogen \
              $(BUILD_DIR)/test_helium \
              $(BUILD_DIR)/test_perturbation \
              $(BUILD_DIR)/test_crank_nicolson \
              $(BUILD_DIR)/test_wkb \
              $(BUILD_DIR)/test_potentials \
              $(BUILD_DIR)/test_central_potential \
              $(BUILD_DIR)/test_central_potential_3d \
              $(BUILD_DIR)/test_rabi \
              $(BUILD_DIR)/test_angular_coupling \
              $(BUILD_DIR)/test_fine_structure \
              $(BUILD_DIR)/test_complex_eigh \
              $(BUILD_DIR)/test_identical \
              $(BUILD_DIR)/test_dirac \
              $(BUILD_DIR)/test_qubits \
              $(BUILD_DIR)/test_scattering \
              $(BUILD_DIR)/test_tridiag \
              $(BUILD_DIR)/test_lindblad \
              $(BUILD_DIR)/test_hartree_fock \
              $(BUILD_DIR)/test_driven \
              $(BUILD_DIR)/test_soft \
              $(BUILD_DIR)/test_fermi_golden_rule \
              $(BUILD_DIR)/test_boson_sampling \
              $(BUILD_DIR)/test_zeeman \
              $(BUILD_DIR)/test_klein_gordon \
              $(BUILD_DIR)/test_tridiag_eigvals \
              $(BUILD_DIR)/test_random \
              $(BUILD_DIR)/test_lanczos \
              $(BUILD_DIR)/test_vmc \
              $(BUILD_DIR)/test_dmc \
              $(BUILD_DIR)/test_pimc \
              $(BUILD_DIR)/test_vqe \
              $(BUILD_DIR)/test_mp2 \
              $(BUILD_DIR)/test_lattice \
              $(BUILD_DIR)/test_quantum_info \
              $(BUILD_DIR)/test_qec \
              $(BUILD_DIR)/test_quantum_algorithms \
              $(BUILD_DIR)/test_second_quant

ifeq ($(PLOT_BACKEND),GR)
    TESTS += $(BUILD_DIR)/test_grplot
endif

# Demo driver
# $(BUILD_DIR)/main: main.c config.h $(BACKEND_SENTINEL) | directories
$(BUILD_DIR)/main: main.c config.h | directories
	$(CC) $(CFLAGS) -I. main.c -o $@

demo: $(BUILD_DIR)/main
	@$(BUILD_DIR)/main

.PHONY: all clean examples tests run-examples run-tests info

# check-backend

check-backend:
	@mkdir -p $(BUILD_DIR)
	@echo "$(PLOT_BACKEND)" > $(BUILD_DIR)/.plot_backend.new
	@if ! cmp -s $(BUILD_DIR)/.plot_backend.new $(BACKEND_SENTINEL) 2>/dev/null; then \
		mv $(BUILD_DIR)/.plot_backend.new $(BACKEND_SENTINEL); \
	else \
		rm -f $(BUILD_DIR)/.plot_backend.new; \
	fi

$(BACKEND_SENTINEL): check-backend
	@:

all: directories $(OUTPUT_DIR) $(EXAMPLES) $(TESTS) $(BUILD_DIR)/main

directories:
	@mkdir -p $(BUILD_DIR)
	@mkdir -p $(BUILD_DIR)/$(CORE_DIR)/linalg
	@mkdir -p $(BUILD_DIR)/$(CORE_DIR)/ode
	@mkdir -p $(BUILD_DIR)/$(CORE_DIR)/special
	@mkdir -p $(BUILD_DIR)/$(CORE_DIR)/fft
	@mkdir -p $(BUILD_DIR)/$(PHYSICS_DIR)
	@mkdir -p $(BUILD_DIR)/$(EXPORT_DIR)
	@mkdir -p $(BUILD_DIR)/$(EXPORT_DIR)/gr/formats
	@mkdir -p $(BUILD_DIR)/$(EXPORT_DIR)/gnuplot
	@mkdir -p $(BUILD_DIR)/$(EXPORT_DIR)/matplotlib
	@mkdir -p $(BUILD_DIR)/$(LATEX_DIR)
# TODO: generic auto create target subdir for all .o instead of using `directories`
# $(BUILD_DIR)/%.o: %.c
# 	@mkdir -p $(dir $@)
# 	$(CC) $(CFLAGS) -Icore -Iexport -I. -c $< -o $@

$(OUTPUT_DIR):
	@mkdir -p $(OUTPUT_DIR)

$(BUILD_DIR)/%.o: %.c
	$(CC) $(CFLAGS) -Icore -Iexport -I. -c $< -o $@

# Examples
$(BUILD_DIR)/eg_%: $(EXAMPLES_DIR)/eg_%.c $(ALL_OBJS) | $(OUTPUT_DIR)
	$(CC) $(CFLAGS) -Icore -Iexport -I. $^ -o $@ $(LDFLAGS)

# Tests
$(BUILD_DIR)/test_%: $(TESTS_DIR)/test_%.c $(ALL_OBJS)
	$(CC) $(CFLAGS) -Icore -Iexport -I. $^ -o $@ $(LDFLAGS)

# Run targets
examples: directories $(OUTPUT_DIR) $(EXAMPLES)
tests:    directories $(TESTS)

# Run a all tests
run-tests: tests
	@echo "Running all tests..."
	@passed=0; failed=0; \
	for t in $(TESTS); do \
		name=$$(basename $$t); \
		printf "Running $$name... "; \
		if $$t > /tmp/$$name.out 2>&1; then \
			printf "\033[32mPASS\033[0m\n"; \
			passed=$$((passed+1)); \
		else \
			printf "\033[31mFAIL\033[0m\n"; \
			cat /tmp/$$name.out; \
			failed=$$((failed+1)); \
		fi; \
	done; \
	echo ""; \
	echo "Results: $$passed passed, $$failed failed"; \
	[ $$failed -eq 0 ]

# Run a specific test
test-%: $(BUILD_DIR)/test_%
	@$

# Run a all examples
run-examples: examples $(OUTPUT_DIR)
	@echo "    Running examples"
	@for ex in $(EXAMPLES); do \
		name=$$(basename $$ex); \
		echo ""; \
		echo " -> $$name ";\
		$$ex || true; \
	done
	@echo ""
	@echo "Output files in: $(OUTPUT_DIR)/"

# Info target
info:
	@echo "Plot backend : $(PLOT_BACKEND)"
	@echo "GR prefix    : $(GR_PREFIX)"
	@echo "Output dir   : $(OUTPUT_DIR)"

# Clean build,test and examples outputs
clean:
	rm -rf $(BUILD_DIR) $(OUTPUT_DIR)
	rm -f *.png *.pdf *.svg *.jpg *.dat
	rm -f *.aux *.log *.out *.toc *.tex
