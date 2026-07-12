CC     = gcc
CFLAGS = -Wall -Wextra -O2 -fPIC
CFLAGS += -fsanitize=address -g

OUTPUT_DIR  = output
CFLAGS     += -DQMC_OUTPUT_DIR=\"$(OUTPUT_DIR)\"

# Plot backend selection: GR (default) -> GNUPLOT -> MATPLOTLIB -> NONE
PLOT_BACKEND ?= GR
# GR_PREFIX    ?= $(PWD)/third_party/gr
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

ifeq ($(PLOT_BACKEND),GR)
    PLOT_SRC     = export/plot_gr.c export/gr/gr_plot.c
    PLOT_LIBS    = -lGR
    PLOT_CFLAGS  = -DUSE_GR -I$(GR_INC) -Iexport
    PLOT_LDFLAGS = -L$(GR_LIB) -Wl,-rpath,$(GR_LIB)
    $(info Plot backend: GR ($(GR_LIB)))
else ifeq ($(PLOT_BACKEND),GNUPLOT)
    PLOT_SRC     = export/plot_gnuplot.c export/gnuplot/gnuplot_pipe.c
    PLOT_CFLAGS  = -DUSE_GNUPLOT -Iexport
    PLOT_LDFLAGS =
    $(info Plot backend: GNUPLOT)
else ifeq ($(PLOT_BACKEND),MATPLOTLIB)
    PLOT_SRC     = export/plot_matplotlib.c export/matplotlib/matplotlib_pipe.c
    PLOT_CFLAGS  = -DUSE_MATPLOTLIB -Iexport
    PLOT_LDFLAGS =
    $(info Plot backend: MATPLOTLIB)
else
    PLOT_BACKEND := NONE
    PLOT_SRC     = export/plot_none.c
    PLOT_CFLAGS  = -Iexport
    PLOT_LDFLAGS =
    $(warning No usable plotting backend found (GR/gnuplot/matplotlib all unavailable) — building NO-OP backend; no plots will be generated)
endif

CFLAGS  += $(PLOT_CFLAGS)
LDFLAGS  = -lm $(PLOT_LDFLAGS) $(PLOT_LIBS)

# Directories
CORE_DIR     = core
PHYSICS_DIR  = physics
EXPORT_DIR   = export
LATEX_DIR    = latex
EXAMPLES_DIR = examples
TESTS_DIR    = tests
BUILD_DIR    = build

# Source files
# CORE_SRCS = $(CORE_DIR)/complex.c
CORE_SRCS    = $(CORE_DIR)/vector.c \
               $(CORE_DIR)/matrix.c \
               $(CORE_DIR)/utils.c \
               $(CORE_DIR)/sparse.c \
               $(CORE_DIR)/fft/fft.c \
               $(CORE_DIR)/fft/fft_wrapper.c \
               $(CORE_DIR)/linalg/eigen_generic.c \
               $(CORE_DIR)/linalg/tridiag_eigen.c \
               $(CORE_DIR)/linalg/tridiag_eigh.c \
               $(CORE_DIR)/linalg/qr.c \
               $(CORE_DIR)/linalg/lu.c \
               $(CORE_DIR)/linalg/svd.c \
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
               $(PHYSICS_DIR)/hydrogen.c \
               $(PHYSICS_DIR)/perturbation.c \
               $(PHYSICS_DIR)/variational.c \
               $(PHYSICS_DIR)/wkb.c \
               $(PHYSICS_DIR)/scattering.c \
               $(PHYSICS_DIR)/identical.c \
               $(PHYSICS_DIR)/relativistic.c

LATEX_SRCS   = $(LATEX_DIR)/latex_gen.c

PLOT_SRCS    = $(PLOT_SRC)

# Object files
# CORE_OBJS    = $(patsubst %.c,$(BUILD_DIR)/%.o,$(CORE_SRCS))
# PHYSICS_OBJS = $(patsubst %.c,$(BUILD_DIR)/%.o,$(PHYSICS_SRCS))
# LATEX_OBJS   = $(patsubst %.c,$(BUILD_DIR)/%.o,$(LATEX_SRCS))
# PLOT_OBJS    = $(patsubst %.c,$(BUILD_DIR)/%.o,$(PLOT_SRC))
#
# ALL_OBJS  = $(CORE_OBJS) $(PHYSICS_OBJS) $(LATEX_OBJS) $(PLOT_OBJS)
ALL_SRCS    = $(CORE_SRCS) $(PHYSICS_SRCS) $(LATEX_SRCS) $(PLOT_SRC)
ALL_OBJS    = $(patsubst %.c,$(BUILD_DIR)/%.o,$(ALL_SRCS))

# Targets
EXAMPLES    = $(BUILD_DIR)/eg_01_particle_box \
              $(BUILD_DIR)/eg_02_harmonic \
              $(BUILD_DIR)/eg_03_hydrogen \
              $(BUILD_DIR)/eg_04_perturbation \
              $(BUILD_DIR)/eg_05_tunnelling \
              $(BUILD_DIR)/eg_06_finite_well \
              $(BUILD_DIR)/eg_07_infinite_well

TESTS       = $(BUILD_DIR)/test_complex \
              $(BUILD_DIR)/test_matrix \
              $(BUILD_DIR)/test_numerov \
              $(BUILD_DIR)/test_rk4 \
              $(BUILD_DIR)/test_fft \
              $(BUILD_DIR)/test_hydrogen \
              $(BUILD_DIR)/test_perturbation \
              $(BUILD_DIR)/test_crank_nicolson \
              $(BUILD_DIR)/test_wkb \
              $(BUILD_DIR)/test_potentials

# ifeq ($(PLOT_BACKEND),GR)
#     TESTS += $(BUILD_DIR)/test_grplot
# endif

.PHONY: all clean examples tests test run-examples info

all: directories $(OUTPUT_DIR) $(EXAMPLES) $(TESTS)

directories:
	@mkdir -p $(BUILD_DIR)
	@mkdir -p $(BUILD_DIR)/$(CORE_DIR)/linalg
	@mkdir -p $(BUILD_DIR)/$(CORE_DIR)/ode
	@mkdir -p $(BUILD_DIR)/$(CORE_DIR)/special
	@mkdir -p $(BUILD_DIR)/$(CORE_DIR)/fft
	@mkdir -p $(BUILD_DIR)/$(PHYSICS_DIR)
	@mkdir -p $(BUILD_DIR)/$(EXPORT_DIR)/gr
	@mkdir -p $(BUILD_DIR)/$(EXPORT_DIR)/gnuplot
	@mkdir -p $(BUILD_DIR)/$(LATEX_DIR)

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

test: tests
	@echo "Running all tests..."
	@failed=0; \
	for t in $(TESTS); do \
		name=$$(basename $$t); \
		printf "Running $$name... "; \
		if $$t > /tmp/$$name.out 2>&1; then \
			printf "\033[32mPASS\033[0m\n"; \
		else \
			printf "\033[31mFAIL\033[0m\n"; \
			cat /tmp/$$name.out; \
			failed=$$((failed+1)); \
		fi; \
	done; \
	echo ""; \
	echo "Results: $$passed passed, $$failed failed"; \
	[ $$failed -eq 0 ]
	# if [ $$failed -eq 0 ]; then \
	# 	echo "All tests passed."; \
	# else \
	# 	echo "$$failed test(s) failed."; \
	# 	exit 1; \
	# fi

# Run a specific test
test-%: $(BUILD_DIR)/test_%
	@$<

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

clean:
	rm -rf $(BUILD_DIR) $(OUTPUT_DIR)
