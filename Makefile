CC = gcc
CFLAGS = -Wall -Wextra -O2 -fPIC
CFLAGS += -fsanitize=address -g
LDFLAGS = -lm

# Directories
CORE_DIR = core
PHYSICS_DIR = physics
EXPORT_DIR = export
LATEX_DIR = latex
EXAMPLES_DIR = examples
TESTS_DIR = tests
BUILD_DIR = build

# Source files
# CORE_SRCS = $(CORE_DIR)/complex.c
CORE_SRCS = $(CORE_DIR)/vector.c \
            $(CORE_DIR)/matrix.c \
            $(CORE_DIR)/utils.c \
            $(CORE_DIR)/sparse.c \
            $(CORE_DIR)/fft/fft.c \
            $(CORE_DIR)/fft/fft_wrapper.c \
            $(CORE_DIR)/linalg/eigen_generic.c \
            $(CORE_DIR)/linalg/tridiag_eigen.c \
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

EXPORT_SRCS = $(EXPORT_DIR)/gnuplot_pipe.c \
            $(EXPORT_DIR)/matplotlib_pipe.c \
            $(EXPORT_DIR)/svg_writer.c \
            $(EXPORT_DIR)/terminal_plot.c \
            $(EXPORT_DIR)/csv_writer.c
            # TODO:
            # $(EXPORT_DIR)/hdf5_writer.c \
            # $(EXPORT_DIR)/netcdf_writer.c

LATEX_SRCS = $(LATEX_DIR)/latex_gen.c

# Object files
CORE_OBJS   = $(patsubst %.c,$(BUILD_DIR)/%.o,$(CORE_SRCS))
PHYSICS_OBJS = $(patsubst %.c,$(BUILD_DIR)/%.o,$(PHYSICS_SRCS))
EXPORT_OBJS   = $(patsubst %.c,$(BUILD_DIR)/%.o,$(EXPORT_SRCS))
LATEX_OBJS  = $(patsubst %.c,$(BUILD_DIR)/%.o,$(LATEX_SRCS))

ALL_OBJS = $(CORE_OBJS) $(PHYSICS_OBJS) $(EXPORT_OBJS) $(LATEX_OBJS)

# Targets
EXAMPLES = $(BUILD_DIR)/eg_01_particle_box \
           $(BUILD_DIR)/eg_02_harmonic \
           $(BUILD_DIR)/eg_03_hydrogen \
           $(BUILD_DIR)/eg_04_perturbation

TESTS = $(BUILD_DIR)/test_complex \
        $(BUILD_DIR)/test_matrix \
        $(BUILD_DIR)/test_numerov \
        $(BUILD_DIR)/test_fft \
        $(BUILD_DIR)/test_hydrogen \
        $(BUILD_DIR)/test_perturbation

.PHONY: all clean examples tests test run-examples

all: directories $(EXAMPLES) $(TESTS)

directories:
	@mkdir -p $(BUILD_DIR)
	@mkdir -p $(BUILD_DIR)/$(CORE_DIR)/linalg
	@mkdir -p $(BUILD_DIR)/$(CORE_DIR)/ode
	@mkdir -p $(BUILD_DIR)/$(CORE_DIR)/special
	@mkdir -p $(BUILD_DIR)/$(CORE_DIR)/fft
	@mkdir -p $(BUILD_DIR)/$(PHYSICS_DIR)
	@mkdir -p $(BUILD_DIR)/$(EXPORT_DIR)
	@mkdir -p $(BUILD_DIR)/$(LATEX_DIR)

$(BUILD_DIR)/%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Examples
$(BUILD_DIR)/eg_01_particle_box: $(EXAMPLES_DIR)/eg_01_particle_box.c $(ALL_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

$(BUILD_DIR)/eg_02_harmonic: $(EXAMPLES_DIR)/eg_02_harmonic.c $(ALL_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

$(BUILD_DIR)/eg_03_hydrogen: $(EXAMPLES_DIR)/eg_03_hydrogen.c $(ALL_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

$(BUILD_DIR)/eg_04_perturbation: $(EXAMPLES_DIR)/eg_04_perturbation.c $(ALL_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

# Tests
$(BUILD_DIR)/test_complex: $(TESTS_DIR)/test_complex.c $(ALL_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

$(BUILD_DIR)/test_matrix: $(TESTS_DIR)/test_matrix.c $(ALL_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

$(BUILD_DIR)/test_numerov: $(TESTS_DIR)/test_numerov.c $(ALL_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

$(BUILD_DIR)/test_fft: $(TESTS_DIR)/test_fft.c $(ALL_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

$(BUILD_DIR)/test_hydrogen: $(TESTS_DIR)/test_hydrogen.c $(ALL_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

$(BUILD_DIR)/test_perturbation: $(TESTS_DIR)/test_perturbation.c $(ALL_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

examples: $(EXAMPLES)

tests: $(TESTS)

# Run all tests
test: tests
	@echo "Running all tests..."
	@for t in $(TESTS); do \
		echo "Running $$t"; \
		$$t; \
	done

# run a specific test (e.g., `make test-complex`)
test-%: $(BUILD_DIR)/test_%
	@echo "Running $<..."
	@$<

run-examples: examples
	@echo "Running particle box example..."
	@$(BUILD_DIR)/eg_01_particle_box
	@echo ""
	@echo "Running harmonic oscillator example..."
	@$(BUILD_DIR)/eg_02_harmonic

clean:
	rm -rf $(BUILD_DIR)
