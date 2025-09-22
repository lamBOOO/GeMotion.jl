# GeMotion.jl Unit Tests

This directory contains comprehensive unit tests for the GeMotion.jl package, focusing on fluid dynamics simulations using finite element methods.

## Test Structure

### Core Test Files

1. **`runtests.jl`** - Main test runner
   - Orchestrates all test execution
   - Handles graceful degradation when full module cannot be loaded
   - Includes comprehensive integration tests when GeMotion is available

2. **`test_unit_solver.jl`** - Solver component tests (28 tests)
   - Parameter validation for physical quantities (Pr, Ra, n)
   - Gridap finite element integration
   - Mathematical utility function validation
   - Boundary condition setup validation

3. **`test_unit_postprocessing.jl`** - Postprocessing tests (26 tests)
   - Color map and visualization parameter validation
   - Coordinate transformation testing
   - File I/O validation
   - Statistical and mathematical functions for analysis

4. **`test_integration.jl`** - Integration tests (25 tests)
   - End-to-end simulation setup without full execution
   - Finite element space construction
   - Weak form component validation
   - Nonlinear solver configuration testing

## Test Categories

### Mathematical Foundation Tests
- Basic mathematical operations (abs, sqrt, exp, log)
- Linear algebra operations (vectors, matrices, dot products, norms)
- Coordinate transformations (polar to Cartesian)
- Interpolation mathematics (linear, bilinear)

### Finite Element Method Tests
- Gridap model creation and validation
- Function space construction (scalar, vector, mixed)
- Triangulation and boundary triangulation setup
- Test and trial function space creation
- Integration measure validation
- Bilinear form assembly testing

### Physical Parameter Validation
- Prandtl number (Pr) range and type validation
- Rayleigh number (Ra) range and type validation  
- Power law exponent (n) validation
- Boundary condition consistency checking

### Solver Configuration Tests
- Newton method parameter validation
- Line search algorithm setup
- Initial guess strategy validation
- Convergence criteria validation

### Postprocessing Tests
- Color map creation and validation
- Plot parameter range checking
- File output validation (PDF, VTU, CSV)
- Level set computation for contour plots

### Error Handling Tests
- Mismatched boundary condition array lengths
- Invalid parameter ranges
- Type validation for inputs

## Running Tests

### Run all tests:
```bash
julia --project=. -e "using Pkg; Pkg.test()"
```

### Run individual test files:
```bash
julia --project=. test/test_unit_solver.jl
julia --project=. test/test_unit_postprocessing.jl  
julia --project=. test/test_integration.jl
```

### Quick validation test:
```bash
julia --project=. -e "using Test; include(\"test/test_unit_solver.jl\")"
```

## Test Design Philosophy

1. **Robust Testing**: Tests are designed to pass even when full dependencies are not available
2. **Mathematical Validation**: Focus on validating mathematical correctness and physical consistency
3. **Component Isolation**: Tests can be run independently to isolate issues
4. **Comprehensive Coverage**: Tests cover input validation, mathematical operations, and output structure
5. **Performance Awareness**: Tests use minimal problem sizes to ensure fast execution

## Test Coverage Summary

- **Total Tests**: 79 individual test cases
- **Solver Components**: 28 tests
- **Postprocessing**: 26 tests  
- **Integration**: 25 tests
- **Coverage Areas**: Parameter validation, mathematical operations, FE method setup, solver configuration, postprocessing, error handling

The test suite provides comprehensive validation of the GeMotion.jl package while maintaining fast execution times and robust error handling.