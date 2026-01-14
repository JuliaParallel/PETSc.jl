# TS (Time Stepping) - Low-level Interface

The TS (Time Stepping) component provides methods for solving time-dependent differential equations, including ordinary differential equations (ODEs), differential-algebraic equations (DAEs), and time-dependent partial differential equations (PDEs).

## Overview

TS supports various time integration methods including:
- **Explicit methods**: Forward Euler, Runge-Kutta methods (RK2, RK3, RK4, etc.)
- **Implicit methods**: Backward Euler, Crank-Nicolson, BDF methods
- **IMEX methods**: Implicit-Explicit combinations for stiff-nonstiff problems
- **Adaptive methods**: Adaptive time stepping with error control
- **Specialized methods**: Theta methods, Rosenbrock methods, SSP methods

The TS object manages the time integration loop, adaptive time stepping, event detection, and provides interfaces to various time stepping schemes.

## Basic Usage

```julia
using PETSc

# Initialize PETSc
petsclib = PETSc.getlib()

# Create a TS object
ts = LibPETSc.TSCreate(petsclib, LibPETSc.PETSC_COMM_SELF)

# Set the problem type (ODE or DAE)
LibPETSc.TSSetProblemType(petsclib, ts, LibPETSc.TS_NONLINEAR)

# Set the time stepping method (e.g., BDF, RK, Theta)
LibPETSc.TSSetType(petsclib, ts, "bdf")  # String convenience wrapper

# Set time span
LibPETSc.TSSetTime(petsclib, ts, 0.0)  # Initial time
LibPETSc.TSSetMaxTime(petsclib, ts, 1.0)  # Final time
LibPETSc.TSSetExactFinalTime(petsclib, ts, LibPETSc.TS_EXACTFINALTIME_STEPOVER)

# Set initial time step
LibPETSc.TSSetTimeStep(petsclib, ts, 0.01)

# Set the right-hand side function (for ODE: du/dt = f(t,u))
# LibPETSc.TSSetRHSFunction(petsclib, ts, C_NULL, rhs_function_ptr, C_NULL)

# Set options from command line/options database
LibPETSc.TSSetFromOptions(petsclib, ts)

# Set initial condition
# LibPETSc.TSSetSolution(petsclib, ts, initial_vec)

# Solve
# LibPETSc.TSSolve(petsclib, ts, solution_vec)

# Get solution time (returns value directly)
final_time = LibPETSc.TSGetSolveTime(petsclib, ts)

# Cleanup
LibPETSc.TSDestroy(petsclib, ts)
```

## Common Workflow

1. **Create and configure TS**: Use `TSCreate`, `TSSetType`, `TSSetProblemType`
2. **Set time parameters**: `TSSetTime`, `TSSetMaxTime`, `TSSetTimeStep`, `TSSetMaxSteps`
3. **Define equations**: `TSSetRHSFunction`, `TSSetIFunction`, `TSSetIJacobian`
4. **Configure adaptivity**: `TSAdaptSetType`, `TSSetTolerances`
5. **Set initial condition**: `TSSetSolution`
6. **Solve**: `TSSolve`
7. **Retrieve solution**: `TSGetSolution`, `TSGetTime`, `TSGetStepNumber`

## Time Integration Schemes

Available through `TSSetType`:
- `TSEULER`: Forward/Backward Euler
- `TSRK`: Runge-Kutta methods (various orders)
- `TSBDF`: Backward Differentiation Formulas
- `TSTHETA`: Theta method (generalizes Euler, Crank-Nicolson)
- `TSROSW`: Rosenbrock-W methods for stiff problems
- `TSARKIMEX`: Additive Runge-Kutta IMEX methods
- `TSGLEE`: Generalized Linear Evolution Equations
- `TSALPHA`: Alpha methods for second-order systems
- `TSSSP`: Strong Stability Preserving methods

## Event Detection

TS supports event detection (zero-crossing of event functions) during time integration:
- `TSSetEventHandler`: Configure events to detect
- Events can trigger actions like termination, adaptation, or post-event processing

## Adaptive Time Stepping

TS includes sophisticated adaptive time stepping:
- `TSAdaptChoose`: Select time step based on error estimates
- Various adapt types: `TSADAPTBASIC`, `TSADAPTNONE`, `TSADAPTDSP`
- Tolerance control: `TSSetTolerances` for absolute/relative error

## Function Reference

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/TS_wrappers.jl"]
Order   = [:function]
```

## TS Add-ons

Additional TS utilities and helper functions:

```@autodocs
Modules = [PETSc.LibPETSc]
Pages   = ["autowrapped/TSaddons_wrappers.jl"]
Order   = [:function]
```
