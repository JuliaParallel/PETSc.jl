# autodefined type arguments for class ------
mutable struct _n_TSTrajectory end
const TSTrajectory = Ptr{_n_TSTrajectory}

mutable struct _n_TSAdapt end
const TSAdapt = Ptr{_n_TSAdapt}

mutable struct _n_TSMonitorVTKCtx end
const TSMonitorVTKCtx = Ptr{_n_TSMonitorVTKCtx}

mutable struct TSRHSJacobianPFn end

mutable struct _n_TSGLLEAdapt end
const TSGLLEAdapt = Ptr{_n_TSGLLEAdapt}

mutable struct TSGLLEAcceptFn end

mutable struct TSAlpha2PredictorFn end

#mutable struct TSRHSFunctionFn end
#mutable struct TSSolutionFn end
#mutable struct TSForcingFn end
#mutable struct TSRHSJacobianFn end
#mutable struct TSIFunctionFn end
#mutable struct TSIJacobianFn end
#mutable struct TSI2FunctionFn end
#mutable struct TSI2JacobianFn end
#mutable struct TSTransientVariableFn end
# -------------------------------------------------------
"""
	TSSetFromOptions(petsclib::PetscLibType,ts::TS) 
Sets various `TS` parameters from the options database

Collective

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Options Database Keys:
- `-ts_type <type>`                                                    - EULER, BEULER, SUNDIALS, PSEUDO, CN, RK, THETA, ALPHA, GLLE,  SSP, GLEE, BSYMP, IRK, see `TSType`
- `-ts_save_trajectory`                                                - checkpoint the solution at each time-step
- `-ts_max_time <time>`                                                - maximum time to compute to
- `-ts_time_span <t0,...tf>`                                           - sets the time span, solutions are computed and stored for each indicated time, init_time and max_time are set
- `-ts_eval_times <t0,...tn>`                                          - time points where solutions are computed and stored for each indicated time
- `-ts_max_steps <steps>`                                              - maximum time-step number to execute until (possibly with nonzero starting value)
- `-ts_run_steps <steps>`                                              - maximum number of time steps for TSSolve to take on each call
- `-ts_init_time <time>`                                               - initial time to start computation
- `-ts_final_time <time>`                                              - final time to compute to (deprecated: use `-ts_max_time`)
- `-ts_dt <dt>`                                                        - initial time step
- `-ts_exact_final_time <stepover,interpolate,matchstep>`              - whether to stop at the exact given final time and how to compute the solution at that time
- `-ts_max_snes_failures <maxfailures>`                                - Maximum number of nonlinear solve failures allowed
- `-ts_max_reject <maxrejects>`                                        - Maximum number of step rejections before step fails
- `-ts_error_if_step_fails <true,false>`                               - Error if no step succeeds
- `-ts_rtol <rtol>`                                                    - relative tolerance for local truncation error
- `-ts_atol <atol>`                                                    - Absolute tolerance for local truncation error
- `-ts_rhs_jacobian_test_mult -mat_shell_test_mult_view`               - test the Jacobian at each iteration against finite difference with RHS function
- `-ts_rhs_jacobian_test_mult_transpose`                               - test the Jacobian at each iteration against finite difference with RHS function
- `-ts_adjoint_solve <yes,no>`                                         - After solving the ODE/DAE solve the adjoint problem (requires `-ts_save_trajectory`)
- `-ts_fd_color`                                                       - Use finite differences with coloring to compute IJacobian
- `-ts_monitor`                                                        - print information at each timestep
- `-ts_monitor_cancel`                                                 - Cancel all monitors
- `-ts_monitor_wall_clock_time`                                        - Monitor wall-clock time, KSP iterations, and SNES iterations per step
- `-ts_monitor_lg_solution`                                            - Monitor solution graphically
- `-ts_monitor_lg_error`                                               - Monitor error graphically
- `-ts_monitor_error`                                                  - Monitors norm of error
- `-ts_monitor_lg_timestep`                                            - Monitor timestep size graphically
- `-ts_monitor_lg_timestep_log`                                        - Monitor log timestep size graphically
- `-ts_monitor_lg_snes_iterations`                                     - Monitor number nonlinear iterations for each timestep graphically
- `-ts_monitor_lg_ksp_iterations`                                      - Monitor number nonlinear iterations for each timestep graphically
- `-ts_monitor_sp_eig`                                                 - Monitor eigenvalues of linearized operator graphically
- `-ts_monitor_draw_solution`                                          - Monitor solution graphically
- `-ts_monitor_draw_solution_phase  <xleft,yleft,xright,yright>`       - Monitor solution graphically with phase diagram, requires problem with exactly 2 degrees of freedom
- `-ts_monitor_draw_error`                                             - Monitor error graphically, requires use to have provided TSSetSolutionFunction()
- `-ts_monitor_solution [ascii binary draw][:filename][:viewerformat]` - monitors the solution at each timestep
- `-ts_monitor_solution_interval <interval>`                           - output once every interval (default=1) time steps. Use -1 to only output at the end of the simulation
- `-ts_monitor_solution_skip_initial`                                  - skip writing of initial condition
- `-ts_monitor_solution_vtk <filename.vts,filename.vtu>`               - Save each time step to a binary file, use filename-%%03" PetscInt_FMT ".vts (filename-%%03" PetscInt_FMT ".vtu)
- `-ts_monitor_solution_vtk_interval <interval>`                       - output once every interval (default=1) time steps. Use -1 to only output at the end of the simulation
- `-ts_monitor_envelope`                                               - determine maximum and minimum value of each component of the solution over the solution time

Level: beginner

-seealso: [](ch_ts), `TS`, `TSGetType()`

# External Links
$(_doc_external("Ts/TSSetFromOptions"))
"""
function TSSetFromOptions(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSSetFromOptions(petsclib::$UnionPetscLib, ts::TS )

    @chk ccall(
               (:TSSetFromOptions, $petsc_library),
               PetscErrorCode,
               (CTS,),
               ts,
              )


	return nothing
end 

"""
	TSGetTrajectory(petsclib::PetscLibType,ts::TS, tr::TSTrajectory) 
Gets the trajectory from a `TS` if it exists

Collective

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Output Parameter:
- `tr` - the `TSTrajectory` object, if it exists

Level: advanced

-seealso: [](ch_ts), `TS`, `TSTrajectory`, `TSAdjointSolve()`, `TSTrajectoryCreate()`

# External Links
$(_doc_external("Ts/TSGetTrajectory"))
"""
function TSGetTrajectory(petsclib::PetscLibType, ts::TS, tr::TSTrajectory) end

@for_petsc function TSGetTrajectory(petsclib::$UnionPetscLib, ts::TS, tr::TSTrajectory )

    @chk ccall(
               (:TSGetTrajectory, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{TSTrajectory}),
               ts, tr,
              )


	return nothing
end 

"""
	TSSetSaveTrajectory(petsclib::PetscLibType,ts::TS) 
Causes the `TS` to save its solutions as it iterates forward in time in a `TSTrajectory` object

Collective

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Options Database Keys:
- `-ts_save_trajectory`      - saves the trajectory to a file
- `-ts_trajectory_type type` - set trajectory type

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSTrajectory`, `TSGetTrajectory()`, `TSAdjointSolve()`

# External Links
$(_doc_external("Ts/TSSetSaveTrajectory"))
"""
function TSSetSaveTrajectory(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSSetSaveTrajectory(petsclib::$UnionPetscLib, ts::TS )

    @chk ccall(
               (:TSSetSaveTrajectory, $petsc_library),
               PetscErrorCode,
               (CTS,),
               ts,
              )


	return nothing
end 

"""
	TSResetTrajectory(petsclib::PetscLibType,ts::TS) 
Destroys and recreates the internal `TSTrajectory` object

Collective

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Level: intermediate

-seealso: [](ch_ts), `TSTrajectory`, `TSGetTrajectory()`, `TSAdjointSolve()`, `TSRemoveTrajectory()`

# External Links
$(_doc_external("Ts/TSResetTrajectory"))
"""
function TSResetTrajectory(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSResetTrajectory(petsclib::$UnionPetscLib, ts::TS )

    @chk ccall(
               (:TSResetTrajectory, $petsc_library),
               PetscErrorCode,
               (CTS,),
               ts,
              )


	return nothing
end 

"""
	TSRemoveTrajectory(petsclib::PetscLibType,ts::TS) 
Destroys and removes the internal `TSTrajectory` object from a `TS`

Collective

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Level: intermediate

-seealso: [](ch_ts), `TSTrajectory`, `TSResetTrajectory()`, `TSAdjointSolve()`

# External Links
$(_doc_external("Ts/TSRemoveTrajectory"))
"""
function TSRemoveTrajectory(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSRemoveTrajectory(petsclib::$UnionPetscLib, ts::TS )

    @chk ccall(
               (:TSRemoveTrajectory, $petsc_library),
               PetscErrorCode,
               (CTS,),
               ts,
              )


	return nothing
end 

"""
	TSComputeRHSJacobian(petsclib::PetscLibType,ts::TS, t::PetscReal, U::PetscVec, A::PetscMat, B::PetscMat) 
Computes the Jacobian matrix that has been
set with `TSSetRHSJacobian()`.

Collective

Input Parameters:
- `ts` - the `TS` context
- `t`  - current timestep
- `U`  - input vector

Output Parameters:
- `A` - Jacobian matrix
- `B` - optional matrix used to compute the preconditioner, often the same as `A`

Level: developer

-seealso: [](ch_ts), `TS`, `TSSetRHSJacobian()`, `KSPSetOperators()`

# External Links
$(_doc_external("Ts/TSComputeRHSJacobian"))
"""
function TSComputeRHSJacobian(petsclib::PetscLibType, ts::TS, t::PetscReal, U::PetscVec, A::PetscMat, B::PetscMat) end

@for_petsc function TSComputeRHSJacobian(petsclib::$UnionPetscLib, ts::TS, t::$PetscReal, U::PetscVec, A::PetscMat, B::PetscMat )

    @chk ccall(
               (:TSComputeRHSJacobian, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal, CVec, CMat, CMat),
               ts, t, U, A, B,
              )


	return nothing
end 

"""
	TSComputeRHSFunction(petsclib::PetscLibType,ts::TS, t::PetscReal, U::PetscVec, y::PetscVec) 
Evaluates the right

Collective

Input Parameters:
- `ts` - the `TS` context
- `t`  - current time
- `U`  - state vector

Output Parameter:
- `y` - right-hand side

Level: developer

-seealso: [](ch_ts), `TS`, `TSSetRHSFunction()`, `TSComputeIFunction()`

# External Links
$(_doc_external("Ts/TSComputeRHSFunction"))
"""
function TSComputeRHSFunction(petsclib::PetscLibType, ts::TS, t::PetscReal, U::PetscVec, y::PetscVec) end

@for_petsc function TSComputeRHSFunction(petsclib::$UnionPetscLib, ts::TS, t::$PetscReal, U::PetscVec, y::PetscVec )

    @chk ccall(
               (:TSComputeRHSFunction, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal, CVec, CVec),
               ts, t, U, y,
              )


	return nothing
end 

"""
	TSComputeSolutionFunction(petsclib::PetscLibType,ts::TS, t::PetscReal, U::PetscVec) 
Evaluates the solution function.

Collective

Input Parameters:
- `ts` - the `TS` context
- `t`  - current time

Output Parameter:
- `U` - the solution

Level: developer

-seealso: [](ch_ts), `TS`, `TSSetSolutionFunction()`, `TSSetRHSFunction()`, `TSComputeIFunction()`

# External Links
$(_doc_external("Ts/TSComputeSolutionFunction"))
"""
function TSComputeSolutionFunction(petsclib::PetscLibType, ts::TS, t::PetscReal, U::PetscVec) end

@for_petsc function TSComputeSolutionFunction(petsclib::$UnionPetscLib, ts::TS, t::$PetscReal, U::PetscVec )

    @chk ccall(
               (:TSComputeSolutionFunction, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal, CVec),
               ts, t, U,
              )


	return nothing
end 

"""
	TSComputeForcingFunction(petsclib::PetscLibType,ts::TS, t::PetscReal, U::PetscVec) 
Evaluates the forcing function.

Collective

Input Parameters:
- `ts` - the `TS` context
- `t`  - current time

Output Parameter:
- `U` - the function value

Level: developer

-seealso: [](ch_ts), `TS`, `TSSetSolutionFunction()`, `TSSetRHSFunction()`, `TSComputeIFunction()`

# External Links
$(_doc_external("Ts/TSComputeForcingFunction"))
"""
function TSComputeForcingFunction(petsclib::PetscLibType, ts::TS, t::PetscReal, U::PetscVec) end

@for_petsc function TSComputeForcingFunction(petsclib::$UnionPetscLib, ts::TS, t::$PetscReal, U::PetscVec )

    @chk ccall(
               (:TSComputeForcingFunction, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal, CVec),
               ts, t, U,
              )


	return nothing
end 

"""
	TSComputeIFunction(petsclib::PetscLibType,ts::TS, t::PetscReal, U::PetscVec, Udot::PetscVec, Y::PetscVec, imex::PetscBool) 
Evaluates the DAE residual written in the implicit form F(t,U,Udot)=0

Collective

Input Parameters:
- `ts`   - the `TS` context
- `t`    - current time
- `U`    - state vector
- `Udot` - time derivative of state vector
- `imex` - flag indicates if the method is `TSARKIMEX` so that the RHSFunction should be kept separate

Output Parameter:
- `Y` - right-hand side

Level: developer

-seealso: [](ch_ts), `TS`, `TSSetIFunction()`, `TSComputeRHSFunction()`

# External Links
$(_doc_external("Ts/TSComputeIFunction"))
"""
function TSComputeIFunction(petsclib::PetscLibType, ts::TS, t::PetscReal, U::PetscVec, Udot::PetscVec, Y::PetscVec, imex::PetscBool) end

@for_petsc function TSComputeIFunction(petsclib::$UnionPetscLib, ts::TS, t::$PetscReal, U::PetscVec, Udot::PetscVec, Y::PetscVec, imex::PetscBool )

    @chk ccall(
               (:TSComputeIFunction, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal, CVec, CVec, CVec, PetscBool),
               ts, t, U, Udot, Y, imex,
              )


	return nothing
end 

"""
	TSComputeIJacobian(petsclib::PetscLibType,ts::TS, t::PetscReal, U::PetscVec, Udot::PetscVec, shift::PetscReal, A::PetscMat, B::PetscMat, imex::PetscBool) 
Evaluates the Jacobian of the DAE

Collective

Input Parameters:
- `ts`    - the `TS` context
- `t`     - current timestep
- `U`     - state vector
- `Udot`  - time derivative of state vector
- `shift` - shift to apply, see note below
- `imex`  - flag indicates if the method is `TSARKIMEX` so that the RHSJacobian should be kept separate

Output Parameters:
- `A` - Jacobian matrix
- `B` - matrix from which the preconditioner is constructed; often the same as `A`

Level: developer

-seealso: [](ch_ts), `TS`, `TSSetIJacobian()`

# External Links
$(_doc_external("Ts/TSComputeIJacobian"))
"""
function TSComputeIJacobian(petsclib::PetscLibType, ts::TS, t::PetscReal, U::PetscVec, Udot::PetscVec, shift::PetscReal, A::PetscMat, B::PetscMat, imex::PetscBool) end

@for_petsc function TSComputeIJacobian(petsclib::$UnionPetscLib, ts::TS, t::$PetscReal, U::PetscVec, Udot::PetscVec, shift::$PetscReal, A::PetscMat, B::PetscMat, imex::PetscBool )

    @chk ccall(
               (:TSComputeIJacobian, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal, CVec, CVec, $PetscReal, CMat, CMat, PetscBool),
               ts, t, U, Udot, shift, A, B, imex,
              )


	return nothing
end 

"""
	TSSetRHSFunction(petsclib::PetscLibType,ts::TS, r::PetscVec, f::TSRHSFunctionFn, ctx::Cvoid) 
Sets the routine for evaluating the function,
where U_t = G(t,u).

Logically Collective

Input Parameters:
- `ts`  - the `TS` context obtained from `TSCreate()`
- `r`   - vector to put the computed right-hand side (or `NULL` to have it created)
- `f`   - routine for evaluating the right-hand-side function
- `ctx` - [optional] user-defined context for private data for the function evaluation routine (may be `NULL`)

Level: beginner

-seealso: [](ch_ts), `TS`, `TSRHSFunctionFn`, `TSSetRHSJacobian()`, `TSSetIJacobian()`, `TSSetIFunction()`

# External Links
$(_doc_external("Ts/TSSetRHSFunction"))
"""
function TSSetRHSFunction(petsclib::PetscLibType, ts::TS, r::PetscVec, f::TSRHSFunctionFn, ctx::Cvoid) end

@for_petsc function TSSetRHSFunction(petsclib::$UnionPetscLib, ts::TS, r::PetscVec, f::TSRHSFunctionFn, ctx::Cvoid )

    @chk ccall(
               (:TSSetRHSFunction, $petsc_library),
               PetscErrorCode,
               (CTS, CVec, Ptr{TSRHSFunctionFn}, Ptr{Cvoid}),
               ts, r, f, ctx,
              )


	return nothing
end 

"""
	TSSetSolutionFunction(petsclib::PetscLibType,ts::TS, f::TSSolutionFn, ctx::Cvoid) 
Provide a function that computes the solution of the ODE or DAE

Logically Collective

Input Parameters:
- `ts`  - the `TS` context obtained from `TSCreate()`
- `f`   - routine for evaluating the solution
- `ctx` - [optional] user-defined context for private data for the
function evaluation routine (may be `NULL`)

Options Database Keys:
- `-ts_monitor_lg_error`   - create a graphical monitor of error history, requires user to have provided `TSSetSolutionFunction()`
- `-ts_monitor_draw_error` - Monitor error graphically, requires user to have provided `TSSetSolutionFunction()`

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSSolutionFn`, `TSSetRHSJacobian()`, `TSSetIJacobian()`, `TSComputeSolutionFunction()`, `TSSetForcingFunction()`, `TSSetSolution()`, `TSGetSolution()`, `TSMonitorLGError()`, `TSMonitorDrawError()`

# External Links
$(_doc_external("Ts/TSSetSolutionFunction"))
"""
function TSSetSolutionFunction(petsclib::PetscLibType, ts::TS, f::TSSolutionFn, ctx::Cvoid) end

@for_petsc function TSSetSolutionFunction(petsclib::$UnionPetscLib, ts::TS, f::TSSolutionFn, ctx::Cvoid )

    @chk ccall(
               (:TSSetSolutionFunction, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{TSSolutionFn}, Ptr{Cvoid}),
               ts, f, ctx,
              )


	return nothing
end 

"""
	TSSetForcingFunction(petsclib::PetscLibType,ts::TS, func::TSForcingFn, ctx::Cvoid) 
Provide a function that computes a forcing term for a ODE or PDE

Logically Collective

Input Parameters:
- `ts`   - the `TS` context obtained from `TSCreate()`
- `func` - routine for evaluating the forcing function
- `ctx`  - [optional] user-defined context for private data for the function evaluation routine
(may be `NULL`)

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSForcingFn`, `TSSetRHSJacobian()`, `TSSetIJacobian()`,
`TSComputeSolutionFunction()`, `TSSetSolutionFunction()`

# External Links
$(_doc_external("Ts/TSSetForcingFunction"))
"""
function TSSetForcingFunction(petsclib::PetscLibType, ts::TS, func::TSForcingFn, ctx::Cvoid) end

@for_petsc function TSSetForcingFunction(petsclib::$UnionPetscLib, ts::TS, func::TSForcingFn, ctx::Cvoid )

    @chk ccall(
               (:TSSetForcingFunction, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{TSForcingFn}, Ptr{Cvoid}),
               ts, func, ctx,
              )


	return nothing
end 

"""
	TSSetRHSJacobian(petsclib::PetscLibType,ts::TS, Amat::PetscMat, Pmat::PetscMat, f::TSRHSJacobianFn, ctx::Cvoid) 
Sets the function to compute the Jacobian of G,
where U_t = G(U,t), as well as the location to store the matrix.

Logically Collective

Input Parameters:
- `ts`   - the `TS` context obtained from `TSCreate()`
- `Amat` - (approximate) location to store Jacobian matrix entries computed by `f`
- `Pmat` - matrix from which preconditioner is to be constructed (usually the same as `Amat`)
- `f`    - the Jacobian evaluation routine
- `ctx`  - [optional] user-defined context for private data for the Jacobian evaluation routine (may be `NULL`)

Level: beginner

-seealso: [](ch_ts), `TS`, `TSRHSJacobianFn`, `SNESComputeJacobianDefaultColor()`,
`TSSetRHSFunction()`, `TSRHSJacobianSetReuse()`, `TSSetIJacobian()`, `TSRHSFunctionFn`, `TSIFunctionFn`

# External Links
$(_doc_external("Ts/TSSetRHSJacobian"))
"""
function TSSetRHSJacobian(petsclib::PetscLibType, ts::TS, Amat::PetscMat, Pmat::PetscMat, f::TSRHSJacobianFn, ctx::Cvoid) end

@for_petsc function TSSetRHSJacobian(petsclib::$UnionPetscLib, ts::TS, Amat::PetscMat, Pmat::PetscMat, f::TSRHSJacobianFn, ctx::Cvoid )

    @chk ccall(
               (:TSSetRHSJacobian, $petsc_library),
               PetscErrorCode,
               (CTS, CMat, CMat, Ptr{TSRHSJacobianFn}, Ptr{Cvoid}),
               ts, Amat, Pmat, f, ctx,
              )


	return nothing
end 

"""
	TSSetIFunction(petsclib::PetscLibType,ts::TS, r::PetscVec, f::TSIFunctionFn, ctx::Cvoid) 
Set the function to compute F(t,U,U_t) where F() = 0 is the DAE to be solved.

Logically Collective

Input Parameters:
- `ts`  - the `TS` context obtained from `TSCreate()`
- `r`   - vector to hold the residual (or `NULL` to have it created internally)
- `f`   - the function evaluation routine
- `ctx` - user-defined context for private data for the function evaluation routine (may be `NULL`)

Level: beginner

-seealso: [](ch_ts), `TS`, `TSIFunctionFn`, `TSSetRHSJacobian()`, `TSSetRHSFunction()`,
`TSSetIJacobian()`

# External Links
$(_doc_external("Ts/TSSetIFunction"))
"""
function TSSetIFunction(petsclib::PetscLibType, ts::TS, r::PetscVec, f::TSIFunctionFn, ctx::Cvoid) end

@for_petsc function TSSetIFunction(petsclib::$UnionPetscLib, ts::TS, r::PetscVec, f::TSIFunctionFn, ctx::Cvoid )

    @chk ccall(
               (:TSSetIFunction, $petsc_library),
               PetscErrorCode,
               (CTS, CVec, Ptr{TSIFunctionFn}, Ptr{Cvoid}),
               ts, r, f, ctx,
              )


	return nothing
end 

"""
	TSGetIFunction(petsclib::PetscLibType,ts::TS, r::PetscVec, func::TSIFunctionFn, ctx::Cvoid) 
Returns the vector where the implicit residual is stored and the function/context to compute it.

Not Collective

Input Parameter:
- `ts` - the `TS` context

Output Parameters:
- `r`    - vector to hold residual (or `NULL`)
- `func` - the function to compute residual (or `NULL`)
- `ctx`  - the function context (or `NULL`)

Level: advanced

-seealso: [](ch_ts), `TS`, `TSSetIFunction()`, `SNESGetFunction()`

# External Links
$(_doc_external("Ts/TSGetIFunction"))
"""
function TSGetIFunction(petsclib::PetscLibType, ts::TS, r::PetscVec, func::TSIFunctionFn, ctx::Cvoid) end

@for_petsc function TSGetIFunction(petsclib::$UnionPetscLib, ts::TS, r::PetscVec, func::TSIFunctionFn, ctx::Cvoid )
	r_ = Ref(r.ptr)

    @chk ccall(
               (:TSGetIFunction, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{CVec}, TSIFunctionFn, Cvoid),
               ts, r_, func, ctx,
              )

	r.ptr = C_NULL

	return nothing
end 

"""
	TSGetRHSFunction(petsclib::PetscLibType,ts::TS, r::PetscVec, func::TSRHSFunctionFn, ctx::Cvoid) 
Returns the vector where the right

Not Collective

Input Parameter:
- `ts` - the `TS` context

Output Parameters:
- `r`    - vector to hold computed right-hand side (or `NULL`)
- `func` - the function to compute right-hand side (or `NULL`)
- `ctx`  - the function context (or `NULL`)

Level: advanced

-seealso: [](ch_ts), `TS`, `TSSetRHSFunction()`, `SNESGetFunction()`

# External Links
$(_doc_external("Ts/TSGetRHSFunction"))
"""
function TSGetRHSFunction(petsclib::PetscLibType, ts::TS, r::PetscVec, func::TSRHSFunctionFn, ctx::Cvoid) end

@for_petsc function TSGetRHSFunction(petsclib::$UnionPetscLib, ts::TS, r::PetscVec, func::TSRHSFunctionFn, ctx::Cvoid )
	r_ = Ref(r.ptr)

    @chk ccall(
               (:TSGetRHSFunction, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{CVec}, TSRHSFunctionFn, Cvoid),
               ts, r_, func, ctx,
              )

	r.ptr = C_NULL

	return nothing
end 

"""
	TSSetIJacobian(petsclib::PetscLibType,ts::TS, Amat::PetscMat, Pmat::PetscMat, f::TSIJacobianFn, ctx::Cvoid) 
Set the function to compute the matrix dF/dU + a*dF/dU_t where F(t,U,U_t) is the function
provided with `TSSetIFunction()`.

Logically Collective

Input Parameters:
- `ts`   - the `TS` context obtained from `TSCreate()`
- `Amat` - (approximate) matrix to store Jacobian entries computed by `f`
- `Pmat` - matrix used to compute preconditioner (usually the same as `Amat`)
- `f`    - the Jacobian evaluation routine
- `ctx`  - user-defined context for private data for the Jacobian evaluation routine (may be `NULL`)

Level: beginner

-seealso: [](ch_ts), `TS`, `TSIJacobianFn`, `TSSetIFunction()`, `TSSetRHSJacobian()`,
`SNESComputeJacobianDefaultColor()`, `SNESComputeJacobianDefault()`, `TSSetRHSFunction()`

# External Links
$(_doc_external("Ts/TSSetIJacobian"))
"""
function TSSetIJacobian(petsclib::PetscLibType, ts::TS, Amat::PetscMat, Pmat::PetscMat, f::TSIJacobianFn, ctx::Cvoid) end

@for_petsc function TSSetIJacobian(petsclib::$UnionPetscLib, ts::TS, Amat::PetscMat, Pmat::PetscMat, f::TSIJacobianFn, ctx::Cvoid )

    @chk ccall(
               (:TSSetIJacobian, $petsc_library),
               PetscErrorCode,
               (CTS, CMat, CMat, Ptr{TSIJacobianFn}, Ptr{Cvoid}),
               ts, Amat, Pmat, f, ctx,
              )


	return nothing
end 

"""
	TSRHSJacobianSetReuse(petsclib::PetscLibType,ts::TS, reuse::PetscBool) 
restore the RHS Jacobian before calling the user

Logically Collective

Input Parameters:
- `ts`    - `TS` context obtained from `TSCreate()`
- `reuse` - `PETSC_TRUE` if the RHS Jacobian

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSSetRHSJacobian()`, `TSComputeRHSJacobianConstant()`

# External Links
$(_doc_external("Ts/TSRHSJacobianSetReuse"))
"""
function TSRHSJacobianSetReuse(petsclib::PetscLibType, ts::TS, reuse::PetscBool) end

@for_petsc function TSRHSJacobianSetReuse(petsclib::$UnionPetscLib, ts::TS, reuse::PetscBool )

    @chk ccall(
               (:TSRHSJacobianSetReuse, $petsc_library),
               PetscErrorCode,
               (CTS, PetscBool),
               ts, reuse,
              )


	return nothing
end 

"""
	TSSetI2Function(petsclib::PetscLibType,ts::TS, F::PetscVec, fun::TSI2FunctionFn, ctx::Cvoid) 
Set the function to compute F(t,U,U_t,U_tt) where F = 0 is the DAE to be solved.

Logically Collective

Input Parameters:
- `ts`  - the `TS` context obtained from `TSCreate()`
- `F`   - vector to hold the residual (or `NULL` to have it created internally)
- `fun` - the function evaluation routine
- `ctx` - user-defined context for private data for the function evaluation routine (may be `NULL`)

Level: beginner

-seealso: [](ch_ts), `TS`, `TSI2FunctionFn`, `TSSetI2Jacobian()`, `TSSetIFunction()`,
`TSCreate()`, `TSSetRHSFunction()`

# External Links
$(_doc_external("Ts/TSSetI2Function"))
"""
function TSSetI2Function(petsclib::PetscLibType, ts::TS, F::PetscVec, fun::TSI2FunctionFn, ctx::Cvoid) end

@for_petsc function TSSetI2Function(petsclib::$UnionPetscLib, ts::TS, F::PetscVec, fun::TSI2FunctionFn, ctx::Cvoid )

    @chk ccall(
               (:TSSetI2Function, $petsc_library),
               PetscErrorCode,
               (CTS, CVec, Ptr{TSI2FunctionFn}, Ptr{Cvoid}),
               ts, F, fun, ctx,
              )


	return nothing
end 

"""
	TSGetI2Function(petsclib::PetscLibType,ts::TS, r::PetscVec, fun::TSI2FunctionFn, ctx::Cvoid) 
Returns the vector where the implicit residual is stored and the function/context to compute it.

Not Collective

Input Parameter:
- `ts` - the `TS` context

Output Parameters:
- `r`   - vector to hold residual (or `NULL`)
- `fun` - the function to compute residual (or `NULL`)
- `ctx` - the function context (or `NULL`)

Level: advanced

-seealso: [](ch_ts), `TS`, `TSSetIFunction()`, `SNESGetFunction()`, `TSCreate()`

# External Links
$(_doc_external("Ts/TSGetI2Function"))
"""
function TSGetI2Function(petsclib::PetscLibType, ts::TS, r::PetscVec, fun::TSI2FunctionFn, ctx::Cvoid) end

@for_petsc function TSGetI2Function(petsclib::$UnionPetscLib, ts::TS, r::PetscVec, fun::TSI2FunctionFn, ctx::Cvoid )
	r_ = Ref(r.ptr)

    @chk ccall(
               (:TSGetI2Function, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{CVec}, TSI2FunctionFn, Cvoid),
               ts, r_, fun, ctx,
              )

	r.ptr = C_NULL

	return nothing
end 

"""
	TSSetI2Jacobian(petsclib::PetscLibType,ts::TS, J::PetscMat, P::PetscMat, jac::TSI2JacobianFn, ctx::Cvoid) 
Set the function to compute the matrix dF/dU + v*dF/dU_t  + a*dF/dU_tt
where F(t,U,U_t,U_tt) is the function you provided with `TSSetI2Function()`.

Logically Collective

Input Parameters:
- `ts`  - the `TS` context obtained from `TSCreate()`
- `J`   - matrix to hold the Jacobian values
- `P`   - matrix for constructing the preconditioner (may be same as `J`)
- `jac` - the Jacobian evaluation routine, see `TSI2JacobianFn` for the calling sequence
- `ctx` - user-defined context for private data for the Jacobian evaluation routine (may be `NULL`)

Level: beginner

-seealso: [](ch_ts), `TS`, `TSI2JacobianFn`, `TSSetI2Function()`, `TSGetI2Jacobian()`

# External Links
$(_doc_external("Ts/TSSetI2Jacobian"))
"""
function TSSetI2Jacobian(petsclib::PetscLibType, ts::TS, J::PetscMat, P::PetscMat, jac::TSI2JacobianFn, ctx::Cvoid) end

@for_petsc function TSSetI2Jacobian(petsclib::$UnionPetscLib, ts::TS, J::PetscMat, P::PetscMat, jac::TSI2JacobianFn, ctx::Cvoid )

    @chk ccall(
               (:TSSetI2Jacobian, $petsc_library),
               PetscErrorCode,
               (CTS, CMat, CMat, Ptr{TSI2JacobianFn}, Ptr{Cvoid}),
               ts, J, P, jac, ctx,
              )


	return nothing
end 

"""
	TSGetI2Jacobian(petsclib::PetscLibType,ts::TS, J::PetscMat, P::PetscMat, jac::TSI2JacobianFn, ctx::Cvoid) 
Returns the implicit Jacobian at the present timestep.

Not Collective, but parallel objects are returned if `TS` is parallel

Input Parameter:
- `ts` - The `TS` context obtained from `TSCreate()`

Output Parameters:
- `J`   - The (approximate) Jacobian of F(t,U,U_t,U_tt)
- `P`   - The matrix from which the preconditioner is constructed, often the same as `J`
- `jac` - The function to compute the Jacobian matrices
- `ctx` - User-defined context for Jacobian evaluation routine

Level: advanced

-seealso: [](ch_ts), `TS`, `TSGetTimeStep()`, `TSGetMatrices()`, `TSGetTime()`, `TSGetStepNumber()`, `TSSetI2Jacobian()`, `TSGetI2Function()`, `TSCreate()`

# External Links
$(_doc_external("Ts/TSGetI2Jacobian"))
"""
function TSGetI2Jacobian(petsclib::PetscLibType, ts::TS, J::PetscMat, P::PetscMat, jac::TSI2JacobianFn, ctx::Cvoid) end

@for_petsc function TSGetI2Jacobian(petsclib::$UnionPetscLib, ts::TS, J::PetscMat, P::PetscMat, jac::TSI2JacobianFn, ctx::Cvoid )
	J_ = Ref(J.ptr)
	P_ = Ref(P.ptr)

    @chk ccall(
               (:TSGetI2Jacobian, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{CMat}, Ptr{CMat}, TSI2JacobianFn, Cvoid),
               ts, J_, P_, jac, ctx,
              )

	J.ptr = C_NULL
	P.ptr = C_NULL

	return nothing
end 

"""
	TSComputeI2Function(petsclib::PetscLibType,ts::TS, t::PetscReal, U::PetscVec, V::PetscVec, A::PetscVec, F::PetscVec) 
Evaluates the DAE residual written in implicit form F(t,U,U_t,U_tt) = 0

Collective

Input Parameters:
- `ts` - the `TS` context
- `t`  - current time
- `U`  - state vector
- `V`  - time derivative of state vector (U_t)
- `A`  - second time derivative of state vector (U_tt)

Output Parameter:
- `F` - the residual vector

Level: developer

-seealso: [](ch_ts), `TS`, `TSSetI2Function()`, `TSGetI2Function()`

# External Links
$(_doc_external("Ts/TSComputeI2Function"))
"""
function TSComputeI2Function(petsclib::PetscLibType, ts::TS, t::PetscReal, U::PetscVec, V::PetscVec, A::PetscVec, F::PetscVec) end

@for_petsc function TSComputeI2Function(petsclib::$UnionPetscLib, ts::TS, t::$PetscReal, U::PetscVec, V::PetscVec, A::PetscVec, F::PetscVec )

    @chk ccall(
               (:TSComputeI2Function, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal, CVec, CVec, CVec, CVec),
               ts, t, U, V, A, F,
              )


	return nothing
end 

"""
	TSComputeI2Jacobian(petsclib::PetscLibType,ts::TS, t::PetscReal, U::PetscVec, V::PetscVec, A::PetscVec, shiftV::PetscReal, shiftA::PetscReal, J::PetscMat, P::PetscMat) 
Evaluates the Jacobian of the DAE

Collective

Input Parameters:
- `ts`     - the `TS` context
- `t`      - current timestep
- `U`      - state vector
- `V`      - time derivative of state vector
- `A`      - second time derivative of state vector
- `shiftV` - shift to apply, see note below
- `shiftA` - shift to apply, see note below

Output Parameters:
- `J` - Jacobian matrix
- `P` - optional matrix used to construct the preconditioner

Level: developer

-seealso: [](ch_ts), `TS`, `TSSetI2Jacobian()`

# External Links
$(_doc_external("Ts/TSComputeI2Jacobian"))
"""
function TSComputeI2Jacobian(petsclib::PetscLibType, ts::TS, t::PetscReal, U::PetscVec, V::PetscVec, A::PetscVec, shiftV::PetscReal, shiftA::PetscReal, J::PetscMat, P::PetscMat) end

@for_petsc function TSComputeI2Jacobian(petsclib::$UnionPetscLib, ts::TS, t::$PetscReal, U::PetscVec, V::PetscVec, A::PetscVec, shiftV::$PetscReal, shiftA::$PetscReal, J::PetscMat, P::PetscMat )

    @chk ccall(
               (:TSComputeI2Jacobian, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal, CVec, CVec, CVec, $PetscReal, $PetscReal, CMat, CMat),
               ts, t, U, V, A, shiftV, shiftA, J, P,
              )


	return nothing
end 

"""
	TSSetTransientVariable(petsclib::PetscLibType,ts::TS, tvar::TSTransientVariableFn, ctx::Cvoid) 
sets function to transform from state to transient variables

Logically Collective

Input Parameters:
- `ts`   - time stepping context on which to change the transient variable
- `tvar` - a function that transforms to transient variables, see `TSTransientVariableFn` for the calling sequence
- `ctx`  - a context for tvar

Level: advanced

-seealso: [](ch_ts), `TS`, `TSBDF`, `TSTransientVariableFn`, `DMTSSetTransientVariable()`, `DMTSGetTransientVariable()`, `TSSetIFunction()`, `TSSetIJacobian()`

# External Links
$(_doc_external("Ts/TSSetTransientVariable"))
"""
function TSSetTransientVariable(petsclib::PetscLibType, ts::TS, tvar::TSTransientVariableFn, ctx::Cvoid) end

@for_petsc function TSSetTransientVariable(petsclib::$UnionPetscLib, ts::TS, tvar::TSTransientVariableFn, ctx::Cvoid )

    @chk ccall(
               (:TSSetTransientVariable, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{TSTransientVariableFn}, Ptr{Cvoid}),
               ts, tvar, ctx,
              )


	return nothing
end 

"""
	TSComputeTransientVariable(petsclib::PetscLibType,ts::TS, U::PetscVec, C::PetscVec) 
transforms state (primitive) variables to transient (conservative) variables

Logically Collective

Input Parameters:
- `ts` - TS on which to compute
- `U`  - state vector to be transformed to transient variables

Output Parameter:
- `C` - transient (conservative) variable

Level: developer

-seealso: [](ch_ts), `TS`, `TSBDF`, `DMTSSetTransientVariable()`, `TSComputeIFunction()`, `TSComputeIJacobian()`

# External Links
$(_doc_external("Ts/TSComputeTransientVariable"))
"""
function TSComputeTransientVariable(petsclib::PetscLibType, ts::TS, U::PetscVec, C::PetscVec) end

@for_petsc function TSComputeTransientVariable(petsclib::$UnionPetscLib, ts::TS, U::PetscVec, C::PetscVec )

    @chk ccall(
               (:TSComputeTransientVariable, $petsc_library),
               PetscErrorCode,
               (CTS, CVec, CVec),
               ts, U, C,
              )


	return nothing
end 

"""
	has::PetscBool = TSHasTransientVariable(petsclib::PetscLibType,ts::TS) 
determine whether transient variables have been set

Logically Collective

Input Parameter:
- `ts` - `TS` on which to compute

Output Parameter:
- `has` - `PETSC_TRUE` if transient variables have been set

Level: developer

-seealso: [](ch_ts), `TS`, `TSBDF`, `DMTSSetTransientVariable()`, `TSComputeTransientVariable()`

# External Links
$(_doc_external("Ts/TSHasTransientVariable"))
"""
function TSHasTransientVariable(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSHasTransientVariable(petsclib::$UnionPetscLib, ts::TS )
	has_ = Ref{PetscBool}()

    @chk ccall(
               (:TSHasTransientVariable, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{PetscBool}),
               ts, has_,
              )

	has = has_[]

	return has
end 

"""
	TS2SetSolution(petsclib::PetscLibType,ts::TS, u::PetscVec, v::PetscVec) 
Sets the initial solution and time derivative vectors
for use by the `TS` routines handling second order equations.

Logically Collective

Input Parameters:
- `ts` - the `TS` context obtained from `TSCreate()`
- `u`  - the solution vector
- `v`  - the time derivative vector

Level: beginner

-seealso: [](ch_ts), `TS`

# External Links
$(_doc_external("Ts/TS2SetSolution"))
"""
function TS2SetSolution(petsclib::PetscLibType, ts::TS, u::PetscVec, v::PetscVec) end

@for_petsc function TS2SetSolution(petsclib::$UnionPetscLib, ts::TS, u::PetscVec, v::PetscVec )

    @chk ccall(
               (:TS2SetSolution, $petsc_library),
               PetscErrorCode,
               (CTS, CVec, CVec),
               ts, u, v,
              )


	return nothing
end 

"""
	TS2GetSolution(petsclib::PetscLibType,ts::TS, u::PetscVec, v::PetscVec) 
Returns the solution and time derivative at the present timestep
for second order equations.

Not Collective

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Output Parameters:
- `u` - the vector containing the solution
- `v` - the vector containing the time derivative

Level: intermediate

-seealso: [](ch_ts), `TS`, `TS2SetSolution()`, `TSGetTimeStep()`, `TSGetTime()`

# External Links
$(_doc_external("Ts/TS2GetSolution"))
"""
function TS2GetSolution(petsclib::PetscLibType, ts::TS, u::PetscVec, v::PetscVec) end

@for_petsc function TS2GetSolution(petsclib::$UnionPetscLib, ts::TS, u::PetscVec, v::PetscVec )
	u_ = Ref(u.ptr)
	v_ = Ref(v.ptr)

    @chk ccall(
               (:TS2GetSolution, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{CVec}, Ptr{CVec}),
               ts, u_, v_,
              )

	u.ptr = C_NULL
	v.ptr = C_NULL

	return nothing
end 

"""
	TSLoad(petsclib::PetscLibType,ts::TS, viewer::PetscViewer) 
Loads a `TS` that has been stored in binary  with `TSView()`.

Collective

Input Parameters:
- `ts`     - the newly loaded `TS`, this needs to have been created with `TSCreate()` or
some related function before a call to `TSLoad()`.
- `viewer` - binary file viewer, obtained from `PetscViewerBinaryOpen()`

Level: intermediate

-seealso: [](ch_ts), `TS`, `PetscViewer`, `PetscViewerBinaryOpen()`, `TSView()`, `MatLoad()`, `VecLoad()`

# External Links
$(_doc_external("Ts/TSLoad"))
"""
function TSLoad(petsclib::PetscLibType, ts::TS, viewer::PetscViewer) end

@for_petsc function TSLoad(petsclib::$UnionPetscLib, ts::TS, viewer::PetscViewer )

    @chk ccall(
               (:TSLoad, $petsc_library),
               PetscErrorCode,
               (CTS, PetscViewer),
               ts, viewer,
              )


	return nothing
end 

"""
	TSViewFromOptions(petsclib::PetscLibType,ts::TS, obj::PetscObject, name::String) 
View a `TS` based on values in the options database

Collective

Input Parameters:
- `ts`   - the `TS` context
- `obj`  - Optional object that provides the prefix for the options database keys
- `name` - command line option string to be passed by user

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSView`, `PetscObjectViewFromOptions()`, `TSCreate()`

# External Links
$(_doc_external("Ts/TSViewFromOptions"))
"""
function TSViewFromOptions(petsclib::PetscLibType, ts::TS, obj::PetscObject, name::String) end

@for_petsc function TSViewFromOptions(petsclib::$UnionPetscLib, ts::TS, obj::PetscObject, name::String )

    @chk ccall(
               (:TSViewFromOptions, $petsc_library),
               PetscErrorCode,
               (CTS, PetscObject, Ptr{Cchar}),
               ts, obj, name,
              )


	return nothing
end 

"""
	TSView(petsclib::PetscLibType,ts::TS, viewer::PetscViewer) 
Prints the `TS` data structure.

Collective

Input Parameters:
- `ts`     - the `TS` context obtained from `TSCreate()`
- `viewer` - visualization context

Options Database Key:
- `-ts_view` - calls `TSView()` at end of `TSStep()`

Level: beginner

-seealso: [](ch_ts), `TS`, `PetscViewer`, `PetscViewerASCIIOpen()`

# External Links
$(_doc_external("Ts/TSView"))
"""
function TSView(petsclib::PetscLibType, ts::TS, viewer::PetscViewer) end

@for_petsc function TSView(petsclib::$UnionPetscLib, ts::TS, viewer::PetscViewer )

    @chk ccall(
               (:TSView, $petsc_library),
               PetscErrorCode,
               (CTS, PetscViewer),
               ts, viewer,
              )


	return nothing
end 

"""
	TSSetApplicationContext(petsclib::PetscLibType,ts::TS, ctx::PeCtx) 
Sets an optional user
`TS` callbacks with `TSGetApplicationContext()`

Logically Collective

Input Parameters:
- `ts`  - the `TS` context obtained from `TSCreate()`
- `ctx` - user context

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSGetApplicationContext()`

# External Links
$(_doc_external("Ts/TSSetApplicationContext"))
"""
function TSSetApplicationContext(petsclib::PetscLibType, ts::TS, ctx::PeCtx) end

@for_petsc function TSSetApplicationContext(petsclib::$UnionPetscLib, ts::TS, ctx::PeCtx )

    @chk ccall(
               (:TSSetApplicationContext, $petsc_library),
               PetscErrorCode,
               (CTS, PeCtx),
               ts, ctx,
              )


	return nothing
end 

"""
	TSGetApplicationContext(petsclib::PetscLibType,ts::TS, ctx::Cvoid) 
Gets the user
timestepper that was set with `TSSetApplicationContext()`

Not Collective

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Output Parameter:
- `ctx` - a pointer to the user context

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSSetApplicationContext()`

# External Links
$(_doc_external("Ts/TSGetApplicationContext"))
"""
function TSGetApplicationContext(petsclib::PetscLibType, ts::TS, ctx::Cvoid) end

@for_petsc function TSGetApplicationContext(petsclib::$UnionPetscLib, ts::TS, ctx::Cvoid )

    @chk ccall(
               (:TSGetApplicationContext, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{Cvoid}),
               ts, ctx,
              )


	return nothing
end 

"""
	steps::PetscInt = TSGetStepNumber(petsclib::PetscLibType,ts::TS) 
Gets the number of time steps completed.

Not Collective

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Output Parameter:
- `steps` - number of steps completed so far

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSGetTime()`, `TSGetTimeStep()`, `TSSetPreStep()`, `TSSetPreStage()`, `TSSetPostStage()`, `TSSetPostStep()`

# External Links
$(_doc_external("Ts/TSGetStepNumber"))
"""
function TSGetStepNumber(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSGetStepNumber(petsclib::$UnionPetscLib, ts::TS )
	steps_ = Ref{$PetscInt}()

    @chk ccall(
               (:TSGetStepNumber, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{$PetscInt}),
               ts, steps_,
              )

	steps = steps_[]

	return steps
end 

"""
	TSSetStepNumber(petsclib::PetscLibType,ts::TS, steps::PetscInt) 
Sets the number of steps completed.

Logically Collective

Input Parameters:
- `ts`    - the `TS` context
- `steps` - number of steps completed so far

Level: developer

-seealso: [](ch_ts), `TS`, `TSGetStepNumber()`, `TSSetTime()`, `TSSetTimeStep()`, `TSSetSolution()`

# External Links
$(_doc_external("Ts/TSSetStepNumber"))
"""
function TSSetStepNumber(petsclib::PetscLibType, ts::TS, steps::PetscInt) end

@for_petsc function TSSetStepNumber(petsclib::$UnionPetscLib, ts::TS, steps::$PetscInt )

    @chk ccall(
               (:TSSetStepNumber, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt),
               ts, steps,
              )


	return nothing
end 

"""
	TSSetTimeStep(petsclib::PetscLibType,ts::TS, time_step::PetscReal) 
Allows one to reset the timestep at any time,
useful for simple pseudo-timestepping codes.

Logically Collective

Input Parameters:
- `ts`        - the `TS` context obtained from `TSCreate()`
- `time_step` - the size of the timestep

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSPSEUDO`, `TSGetTimeStep()`, `TSSetTime()`

# External Links
$(_doc_external("Ts/TSSetTimeStep"))
"""
function TSSetTimeStep(petsclib::PetscLibType, ts::TS, time_step::PetscReal) end

@for_petsc function TSSetTimeStep(petsclib::$UnionPetscLib, ts::TS, time_step::$PetscReal )

    @chk ccall(
               (:TSSetTimeStep, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal),
               ts, time_step,
              )


	return nothing
end 

"""
	TSSetExactFinalTime(petsclib::PetscLibType,ts::TS, eftopt::TSExactFinalTimeOption) 
Determines whether to adapt the final time step to
match the exact final time, to interpolate the solution to the exact final time,
or to just return at the final time `TS` computed (which may be slightly larger
than the requested final time).

Logically Collective

Input Parameters:
- `ts`     - the time-step context
- `eftopt` - exact final time option
-seealso: [](ch_ts), `TS`, `TSExactFinalTimeOption`, `TSGetExactFinalTime()`

# External Links
$(_doc_external("Ts/TSSetExactFinalTime"))
"""
function TSSetExactFinalTime(petsclib::PetscLibType, ts::TS, eftopt::TSExactFinalTimeOption) end

@for_petsc function TSSetExactFinalTime(petsclib::$UnionPetscLib, ts::TS, eftopt::TSExactFinalTimeOption )

    @chk ccall(
               (:TSSetExactFinalTime, $petsc_library),
               PetscErrorCode,
               (CTS, TSExactFinalTimeOption),
               ts, eftopt,
              )


	return nothing
end 

"""
	TSGetExactFinalTime(petsclib::PetscLibType,ts::TS, eftopt::TSExactFinalTimeOption) 
Gets the exact final time option set with `TSSetExactFinalTime()`

Not Collective

Input Parameter:
- `ts` - the `TS` context

Output Parameter:
- `eftopt` - exact final time option

Level: beginner

-seealso: [](ch_ts), `TS`, `TSExactFinalTimeOption`, `TSSetExactFinalTime()`

# External Links
$(_doc_external("Ts/TSGetExactFinalTime"))
"""
function TSGetExactFinalTime(petsclib::PetscLibType, ts::TS, eftopt::TSExactFinalTimeOption) end

@for_petsc function TSGetExactFinalTime(petsclib::$UnionPetscLib, ts::TS, eftopt::TSExactFinalTimeOption )

    @chk ccall(
               (:TSGetExactFinalTime, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{TSExactFinalTimeOption}),
               ts, eftopt,
              )


	return nothing
end 

"""
	dt::PetscReal = TSGetTimeStep(petsclib::PetscLibType,ts::TS) 
Gets the current timestep size.

Not Collective

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Output Parameter:
- `dt` - the current timestep size

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSSetTimeStep()`, `TSGetTime()`

# External Links
$(_doc_external("Ts/TSGetTimeStep"))
"""
function TSGetTimeStep(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSGetTimeStep(petsclib::$UnionPetscLib, ts::TS )
	dt_ = Ref{$PetscReal}()

    @chk ccall(
               (:TSGetTimeStep, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{$PetscReal}),
               ts, dt_,
              )

	dt = dt_[]

	return dt
end 

"""
	TSGetSolution(petsclib::PetscLibType,ts::TS, v::PetscVec) 
Returns the solution at the present timestep. It
is valid to call this routine inside the function that you are evaluating
in order to move to the new timestep. This vector not changed until
the solution at the next timestep has been calculated.

Not Collective, but v returned is parallel if ts is parallel

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Output Parameter:
- `v` - the vector containing the solution

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSGetTimeStep()`, `TSGetTime()`, `TSGetSolveTime()`, `TSGetSolutionComponents()`, `TSSetSolutionFunction()`

# External Links
$(_doc_external("Ts/TSGetSolution"))
"""
function TSGetSolution(petsclib::PetscLibType, ts::TS, v::PetscVec) end

@for_petsc function TSGetSolution(petsclib::$UnionPetscLib, ts::TS, v::PetscVec )
	v_ = Ref(v.ptr)

    @chk ccall(
               (:TSGetSolution, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{CVec}),
               ts, v_,
              )

	v.ptr = C_NULL

	return nothing
end 

"""
	TSGetSolutionComponents(petsclib::PetscLibType,ts::TS, n::PetscInt, v::PetscVec) 
Returns any solution components at the present
timestep, if available for the time integration method being used.
Solution components are quantities that share the same size and
structure as the solution vector.

Not Collective, but v returned is parallel if ts is parallel

Input Parameters:
- `ts` - the `TS` context obtained from `TSCreate()` (input parameter).
- `n`  - If v is `NULL`, then the number of solution components is
returned through n, else the n-th solution component is
returned in v.
- `v`  - the vector containing the n-th solution component
(may be `NULL` to use this function to find out
the number of solutions components).

Level: advanced

-seealso: [](ch_ts), `TS`, `TSGetSolution()`

# External Links
$(_doc_external("Ts/TSGetSolutionComponents"))
"""
function TSGetSolutionComponents(petsclib::PetscLibType, ts::TS, n::PetscInt, v::PetscVec) end

@for_petsc function TSGetSolutionComponents(petsclib::$UnionPetscLib, ts::TS, n::$PetscInt, v::PetscVec )
	v_ = Ref(v.ptr)

    @chk ccall(
               (:TSGetSolutionComponents, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{$PetscInt}, Ptr{CVec}),
               ts, n, v_,
              )

	v.ptr = C_NULL

	return nothing
end 

"""
	TSGetAuxSolution(petsclib::PetscLibType,ts::TS, v::PetscVec) 
Returns an auxiliary solution at the present
timestep, if available for the time integration method being used.

Not Collective, but v returned is parallel if ts is parallel

Input Parameters:
- `ts` - the `TS` context obtained from `TSCreate()` (input parameter).
- `v`  - the vector containing the auxiliary solution

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSGetSolution()`

# External Links
$(_doc_external("Ts/TSGetAuxSolution"))
"""
function TSGetAuxSolution(petsclib::PetscLibType, ts::TS, v::PetscVec) end

@for_petsc function TSGetAuxSolution(petsclib::$UnionPetscLib, ts::TS, v::PetscVec )
	v_ = Ref(v.ptr)

    @chk ccall(
               (:TSGetAuxSolution, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{CVec}),
               ts, v_,
              )

	v.ptr = C_NULL

	return nothing
end 

"""
	TSGetTimeError(petsclib::PetscLibType,ts::TS, n::PetscInt, v::PetscVec) 
Returns the estimated error vector, if the chosen
`TSType` has an error estimation functionality and `TSSetTimeError()` was called

Not Collective, but v returned is parallel if ts is parallel

Input Parameters:
- `ts` - the `TS` context obtained from `TSCreate()` (input parameter).
- `n`  - current estimate (n=0) or previous one (n=-1)
- `v`  - the vector containing the error (same size as the solution).

Level: intermediate

-seealso: [](ch_ts), `TSGetSolution()`, `TSSetTimeError()`

# External Links
$(_doc_external("Ts/TSGetTimeError"))
"""
function TSGetTimeError(petsclib::PetscLibType, ts::TS, n::PetscInt, v::PetscVec) end

@for_petsc function TSGetTimeError(petsclib::$UnionPetscLib, ts::TS, n::$PetscInt, v::PetscVec )
	v_ = Ref(v.ptr)

    @chk ccall(
               (:TSGetTimeError, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt, Ptr{CVec}),
               ts, n, v_,
              )

	v.ptr = C_NULL

	return nothing
end 

"""
	TSSetTimeError(petsclib::PetscLibType,ts::TS, v::PetscVec) 
Sets the estimated error vector, if the chosen
`TSType` has an error estimation functionality. This can be used
to restart such a time integrator with a given error vector.

Not Collective, but v returned is parallel if ts is parallel

Input Parameters:
- `ts` - the `TS` context obtained from `TSCreate()` (input parameter).
- `v`  - the vector containing the error (same size as the solution).

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSSetSolution()`, `TSGetTimeError()`

# External Links
$(_doc_external("Ts/TSSetTimeError"))
"""
function TSSetTimeError(petsclib::PetscLibType, ts::TS, v::PetscVec) end

@for_petsc function TSSetTimeError(petsclib::$UnionPetscLib, ts::TS, v::PetscVec )

    @chk ccall(
               (:TSSetTimeError, $petsc_library),
               PetscErrorCode,
               (CTS, CVec),
               ts, v,
              )


	return nothing
end 

"""
	TSSetProblemType(petsclib::PetscLibType,ts::TS, type::TSProblemType) 
Sets the type of problem to be solved.

Not collective

Input Parameters:
- `ts`   - The `TS`
- `type` - One of `TS_LINEAR`, `TS_NONLINEAR` where these types refer to problems of the forms
-seealso: [](ch_ts), `TSSetUp()`, `TSProblemType`, `TS`

# External Links
$(_doc_external("Ts/TSSetProblemType"))
"""
function TSSetProblemType(petsclib::PetscLibType, ts::TS, type::TSProblemType) end

@for_petsc function TSSetProblemType(petsclib::$UnionPetscLib, ts::TS, type::TSProblemType )

    @chk ccall(
               (:TSSetProblemType, $petsc_library),
               PetscErrorCode,
               (CTS, TSProblemType),
               ts, type,
              )


	return nothing
end 

"""
	type::TSProblemType = TSGetProblemType(petsclib::PetscLibType,ts::TS) 
Gets the type of problem to be solved.

Not collective

Input Parameter:
- `ts` - The `TS`

Output Parameter:
- `type` - One of `TS_LINEAR`, `TS_NONLINEAR` where these types refer to problems of the forms
-seealso: [](ch_ts), `TSSetUp()`, `TSProblemType`, `TS`

# External Links
$(_doc_external("Ts/TSGetProblemType"))
"""
function TSGetProblemType(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSGetProblemType(petsclib::$UnionPetscLib, ts::TS )
	type_ = Ref{TSProblemType}()

    @chk ccall(
               (:TSGetProblemType, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{TSProblemType}),
               ts, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	TSSetUp(petsclib::PetscLibType,ts::TS) 
Sets up the internal data structures for the later use of a timestepper.

Collective

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Level: advanced

-seealso: [](ch_ts), `TSCreate()`, `TS`, `TSStep()`, `TSDestroy()`, `TSSolve()`

# External Links
$(_doc_external("Ts/TSSetUp"))
"""
function TSSetUp(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSSetUp(petsclib::$UnionPetscLib, ts::TS )

    @chk ccall(
               (:TSSetUp, $petsc_library),
               PetscErrorCode,
               (CTS,),
               ts,
              )


	return nothing
end 

"""
	TSReset(petsclib::PetscLibType,ts::TS) 
Resets a `TS` context to the state it was in before `TSSetUp()` was called and removes any allocated `Vec` and `Mat` from its data structures

Collective

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Level: developer

-seealso: [](ch_ts), `TS`, `TSCreate()`, `TSSetUp()`, `TSDestroy()`, `TSSetResize()`

# External Links
$(_doc_external("Ts/TSReset"))
"""
function TSReset(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSReset(petsclib::$UnionPetscLib, ts::TS )

    @chk ccall(
               (:TSReset, $petsc_library),
               PetscErrorCode,
               (CTS,),
               ts,
              )


	return nothing
end 

"""
	TSDestroy(petsclib::PetscLibType,ts::TS) 
Destroys the timestepper context that was created
with `TSCreate()`.

Collective

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Level: beginner

-seealso: [](ch_ts), `TS`, `TSCreate()`, `TSSetUp()`, `TSSolve()`

# External Links
$(_doc_external("Ts/TSDestroy"))
"""
function TSDestroy(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSDestroy(petsclib::$UnionPetscLib, ts::TS )
	ts_ = Ref(ts.ptr)

    @chk ccall(
               (:TSDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{CTS},),
               ts_,
              )

	ts.ptr = C_NULL

	return nothing
end 

"""
	TSGetSNES(petsclib::PetscLibType,ts::TS, snes::PetscSNES) 
Returns the `SNES` (nonlinear solver) associated with
a `TS` (timestepper) context. Valid only for nonlinear problems.

Not Collective, but snes is parallel if ts is parallel

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Output Parameter:
- `snes` - the nonlinear solver context

Level: beginner

-seealso: [](ch_ts), `TS`, `SNES`, `TSCreate()`, `TSSetUp()`, `TSSolve()`

# External Links
$(_doc_external("Ts/TSGetSNES"))
"""
function TSGetSNES(petsclib::PetscLibType, ts::TS, snes::PetscSNES) end

@for_petsc function TSGetSNES(petsclib::$UnionPetscLib, ts::TS, snes::PetscSNES )
	snes_ = Ref(snes.ptr)

    @chk ccall(
               (:TSGetSNES, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{CSNES}),
               ts, snes_,
              )

	snes.ptr = C_NULL

	return nothing
end 

"""
	TSSetSNES(petsclib::PetscLibType,ts::TS, snes::PetscSNES) 
Set the `SNES` (nonlinear solver) to be used by the `TS` timestepping context

Collective

Input Parameters:
- `ts`   - the `TS` context obtained from `TSCreate()`
- `snes` - the nonlinear solver context

Level: developer

-seealso: [](ch_ts), `TS`, `SNES`, `TSCreate()`, `TSSetUp()`, `TSSolve()`, `TSGetSNES()`

# External Links
$(_doc_external("Ts/TSSetSNES"))
"""
function TSSetSNES(petsclib::PetscLibType, ts::TS, snes::PetscSNES) end

@for_petsc function TSSetSNES(petsclib::$UnionPetscLib, ts::TS, snes::PetscSNES )

    @chk ccall(
               (:TSSetSNES, $petsc_library),
               PetscErrorCode,
               (CTS, CSNES),
               ts, snes,
              )


	return nothing
end 

"""
	TSGetKSP(petsclib::PetscLibType,ts::TS, ksp::PetscKSP) 
Returns the `KSP` (linear solver) associated with
a `TS` (timestepper) context.

Not Collective, but `ksp` is parallel if `ts` is parallel

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Output Parameter:
- `ksp` - the nonlinear solver context

Level: beginner

-seealso: [](ch_ts), `TS`, `SNES`, `KSP`, `TSCreate()`, `TSSetUp()`, `TSSolve()`, `TSGetSNES()`

# External Links
$(_doc_external("Ts/TSGetKSP"))
"""
function TSGetKSP(petsclib::PetscLibType, ts::TS, ksp::PetscKSP) end

@for_petsc function TSGetKSP(petsclib::$UnionPetscLib, ts::TS, ksp::PetscKSP )
	ksp_ = Ref(ksp.ptr)

    @chk ccall(
               (:TSGetKSP, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{CKSP}),
               ts, ksp_,
              )

	ksp.ptr = C_NULL

	return nothing
end 

"""
	TSSetMaxSteps(petsclib::PetscLibType,ts::TS, maxsteps::PetscInt) 
Sets the maximum number of steps to use.

Logically Collective

Input Parameters:
- `ts`       - the `TS` context obtained from `TSCreate()`
- `maxsteps` - maximum number of steps to use

Options Database Key:
- `-ts_max_steps <maxsteps>` - Sets maxsteps

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSGetMaxSteps()`, `TSSetMaxTime()`, `TSSetExactFinalTime()`

# External Links
$(_doc_external("Ts/TSSetMaxSteps"))
"""
function TSSetMaxSteps(petsclib::PetscLibType, ts::TS, maxsteps::PetscInt) end

@for_petsc function TSSetMaxSteps(petsclib::$UnionPetscLib, ts::TS, maxsteps::$PetscInt )

    @chk ccall(
               (:TSSetMaxSteps, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt),
               ts, maxsteps,
              )


	return nothing
end 

"""
	maxsteps::PetscInt = TSGetMaxSteps(petsclib::PetscLibType,ts::TS) 
Gets the maximum number of steps to use.

Not Collective

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Output Parameter:
- `maxsteps` - maximum number of steps to use

Level: advanced

-seealso: [](ch_ts), `TS`, `TSSetMaxSteps()`, `TSGetMaxTime()`, `TSSetMaxTime()`

# External Links
$(_doc_external("Ts/TSGetMaxSteps"))
"""
function TSGetMaxSteps(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSGetMaxSteps(petsclib::$UnionPetscLib, ts::TS )
	maxsteps_ = Ref{$PetscInt}()

    @chk ccall(
               (:TSGetMaxSteps, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{$PetscInt}),
               ts, maxsteps_,
              )

	maxsteps = maxsteps_[]

	return maxsteps
end 

"""
	TSSetRunSteps(petsclib::PetscLibType,ts::TS, runsteps::PetscInt) 
Sets the maximum number of steps to take in each call to `TSSolve()`.

If the step count when `TSSolve()` is `start_step`, this will stop the simulation once `current_step - start_step >= run_steps`.
Comparatively, `TSSetMaxSteps()` will stop if `current_step >= max_steps`.
The simulation will stop when either condition is reached.

Logically Collective

Input Parameters:
- `ts`       - the `TS` context obtained from `TSCreate()`
- `runsteps` - maximum number of steps to take in each call to `TSSolve()`;

Options Database Key:
- `-ts_run_steps <runsteps>` - Sets runsteps

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSGetRunSteps()`, `TSSetMaxTime()`, `TSSetExactFinalTime()`, `TSSetMaxSteps()`

# External Links
$(_doc_external("Ts/TSSetRunSteps"))
"""
function TSSetRunSteps(petsclib::PetscLibType, ts::TS, runsteps::PetscInt) end

@for_petsc function TSSetRunSteps(petsclib::$UnionPetscLib, ts::TS, runsteps::$PetscInt )

    @chk ccall(
               (:TSSetRunSteps, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt),
               ts, runsteps,
              )


	return nothing
end 

"""
	runsteps::PetscInt = TSGetRunSteps(petsclib::PetscLibType,ts::TS) 
Gets the maximum number of steps to take in each call to `TSSolve()`.

Not Collective

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Output Parameter:
- `runsteps` - maximum number of steps to take in each call to `TSSolve`.

Level: advanced

-seealso: [](ch_ts), `TS`, `TSSetRunSteps()`, `TSGetMaxTime()`, `TSSetMaxTime()`, `TSGetMaxSteps()`

# External Links
$(_doc_external("Ts/TSGetRunSteps"))
"""
function TSGetRunSteps(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSGetRunSteps(petsclib::$UnionPetscLib, ts::TS )
	runsteps_ = Ref{$PetscInt}()

    @chk ccall(
               (:TSGetRunSteps, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{$PetscInt}),
               ts, runsteps_,
              )

	runsteps = runsteps_[]

	return runsteps
end 

"""
	TSSetMaxTime(petsclib::PetscLibType,ts::TS, maxtime::PetscReal) 
Sets the maximum (or final) time for timestepping.

Logically Collective

Input Parameters:
- `ts`      - the `TS` context obtained from `TSCreate()`
- `maxtime` - final time to step to

Options Database Key:
- `-ts_max_time <maxtime>` - Sets maxtime

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSGetMaxTime()`, `TSSetMaxSteps()`, `TSSetExactFinalTime()`

# External Links
$(_doc_external("Ts/TSSetMaxTime"))
"""
function TSSetMaxTime(petsclib::PetscLibType, ts::TS, maxtime::PetscReal) end

@for_petsc function TSSetMaxTime(petsclib::$UnionPetscLib, ts::TS, maxtime::$PetscReal )

    @chk ccall(
               (:TSSetMaxTime, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal),
               ts, maxtime,
              )


	return nothing
end 

"""
	maxtime::PetscReal = TSGetMaxTime(petsclib::PetscLibType,ts::TS) 
Gets the maximum (or final) time for timestepping.

Not Collective

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Output Parameter:
- `maxtime` - final time to step to

Level: advanced

-seealso: [](ch_ts), `TS`, `TSSetMaxTime()`, `TSGetMaxSteps()`, `TSSetMaxSteps()`

# External Links
$(_doc_external("Ts/TSGetMaxTime"))
"""
function TSGetMaxTime(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSGetMaxTime(petsclib::$UnionPetscLib, ts::TS )
	maxtime_ = Ref{$PetscReal}()

    @chk ccall(
               (:TSGetMaxTime, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{$PetscReal}),
               ts, maxtime_,
              )

	maxtime = maxtime_[]

	return maxtime
end 

"""
	TSSetSolution(petsclib::PetscLibType,ts::TS, u::PetscVec) 
Sets the initial solution vector
for use by the `TS` routines.

Logically Collective

Input Parameters:
- `ts` - the `TS` context obtained from `TSCreate()`
- `u`  - the solution vector

Level: beginner

-seealso: [](ch_ts), `TS`, `TSSetSolutionFunction()`, `TSGetSolution()`, `TSCreate()`

# External Links
$(_doc_external("Ts/TSSetSolution"))
"""
function TSSetSolution(petsclib::PetscLibType, ts::TS, u::PetscVec) end

@for_petsc function TSSetSolution(petsclib::$UnionPetscLib, ts::TS, u::PetscVec )

    @chk ccall(
               (:TSSetSolution, $petsc_library),
               PetscErrorCode,
               (CTS, CVec),
               ts, u,
              )


	return nothing
end 

"""
	TSSetPreStep(petsclib::PetscLibType,ts::TS, func::external) 
Sets the general
called once at the beginning of each time step.

Logically Collective

Input Parameters:
- `ts`   - The `TS` context obtained from `TSCreate()`
- `func` - The function

Calling sequence of `func`:
- `ts` - the `TS` context

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSSetPreStage()`, `TSSetPostStage()`, `TSSetPostStep()`, `TSStep()`, `TSRestartStep()`

# External Links
$(_doc_external("Ts/TSSetPreStep"))
"""
function TSSetPreStep(petsclib::PetscLibType, ts::TS, func::external) end

@for_petsc function TSSetPreStep(petsclib::$UnionPetscLib, ts::TS, func::external )

    @chk ccall(
               (:TSSetPreStep, $petsc_library),
               PetscErrorCode,
               (CTS, external),
               ts, func,
              )


	return nothing
end 

"""
	TSPreStep(petsclib::PetscLibType,ts::TS) 
Runs the user

Collective

Input Parameter:
- `ts` - The `TS` context obtained from `TSCreate()`

Level: developer

-seealso: [](ch_ts), `TS`, `TSSetPreStep()`, `TSPreStage()`, `TSPostStage()`, `TSPostStep()`

# External Links
$(_doc_external("Ts/TSPreStep"))
"""
function TSPreStep(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSPreStep(petsclib::$UnionPetscLib, ts::TS )

    @chk ccall(
               (:TSPreStep, $petsc_library),
               PetscErrorCode,
               (CTS,),
               ts,
              )


	return nothing
end 

"""
	TSSetPreStage(petsclib::PetscLibType,ts::TS, func::external) 
Sets the general
called once at the beginning of each stage.

Logically Collective

Input Parameters:
- `ts`   - The `TS` context obtained from `TSCreate()`
- `func` - The function

Calling sequence of `func`:
- `ts`        - the `TS` context
- `stagetime` - the stage time

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSSetPostStage()`, `TSSetPreStep()`, `TSSetPostStep()`, `TSGetApplicationContext()`

# External Links
$(_doc_external("Ts/TSSetPreStage"))
"""
function TSSetPreStage(petsclib::PetscLibType, ts::TS, func::external) end

@for_petsc function TSSetPreStage(petsclib::$UnionPetscLib, ts::TS, func::external )

    @chk ccall(
               (:TSSetPreStage, $petsc_library),
               PetscErrorCode,
               (CTS, external),
               ts, func,
              )


	return nothing
end 

"""
	TSSetPostStage(petsclib::PetscLibType,ts::TS, func::external) 
Sets the general
called once at the end of each stage.

Logically Collective

Input Parameters:
- `ts`   - The `TS` context obtained from `TSCreate()`
- `func` - The function

Calling sequence of `func`:
- `ts`         - the `TS` context
- `stagetime`  - the stage time
- `stageindex` - the stage index
- `Y`          - Array of vectors (of size = total number of stages) with the stage solutions

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSSetPreStage()`, `TSSetPreStep()`, `TSSetPostStep()`, `TSGetApplicationContext()`

# External Links
$(_doc_external("Ts/TSSetPostStage"))
"""
function TSSetPostStage(petsclib::PetscLibType, ts::TS, func::external) end

@for_petsc function TSSetPostStage(petsclib::$UnionPetscLib, ts::TS, func::external )

    @chk ccall(
               (:TSSetPostStage, $petsc_library),
               PetscErrorCode,
               (CTS, external),
               ts, func,
              )


	return nothing
end 

"""
	TSSetPostEvaluate(petsclib::PetscLibType,ts::TS, func::external) 
Sets the general
called at the end of each step evaluation.

Logically Collective

Input Parameters:
- `ts`   - The `TS` context obtained from `TSCreate()`
- `func` - The function

Calling sequence of `func`:
- `ts` - the `TS` context

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSSetPreStage()`, `TSSetPreStep()`, `TSSetPostStep()`, `TSGetApplicationContext()`

# External Links
$(_doc_external("Ts/TSSetPostEvaluate"))
"""
function TSSetPostEvaluate(petsclib::PetscLibType, ts::TS, func::external) end

@for_petsc function TSSetPostEvaluate(petsclib::$UnionPetscLib, ts::TS, func::external )

    @chk ccall(
               (:TSSetPostEvaluate, $petsc_library),
               PetscErrorCode,
               (CTS, external),
               ts, func,
              )


	return nothing
end 

"""
	TSPreStage(petsclib::PetscLibType,ts::TS, stagetime::PetscReal) 
Runs the user

Collective

Input Parameters:
- `ts`        - The `TS` context obtained from `TSCreate()`
- `stagetime` - The absolute time of the current stage

Level: developer

-seealso: [](ch_ts), `TS`, `TSPostStage()`, `TSSetPreStep()`, `TSPreStep()`, `TSPostStep()`

# External Links
$(_doc_external("Ts/TSPreStage"))
"""
function TSPreStage(petsclib::PetscLibType, ts::TS, stagetime::PetscReal) end

@for_petsc function TSPreStage(petsclib::$UnionPetscLib, ts::TS, stagetime::$PetscReal )

    @chk ccall(
               (:TSPreStage, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal),
               ts, stagetime,
              )


	return nothing
end 

"""
	TSPostStage(petsclib::PetscLibType,ts::TS, stagetime::PetscReal, stageindex::PetscInt, Y::Vector{PetscVec}) 
Runs the user

Collective

Input Parameters:
- `ts`         - The `TS` context obtained from `TSCreate()`
- `stagetime`  - The absolute time of the current stage
- `stageindex` - Stage number
- `Y`          - Array of vectors (of size = total number of stages) with the stage solutions

Level: developer

-seealso: [](ch_ts), `TS`, `TSPreStage()`, `TSSetPreStep()`, `TSPreStep()`, `TSPostStep()`

# External Links
$(_doc_external("Ts/TSPostStage"))
"""
function TSPostStage(petsclib::PetscLibType, ts::TS, stagetime::PetscReal, stageindex::PetscInt, Y::Vector{PetscVec}) end

@for_petsc function TSPostStage(petsclib::$UnionPetscLib, ts::TS, stagetime::$PetscReal, stageindex::$PetscInt, Y::Vector{PetscVec} )

    @chk ccall(
               (:TSPostStage, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal, $PetscInt, Ptr{CVec}),
               ts, stagetime, stageindex, Y,
              )


	return nothing
end 

"""
	TSPostEvaluate(petsclib::PetscLibType,ts::TS) 
Runs the user

Collective

Input Parameter:
- `ts` - The `TS` context obtained from `TSCreate()`

Level: developer

-seealso: [](ch_ts), `TS`, `TSSetPostEvaluate()`, `TSSetPreStep()`, `TSPreStep()`, `TSPostStep()`

# External Links
$(_doc_external("Ts/TSPostEvaluate"))
"""
function TSPostEvaluate(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSPostEvaluate(petsclib::$UnionPetscLib, ts::TS )

    @chk ccall(
               (:TSPostEvaluate, $petsc_library),
               PetscErrorCode,
               (CTS,),
               ts,
              )


	return nothing
end 

"""
	TSSetPostStep(petsclib::PetscLibType,ts::TS, func::external) 
Sets the general
called once at the end of each successful time step.

Logically Collective

Input Parameters:
- `ts`   - The `TS` context obtained from `TSCreate()`
- `func` - The function

Calling sequence of `func`:
- `ts` - the `TS` context

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSSetPreStep()`, `TSSetPreStage()`, `TSSetPostEvaluate()`, `TSGetTimeStep()`, `TSGetStepNumber()`, `TSGetTime()`, `TSRestartStep()`

# External Links
$(_doc_external("Ts/TSSetPostStep"))
"""
function TSSetPostStep(petsclib::PetscLibType, ts::TS, func::external) end

@for_petsc function TSSetPostStep(petsclib::$UnionPetscLib, ts::TS, func::external )

    @chk ccall(
               (:TSSetPostStep, $petsc_library),
               PetscErrorCode,
               (CTS, external),
               ts, func,
              )


	return nothing
end 

"""
	TSPostStep(petsclib::PetscLibType,ts::TS) 
Runs the user

Collective

Input Parameter:
- `ts` - The `TS` context obtained from `TSCreate()`

-seealso: [](ch_ts), `TS`, `TSSetPreStep()`, `TSSetPreStage()`, `TSSetPostEvaluate()`, `TSGetTimeStep()`, `TSGetStepNumber()`, `TSGetTime()`, `TSSetPostStep()`

# External Links
$(_doc_external("Ts/TSPostStep"))
"""
function TSPostStep(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSPostStep(petsclib::$UnionPetscLib, ts::TS )

    @chk ccall(
               (:TSPostStep, $petsc_library),
               PetscErrorCode,
               (CTS,),
               ts,
              )


	return nothing
end 

"""
	TSInterpolate(petsclib::PetscLibType,ts::TS, t::PetscReal, U::PetscVec) 
Interpolate the solution computed during the previous step to an arbitrary location in the interval

Collective

Input Parameters:
- `ts` - time stepping context
- `t`  - time to interpolate to

Output Parameter:
- `U` - state at given time

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSSetExactFinalTime()`, `TSSolve()`

# External Links
$(_doc_external("Ts/TSInterpolate"))
"""
function TSInterpolate(petsclib::PetscLibType, ts::TS, t::PetscReal, U::PetscVec) end

@for_petsc function TSInterpolate(petsclib::$UnionPetscLib, ts::TS, t::$PetscReal, U::PetscVec )

    @chk ccall(
               (:TSInterpolate, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal, CVec),
               ts, t, U,
              )


	return nothing
end 

"""
	TSStep(petsclib::PetscLibType,ts::TS) 
Steps one time step

Collective

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Level: developer

-seealso: [](ch_ts), `TS`, `TSCreate()`, `TSSetUp()`, `TSDestroy()`, `TSSolve()`, `TSSetPreStep()`, `TSSetPreStage()`, `TSSetPostStage()`, `TSInterpolate()`

# External Links
$(_doc_external("Ts/TSStep"))
"""
function TSStep(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSStep(petsclib::$UnionPetscLib, ts::TS )

    @chk ccall(
               (:TSStep, $petsc_library),
               PetscErrorCode,
               (CTS,),
               ts,
              )


	return nothing
end 

"""
	order::PetscInt,wlte::PetscReal = TSEvaluateWLTE(petsclib::PetscLibType,ts::TS, wnormtype::NormType) 
Evaluate the weighted local truncation error norm
at the end of a time step with a given order of accuracy.

Collective

Input Parameters:
- `ts`        - time stepping context
- `wnormtype` - norm type, either `NORM_2` or `NORM_INFINITY`

Input/Output Parameter:
- `order` - optional, desired order for the error evaluation or `PETSC_DECIDE`;
on output, the actual order of the error evaluation

Output Parameter:
- `wlte` - the weighted local truncation error norm

Level: advanced

-seealso: [](ch_ts), `TS`, `TSStep()`, `TSAdapt`, `TSErrorWeightedNorm()`

# External Links
$(_doc_external("Ts/TSEvaluateWLTE"))
"""
function TSEvaluateWLTE(petsclib::PetscLibType, ts::TS, wnormtype::NormType) end

@for_petsc function TSEvaluateWLTE(petsclib::$UnionPetscLib, ts::TS, wnormtype::NormType )
	order_ = Ref{$PetscInt}()
	wlte_ = Ref{$PetscReal}()

    @chk ccall(
               (:TSEvaluateWLTE, $petsc_library),
               PetscErrorCode,
               (CTS, NormType, Ptr{$PetscInt}, Ptr{$PetscReal}),
               ts, wnormtype, order_, wlte_,
              )

	order = order_[]
	wlte = wlte_[]

	return order,wlte
end 

"""
	TSEvaluateStep(petsclib::PetscLibType,ts::TS, order::PetscInt, U::PetscVec, done::PetscBool) 
Evaluate the solution at the end of a time step with a given order of accuracy.

Collective

Input Parameters:
- `ts`    - time stepping context
- `order` - desired order of accuracy
- `done`  - whether the step was evaluated at this order (pass `NULL` to generate an error if not available)

Output Parameter:
- `U` - state at the end of the current step

Level: advanced

-seealso: [](ch_ts), `TS`, `TSStep()`, `TSAdapt`

# External Links
$(_doc_external("Ts/TSEvaluateStep"))
"""
function TSEvaluateStep(petsclib::PetscLibType, ts::TS, order::PetscInt, U::PetscVec, done::PetscBool) end

@for_petsc function TSEvaluateStep(petsclib::$UnionPetscLib, ts::TS, order::$PetscInt, U::PetscVec, done::PetscBool )

    @chk ccall(
               (:TSEvaluateStep, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt, CVec, Ptr{PetscBool}),
               ts, order, U, done,
              )


	return nothing
end 

"""
	TSSetComputeInitialCondition(petsclib::PetscLibType,ts::TS, initCondition::external) 
Set the function used to automatically compute an initial condition for the timestepping.

Logically collective

Input Parameters:
- `ts`            - time stepping context
- `initCondition` - The function which computes an initial condition

Calling sequence of `initCondition`:
- `ts` - The timestepping context
- `e`  - The input vector in which the initial condition is to be stored

Level: advanced

-seealso: [](ch_ts), `TS`, `TSGetComputeInitialCondition()`, `TSComputeInitialCondition()`

# External Links
$(_doc_external("Ts/TSSetComputeInitialCondition"))
"""
function TSSetComputeInitialCondition(petsclib::PetscLibType, ts::TS, initCondition::external) end

@for_petsc function TSSetComputeInitialCondition(petsclib::$UnionPetscLib, ts::TS, initCondition::external )

    @chk ccall(
               (:TSSetComputeInitialCondition, $petsc_library),
               PetscErrorCode,
               (CTS, external),
               ts, initCondition,
              )


	return nothing
end 

"""
	TSComputeInitialCondition(petsclib::PetscLibType,ts::TS, u::PetscVec) 
Compute an initial condition for the timestepping using the function previously set with `TSSetComputeInitialCondition()`

Collective

Input Parameters:
- `ts` - time stepping context
- `u`  - The `Vec` to store the condition in which will be used in `TSSolve()`

Level: advanced

-seealso: [](ch_ts), `TS`, `TSGetComputeInitialCondition()`, `TSSetComputeInitialCondition()`, `TSSolve()`

# External Links
$(_doc_external("Ts/TSComputeInitialCondition"))
"""
function TSComputeInitialCondition(petsclib::PetscLibType, ts::TS, u::PetscVec) end

@for_petsc function TSComputeInitialCondition(petsclib::$UnionPetscLib, ts::TS, u::PetscVec )

    @chk ccall(
               (:TSComputeInitialCondition, $petsc_library),
               PetscErrorCode,
               (CTS, CVec),
               ts, u,
              )


	return nothing
end 

"""
	TSSetComputeExactError(petsclib::PetscLibType,ts::TS, exactError::external) 
Set the function used to automatically compute the exact error for the timestepping.

Logically collective

Input Parameters:
- `ts`         - time stepping context
- `exactError` - The function which computes the solution error

Calling sequence of `exactError`:
- `ts` - The timestepping context
- `u`  - The approximate solution vector
- `e`  - The  vector in which the error is stored

Level: advanced

-seealso: [](ch_ts), `TS`, `TSGetComputeExactError()`, `TSComputeExactError()`

# External Links
$(_doc_external("Ts/TSSetComputeExactError"))
"""
function TSSetComputeExactError(petsclib::PetscLibType, ts::TS, exactError::external) end

@for_petsc function TSSetComputeExactError(petsclib::$UnionPetscLib, ts::TS, exactError::external )

    @chk ccall(
               (:TSSetComputeExactError, $petsc_library),
               PetscErrorCode,
               (CTS, external),
               ts, exactError,
              )


	return nothing
end 

"""
	TSComputeExactError(petsclib::PetscLibType,ts::TS, u::PetscVec, e::PetscVec) 
Compute the solution error for the timestepping using the function previously set with `TSSetComputeExactError()`

Collective

Input Parameters:
- `ts` - time stepping context
- `u`  - The approximate solution
- `e`  - The `Vec` used to store the error

Level: advanced

-seealso: [](ch_ts), `TS`, `TSGetComputeInitialCondition()`, `TSSetComputeInitialCondition()`, `TSSolve()`

# External Links
$(_doc_external("Ts/TSComputeExactError"))
"""
function TSComputeExactError(petsclib::PetscLibType, ts::TS, u::PetscVec, e::PetscVec) end

@for_petsc function TSComputeExactError(petsclib::$UnionPetscLib, ts::TS, u::PetscVec, e::PetscVec )

    @chk ccall(
               (:TSComputeExactError, $petsc_library),
               PetscErrorCode,
               (CTS, CVec, CVec),
               ts, u, e,
              )


	return nothing
end 

"""
	TSSetResize(petsclib::PetscLibType,ts::TS, rollback::PetscBool, setup::external, transfer::external, ctx::Cvoid) 
Sets the resize callbacks.

Logically Collective

Input Parameters:
- `ts`       - The `TS` context obtained from `TSCreate()`
- `rollback` - Whether a resize will restart the step
- `setup`    - The setup function
- `transfer` - The transfer function
- `ctx`      - [optional] The user-defined context

Calling sequence of `setup`:
- `ts`     - the `TS` context
- `step`   - the current step
- `time`   - the current time
- `state`  - the current vector of state
- `resize` - (output parameter) `PETSC_TRUE` if need resizing, `PETSC_FALSE` otherwise
- `ctx`    - user defined context

Calling sequence of `transfer`:
- `ts`      - the `TS` context
- `nv`      - the number of vectors to be transferred
- `vecsin`  - array of vectors to be transferred
- `vecsout` - array of transferred vectors
- `ctx`     - user defined context

-seealso: [](ch_ts), `TS`, `TSSetDM()`, `TSSetIJacobian()`, `TSSetRHSJacobian()`

# External Links
$(_doc_external("Ts/TSSetResize"))
"""
function TSSetResize(petsclib::PetscLibType, ts::TS, rollback::PetscBool, setup::external, transfer::external, ctx::Cvoid) end

@for_petsc function TSSetResize(petsclib::$UnionPetscLib, ts::TS, rollback::PetscBool, setup::external, transfer::external, ctx::Cvoid )

    @chk ccall(
               (:TSSetResize, $petsc_library),
               PetscErrorCode,
               (CTS, PetscBool, external, external, Ptr{Cvoid}),
               ts, rollback, setup, transfer, ctx,
              )


	return nothing
end 

"""
	TSResizeRegisterVec(petsclib::PetscLibType,ts::TS, name::String, vec::PetscVec) 
Register a vector to be transferred with `TSResize()`.

Collective

Input Parameters:
- `ts`   - The `TS` context obtained from `TSCreate()`
- `name` - A string identifying the vector
- `vec`  - The vector

Level: developer

-seealso: [](ch_ts), `TS`, `TSSetResize()`, `TSResize()`, `TSResizeRetrieveVec()`

# External Links
$(_doc_external("Ts/TSResizeRegisterVec"))
"""
function TSResizeRegisterVec(petsclib::PetscLibType, ts::TS, name::String, vec::PetscVec) end

@for_petsc function TSResizeRegisterVec(petsclib::$UnionPetscLib, ts::TS, name::String, vec::PetscVec )

    @chk ccall(
               (:TSResizeRegisterVec, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{Cchar}, CVec),
               ts, name, vec,
              )


	return nothing
end 

"""
	TSResizeRetrieveVec(petsclib::PetscLibType,ts::TS, name::String, vec::PetscVec) 
Retrieve a vector registered with `TSResizeRegisterVec()`.

Collective

Input Parameters:
- `ts`   - The `TS` context obtained from `TSCreate()`
- `name` - A string identifying the vector
- `vec`  - The vector

Level: developer

-seealso: [](ch_ts), `TS`, `TSSetResize()`, `TSResize()`, `TSResizeRegisterVec()`

# External Links
$(_doc_external("Ts/TSResizeRetrieveVec"))
"""
function TSResizeRetrieveVec(petsclib::PetscLibType, ts::TS, name::String, vec::PetscVec) end

@for_petsc function TSResizeRetrieveVec(petsclib::$UnionPetscLib, ts::TS, name::String, vec::PetscVec )
	vec_ = Ref(vec.ptr)

    @chk ccall(
               (:TSResizeRetrieveVec, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{Cchar}, Ptr{CVec}),
               ts, name, vec_,
              )

	vec.ptr = C_NULL

	return nothing
end 

"""
	TSResize(petsclib::PetscLibType,ts::TS) 
Runs the user

Collective

Input Parameter:
- `ts` - The `TS` context obtained from `TSCreate()`

Level: developer

-seealso: [](ch_ts), `TS`, `TSSetResize()`

# External Links
$(_doc_external("Ts/TSResize"))
"""
function TSResize(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSResize(petsclib::$UnionPetscLib, ts::TS )

    @chk ccall(
               (:TSResize, $petsc_library),
               PetscErrorCode,
               (CTS,),
               ts,
              )


	return nothing
end 

"""
	TSSolve(petsclib::PetscLibType,ts::TS, u::PetscVec) 
Steps the requested number of timesteps.

Collective

Input Parameters:
- `ts` - the `TS` context obtained from `TSCreate()`
- `u`  - the solution vector  (can be null if `TSSetSolution()` was used and `TSSetExactFinalTime`(ts,`TS_EXACTFINALTIME_MATCHSTEP`) was not used,
otherwise it must contain the initial conditions and will contain the solution at the final requested time

Level: beginner

-seealso: [](ch_ts), `TS`, `TSCreate()`, `TSSetSolution()`, `TSStep()`, `TSGetTime()`, `TSGetSolveTime()`

# External Links
$(_doc_external("Ts/TSSolve"))
"""
function TSSolve(petsclib::PetscLibType, ts::TS, u::PetscVec) end

@for_petsc function TSSolve(petsclib::$UnionPetscLib, ts::TS, u::PetscVec )

    @chk ccall(
               (:TSSolve, $petsc_library),
               PetscErrorCode,
               (CTS, CVec),
               ts, u,
              )


	return nothing
end 

"""
	t::PetscReal = TSGetTime(petsclib::PetscLibType,ts::TS) 
Gets the time of the most recently completed step.

Not Collective

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Output Parameter:
- `t` - the current time. This time may not corresponds to the final time set with `TSSetMaxTime()`, use `TSGetSolveTime()`.

Level: beginner

-seealso: [](ch_ts), `TS`, `TSGetSolveTime()`, `TSSetTime()`, `TSGetTimeStep()`, `TSGetStepNumber()`

# External Links
$(_doc_external("Ts/TSGetTime"))
"""
function TSGetTime(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSGetTime(petsclib::$UnionPetscLib, ts::TS )
	t_ = Ref{$PetscReal}()

    @chk ccall(
               (:TSGetTime, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{$PetscReal}),
               ts, t_,
              )

	t = t_[]

	return t
end 

"""
	t::PetscReal = TSGetPrevTime(petsclib::PetscLibType,ts::TS) 
Gets the starting time of the previously completed step.

Not Collective

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Output Parameter:
- `t` - the previous time

Level: beginner

-seealso: [](ch_ts), `TS`, `TSGetTime()`, `TSGetSolveTime()`, `TSGetTimeStep()`

# External Links
$(_doc_external("Ts/TSGetPrevTime"))
"""
function TSGetPrevTime(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSGetPrevTime(petsclib::$UnionPetscLib, ts::TS )
	t_ = Ref{$PetscReal}()

    @chk ccall(
               (:TSGetPrevTime, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{$PetscReal}),
               ts, t_,
              )

	t = t_[]

	return t
end 

"""
	TSSetTime(petsclib::PetscLibType,ts::TS, t::PetscReal) 
Allows one to reset the time.

Logically Collective

Input Parameters:
- `ts` - the `TS` context obtained from `TSCreate()`
- `t`  - the time

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSGetTime()`, `TSSetMaxSteps()`

# External Links
$(_doc_external("Ts/TSSetTime"))
"""
function TSSetTime(petsclib::PetscLibType, ts::TS, t::PetscReal) end

@for_petsc function TSSetTime(petsclib::$UnionPetscLib, ts::TS, t::$PetscReal )

    @chk ccall(
               (:TSSetTime, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal),
               ts, t,
              )


	return nothing
end 

"""
	TSSetOptionsPrefix(petsclib::PetscLibType,ts::TS, prefix::String) 
Sets the prefix used for searching for all
TS options in the database.

Logically Collective

Input Parameters:
- `ts`     - The `TS` context
- `prefix` - The prefix to prepend to all option names

Level: advanced

-seealso: [](ch_ts), `TS`, `TSSetFromOptions()`, `TSAppendOptionsPrefix()`

# External Links
$(_doc_external("Ts/TSSetOptionsPrefix"))
"""
function TSSetOptionsPrefix(petsclib::PetscLibType, ts::TS, prefix::String) end

@for_petsc function TSSetOptionsPrefix(petsclib::$UnionPetscLib, ts::TS, prefix::String )

    @chk ccall(
               (:TSSetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{Cchar}),
               ts, prefix,
              )


	return nothing
end 

"""
	TSAppendOptionsPrefix(petsclib::PetscLibType,ts::TS, prefix::String) 
Appends to the prefix used for searching for all
TS options in the database.

Logically Collective

Input Parameters:
- `ts`     - The `TS` context
- `prefix` - The prefix to prepend to all option names

Level: advanced

-seealso: [](ch_ts), `TS`, `TSGetOptionsPrefix()`, `TSSetOptionsPrefix()`, `TSSetFromOptions()`

# External Links
$(_doc_external("Ts/TSAppendOptionsPrefix"))
"""
function TSAppendOptionsPrefix(petsclib::PetscLibType, ts::TS, prefix::String) end

@for_petsc function TSAppendOptionsPrefix(petsclib::$UnionPetscLib, ts::TS, prefix::String )

    @chk ccall(
               (:TSAppendOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{Cchar}),
               ts, prefix,
              )


	return nothing
end 

"""
	TSGetOptionsPrefix(petsclib::PetscLibType,ts::TS, prefix::String) 
Sets the prefix used for searching for all
`TS` options in the database.

Not Collective

Input Parameter:
- `ts` - The `TS` context

Output Parameter:
- `prefix` - A pointer to the prefix string used

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSAppendOptionsPrefix()`, `TSSetFromOptions()`

# External Links
$(_doc_external("Ts/TSGetOptionsPrefix"))
"""
function TSGetOptionsPrefix(petsclib::PetscLibType, ts::TS, prefix::String) end

@for_petsc function TSGetOptionsPrefix(petsclib::$UnionPetscLib, ts::TS, prefix::String )
	prefix_ = Ref(pointer(prefix))

    @chk ccall(
               (:TSGetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{Ptr{Cchar}}),
               ts, prefix_,
              )


	return nothing
end 

"""
	TSGetRHSJacobian(petsclib::PetscLibType,ts::TS, Amat::PetscMat, Pmat::PetscMat, func::TSRHSJacobianFn, ctx::Cvoid) 
Returns the Jacobian J at the present timestep.

Not Collective, but parallel objects are returned if ts is parallel

Input Parameter:
- `ts` - The `TS` context obtained from `TSCreate()`

Output Parameters:
- `Amat` - The (approximate) Jacobian J of G, where U_t = G(U,t)  (or `NULL`)
- `Pmat` - The matrix from which the preconditioner is constructed, usually the same as `Amat`  (or `NULL`)
- `func` - Function to compute the Jacobian of the RHS  (or `NULL`)
- `ctx`  - User-defined context for Jacobian evaluation routine  (or `NULL`)

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSGetTimeStep()`, `TSGetMatrices()`, `TSGetTime()`, `TSGetStepNumber()`


# External Links
$(_doc_external("Ts/TSGetRHSJacobian"))
"""
function TSGetRHSJacobian(petsclib::PetscLibType, ts::TS, Amat::PetscMat, Pmat::PetscMat, func::TSRHSJacobianFn, ctx::Cvoid) end

@for_petsc function TSGetRHSJacobian(petsclib::$UnionPetscLib, ts::TS, Amat::PetscMat, Pmat::PetscMat, func::TSRHSJacobianFn, ctx::Cvoid )
	Amat_ = Ref(Amat.ptr)
	Pmat_ = Ref(Pmat.ptr)

    @chk ccall(
               (:TSGetRHSJacobian, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{CMat}, Ptr{CMat}, TSRHSJacobianFn, Cvoid),
               ts, Amat_, Pmat_, func, ctx,
              )

	Amat.ptr = C_NULL
	Pmat.ptr = C_NULL

	return nothing
end 

"""
	TSGetIJacobian(petsclib::PetscLibType,ts::TS, Amat::PetscMat, Pmat::PetscMat, f::TSIJacobianFn, ctx::Cvoid) 
Returns the implicit Jacobian at the present timestep.

Not Collective, but parallel objects are returned if ts is parallel

Input Parameter:
- `ts` - The `TS` context obtained from `TSCreate()`

Output Parameters:
- `Amat` - The (approximate) Jacobian of F(t,U,U_t)
- `Pmat` - The matrix from which the preconditioner is constructed, often the same as `Amat`
- `f`    - The function to compute the matrices
- `ctx`  - User-defined context for Jacobian evaluation routine

Level: advanced

-seealso: [](ch_ts), `TS`, `TSGetTimeStep()`, `TSGetRHSJacobian()`, `TSGetMatrices()`, `TSGetTime()`, `TSGetStepNumber()`

# External Links
$(_doc_external("Ts/TSGetIJacobian"))
"""
function TSGetIJacobian(petsclib::PetscLibType, ts::TS, Amat::PetscMat, Pmat::PetscMat, f::TSIJacobianFn, ctx::Cvoid) end

@for_petsc function TSGetIJacobian(petsclib::$UnionPetscLib, ts::TS, Amat::PetscMat, Pmat::PetscMat, f::TSIJacobianFn, ctx::Cvoid )
	Amat_ = Ref(Amat.ptr)
	Pmat_ = Ref(Pmat.ptr)

    @chk ccall(
               (:TSGetIJacobian, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{CMat}, Ptr{CMat}, TSIJacobianFn, Cvoid),
               ts, Amat_, Pmat_, f, ctx,
              )

	Amat.ptr = C_NULL
	Pmat.ptr = C_NULL

	return nothing
end 

"""
	TSSetDM(petsclib::PetscLibType,ts::TS, dm::PetscDM) 
Sets the `DM` that may be used by some nonlinear solvers or preconditioners under the `TS`

Logically Collective

Input Parameters:
- `ts` - the `TS` integrator object
- `dm` - the dm, cannot be `NULL`

Level: intermediate

-seealso: [](ch_ts), `TS`, `DM`, `TSGetDM()`, `SNESSetDM()`, `SNESGetDM()`

# External Links
$(_doc_external("Ts/TSSetDM"))
"""
function TSSetDM(petsclib::PetscLibType, ts::TS, dm::PetscDM) end

@for_petsc function TSSetDM(petsclib::$UnionPetscLib, ts::TS, dm::PetscDM )

    @chk ccall(
               (:TSSetDM, $petsc_library),
               PetscErrorCode,
               (CTS, CDM),
               ts, dm,
              )


	return nothing
end 

"""
	dm::PetscDM = TSGetDM(petsclib::PetscLibType,ts::TS, dm::PetscDM) 
Gets the `DM` that may be used by some preconditioners

Not Collective

Input Parameter:
- `ts` - the `TS`

Output Parameter:
- `dm` - the `DM`

Level: intermediate

-seealso: [](ch_ts), `TS`, `DM`, `TSSetDM()`, `SNESSetDM()`, `SNESGetDM()`

# External Links
$(_doc_external("Ts/TSGetDM"))
"""
function TSGetDM(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSGetDM(petsclib::$UnionPetscLib, ts::TS )
	dm_ = Ref{CDM}()

    @chk ccall(
               (:TSGetDM, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{CDM}),
               ts, dm_,
              )

    dm = PetscDM(dm_[], petsclib)
	return dm
end 

"""
	TSComputeRHSFunctionLinear(petsclib::PetscLibType,ts::TS, t::PetscReal, U::PetscVec, F::PetscVec, ctx::Cvoid) 
Evaluate the right

Collective

Input Parameters:
- `ts`  - time stepping context
- `t`   - time at which to evaluate
- `U`   - state at which to evaluate
- `ctx` - context

Output Parameter:
- `F` - right-hand side

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSSetRHSFunction()`, `TSSetRHSJacobian()`, `TSComputeRHSJacobianConstant()`

# External Links
$(_doc_external("Ts/TSComputeRHSFunctionLinear"))
"""
function TSComputeRHSFunctionLinear(petsclib::PetscLibType, ts::TS, t::PetscReal, U::PetscVec, F::PetscVec, ctx::Cvoid) end

@for_petsc function TSComputeRHSFunctionLinear(petsclib::$UnionPetscLib, ts::TS, t::$PetscReal, U::PetscVec, F::PetscVec, ctx::Cvoid )

    @chk ccall(
               (:TSComputeRHSFunctionLinear, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal, CVec, CVec, Ptr{Cvoid}),
               ts, t, U, F, ctx,
              )


	return nothing
end 

"""
	TSComputeRHSJacobianConstant(petsclib::PetscLibType,ts::TS, t::PetscReal, U::PetscVec, A::PetscMat, B::PetscMat, ctx::Cvoid) 
Reuses a Jacobian that is time

Collective

Input Parameters:
- `ts`  - time stepping context
- `t`   - time at which to evaluate
- `U`   - state at which to evaluate
- `ctx` - context

Output Parameters:
- `A` - Jacobian
- `B` - matrix used to construct the preconditioner, often the same as `A`

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSSetRHSFunction()`, `TSSetRHSJacobian()`, `TSComputeRHSFunctionLinear()`

# External Links
$(_doc_external("Ts/TSComputeRHSJacobianConstant"))
"""
function TSComputeRHSJacobianConstant(petsclib::PetscLibType, ts::TS, t::PetscReal, U::PetscVec, A::PetscMat, B::PetscMat, ctx::Cvoid) end

@for_petsc function TSComputeRHSJacobianConstant(petsclib::$UnionPetscLib, ts::TS, t::$PetscReal, U::PetscVec, A::PetscMat, B::PetscMat, ctx::Cvoid )

    @chk ccall(
               (:TSComputeRHSJacobianConstant, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal, CVec, CMat, CMat, Ptr{Cvoid}),
               ts, t, U, A, B, ctx,
              )


	return nothing
end 

"""
	TSComputeIFunctionLinear(petsclib::PetscLibType,ts::TS, t::PetscReal, U::PetscVec, Udot::PetscVec, F::PetscVec, ctx::Cvoid) 
Evaluate the left hand side via the user

Collective

Input Parameters:
- `ts`   - time stepping context
- `t`    - time at which to evaluate
- `U`    - state at which to evaluate
- `Udot` - time derivative of state vector
- `ctx`  - context

Output Parameter:
- `F` - left hand side

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSSetIFunction()`, `TSSetIJacobian()`, `TSComputeIJacobianConstant()`, `TSComputeRHSFunctionLinear()`

# External Links
$(_doc_external("Ts/TSComputeIFunctionLinear"))
"""
function TSComputeIFunctionLinear(petsclib::PetscLibType, ts::TS, t::PetscReal, U::PetscVec, Udot::PetscVec, F::PetscVec, ctx::Cvoid) end

@for_petsc function TSComputeIFunctionLinear(petsclib::$UnionPetscLib, ts::TS, t::$PetscReal, U::PetscVec, Udot::PetscVec, F::PetscVec, ctx::Cvoid )

    @chk ccall(
               (:TSComputeIFunctionLinear, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal, CVec, CVec, CVec, Ptr{Cvoid}),
               ts, t, U, Udot, F, ctx,
              )


	return nothing
end 

"""
	TSComputeIJacobianConstant(petsclib::PetscLibType,ts::TS, t::PetscReal, U::PetscVec, Udot::PetscVec, shift::PetscReal, A::PetscMat, B::PetscMat, ctx::Cvoid) 
Reuses the matrix previously computed with the provided `TSIJacobianFn` for a semi

Collective

Input Parameters:
- `ts`    - time stepping context
- `t`     - time at which to evaluate
- `U`     - state at which to evaluate
- `Udot`  - time derivative of state vector
- `shift` - shift to apply
- `ctx`   - context

Output Parameters:
- `A` - pointer to operator
- `B` - pointer to matrix from which the preconditioner is built (often `A`)

Level: advanced

-seealso: [](ch_ts), `TS`, `TSROSW`, `TSARKIMEX`, `TSSetIFunction()`, `TSSetIJacobian()`, `TSComputeIFunctionLinear()`

# External Links
$(_doc_external("Ts/TSComputeIJacobianConstant"))
"""
function TSComputeIJacobianConstant(petsclib::PetscLibType, ts::TS, t::PetscReal, U::PetscVec, Udot::PetscVec, shift::PetscReal, A::PetscMat, B::PetscMat, ctx::Cvoid) end

@for_petsc function TSComputeIJacobianConstant(petsclib::$UnionPetscLib, ts::TS, t::$PetscReal, U::PetscVec, Udot::PetscVec, shift::$PetscReal, A::PetscMat, B::PetscMat, ctx::Cvoid )

    @chk ccall(
               (:TSComputeIJacobianConstant, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal, CVec, CVec, $PetscReal, CMat, CMat, Ptr{Cvoid}),
               ts, t, U, Udot, shift, A, B, ctx,
              )


	return nothing
end 

"""
	equation_type::TSEquationType = TSGetEquationType(petsclib::PetscLibType,ts::TS) 
Gets the type of the equation that `TS` is solving.

Not Collective

Input Parameter:
- `ts` - the `TS` context

Output Parameter:
- `equation_type` - see `TSEquationType`

Level: beginner

-seealso: [](ch_ts), `TS`, `TSSetEquationType()`, `TSEquationType`

# External Links
$(_doc_external("Ts/TSGetEquationType"))
"""
function TSGetEquationType(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSGetEquationType(petsclib::$UnionPetscLib, ts::TS )
	equation_type_ = Ref{TSEquationType}()

    @chk ccall(
               (:TSGetEquationType, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{TSEquationType}),
               ts, equation_type_,
              )

	equation_type = unsafe_string(equation_type_[])

	return equation_type
end 

"""
	TSSetEquationType(petsclib::PetscLibType,ts::TS, equation_type::TSEquationType) 
Sets the type of the equation that `TS` is solving.

Not Collective

Input Parameters:
- `ts`            - the `TS` context
- `equation_type` - see `TSEquationType`

Level: advanced

-seealso: [](ch_ts), `TS`, `TSGetEquationType()`, `TSEquationType`

# External Links
$(_doc_external("Ts/TSSetEquationType"))
"""
function TSSetEquationType(petsclib::PetscLibType, ts::TS, equation_type::TSEquationType) end

@for_petsc function TSSetEquationType(petsclib::$UnionPetscLib, ts::TS, equation_type::TSEquationType )

    @chk ccall(
               (:TSSetEquationType, $petsc_library),
               PetscErrorCode,
               (CTS, TSEquationType),
               ts, equation_type,
              )


	return nothing
end 

"""
	TSGetConvergedReason(petsclib::PetscLibType,ts::TS, reason::TSConvergedReason) 
Gets the reason the `TS` iteration was stopped.

Not Collective

Input Parameter:
- `ts` - the `TS` context

Output Parameter:
- `reason` - negative value indicates diverged, positive value converged, see `TSConvergedReason` or the
manual pages for the individual convergence tests for complete lists

Level: beginner

-seealso: [](ch_ts), `TS`, `TSSolve()`, `TSConvergedReason`

# External Links
$(_doc_external("Ts/TSGetConvergedReason"))
"""
function TSGetConvergedReason(petsclib::PetscLibType, ts::TS, reason::TSConvergedReason) end

@for_petsc function TSGetConvergedReason(petsclib::$UnionPetscLib, ts::TS, reason::TSConvergedReason )

    @chk ccall(
               (:TSGetConvergedReason, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{TSConvergedReason}),
               ts, reason,
              )


	return nothing
end 

"""
	TSSetConvergedReason(petsclib::PetscLibType,ts::TS, reason::TSConvergedReason) 
Sets the reason for handling the convergence of `TSSolve()`.

Logically Collective; reason must contain common value

Input Parameters:
- `ts`     - the `TS` context
- `reason` - negative value indicates diverged, positive value converged, see `TSConvergedReason` or the
manual pages for the individual convergence tests for complete lists

Level: advanced

-seealso: [](ch_ts), `TS`, `TSSolve()`, `TSConvergedReason`

# External Links
$(_doc_external("Ts/TSSetConvergedReason"))
"""
function TSSetConvergedReason(petsclib::PetscLibType, ts::TS, reason::TSConvergedReason) end

@for_petsc function TSSetConvergedReason(petsclib::$UnionPetscLib, ts::TS, reason::TSConvergedReason )

    @chk ccall(
               (:TSSetConvergedReason, $petsc_library),
               PetscErrorCode,
               (CTS, TSConvergedReason),
               ts, reason,
              )


	return nothing
end 

"""
	ftime::PetscReal = TSGetSolveTime(petsclib::PetscLibType,ts::TS) 
Gets the time after a call to `TSSolve()`

Not Collective

Input Parameter:
- `ts` - the `TS` context

Output Parameter:
- `ftime` - the final time. This time corresponds to the final time set with `TSSetMaxTime()`

Level: beginner

-seealso: [](ch_ts), `TS`, `TSSolve()`, `TSConvergedReason`

# External Links
$(_doc_external("Ts/TSGetSolveTime"))
"""
function TSGetSolveTime(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSGetSolveTime(petsclib::$UnionPetscLib, ts::TS )
	ftime_ = Ref{$PetscReal}()

    @chk ccall(
               (:TSGetSolveTime, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{$PetscReal}),
               ts, ftime_,
              )

	ftime = ftime_[]

	return ftime
end 

"""
	nits::PetscInt = TSGetSNESIterations(petsclib::PetscLibType,ts::TS) 
Gets the total number of nonlinear iterations
used by the time integrator.

Not Collective

Input Parameter:
- `ts` - `TS` context

Output Parameter:
- `nits` - number of nonlinear iterations

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSSolve()`, `TSGetKSPIterations()`

# External Links
$(_doc_external("Ts/TSGetSNESIterations"))
"""
function TSGetSNESIterations(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSGetSNESIterations(petsclib::$UnionPetscLib, ts::TS )
	nits_ = Ref{$PetscInt}()

    @chk ccall(
               (:TSGetSNESIterations, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{$PetscInt}),
               ts, nits_,
              )

	nits = nits_[]

	return nits
end 

"""
	lits::PetscInt = TSGetKSPIterations(petsclib::PetscLibType,ts::TS) 
Gets the total number of linear iterations
used by the time integrator.

Not Collective

Input Parameter:
- `ts` - `TS` context

Output Parameter:
- `lits` - number of linear iterations

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSSolve()`, `TSGetSNESIterations()`

# External Links
$(_doc_external("Ts/TSGetKSPIterations"))
"""
function TSGetKSPIterations(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSGetKSPIterations(petsclib::$UnionPetscLib, ts::TS )
	lits_ = Ref{$PetscInt}()

    @chk ccall(
               (:TSGetKSPIterations, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{$PetscInt}),
               ts, lits_,
              )

	lits = lits_[]

	return lits
end 

"""
	rejects::PetscInt = TSGetStepRejections(petsclib::PetscLibType,ts::TS) 
Gets the total number of rejected steps.

Not Collective

Input Parameter:
- `ts` - `TS` context

Output Parameter:
- `rejects` - number of steps rejected

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSSolve()`, `TSGetSNESIterations()`, `TSGetKSPIterations()`, `TSSetMaxStepRejections()`, `TSGetSNESFailures()`, `TSSetMaxSNESFailures()`, `TSSetErrorIfStepFails()`

# External Links
$(_doc_external("Ts/TSGetStepRejections"))
"""
function TSGetStepRejections(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSGetStepRejections(petsclib::$UnionPetscLib, ts::TS )
	rejects_ = Ref{$PetscInt}()

    @chk ccall(
               (:TSGetStepRejections, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{$PetscInt}),
               ts, rejects_,
              )

	rejects = rejects_[]

	return rejects
end 

"""
	fails::PetscInt = TSGetSNESFailures(petsclib::PetscLibType,ts::TS) 
Gets the total number of failed `SNES` solves in a `TS`

Not Collective

Input Parameter:
- `ts` - `TS` context

Output Parameter:
- `fails` - number of failed nonlinear solves

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSSolve()`, `TSGetSNESIterations()`, `TSGetKSPIterations()`, `TSSetMaxStepRejections()`, `TSGetStepRejections()`, `TSSetMaxSNESFailures()`

# External Links
$(_doc_external("Ts/TSGetSNESFailures"))
"""
function TSGetSNESFailures(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSGetSNESFailures(petsclib::$UnionPetscLib, ts::TS )
	fails_ = Ref{$PetscInt}()

    @chk ccall(
               (:TSGetSNESFailures, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{$PetscInt}),
               ts, fails_,
              )

	fails = fails_[]

	return fails
end 

"""
	TSSetMaxStepRejections(petsclib::PetscLibType,ts::TS, rejects::PetscInt) 
Sets the maximum number of step rejections before a time step fails

Not Collective

Input Parameters:
- `ts`      - `TS` context
- `rejects` - maximum number of rejected steps, pass `PETSC_UNLIMITED` for unlimited

Options Database Key:
- `-ts_max_reject` - Maximum number of step rejections before a step fails

Level: intermediate

-seealso: [](ch_ts), `TS`, `SNES`, `TSGetSNESIterations()`, `TSGetKSPIterations()`, `TSSetMaxSNESFailures()`, `TSGetStepRejections()`, `TSGetSNESFailures()`, `TSSetErrorIfStepFails()`, `TSGetConvergedReason()`

# External Links
$(_doc_external("Ts/TSSetMaxStepRejections"))
"""
function TSSetMaxStepRejections(petsclib::PetscLibType, ts::TS, rejects::PetscInt) end

@for_petsc function TSSetMaxStepRejections(petsclib::$UnionPetscLib, ts::TS, rejects::$PetscInt )

    @chk ccall(
               (:TSSetMaxStepRejections, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt),
               ts, rejects,
              )


	return nothing
end 

"""
	TSSetMaxSNESFailures(petsclib::PetscLibType,ts::TS, fails::PetscInt) 
Sets the maximum number of failed `SNES` solves

Not Collective

Input Parameters:
- `ts`    - `TS` context
- `fails` - maximum number of failed nonlinear solves, pass `PETSC_UNLIMITED` to allow any number of failures.

Options Database Key:
- `-ts_max_snes_failures` - Maximum number of nonlinear solve failures

Level: intermediate

-seealso: [](ch_ts), `TS`, `SNES`, `TSGetSNESIterations()`, `TSGetKSPIterations()`, `TSSetMaxStepRejections()`, `TSGetStepRejections()`, `TSGetSNESFailures()`, `SNESGetConvergedReason()`, `TSGetConvergedReason()`

# External Links
$(_doc_external("Ts/TSSetMaxSNESFailures"))
"""
function TSSetMaxSNESFailures(petsclib::PetscLibType, ts::TS, fails::PetscInt) end

@for_petsc function TSSetMaxSNESFailures(petsclib::$UnionPetscLib, ts::TS, fails::$PetscInt )

    @chk ccall(
               (:TSSetMaxSNESFailures, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt),
               ts, fails,
              )


	return nothing
end 

"""
	TSSetErrorIfStepFails(petsclib::PetscLibType,ts::TS, err::PetscBool) 
Immediately error if no step succeeds during `TSSolve()`

Not Collective

Input Parameters:
- `ts`  - `TS` context
- `err` - `PETSC_TRUE` to error if no step succeeds, `PETSC_FALSE` to return without failure

Options Database Key:
- `-ts_error_if_step_fails` - Error if no step succeeds

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSGetSNESIterations()`, `TSGetKSPIterations()`, `TSSetMaxStepRejections()`, `TSGetStepRejections()`, `TSGetSNESFailures()`, `TSGetConvergedReason()`

# External Links
$(_doc_external("Ts/TSSetErrorIfStepFails"))
"""
function TSSetErrorIfStepFails(petsclib::PetscLibType, ts::TS, err::PetscBool) end

@for_petsc function TSSetErrorIfStepFails(petsclib::$UnionPetscLib, ts::TS, err::PetscBool )

    @chk ccall(
               (:TSSetErrorIfStepFails, $petsc_library),
               PetscErrorCode,
               (CTS, PetscBool),
               ts, err,
              )


	return nothing
end 

"""
	TSGetAdapt(petsclib::PetscLibType,ts::TS, adapt::TSAdapt) 
Get the adaptive controller context for the current method

Collective if controller has not yet been created

Input Parameter:
- `ts` - time stepping context

Output Parameter:
- `adapt` - adaptive controller

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSAdapt`, `TSAdaptSetType()`, `TSAdaptChoose()`

# External Links
$(_doc_external("Ts/TSGetAdapt"))
"""
function TSGetAdapt(petsclib::PetscLibType, ts::TS, adapt::TSAdapt) end

@for_petsc function TSGetAdapt(petsclib::$UnionPetscLib, ts::TS, adapt::TSAdapt )

    @chk ccall(
               (:TSGetAdapt, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{TSAdapt}),
               ts, adapt,
              )


	return nothing
end 

"""
	TSSetTolerances(petsclib::PetscLibType,ts::TS, atol::PetscReal, vatol::PetscVec, rtol::PetscReal, vrtol::PetscVec) 
Set tolerances for local truncation error when using an adaptive controller

Logically Collective

Input Parameters:
- `ts`    - time integration context
- `atol`  - scalar absolute tolerances
- `vatol` - vector of absolute tolerances or `NULL`, used in preference to `atol` if present
- `rtol`  - scalar relative tolerances
- `vrtol` - vector of relative tolerances or `NULL`, used in preference to `rtol` if present

Options Database Keys:
- `-ts_rtol <rtol>` - relative tolerance for local truncation error
- `-ts_atol <atol>` - Absolute tolerance for local truncation error

Level: beginner

-seealso: [](ch_ts), `TS`, `TSAdapt`, `TSErrorWeightedNorm()`, `TSGetTolerances()`

# External Links
$(_doc_external("Ts/TSSetTolerances"))
"""
function TSSetTolerances(petsclib::PetscLibType, ts::TS, atol::PetscReal, vatol::PetscVec, rtol::PetscReal, vrtol::PetscVec) end

@for_petsc function TSSetTolerances(petsclib::$UnionPetscLib, ts::TS, atol::$PetscReal, vatol::PetscVec, rtol::$PetscReal, vrtol::PetscVec )

    @chk ccall(
               (:TSSetTolerances, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal, CVec, $PetscReal, CVec),
               ts, atol, vatol, rtol, vrtol,
              )


	return nothing
end 

"""
	atol::PetscReal,rtol::PetscReal = TSGetTolerances(petsclib::PetscLibType,ts::TS, vatol::PetscVec, vrtol::PetscVec) 
Get tolerances for local truncation error when using adaptive controller

Logically Collective

Input Parameter:
- `ts` - time integration context

Output Parameters:
- `atol`  - scalar absolute tolerances, `NULL` to ignore
- `vatol` - vector of absolute tolerances, `NULL` to ignore
- `rtol`  - scalar relative tolerances, `NULL` to ignore
- `vrtol` - vector of relative tolerances, `NULL` to ignore

Level: beginner

-seealso: [](ch_ts), `TS`, `TSAdapt`, `TSErrorWeightedNorm()`, `TSSetTolerances()`

# External Links
$(_doc_external("Ts/TSGetTolerances"))
"""
function TSGetTolerances(petsclib::PetscLibType, ts::TS, vatol::PetscVec, vrtol::PetscVec) end

@for_petsc function TSGetTolerances(petsclib::$UnionPetscLib, ts::TS, vatol::PetscVec, vrtol::PetscVec )
	atol_ = Ref{$PetscReal}()
	vatol_ = Ref(vatol.ptr)
	rtol_ = Ref{$PetscReal}()
	vrtol_ = Ref(vrtol.ptr)

    @chk ccall(
               (:TSGetTolerances, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{$PetscReal}, Ptr{CVec}, Ptr{$PetscReal}, Ptr{CVec}),
               ts, atol_, vatol_, rtol_, vrtol_,
              )

	atol = atol_[]
	vatol.ptr = C_NULL
	rtol = rtol_[]
	vrtol.ptr = C_NULL

	return atol,rtol
end 

"""
	norm::PetscReal,norma::PetscReal,normr::PetscReal = TSErrorWeightedNorm(petsclib::PetscLibType,ts::TS, U::PetscVec, Y::PetscVec, wnormtype::NormType) 
compute a weighted norm of the difference between two state vectors based on supplied absolute and relative tolerances

Collective

Input Parameters:
- `ts`        - time stepping context
- `U`         - state vector, usually ts->vec_sol
- `Y`         - state vector to be compared to U
- `wnormtype` - norm type, either `NORM_2` or `NORM_INFINITY`

Output Parameters:
- `norm`  - weighted norm, a value of 1.0 achieves a balance between absolute and relative tolerances
- `norma` - weighted norm, a value of 1.0 means that the error meets the absolute tolerance set by the user
- `normr` - weighted norm, a value of 1.0 means that the error meets the relative tolerance set by the user

Options Database Key:
- `-ts_adapt_wnormtype <wnormtype>` - 2, INFINITY

Level: developer

-seealso: [](ch_ts), `TS`, `VecErrorWeightedNorms()`, `TSErrorWeightedENorm()`

# External Links
$(_doc_external("Ts/TSErrorWeightedNorm"))
"""
function TSErrorWeightedNorm(petsclib::PetscLibType, ts::TS, U::PetscVec, Y::PetscVec, wnormtype::NormType) end

@for_petsc function TSErrorWeightedNorm(petsclib::$UnionPetscLib, ts::TS, U::PetscVec, Y::PetscVec, wnormtype::NormType )
	norm_ = Ref{$PetscReal}()
	norma_ = Ref{$PetscReal}()
	normr_ = Ref{$PetscReal}()

    @chk ccall(
               (:TSErrorWeightedNorm, $petsc_library),
               PetscErrorCode,
               (CTS, CVec, CVec, NormType, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
               ts, U, Y, wnormtype, norm_, norma_, normr_,
              )

	norm = norm_[]
	norma = norma_[]
	normr = normr_[]

	return norm,norma,normr
end 

"""
	norm::PetscReal,norma::PetscReal,normr::PetscReal = TSErrorWeightedENorm(petsclib::PetscLibType,ts::TS, E::PetscVec, U::PetscVec, Y::PetscVec, wnormtype::NormType) 
compute a weighted error norm based on supplied absolute and relative tolerances

Collective

Input Parameters:
- `ts`        - time stepping context
- `E`         - error vector
- `U`         - state vector, usually ts->vec_sol
- `Y`         - state vector, previous time step
- `wnormtype` - norm type, either `NORM_2` or `NORM_INFINITY`

Output Parameters:
- `norm`  - weighted norm, a value of 1.0 achieves a balance between absolute and relative tolerances
- `norma` - weighted norm, a value of 1.0 means that the error meets the absolute tolerance set by the user
- `normr` - weighted norm, a value of 1.0 means that the error meets the relative tolerance set by the user

Options Database Key:
- `-ts_adapt_wnormtype <wnormtype>` - 2, INFINITY

Level: developer

-seealso: [](ch_ts), `TS`, `VecErrorWeightedNorms()`, `TSErrorWeightedNorm()`

# External Links
$(_doc_external("Ts/TSErrorWeightedENorm"))
"""
function TSErrorWeightedENorm(petsclib::PetscLibType, ts::TS, E::PetscVec, U::PetscVec, Y::PetscVec, wnormtype::NormType) end

@for_petsc function TSErrorWeightedENorm(petsclib::$UnionPetscLib, ts::TS, E::PetscVec, U::PetscVec, Y::PetscVec, wnormtype::NormType )
	norm_ = Ref{$PetscReal}()
	norma_ = Ref{$PetscReal}()
	normr_ = Ref{$PetscReal}()

    @chk ccall(
               (:TSErrorWeightedENorm, $petsc_library),
               PetscErrorCode,
               (CTS, CVec, CVec, CVec, NormType, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
               ts, E, U, Y, wnormtype, norm_, norma_, normr_,
              )

	norm = norm_[]
	norma = norma_[]
	normr = normr_[]

	return norm,norma,normr
end 

"""
	TSSetCFLTimeLocal(petsclib::PetscLibType,ts::TS, cfltime::PetscReal) 
Set the local CFL constraint relative to forward Euler

Logically Collective

Input Parameters:
- `ts`      - time stepping context
- `cfltime` - maximum stable time step if using forward Euler (value can be different on each process)

-seealso: [](ch_ts), `TSGetCFLTime()`, `TSADAPTCFL`

# External Links
$(_doc_external("Ts/TSSetCFLTimeLocal"))
"""
function TSSetCFLTimeLocal(petsclib::PetscLibType, ts::TS, cfltime::PetscReal) end

@for_petsc function TSSetCFLTimeLocal(petsclib::$UnionPetscLib, ts::TS, cfltime::$PetscReal )

    @chk ccall(
               (:TSSetCFLTimeLocal, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal),
               ts, cfltime,
              )


	return nothing
end 

"""
	cfltime::PetscReal = TSGetCFLTime(petsclib::PetscLibType,ts::TS) 
Get the maximum stable time step according to CFL criteria applied to forward Euler

Collective

Input Parameter:
- `ts` - time stepping context

Output Parameter:
- `cfltime` - maximum stable time step for forward Euler

Level: advanced

-seealso: [](ch_ts), `TSSetCFLTimeLocal()`

# External Links
$(_doc_external("Ts/TSGetCFLTime"))
"""
function TSGetCFLTime(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSGetCFLTime(petsclib::$UnionPetscLib, ts::TS )
	cfltime_ = Ref{$PetscReal}()

    @chk ccall(
               (:TSGetCFLTime, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{$PetscReal}),
               ts, cfltime_,
              )

	cfltime = cfltime_[]

	return cfltime
end 

"""
	TSVISetVariableBounds(petsclib::PetscLibType,ts::TS, xl::PetscVec, xu::PetscVec) 
Sets the lower and upper bounds for the solution vector. xl <= x <= xu

Input Parameters:
- `ts` - the `TS` context.
- `xl` - lower bound.
- `xu` - upper bound.

Level: advanced

-seealso: [](ch_ts), `TS`

# External Links
$(_doc_external("Ts/TSVISetVariableBounds"))
"""
function TSVISetVariableBounds(petsclib::PetscLibType, ts::TS, xl::PetscVec, xu::PetscVec) end

@for_petsc function TSVISetVariableBounds(petsclib::$UnionPetscLib, ts::TS, xl::PetscVec, xu::PetscVec )

    @chk ccall(
               (:TSVISetVariableBounds, $petsc_library),
               PetscErrorCode,
               (CTS, CVec, CVec),
               ts, xl, xu,
              )


	return nothing
end 

"""
	yr::PetscReal,yi::PetscReal = TSComputeLinearStability(petsclib::PetscLibType,ts::TS, xr::PetscReal, xi::PetscReal) 
computes the linear stability function at a point

Collective

Input Parameters:
- `ts` - the `TS` context
- `xr` - real part of input argument
- `xi` - imaginary part of input argument

Output Parameters:
- `yr` - real part of function value
- `yi` - imaginary part of function value

Level: developer

-seealso: [](ch_ts), `TS`, `TSSetRHSFunction()`, `TSComputeIFunction()`

# External Links
$(_doc_external("Ts/TSComputeLinearStability"))
"""
function TSComputeLinearStability(petsclib::PetscLibType, ts::TS, xr::PetscReal, xi::PetscReal) end

@for_petsc function TSComputeLinearStability(petsclib::$UnionPetscLib, ts::TS, xr::$PetscReal, xi::$PetscReal )
	yr_ = Ref{$PetscReal}()
	yi_ = Ref{$PetscReal}()

    @chk ccall(
               (:TSComputeLinearStability, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal, $PetscReal, Ptr{$PetscReal}, Ptr{$PetscReal}),
               ts, xr, xi, yr_, yi_,
              )

	yr = yr_[]
	yi = yi_[]

	return yr,yi
end 

"""
	TSRestartStep(petsclib::PetscLibType,ts::TS) 
Flags the solver to restart the next step

Collective

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Level: advanced

-seealso: [](ch_ts), `TS`, `TSBDF`, `TSSolve()`, `TSSetPreStep()`, `TSSetPostStep()`

# External Links
$(_doc_external("Ts/TSRestartStep"))
"""
function TSRestartStep(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSRestartStep(petsclib::$UnionPetscLib, ts::TS )

    @chk ccall(
               (:TSRestartStep, $petsc_library),
               PetscErrorCode,
               (CTS,),
               ts,
              )


	return nothing
end 

"""
	TSRollBack(petsclib::PetscLibType,ts::TS) 
Rolls back one time step

Collective

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Level: advanced

-seealso: [](ch_ts), `TS`, `TSGetStepRollBack()`, `TSCreate()`, `TSSetUp()`, `TSDestroy()`, `TSSolve()`, `TSSetPreStep()`, `TSSetPreStage()`, `TSInterpolate()`

# External Links
$(_doc_external("Ts/TSRollBack"))
"""
function TSRollBack(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSRollBack(petsclib::$UnionPetscLib, ts::TS )

    @chk ccall(
               (:TSRollBack, $petsc_library),
               PetscErrorCode,
               (CTS,),
               ts,
              )


	return nothing
end 

"""
	flg::PetscBool = TSGetStepRollBack(petsclib::PetscLibType,ts::TS) 
Get the internal flag indicating if you are rolling back a step

Not collective

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Output Parameter:
- `flg` - the rollback flag

Level: advanced

-seealso: [](ch_ts), `TS`, `TSCreate()`, `TSRollBack()`

# External Links
$(_doc_external("Ts/TSGetStepRollBack"))
"""
function TSGetStepRollBack(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSGetStepRollBack(petsclib::$UnionPetscLib, ts::TS )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:TSGetStepRollBack, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{PetscBool}),
               ts, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	flg::PetscBool = TSGetStepResize(petsclib::PetscLibType,ts::TS) 
Get the internal flag indicating if the current step is after a resize.

Not collective

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Output Parameter:
- `flg` - the resize flag

Level: advanced

-seealso: [](ch_ts), `TS`, `TSCreate()`, `TSSetResize()`

# External Links
$(_doc_external("Ts/TSGetStepResize"))
"""
function TSGetStepResize(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSGetStepResize(petsclib::$UnionPetscLib, ts::TS )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:TSGetStepResize, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{PetscBool}),
               ts, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	ns::PetscInt = TSGetStages(petsclib::PetscLibType,ts::TS, Y::PetscVec) 
Get the number of stages and stage values

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Output Parameters:
- `ns` - the number of stages
- `Y`  - the current stage vectors

Level: advanced

-seealso: [](ch_ts), `TS`, `TSCreate()`

# External Links
$(_doc_external("Ts/TSGetStages"))
"""
function TSGetStages(petsclib::PetscLibType, ts::TS, Y::PetscVec) end

@for_petsc function TSGetStages(petsclib::$UnionPetscLib, ts::TS, Y::PetscVec )
	ns_ = Ref{$PetscInt}()

    @chk ccall(
               (:TSGetStages, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{$PetscInt}, CVec),
               ts, ns_, Y,
              )

	ns = ns_[]

	return ns
end 

"""
	TSComputeIJacobianDefaultColor(petsclib::PetscLibType,ts::TS, t::PetscReal, U::PetscVec, Udot::PetscVec, shift::PetscReal, J::PetscMat, B::PetscMat, ctx::Cvoid) 
Computes the Jacobian using finite differences and coloring to exploit matrix sparsity.

Collective

Input Parameters:
- `ts`    - the `TS` context
- `t`     - current timestep
- `U`     - state vector
- `Udot`  - time derivative of state vector
- `shift` - shift to apply, see note below
- `ctx`   - an optional user context

Output Parameters:
- `J` - Jacobian matrix (not altered in this routine)
- `B` - newly computed Jacobian matrix to use with preconditioner (generally the same as `J`)

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSSetIJacobian()`, `MatFDColoringCreate()`, `MatFDColoringSetFunction()`

# External Links
$(_doc_external("Ts/TSComputeIJacobianDefaultColor"))
"""
function TSComputeIJacobianDefaultColor(petsclib::PetscLibType, ts::TS, t::PetscReal, U::PetscVec, Udot::PetscVec, shift::PetscReal, J::PetscMat, B::PetscMat, ctx::Cvoid) end

@for_petsc function TSComputeIJacobianDefaultColor(petsclib::$UnionPetscLib, ts::TS, t::$PetscReal, U::PetscVec, Udot::PetscVec, shift::$PetscReal, J::PetscMat, B::PetscMat, ctx::Cvoid )

    @chk ccall(
               (:TSComputeIJacobianDefaultColor, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal, CVec, CVec, $PetscReal, CMat, CMat, Ptr{Cvoid}),
               ts, t, U, Udot, shift, J, B, ctx,
              )


	return nothing
end 

"""
	TSSetFunctionDomainError(petsclib::PetscLibType,ts::TS, func::external) 
Set a function that tests if the current state vector is valid

Logically collective

Input Parameters:
- `ts`   - the `TS` context
- `func` - function called within `TSFunctionDomainError()`

Calling sequence of `func`:
- `ts`     - the `TS` context
- `time`   - the current time (of the stage)
- `state`  - the state to check if it is valid
- `accept` - (output parameter) `PETSC_FALSE` if the state is not acceptable, `PETSC_TRUE` if acceptable

Level: intermediate

-seealso: [](ch_ts), `TSAdaptCheckStage()`, `TSFunctionDomainError()`, `SNESSetFunctionDomainError()`, `TSGetSNES()`

# External Links
$(_doc_external("Ts/TSSetFunctionDomainError"))
"""
function TSSetFunctionDomainError(petsclib::PetscLibType, ts::TS, func::external) end

@for_petsc function TSSetFunctionDomainError(petsclib::$UnionPetscLib, ts::TS, func::external )

    @chk ccall(
               (:TSSetFunctionDomainError, $petsc_library),
               PetscErrorCode,
               (CTS, external),
               ts, func,
              )


	return nothing
end 

"""
	accept::PetscBool = TSFunctionDomainError(petsclib::PetscLibType,ts::TS, stagetime::PetscReal, Y::PetscVec) 
Checks if the current state is valid

Collective

Input Parameters:
- `ts`        - the `TS` context
- `stagetime` - time of the simulation
- `Y`         - state vector to check.

Output Parameter:
- `accept` - Set to `PETSC_FALSE` if the current state vector is valid.

Level: developer

-seealso: [](ch_ts), `TS`, `TSSetFunctionDomainError()`

# External Links
$(_doc_external("Ts/TSFunctionDomainError"))
"""
function TSFunctionDomainError(petsclib::PetscLibType, ts::TS, stagetime::PetscReal, Y::PetscVec) end

@for_petsc function TSFunctionDomainError(petsclib::$UnionPetscLib, ts::TS, stagetime::$PetscReal, Y::PetscVec )
	accept_ = Ref{PetscBool}()

    @chk ccall(
               (:TSFunctionDomainError, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal, CVec, Ptr{PetscBool}),
               ts, stagetime, Y, accept_,
              )

	accept = accept_[]

	return accept
end 

"""
	TSClone(petsclib::PetscLibType,tsin::TS, tsout::TS) 
This function clones a time step `TS` object.

Collective

Input Parameter:
- `tsin` - The input `TS`

Output Parameter:
- `tsout` - The output `TS` (cloned)

Level: developer

-seealso: [](ch_ts), `TS`, `SNES`, `TSCreate()`, `TSSetType()`, `TSSetUp()`, `TSDestroy()`, `TSSetProblemType()`

# External Links
$(_doc_external("Ts/TSClone"))
"""
function TSClone(petsclib::PetscLibType, tsin::TS, tsout::TS) end

@for_petsc function TSClone(petsclib::$UnionPetscLib, tsin::TS, tsout::TS )
	tsout_ = Ref(tsout.ptr)

    @chk ccall(
               (:TSClone, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{CTS}),
               tsin, tsout_,
              )

	tsout.ptr = C_NULL

	return nothing
end 

"""
	flg::PetscBool = TSRHSJacobianTest(petsclib::PetscLibType,ts::TS) 
Compares the multiply routine provided to the `MATSHELL` with differencing on the `TS` given RHS function.

Logically Collective

Input Parameter:
- `ts` - the time stepping routine

Output Parameter:
- `flg` - `PETSC_TRUE` if the multiply is likely correct

Options Database Key:
- `-ts_rhs_jacobian_test_mult -mat_shell_test_mult_view` - run the test at each timestep of the integrator

Level: advanced

-seealso: [](ch_ts), `TS`, `Mat`, `MATSHELL`, `MatCreateShell()`, `MatShellGetContext()`, `MatShellGetOperation()`, `MatShellTestMultTranspose()`, `TSRHSJacobianTestTranspose()`

# External Links
$(_doc_external("Ts/TSRHSJacobianTest"))
"""
function TSRHSJacobianTest(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSRHSJacobianTest(petsclib::$UnionPetscLib, ts::TS )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:TSRHSJacobianTest, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{PetscBool}),
               ts, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	flg::PetscBool = TSRHSJacobianTestTranspose(petsclib::PetscLibType,ts::TS) 
Compares the multiply transpose routine provided to the `MATSHELL` with differencing on the `TS` given RHS function.

Logically Collective

Input Parameter:
- `ts` - the time stepping routine

Output Parameter:
- `flg` - `PETSC_TRUE` if the multiply is likely correct

Options Database Key:
- `-ts_rhs_jacobian_test_mult_transpose -mat_shell_test_mult_transpose_view` - run the test at each timestep of the integrator

Level: advanced

-seealso: [](ch_ts), `TS`, `Mat`, `MatCreateShell()`, `MatShellGetContext()`, `MatShellGetOperation()`, `MatShellTestMultTranspose()`, `TSRHSJacobianTest()`

# External Links
$(_doc_external("Ts/TSRHSJacobianTestTranspose"))
"""
function TSRHSJacobianTestTranspose(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSRHSJacobianTestTranspose(petsclib::$UnionPetscLib, ts::TS )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:TSRHSJacobianTestTranspose, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{PetscBool}),
               ts, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	TSSetUseSplitRHSFunction(petsclib::PetscLibType,ts::TS, use_splitrhsfnc::PetscBool) 
Use the split RHSFunction when a multirate method is used.

Logically Collective

Input Parameters:
- `ts`                   - timestepping context
- `use_splitrhsfunction` - `PETSC_TRUE` indicates that the split RHSFunction will be used

Options Database Key:
- `-ts_use_splitrhsfunction` - <true,false>

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSGetUseSplitRHSFunction()`

# External Links
$(_doc_external("Ts/TSSetUseSplitRHSFunction"))
"""
function TSSetUseSplitRHSFunction(petsclib::PetscLibType, ts::TS, use_splitrhsfnc::PetscBool) end

@for_petsc function TSSetUseSplitRHSFunction(petsclib::$UnionPetscLib, ts::TS, use_splitrhsfnc::PetscBool )

    @chk ccall(
               (:TSSetUseSplitRHSFunction, $petsc_library),
               PetscErrorCode,
               (CTS, PetscBool),
               ts, use_splitrhsfnc,
              )


	return nothing
end 

"""
	use_splitrhsfnc::PetscBool = TSGetUseSplitRHSFunction(petsclib::PetscLibType,ts::TS) 
Gets whether to use the split RHSFunction when a multirate method is used.

Not Collective

Input Parameter:
- `ts` - timestepping context

Output Parameter:
- `use_splitrhsfunction` - `PETSC_TRUE` indicates that the split RHSFunction will be used

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSSetUseSplitRHSFunction()`

# External Links
$(_doc_external("Ts/TSGetUseSplitRHSFunction"))
"""
function TSGetUseSplitRHSFunction(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSGetUseSplitRHSFunction(petsclib::$UnionPetscLib, ts::TS )
	use_splitrhsfnc_ = Ref{PetscBool}()

    @chk ccall(
               (:TSGetUseSplitRHSFunction, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{PetscBool}),
               ts, use_splitrhsfnc_,
              )

	use_splitrhsfnc = use_splitrhsfnc_[]

	return use_splitrhsfnc
end 

"""
	TSSetMatStructure(petsclib::PetscLibType,ts::TS, str::MatStructure) 
sets the relationship between the nonzero structure of the RHS Jacobian matrix to the IJacobian matrix.

Logically  Collective

Input Parameters:
- `ts`  - the time-stepper
- `str` - the structure (the default is `UNKNOWN_NONZERO_PATTERN`)

Level: intermediate

-seealso: [](ch_ts), `TS`, `MatAXPY()`, `MatStructure`

# External Links
$(_doc_external("Ts/TSSetMatStructure"))
"""
function TSSetMatStructure(petsclib::PetscLibType, ts::TS, str::MatStructure) end

@for_petsc function TSSetMatStructure(petsclib::$UnionPetscLib, ts::TS, str::MatStructure )

    @chk ccall(
               (:TSSetMatStructure, $petsc_library),
               PetscErrorCode,
               (CTS, MatStructure),
               ts, str,
              )


	return nothing
end 

"""
	TSSetEvaluationTimes(petsclib::PetscLibType,ts::TS, n::PetscInt, time_points::Vector{PetscReal}) 
sets the evaluation points. The solution will be computed and stored for each time requested

Collective

Input Parameters:
- `ts`          - the time-stepper
- `n`           - number of the time points
- `time_points` - array of the time points, must be increasing

Options Database Key:
- `-ts_eval_times <t0,...tn>` - Sets the evaluation times

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSGetEvaluationTimes()`, `TSGetEvaluationSolutions()`, `TSSetTimeSpan()`

# External Links
$(_doc_external("Ts/TSSetEvaluationTimes"))
"""
function TSSetEvaluationTimes(petsclib::PetscLibType, ts::TS, n::PetscInt, time_points::Vector{PetscReal}) end

@for_petsc function TSSetEvaluationTimes(petsclib::$UnionPetscLib, ts::TS, n::$PetscInt, time_points::Vector{$PetscReal} )

    @chk ccall(
               (:TSSetEvaluationTimes, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt, Ptr{$PetscReal}),
               ts, n, time_points,
              )


	return nothing
end 

"""
	n::PetscInt,time_points::Vector{PetscReal} = TSGetEvaluationTimes(petsclib::PetscLibType,ts::TS) 
gets the evaluation times set with `TSSetEvaluationTimes()`

Not Collective

Input Parameter:
- `ts` - the time-stepper

Output Parameters:
- `n`           - number of the time points
- `time_points` - array of the time points

Level: beginner

-seealso: [](ch_ts), `TS`, `TSSetEvaluationTimes()`, `TSGetEvaluationSolutions()`

# External Links
$(_doc_external("Ts/TSGetEvaluationTimes"))
"""
function TSGetEvaluationTimes(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSGetEvaluationTimes(petsclib::$UnionPetscLib, ts::TS )
	n_ = Ref{$PetscInt}()
	time_points_ = Ref{Ptr{$PetscReal}}()

    @chk ccall(
               (:TSGetEvaluationTimes, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{$PetscInt}, Ptr{Ptr{$PetscReal}}),
               ts, n_, time_points_,
              )

	n = n_[]
	time_points = unsafe_wrap(Array, time_points_[], VecGetLocalSize(petsclib, x); own = false)

	return n,time_points
end 

"""
	nsol::PetscInt,sol_times::Vector{PetscReal} = TSGetEvaluationSolutions(petsclib::PetscLibType,ts::TS, Sols::Vector{PetscVec}) 
Get the number of solutions and the solutions at the evaluation time points specified

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Output Parameters:
- `nsol`      - the number of solutions
- `sol_times` - array of solution times corresponding to the solution vectors. See note below
- `Sols`      - the solution vectors

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSSetEvaluationTimes()`, `TSGetEvaluationTimes()`

# External Links
$(_doc_external("Ts/TSGetEvaluationSolutions"))
"""
function TSGetEvaluationSolutions(petsclib::PetscLibType, ts::TS, Sols::Vector{PetscVec}) end

@for_petsc function TSGetEvaluationSolutions(petsclib::$UnionPetscLib, ts::TS, Sols::Vector{PetscVec} )
	nsol_ = Ref{$PetscInt}()
	sol_times_ = Ref{Ptr{$PetscReal}}()
	Sols_ = Ref(pointer(Sols))

    @chk ccall(
               (:TSGetEvaluationSolutions, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{$PetscInt}, Ptr{Ptr{$PetscReal}}, Ptr{Ptr{CVec}}),
               ts, nsol_, sol_times_, Sols_,
              )

	nsol = nsol_[]
	sol_times = unsafe_wrap(Array, sol_times_[], VecGetLocalSize(petsclib, x); own = false)

	return nsol,sol_times
end 

"""
	TSSetTimeSpan(petsclib::PetscLibType,ts::TS, n::PetscInt, span_times::Vector{PetscReal}) 
sets the time span. The solution will be computed and stored for each time requested in the span

Collective

Input Parameters:
- `ts`         - the time-stepper
- `n`          - number of the time points (>=2)
- `span_times` - array of the time points, must be increasing. The first element and the last element are the initial time and the final time respectively.

Options Database Key:
- `-ts_time_span <t0,...tf>` - Sets the time span

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSSetEvaluationTimes()`, `TSGetEvaluationTimes()`, `TSGetEvaluationSolutions()`

# External Links
$(_doc_external("Ts/TSSetTimeSpan"))
"""
function TSSetTimeSpan(petsclib::PetscLibType, ts::TS, n::PetscInt, span_times::Vector{PetscReal}) end

@for_petsc function TSSetTimeSpan(petsclib::$UnionPetscLib, ts::TS, n::$PetscInt, span_times::Vector{$PetscReal} )

    @chk ccall(
               (:TSSetTimeSpan, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt, Ptr{$PetscReal}),
               ts, n, span_times,
              )


	return nothing
end 

"""
	TSPruneIJacobianColor(petsclib::PetscLibType,ts::TS, J::PetscMat, B::PetscMat) 
Remove nondiagonal zeros in the Jacobian matrix and update the `MatMFFD` coloring information.

Collective

Input Parameters:
- `ts` - the `TS` context
- `J`  - Jacobian matrix (not altered in this routine)
- `B`  - newly computed Jacobian matrix to use with preconditioner

Level: intermediate

-seealso: [](ch_ts), `TS`, `MatFDColoring`, `TSComputeIJacobianDefaultColor()`, `MatEliminateZeros()`, `MatFDColoringCreate()`, `MatFDColoringSetFunction()`

# External Links
$(_doc_external("Ts/TSPruneIJacobianColor"))
"""
function TSPruneIJacobianColor(petsclib::PetscLibType, ts::TS, J::PetscMat, B::PetscMat) end

@for_petsc function TSPruneIJacobianColor(petsclib::$UnionPetscLib, ts::TS, J::PetscMat, B::PetscMat )

    @chk ccall(
               (:TSPruneIJacobianColor, $petsc_library),
               PetscErrorCode,
               (CTS, CMat, CMat),
               ts, J, B,
              )


	return nothing
end 

"""
	TSSetType(petsclib::PetscLibType,ts::TS, type::TSType) 
Sets the algorithm/method to be used for integrating the ODE with the given `TS`.

Collective

Input Parameters:
- `ts`   - The `TS` context
- `type` - A known method

Options Database Key:
- `-ts_type <type>` - Sets the method; use -help for a list of available methods (for instance, euler)

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSSolve()`, `TSCreate()`, `TSSetFromOptions()`, `TSDestroy()`, `TSType`

# External Links
$(_doc_external("Ts/TSSetType"))
"""
function TSSetType(petsclib::PetscLibType, ts::TS, type::TSType) end

@for_petsc function TSSetType(petsclib::$UnionPetscLib, ts::TS, type::TSType )

    @chk ccall(
               (:TSSetType, $petsc_library),
               PetscErrorCode,
               (CTS, TSType),
               ts, type,
              )


	return nothing
end 

"""
	type::TSType = TSGetType(petsclib::PetscLibType,ts::TS) 
Gets the `TS` method type (as a string) that is being used to solve the ODE with the given `TS`

Not Collective

Input Parameter:
- `ts` - The `TS`

Output Parameter:
- `type` - The `TSType`

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSType`, `TSSetType()`

# External Links
$(_doc_external("Ts/TSGetType"))
"""
function TSGetType(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSGetType(petsclib::$UnionPetscLib, ts::TS )
	type_ = Ref{TSType}()

    @chk ccall(
               (:TSGetType, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{TSType}),
               ts, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	TSRegister(petsclib::PetscLibType,sname::String, fnc::external) 
Adds a creation method to the `TS` package.

Not Collective, No Fortran Support

Input Parameters:
- `sname`    - The name of a new user-defined creation routine
- `function` - The creation routine itself

Level: advanced

-seealso: [](ch_ts), `TSSetType()`, `TSType`, `TSRegisterAll()`, `TSRegisterDestroy()`

# External Links
$(_doc_external("Ts/TSRegister"))
"""
function TSRegister(petsclib::PetscLibType, sname::String, fnc::external) end

@for_petsc function TSRegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:TSRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	TSFinalizePackage(petsclib::PetscLibType) 
This function destroys everything in the PETSc interface to `TS`. It is
called from `PetscFinalize()`.

Level: developer

-seealso: [](ch_ts), `TS`, `PetscFinalize()`, `TSInitializePackage()`

# External Links
$(_doc_external("Ts/TSFinalizePackage"))
"""
function TSFinalizePackage(petsclib::PetscLibType) end

@for_petsc function TSFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:TSFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	TSInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the `TS` package. It is called
from `PetscDLLibraryRegister_petscts()` when using dynamic libraries, and on the first call to `TSCreate()`
when using shared or static libraries.

Level: developer

-seealso: [](ch_ts), `TS`, `PetscInitialize()`, `TSFinalizePackage()`

# External Links
$(_doc_external("Ts/TSInitializePackage"))
"""
function TSInitializePackage(petsclib::PetscLibType) end

@for_petsc function TSInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:TSInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	TSMonitor(petsclib::PetscLibType,ts::TS, step::PetscInt, ptime::PetscReal, u::PetscVec) 
Runs all user

Collective

Input Parameters:
- `ts`    - time stepping context obtained from `TSCreate()`
- `step`  - step number that has just completed
- `ptime` - model time of the state
- `u`     - state at the current model time

Level: developer

-seealso: `TS`, `TSMonitorSet()`, `TSMonitorSetFromOptions()`

# External Links
$(_doc_external("Ts/TSMonitor"))
"""
function TSMonitor(petsclib::PetscLibType, ts::TS, step::PetscInt, ptime::PetscReal, u::PetscVec) end

@for_petsc function TSMonitor(petsclib::$UnionPetscLib, ts::TS, step::$PetscInt, ptime::$PetscReal, u::PetscVec )

    @chk ccall(
               (:TSMonitor, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt, $PetscReal, CVec),
               ts, step, ptime, u,
              )


	return nothing
end 

"""
	TSMonitorSetFromOptions(petsclib::PetscLibType,ts::TS, name::String, help::String, manual::String, monitor::external, monitorsetup::external) 
Sets a monitor function and viewer appropriate for the type indicated by the user

Collective

Input Parameters:
- `ts`           - `TS` object you wish to monitor
- `name`         - the monitor type one is seeking
- `help`         - message indicating what monitoring is done
- `manual`       - manual page for the monitor
- `monitor`      - the monitor function, this must use a `PetscViewerFormat` as its context
- `monitorsetup` - a function that is called once ONLY if the user selected this monitor that may set additional features of the `TS` or `PetscViewer` objects

Level: developer

-seealso: [](ch_ts), `TS`, `TSMonitorSet()`, `PetscOptionsCreateViewer()`, `PetscOptionsGetReal()`, `PetscOptionsHasName()`, `PetscOptionsGetString()`,
`PetscOptionsGetIntArray()`, `PetscOptionsGetRealArray()`, `PetscOptionsBool()`
`PetscOptionsInt()`, `PetscOptionsString()`, `PetscOptionsReal()`,
`PetscOptionsName()`, `PetscOptionsBegin()`, `PetscOptionsEnd()`, `PetscOptionsHeadBegin()`,
`PetscOptionsStringArray()`, `PetscOptionsRealArray()`, `PetscOptionsScalar()`,
`PetscOptionsBoolGroupBegin()`, `PetscOptionsBoolGroup()`, `PetscOptionsBoolGroupEnd()`,
`PetscOptionsFList()`, `PetscOptionsEList()`

# External Links
$(_doc_external("Ts/TSMonitorSetFromOptions"))
"""
function TSMonitorSetFromOptions(petsclib::PetscLibType, ts::TS, name::String, help::String, manual::String, monitor::external, monitorsetup::external) end

@for_petsc function TSMonitorSetFromOptions(petsclib::$UnionPetscLib, ts::TS, name::String, help::String, manual::String, monitor::external, monitorsetup::external )

    @chk ccall(
               (:TSMonitorSetFromOptions, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, external, external),
               ts, name, help, manual, monitor, monitorsetup,
              )


	return nothing
end 

"""
	TSMonitorSet(petsclib::PetscLibType,ts::TS, monitor::external, mctx::Cvoid, mdestroy::PetscCtxDestroyFn) 
Sets an ADDITIONAL function that is to be used at every
timestep to display the iteration's  progress.

Logically Collective

Input Parameters:
- `ts`       - the `TS` context obtained from `TSCreate()`
- `monitor`  - monitoring routine
- `mctx`     - [optional] user-defined context for private data for the monitor routine (use `NULL` if no context is desired)
- `mdestroy` - [optional] routine that frees monitor context (may be `NULL`), see `PetscCtxDestroyFn` for the calling sequence

Calling sequence of `monitor`:
- `ts`    - the `TS` context
- `steps` - iteration number (after the final time step the monitor routine may be called with a step of -1, this indicates the solution has been interpolated to this time)
- `time`  - current time
- `u`     - current iterate
- `ctx`   - [optional] monitoring context

Level: intermediate

-seealso: [](ch_ts), `TSMonitorDefault()`, `TSMonitorCancel()`, `TSDMSwarmMonitorMoments()`, `TSMonitorExtreme()`, `TSMonitorDrawSolution()`,
`TSMonitorDrawSolutionPhase()`, `TSMonitorDrawSolutionFunction()`, `TSMonitorDrawError()`, `TSMonitorSolution()`, `TSMonitorSolutionVTK()`,
`TSMonitorLGSolution()`, `TSMonitorLGError()`, `TSMonitorSPSwarmSolution()`, `TSMonitorError()`, `TSMonitorEnvelope()`,  `PetscCtxDestroyFn`

# External Links
$(_doc_external("Ts/TSMonitorSet"))
"""
function TSMonitorSet(petsclib::PetscLibType, ts::TS, monitor::external, mctx::Cvoid, mdestroy::PetscCtxDestroyFn) end

@for_petsc function TSMonitorSet(petsclib::$UnionPetscLib, ts::TS, monitor::external, mctx::Cvoid, mdestroy::PetscCtxDestroyFn )

    @chk ccall(
               (:TSMonitorSet, $petsc_library),
               PetscErrorCode,
               (CTS, external, Ptr{Cvoid}, Ptr{PetscCtxDestroyFn}),
               ts, monitor, mctx, mdestroy,
              )


	return nothing
end 

"""
	TSMonitorCancel(petsclib::PetscLibType,ts::TS) 
Clears all the monitors that have been set on a time

Logically Collective

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSMonitorDefault()`, `TSMonitorSet()`

# External Links
$(_doc_external("Ts/TSMonitorCancel"))
"""
function TSMonitorCancel(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSMonitorCancel(petsclib::$UnionPetscLib, ts::TS )

    @chk ccall(
               (:TSMonitorCancel, $petsc_library),
               PetscErrorCode,
               (CTS,),
               ts,
              )


	return nothing
end 

"""
	TSMonitorDefault(petsclib::PetscLibType,ts::TS, step::PetscInt, ptime::PetscReal, v::PetscVec, vf::PetscViewerAndFormat) 
The default monitor, prints the timestep and time for each step

Input Parameters:
- `ts`    - the `TS` context
- `step`  - iteration number (after the final time step the monitor routine may be called with a step of -1, this indicates the solution has been interpolated to this time)
- `ptime` - current time
- `v`     - current iterate
- `vf`    - the viewer and format

Options Database Key:
- `-ts_monitor` - monitors the time integration

Level: intermediate

-seealso: [](ch_ts), `TSMonitorSet()`, `TSDMSwarmMonitorMoments()`, `TSMonitorWallClockTime()`, `TSMonitorExtreme()`, `TSMonitorDrawSolution()`,
`TSMonitorDrawSolutionPhase()`, `TSMonitorDrawSolutionFunction()`, `TSMonitorDrawError()`, `TSMonitorSolution()`, `TSMonitorSolutionVTK()`,
`TSMonitorLGSolution()`, `TSMonitorLGError()`, `TSMonitorSPSwarmSolution()`, `TSMonitorError()`, `TSMonitorEnvelope()`

# External Links
$(_doc_external("Ts/TSMonitorDefault"))
"""
function TSMonitorDefault(petsclib::PetscLibType, ts::TS, step::PetscInt, ptime::PetscReal, v::PetscVec, vf::PetscViewerAndFormat) end

@for_petsc function TSMonitorDefault(petsclib::$UnionPetscLib, ts::TS, step::$PetscInt, ptime::$PetscReal, v::PetscVec, vf::PetscViewerAndFormat )

    @chk ccall(
               (:TSMonitorDefault, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt, $PetscReal, CVec, Ptr{PetscViewerAndFormat}),
               ts, step, ptime, v, vf,
              )


	return nothing
end 

"""
	TSMonitorWallClockTimeSetUp(petsclib::PetscLibType,ts::TS, vf::PetscViewerAndFormat) 
Setup routine passed to `TSMonitorSetFromOptions()` when using `

Input Parameters:
- `ts` - the `TS` context
- `vf` - the viewer and format

Level: intermediate

-seealso: [](ch_ts), `TSMonitorSet()`

# External Links
$(_doc_external("Ts/TSMonitorWallClockTimeSetUp"))
"""
function TSMonitorWallClockTimeSetUp(petsclib::PetscLibType, ts::TS, vf::PetscViewerAndFormat) end

@for_petsc function TSMonitorWallClockTimeSetUp(petsclib::$UnionPetscLib, ts::TS, vf::PetscViewerAndFormat )

    @chk ccall(
               (:TSMonitorWallClockTimeSetUp, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{PetscViewerAndFormat}),
               ts, vf,
              )


	return nothing
end 

"""
	TSMonitorWallClockTime(petsclib::PetscLibType,ts::TS, step::PetscInt, ptime::PetscReal, v::PetscVec, vf::PetscViewerAndFormat) 
Monitor wall

Input Parameters:
- `ts`    - the `TS` context
- `step`  - iteration number (after the final time step the monitor routine may be called with a step of -1, this indicates the solution has been interpolated to this time)
- `ptime` - current time
- `v`     - current solution
- `vf`    - the viewer and format

Options Database Key:
- `-ts_monitor_wall_clock_time` - Monitor wall-clock time, KSP iterations, and SNES iterations per step.

Level: intermediate

-seealso: [](ch_ts), `TSMonitorSet()`, `TSMonitorDefault()`, `TSMonitorExtreme()`, `TSMonitorDrawSolution()`,
`TSMonitorDrawSolutionPhase()`, `TSMonitorDrawSolutionFunction()`, `TSMonitorDrawError()`, `TSMonitorSolution()`, `TSMonitorSolutionVTK()`,
`TSMonitorLGSolution()`, `TSMonitorLGError()`, `TSMonitorSPSwarmSolution()`, `TSMonitorError()`, `TSMonitorEnvelope()`, `TSDMSwarmMonitorMoments()`

# External Links
$(_doc_external("Ts/TSMonitorWallClockTime"))
"""
function TSMonitorWallClockTime(petsclib::PetscLibType, ts::TS, step::PetscInt, ptime::PetscReal, v::PetscVec, vf::PetscViewerAndFormat) end

@for_petsc function TSMonitorWallClockTime(petsclib::$UnionPetscLib, ts::TS, step::$PetscInt, ptime::$PetscReal, v::PetscVec, vf::PetscViewerAndFormat )

    @chk ccall(
               (:TSMonitorWallClockTime, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt, $PetscReal, CVec, Ptr{PetscViewerAndFormat}),
               ts, step, ptime, v, vf,
              )


	return nothing
end 

"""
	TSMonitorExtreme(petsclib::PetscLibType,ts::TS, step::PetscInt, ptime::PetscReal, v::PetscVec, vf::PetscViewerAndFormat) 
Prints the extreme values of the solution at each timestep

Input Parameters:
- `ts`    - the `TS` context
- `step`  - iteration number (after the final time step the monitor routine may be called with a step of -1, this indicates the solution has been interpolated to this time)
- `ptime` - current time
- `v`     - current iterate
- `vf`    - the viewer and format

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSMonitorSet()`

# External Links
$(_doc_external("Ts/TSMonitorExtreme"))
"""
function TSMonitorExtreme(petsclib::PetscLibType, ts::TS, step::PetscInt, ptime::PetscReal, v::PetscVec, vf::PetscViewerAndFormat) end

@for_petsc function TSMonitorExtreme(petsclib::$UnionPetscLib, ts::TS, step::$PetscInt, ptime::$PetscReal, v::PetscVec, vf::PetscViewerAndFormat )

    @chk ccall(
               (:TSMonitorExtreme, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt, $PetscReal, CVec, Ptr{PetscViewerAndFormat}),
               ts, step, ptime, v, vf,
              )


	return nothing
end 

"""
	TSMonitorLGTimeStep(petsclib::PetscLibType,ts::TS, step::PetscInt, ptime::PetscReal, v::PetscVec, monctx::Cvoid) 
Monitors a `TS` by printing the time

Collective

Input Parameters:
- `ts`     - the time integrator
- `step`   - the current time step
- `ptime`  - the current time
- `v`      - the current state
- `monctx` - the monitor context obtained with `TSMonitorLGCtxCreate()`

Level: advanced

-seealso: [](ch_ts), `TS`, `TSMonitorLGCtxCreate()`, `TSMonitorSet()`, `TSMonitorLGCtxDestroy()`

# External Links
$(_doc_external("Ts/TSMonitorLGTimeStep"))
"""
function TSMonitorLGTimeStep(petsclib::PetscLibType, ts::TS, step::PetscInt, ptime::PetscReal, v::PetscVec, monctx::Cvoid) end

@for_petsc function TSMonitorLGTimeStep(petsclib::$UnionPetscLib, ts::TS, step::$PetscInt, ptime::$PetscReal, v::PetscVec, monctx::Cvoid )

    @chk ccall(
               (:TSMonitorLGTimeStep, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt, $PetscReal, CVec, Ptr{Cvoid}),
               ts, step, ptime, v, monctx,
              )


	return nothing
end 

"""
	TSMonitorDrawSolution(petsclib::PetscLibType,ts::TS, step::PetscInt, ptime::PetscReal, u::PetscVec, dummy::Cvoid) 
Monitors progress of the `TS` solvers by calling
`VecView()` for the solution at each timestep

Collective

Input Parameters:
- `ts`    - the `TS` context
- `step`  - current time-step
- `ptime` - current time
- `u`     - the solution at the current time
- `dummy` - either a viewer or `NULL`

Options Database Keys:
- `-ts_monitor_draw_solution`         - draw the solution at each time-step
- `-ts_monitor_draw_solution_initial` - show initial solution as well as current solution

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSMonitorSet()`, `TSMonitorDefault()`, `VecView()`, `TSMonitorDrawCtxCreate()`, `TSMonitorDrawCtxDestroy()`

# External Links
$(_doc_external("Ts/TSMonitorDrawSolution"))
"""
function TSMonitorDrawSolution(petsclib::PetscLibType, ts::TS, step::PetscInt, ptime::PetscReal, u::PetscVec, dummy::Cvoid) end

@for_petsc function TSMonitorDrawSolution(petsclib::$UnionPetscLib, ts::TS, step::$PetscInt, ptime::$PetscReal, u::PetscVec, dummy::Cvoid )

    @chk ccall(
               (:TSMonitorDrawSolution, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt, $PetscReal, CVec, Ptr{Cvoid}),
               ts, step, ptime, u, dummy,
              )


	return nothing
end 

"""
	TSMonitorDrawSolutionPhase(petsclib::PetscLibType,ts::TS, step::PetscInt, ptime::PetscReal, u::PetscVec, dummy::Cvoid) 
Monitors progress of the `TS` solvers by plotting the solution as a phase diagram

Collective

Input Parameters:
- `ts`    - the `TS` context
- `step`  - current time-step
- `ptime` - current time
- `u`     - the solution at the current time
- `dummy` - either a viewer or `NULL`

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSMonitorSet()`, `TSMonitorDefault()`, `VecView()`

# External Links
$(_doc_external("Ts/TSMonitorDrawSolutionPhase"))
"""
function TSMonitorDrawSolutionPhase(petsclib::PetscLibType, ts::TS, step::PetscInt, ptime::PetscReal, u::PetscVec, dummy::Cvoid) end

@for_petsc function TSMonitorDrawSolutionPhase(petsclib::$UnionPetscLib, ts::TS, step::$PetscInt, ptime::$PetscReal, u::PetscVec, dummy::Cvoid )

    @chk ccall(
               (:TSMonitorDrawSolutionPhase, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt, $PetscReal, CVec, Ptr{Cvoid}),
               ts, step, ptime, u, dummy,
              )


	return nothing
end 

"""
	TSMonitorDrawSolutionFunction(petsclib::PetscLibType,ts::TS, step::PetscInt, ptime::PetscReal, u::PetscVec, dummy::Cvoid) 
Monitors progress of the `TS` solvers by calling
`VecView()` for the solution provided by `TSSetSolutionFunction()` at each timestep

Collective

Input Parameters:
- `ts`    - the `TS` context
- `step`  - current time-step
- `ptime` - current time
- `u`     - solution at current time
- `dummy` - either a viewer or `NULL`

Options Database Key:
- `-ts_monitor_draw_solution_function` - Monitor error graphically, requires user to have provided `TSSetSolutionFunction()`

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSMonitorSet()`, `TSMonitorDefault()`, `VecView()`, `TSSetSolutionFunction()`

# External Links
$(_doc_external("Ts/TSMonitorDrawSolutionFunction"))
"""
function TSMonitorDrawSolutionFunction(petsclib::PetscLibType, ts::TS, step::PetscInt, ptime::PetscReal, u::PetscVec, dummy::Cvoid) end

@for_petsc function TSMonitorDrawSolutionFunction(petsclib::$UnionPetscLib, ts::TS, step::$PetscInt, ptime::$PetscReal, u::PetscVec, dummy::Cvoid )

    @chk ccall(
               (:TSMonitorDrawSolutionFunction, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt, $PetscReal, CVec, Ptr{Cvoid}),
               ts, step, ptime, u, dummy,
              )


	return nothing
end 

"""
	TSMonitorDrawError(petsclib::PetscLibType,ts::TS, step::PetscInt, ptime::PetscReal, u::PetscVec, dummy::Cvoid) 
Monitors progress of the `TS` solvers by calling
`VecView()` for the error at each timestep

Collective

Input Parameters:
- `ts`    - the `TS` context
- `step`  - current time-step
- `ptime` - current time
- `u`     - solution at current time
- `dummy` - either a viewer or `NULL`

Options Database Key:
- `-ts_monitor_draw_error` - Monitor error graphically, requires user to have provided `TSSetSolutionFunction()`

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSMonitorSet()`, `TSMonitorDefault()`, `VecView()`, `TSSetSolutionFunction()`

# External Links
$(_doc_external("Ts/TSMonitorDrawError"))
"""
function TSMonitorDrawError(petsclib::PetscLibType, ts::TS, step::PetscInt, ptime::PetscReal, u::PetscVec, dummy::Cvoid) end

@for_petsc function TSMonitorDrawError(petsclib::$UnionPetscLib, ts::TS, step::$PetscInt, ptime::$PetscReal, u::PetscVec, dummy::Cvoid )

    @chk ccall(
               (:TSMonitorDrawError, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt, $PetscReal, CVec, Ptr{Cvoid}),
               ts, step, ptime, u, dummy,
              )


	return nothing
end 

"""
	TSMonitorSolutionSetup(petsclib::PetscLibType,ts::TS, vf::PetscViewerAndFormat) 
Setups the context for `TSMonitorSolution()`

Collective

Input Parameters:
- `ts` - the `TS` context
- `vf` - viewer and its format

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSMonitorSolution()`, `TSMonitorSet()`, `TSMonitorDefault()`, `VecView()`, `TSMonitorSetFromOptions()`

# External Links
$(_doc_external("Ts/TSMonitorSolutionSetup"))
"""
function TSMonitorSolutionSetup(petsclib::PetscLibType, ts::TS, vf::PetscViewerAndFormat) end

@for_petsc function TSMonitorSolutionSetup(petsclib::$UnionPetscLib, ts::TS, vf::PetscViewerAndFormat )

    @chk ccall(
               (:TSMonitorSolutionSetup, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{PetscViewerAndFormat}),
               ts, vf,
              )


	return nothing
end 

"""
	TSMonitorSolution(petsclib::PetscLibType,ts::TS, step::PetscInt, ptime::PetscReal, u::PetscVec, vf::PetscViewerAndFormat) 
Monitors progress of the `TS` solvers by `VecView()` for the solution at each timestep. Normally the viewer is a binary file or a `PetscDraw` object

Collective

Input Parameters:
- `ts`    - the `TS` context
- `step`  - current time-step
- `ptime` - current time
- `u`     - current state
- `vf`    - viewer and its format

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSMonitorSet()`, `TSMonitorDefault()`, `VecView()`, `TSMonitorSolutionSetup()`,

# External Links
$(_doc_external("Ts/TSMonitorSolution"))
"""
function TSMonitorSolution(petsclib::PetscLibType, ts::TS, step::PetscInt, ptime::PetscReal, u::PetscVec, vf::PetscViewerAndFormat) end

@for_petsc function TSMonitorSolution(petsclib::$UnionPetscLib, ts::TS, step::$PetscInt, ptime::$PetscReal, u::PetscVec, vf::PetscViewerAndFormat )

    @chk ccall(
               (:TSMonitorSolution, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt, $PetscReal, CVec, Ptr{PetscViewerAndFormat}),
               ts, step, ptime, u, vf,
              )


	return nothing
end 

"""
	TSMonitorSolutionVTK(petsclib::PetscLibType,ts::TS, step::PetscInt, ptime::PetscReal, u::PetscVec, ctx::TSMonitorVTKCtx) 
Monitors progress of the `TS` solvers by `VecView()` for the solution at selected timesteps.

Collective

Input Parameters:
- `ts`    - the `TS` context
- `step`  - current time-step
- `ptime` - current time
- `u`     - current state
- `ctx`   - monitor context obtained with `TSMonitorSolutionVTKCtxCreate()`

Level: developer

-seealso: [](ch_ts), `TS`, `TSMonitorSet()`, `TSMonitorDefault()`, `VecView()`

# External Links
$(_doc_external("Ts/TSMonitorSolutionVTK"))
"""
function TSMonitorSolutionVTK(petsclib::PetscLibType, ts::TS, step::PetscInt, ptime::PetscReal, u::PetscVec, ctx::TSMonitorVTKCtx) end

@for_petsc function TSMonitorSolutionVTK(petsclib::$UnionPetscLib, ts::TS, step::$PetscInt, ptime::$PetscReal, u::PetscVec, ctx::TSMonitorVTKCtx )

    @chk ccall(
               (:TSMonitorSolutionVTK, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt, $PetscReal, CVec, TSMonitorVTKCtx),
               ts, step, ptime, u, ctx,
              )


	return nothing
end 

"""
	TSMonitorSolutionVTKDestroy(petsclib::PetscLibType,ctx::TSMonitorVTKCtx) 
Destroy the monitor context created with `TSMonitorSolutionVTKCtxCreate()`

Not Collective

Input Parameter:
- `ctx` - the monitor context

Level: developer

-seealso: [](ch_ts), `TSMonitorSet()`, `TSMonitorSolutionVTK()`

# External Links
$(_doc_external("Ts/TSMonitorSolutionVTKDestroy"))
"""
function TSMonitorSolutionVTKDestroy(petsclib::PetscLibType, ctx::TSMonitorVTKCtx) end

@for_petsc function TSMonitorSolutionVTKDestroy(petsclib::$UnionPetscLib, ctx::TSMonitorVTKCtx )

    @chk ccall(
               (:TSMonitorSolutionVTKDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{TSMonitorVTKCtx},),
               ctx,
              )


	return nothing
end 

"""
	ctx::TSMonitorVTKCtx = TSMonitorSolutionVTKCtxCreate(petsclib::PetscLibType,filenametemplate::String) 
Create the monitor context to be used in `TSMonitorSolutionVTK()`

Not collective

Input Parameter:
- `filenametemplate` - the template file name, e.g. foo-%03d.vts

Output Parameter:
- `ctx` - the monitor context

Level: developer

-seealso: [](ch_ts), `TSMonitorSet()`, `TSMonitorSolutionVTK()`, `TSMonitorSolutionVTKDestroy()`

# External Links
$(_doc_external("Ts/TSMonitorSolutionVTKCtxCreate"))
"""
function TSMonitorSolutionVTKCtxCreate(petsclib::PetscLibType, filenametemplate::String) end

@for_petsc function TSMonitorSolutionVTKCtxCreate(petsclib::$UnionPetscLib, filenametemplate::String )
	ctx_ = Ref{TSMonitorVTKCtx}()

    @chk ccall(
               (:TSMonitorSolutionVTKCtxCreate, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{TSMonitorVTKCtx}),
               filenametemplate, ctx_,
              )

	ctx = ctx_[]

	return ctx
end 

"""
	TSMonitorLGSolution(petsclib::PetscLibType,ts::TS, step::PetscInt, ptime::PetscReal, u::PetscVec, dctx::Cvoid) 
Monitors progress of the `TS` solvers by plotting each component of the solution vector
in a time based line graph

Collective

Input Parameters:
- `ts`    - the `TS` context
- `step`  - current time-step
- `ptime` - current time
- `u`     - current solution
- `dctx`  - the `TSMonitorLGCtx` object that contains all the options for the monitoring, this is created with `TSMonitorLGCtxCreate()`

Options Database Key:
- `-ts_monitor_lg_solution_variables` - enable monitor of lg solution variables

Level: intermediate

-seealso: [](ch_ts), `TSMonitorSet()`, `TSMonitorDefault()`, `VecView()`, `TSMonitorLGCtxCreate()`, `TSMonitorLGCtxSetVariableNames()`, `TSMonitorLGCtxGetVariableNames()`,
`TSMonitorLGSetVariableNames()`, `TSMonitorLGGetVariableNames()`, `TSMonitorLGSetDisplayVariables()`, `TSMonitorLGCtxSetDisplayVariables()`,
`TSMonitorLGCtxSetTransform()`, `TSMonitorLGSetTransform()`, `TSMonitorLGError()`, `TSMonitorLGSNESIterations()`, `TSMonitorLGKSPIterations()`,
`TSMonitorEnvelopeCtxCreate()`, `TSMonitorEnvelopeGetBounds()`, `TSMonitorEnvelopeCtxDestroy()`, `TSMonitorEnvelop()`

# External Links
$(_doc_external("Ts/TSMonitorLGSolution"))
"""
function TSMonitorLGSolution(petsclib::PetscLibType, ts::TS, step::PetscInt, ptime::PetscReal, u::PetscVec, dctx::Cvoid) end

@for_petsc function TSMonitorLGSolution(petsclib::$UnionPetscLib, ts::TS, step::$PetscInt, ptime::$PetscReal, u::PetscVec, dctx::Cvoid )

    @chk ccall(
               (:TSMonitorLGSolution, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt, $PetscReal, CVec, Ptr{Cvoid}),
               ts, step, ptime, u, dctx,
              )


	return nothing
end 

"""
	TSMonitorLGSetVariableNames(petsclib::PetscLibType,ts::TS, names::Cchar) 
Sets the name of each component in the solution vector so that it may be displayed in the plot

Collective

Input Parameters:
- `ts`    - the `TS` context
- `names` - the names of the components, final string must be `NULL`

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSMonitorSet()`, `TSMonitorDefault()`, `VecView()`, `TSMonitorLGSetDisplayVariables()`, `TSMonitorLGCtxSetVariableNames()`

# External Links
$(_doc_external("Ts/TSMonitorLGSetVariableNames"))
"""
function TSMonitorLGSetVariableNames(petsclib::PetscLibType, ts::TS, names::Cchar) end

@for_petsc function TSMonitorLGSetVariableNames(petsclib::$UnionPetscLib, ts::TS, names::Cchar )

    @chk ccall(
               (:TSMonitorLGSetVariableNames, $petsc_library),
               PetscErrorCode,
               (CTS, Cchar),
               ts, names,
              )


	return nothing
end 

"""
	TSMonitorLGGetVariableNames(petsclib::PetscLibType,ts::TS, names::Cchar) 
Gets the name of each component in the solution vector so that it may be displayed in the plot

Collective

Input Parameter:
- `ts` - the `TS` context

Output Parameter:
- `names` - the names of the components, final string must be `NULL`

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSMonitorSet()`, `TSMonitorDefault()`, `VecView()`, `TSMonitorLGSetDisplayVariables()`

# External Links
$(_doc_external("Ts/TSMonitorLGGetVariableNames"))
"""
function TSMonitorLGGetVariableNames(petsclib::PetscLibType, ts::TS, names::Cchar) end

@for_petsc function TSMonitorLGGetVariableNames(petsclib::$UnionPetscLib, ts::TS, names::Cchar )

    @chk ccall(
               (:TSMonitorLGGetVariableNames, $petsc_library),
               PetscErrorCode,
               (CTS, Cchar),
               ts, names,
              )


	return nothing
end 

"""
	TSMonitorLGSetDisplayVariables(petsclib::PetscLibType,ts::TS, displaynames::Cchar) 
Sets the variables that are to be display in the monitor

Collective

Input Parameters:
- `ts`           - the `TS` context
- `displaynames` - the names of the components, final string must be `NULL`

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSMonitorSet()`, `TSMonitorDefault()`, `VecView()`, `TSMonitorLGSetVariableNames()`

# External Links
$(_doc_external("Ts/TSMonitorLGSetDisplayVariables"))
"""
function TSMonitorLGSetDisplayVariables(petsclib::PetscLibType, ts::TS, displaynames::Cchar) end

@for_petsc function TSMonitorLGSetDisplayVariables(petsclib::$UnionPetscLib, ts::TS, displaynames::Cchar )

    @chk ccall(
               (:TSMonitorLGSetDisplayVariables, $petsc_library),
               PetscErrorCode,
               (CTS, Cchar),
               ts, displaynames,
              )


	return nothing
end 

"""
	TSMonitorLGSetTransform(petsclib::PetscLibType,ts::TS, transform::external, destroy::PetscCtxDestroyFn, tctx::Cvoid) 
Solution vector will be transformed by provided function before being displayed

Collective

Input Parameters:
- `ts`        - the `TS` context
- `transform` - the transform function
- `destroy`   - function to destroy the optional context, see `PetscCtxDestroyFn` for its calling sequence
- `tctx`      - optional context used by transform function

Level: intermediate

-seealso: [](ch_ts), `TSMonitorSet()`, `TSMonitorDefault()`, `VecView()`, `TSMonitorLGSetVariableNames()`, `TSMonitorLGCtxSetTransform()`, `PetscCtxDestroyFn`

# External Links
$(_doc_external("Ts/TSMonitorLGSetTransform"))
"""
function TSMonitorLGSetTransform(petsclib::PetscLibType, ts::TS, transform::external, destroy::PetscCtxDestroyFn, tctx::Cvoid) end

@for_petsc function TSMonitorLGSetTransform(petsclib::$UnionPetscLib, ts::TS, transform::external, destroy::PetscCtxDestroyFn, tctx::Cvoid )

    @chk ccall(
               (:TSMonitorLGSetTransform, $petsc_library),
               PetscErrorCode,
               (CTS, external, Ptr{PetscCtxDestroyFn}, Ptr{Cvoid}),
               ts, transform, destroy, tctx,
              )


	return nothing
end 

"""
	TSMonitorLGError(petsclib::PetscLibType,ts::TS, step::PetscInt, ptime::PetscReal, u::PetscVec, dummy::Cvoid) 
Monitors progress of the `TS` solvers by plotting each component of the error
in a time based line graph

Collective

Input Parameters:
- `ts`    - the `TS` context
- `step`  - current time-step
- `ptime` - current time
- `u`     - current solution
- `dummy` - `TSMonitorLGCtx` object created with `TSMonitorLGCtxCreate()`

Options Database Key:
- `-ts_monitor_lg_error` - create a graphical monitor of error history

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSMonitorSet()`, `TSMonitorDefault()`, `VecView()`, `TSSetSolutionFunction()`

# External Links
$(_doc_external("Ts/TSMonitorLGError"))
"""
function TSMonitorLGError(petsclib::PetscLibType, ts::TS, step::PetscInt, ptime::PetscReal, u::PetscVec, dummy::Cvoid) end

@for_petsc function TSMonitorLGError(petsclib::$UnionPetscLib, ts::TS, step::$PetscInt, ptime::$PetscReal, u::PetscVec, dummy::Cvoid )

    @chk ccall(
               (:TSMonitorLGError, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt, $PetscReal, CVec, Ptr{Cvoid}),
               ts, step, ptime, u, dummy,
              )


	return nothing
end 

"""
	TSMonitorSPSwarmSolution(petsclib::PetscLibType,ts::TS, step::PetscInt, ptime::PetscReal, u::PetscVec, dctx::Cvoid) 
Graphically displays phase plots of `DMSWARM` particles on a scatter plot

Input Parameters:
- `ts`    - the `TS` context
- `step`  - current time-step
- `ptime` - current time
- `u`     - current solution
- `dctx`  - the `TSMonitorSPCtx` object that contains all the options for the monitoring, this is created with `TSMonitorSPCtxCreate()`

Options Database Keys:
- `-ts_monitor_sp_swarm <n>`                  - Monitor the solution every n steps, or -1 for plotting only the final solution
- `-ts_monitor_sp_swarm_retain <n>`           - Retain n old points so we can see the history, or -1 for all points
- `-ts_monitor_sp_swarm_multi_species <bool>` - Color each species differently
- `-ts_monitor_sp_swarm_phase <bool>`         - Plot in phase space, as opposed to coordinate space

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSMonitorSet()`, `DMSWARM`, `TSMonitorSPCtxCreate()`

# External Links
$(_doc_external("Ts/TSMonitorSPSwarmSolution"))
"""
function TSMonitorSPSwarmSolution(petsclib::PetscLibType, ts::TS, step::PetscInt, ptime::PetscReal, u::PetscVec, dctx::Cvoid) end

@for_petsc function TSMonitorSPSwarmSolution(petsclib::$UnionPetscLib, ts::TS, step::$PetscInt, ptime::$PetscReal, u::PetscVec, dctx::Cvoid )

    @chk ccall(
               (:TSMonitorSPSwarmSolution, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt, $PetscReal, CVec, Ptr{Cvoid}),
               ts, step, ptime, u, dctx,
              )


	return nothing
end 

"""
	TSMonitorHGSwarmSolution(petsclib::PetscLibType,ts::TS, step::PetscInt, ptime::PetscReal, u::PetscVec, dctx::Cvoid) 
Graphically displays histograms of `DMSWARM` particles

Input Parameters:
- `ts`    - the `TS` context
- `step`  - current time-step
- `ptime` - current time
- `u`     - current solution
- `dctx`  - the `TSMonitorSPCtx` object that contains all the options for the monitoring, this is created with `TSMonitorHGCtxCreate()`

Options Database Keys:
- `-ts_monitor_hg_swarm <n>`             - Monitor the solution every n steps, or -1 for plotting only the final solution
- `-ts_monitor_hg_swarm_species <num>`   - Number of species to histogram
- `-ts_monitor_hg_swarm_bins <num>`      - Number of histogram bins
- `-ts_monitor_hg_swarm_velocity <bool>` - Plot in velocity space, as opposed to coordinate space

Level: intermediate

-seealso: `TSMonitorSet()`

# External Links
$(_doc_external("Ts/TSMonitorHGSwarmSolution"))
"""
function TSMonitorHGSwarmSolution(petsclib::PetscLibType, ts::TS, step::PetscInt, ptime::PetscReal, u::PetscVec, dctx::Cvoid) end

@for_petsc function TSMonitorHGSwarmSolution(petsclib::$UnionPetscLib, ts::TS, step::$PetscInt, ptime::$PetscReal, u::PetscVec, dctx::Cvoid )

    @chk ccall(
               (:TSMonitorHGSwarmSolution, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt, $PetscReal, CVec, Ptr{Cvoid}),
               ts, step, ptime, u, dctx,
              )


	return nothing
end 

"""
	TSMonitorError(petsclib::PetscLibType,ts::TS, step::PetscInt, ptime::PetscReal, u::PetscVec, vf::PetscViewerAndFormat) 
Monitors progress of the `TS` solvers by printing the 2 norm of the error at each timestep

Collective

Input Parameters:
- `ts`    - the `TS` context
- `step`  - current time-step
- `ptime` - current time
- `u`     - current solution
- `vf`    - unused context

Options Database Key:
- `-ts_monitor_error` - create a graphical monitor of error history

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSMonitorSet()`, `TSMonitorDefault()`, `VecView()`, `TSSetSolutionFunction()`

# External Links
$(_doc_external("Ts/TSMonitorError"))
"""
function TSMonitorError(petsclib::PetscLibType, ts::TS, step::PetscInt, ptime::PetscReal, u::PetscVec, vf::PetscViewerAndFormat) end

@for_petsc function TSMonitorError(petsclib::$UnionPetscLib, ts::TS, step::$PetscInt, ptime::$PetscReal, u::PetscVec, vf::PetscViewerAndFormat )

    @chk ccall(
               (:TSMonitorError, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt, $PetscReal, CVec, Ptr{PetscViewerAndFormat}),
               ts, step, ptime, u, vf,
              )


	return nothing
end 

"""
	TSMonitorLGSNESIterations(petsclib::PetscLibType,ts::TS, n::PetscInt, ptime::PetscReal, v::PetscVec, monctx::Cvoid) 

# External Links
$(_doc_external("Ts/TSMonitorLGSNESIterations"))
"""
function TSMonitorLGSNESIterations(petsclib::PetscLibType, ts::TS, n::PetscInt, ptime::PetscReal, v::PetscVec, monctx::Cvoid) end

@for_petsc function TSMonitorLGSNESIterations(petsclib::$UnionPetscLib, ts::TS, n::$PetscInt, ptime::$PetscReal, v::PetscVec, monctx::Cvoid )

    @chk ccall(
               (:TSMonitorLGSNESIterations, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt, $PetscReal, CVec, Ptr{Cvoid}),
               ts, n, ptime, v, monctx,
              )


	return nothing
end 

"""
	TSMonitorLGKSPIterations(petsclib::PetscLibType,ts::TS, n::PetscInt, ptime::PetscReal, v::PetscVec, monctx::Cvoid) 

# External Links
$(_doc_external("Ts/TSMonitorLGKSPIterations"))
"""
function TSMonitorLGKSPIterations(petsclib::PetscLibType, ts::TS, n::PetscInt, ptime::PetscReal, v::PetscVec, monctx::Cvoid) end

@for_petsc function TSMonitorLGKSPIterations(petsclib::$UnionPetscLib, ts::TS, n::$PetscInt, ptime::$PetscReal, v::PetscVec, monctx::Cvoid )

    @chk ccall(
               (:TSMonitorLGKSPIterations, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt, $PetscReal, CVec, Ptr{Cvoid}),
               ts, n, ptime, v, monctx,
              )


	return nothing
end 

"""
	TSMonitorEnvelope(petsclib::PetscLibType,ts::TS, step::PetscInt, ptime::PetscReal, u::PetscVec, dctx::Cvoid) 
Monitors the maximum and minimum value of each component of the solution

Collective

Input Parameters:
- `ts`    - the `TS` context
- `step`  - current time-step
- `ptime` - current time
- `u`     - current solution
- `dctx`  - the envelope context

Options Database Key:
- `-ts_monitor_envelope` - determine maximum and minimum value of each component of the solution over the solution time

Level: intermediate

-seealso: [](ch_ts), `TSMonitorSet()`, `TSMonitorDefault()`, `VecView()`, `TSMonitorEnvelopeGetBounds()`, `TSMonitorEnvelopeCtxCreate()`

# External Links
$(_doc_external("Ts/TSMonitorEnvelope"))
"""
function TSMonitorEnvelope(petsclib::PetscLibType, ts::TS, step::PetscInt, ptime::PetscReal, u::PetscVec, dctx::Cvoid) end

@for_petsc function TSMonitorEnvelope(petsclib::$UnionPetscLib, ts::TS, step::$PetscInt, ptime::$PetscReal, u::PetscVec, dctx::Cvoid )

    @chk ccall(
               (:TSMonitorEnvelope, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt, $PetscReal, CVec, Ptr{Cvoid}),
               ts, step, ptime, u, dctx,
              )


	return nothing
end 

"""
	TSMonitorEnvelopeGetBounds(petsclib::PetscLibType,ts::TS, max::PetscVec, min::PetscVec) 
Gets the bounds for the components of the solution

Collective

Input Parameter:
- `ts` - the `TS` context

Output Parameters:
- `max` - the maximum values
- `min` - the minimum values

Level: intermediate

-seealso: [](ch_ts), `TSMonitorEnvelopeCtx`, `TS`, `TSMonitorSet()`, `TSMonitorDefault()`, `VecView()`, `TSMonitorLGSetDisplayVariables()`

# External Links
$(_doc_external("Ts/TSMonitorEnvelopeGetBounds"))
"""
function TSMonitorEnvelopeGetBounds(petsclib::PetscLibType, ts::TS, max::PetscVec, min::PetscVec) end

@for_petsc function TSMonitorEnvelopeGetBounds(petsclib::$UnionPetscLib, ts::TS, max::PetscVec, min::PetscVec )
	max_ = Ref(max.ptr)
	min_ = Ref(min.ptr)

    @chk ccall(
               (:TSMonitorEnvelopeGetBounds, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{CVec}, Ptr{CVec}),
               ts, max_, min_,
              )

	max.ptr = C_NULL
	min.ptr = C_NULL

	return nothing
end 

"""
	TSDMSwarmMonitorMoments(petsclib::PetscLibType,ts::TS, step::PetscInt, t::PetscReal, U::PetscVec, vf::PetscViewerAndFormat) 
Monitors the first three moments of a `DMSWARM` being evolved by the `TS`

Not Collective

Input Parameters:
- `ts`   - the `TS` context
- `step` - current timestep
- `t`    - current time
- `U`    - current solution
- `vf`   - not used

Options Database Key:
- `-ts_dmswarm_monitor_moments`          - Monitor moments of particle distribution
- `-ts_dmswarm_monitor_moments_interval` - Interval of timesteps between monitor outputs

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSMonitorSet()`, `TSMonitorDefault()`, `DMSWARM`

# External Links
$(_doc_external("Ts/TSDMSwarmMonitorMoments"))
"""
function TSDMSwarmMonitorMoments(petsclib::PetscLibType, ts::TS, step::PetscInt, t::PetscReal, U::PetscVec, vf::PetscViewerAndFormat) end

@for_petsc function TSDMSwarmMonitorMoments(petsclib::$UnionPetscLib, ts::TS, step::$PetscInt, t::$PetscReal, U::PetscVec, vf::PetscViewerAndFormat )

    @chk ccall(
               (:TSDMSwarmMonitorMoments, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt, $PetscReal, CVec, Ptr{PetscViewerAndFormat}),
               ts, step, t, U, vf,
              )


	return nothing
end 

"""
	ts::TS = TSCreate(petsclib::PetscLibType,comm::MPI_Comm) 
This function creates an empty timestepper. The problem type can then be set with `TSSetProblemType()` and the
type of solver can then be set with `TSSetType()`.

Collective

Input Parameter:
- `comm` - The communicator

Output Parameter:
- `ts` - The `TS`

Level: beginner

-seealso: [](ch_ts), `TS`, `SNES`, `TSSetType()`, `TSSetUp()`, `TSDestroy()`, `TSSetProblemType()`

# External Links
$(_doc_external("Ts/TSCreate"))
"""
function TSCreate(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function TSCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	ts_ = Ref{CTS}()

    @chk ccall(
               (:TSCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{CTS}),
               comm, ts_,
              )

	ts = TS(ts_[], petsclib)

	return ts
end 

"""
	TSMonitorSPEig(petsclib::PetscLibType,ts::TS, step::PetscInt, ptime::PetscReal, v::PetscVec, monctx::Cvoid) 

# External Links
$(_doc_external("Ts/TSMonitorSPEig"))
"""
function TSMonitorSPEig(petsclib::PetscLibType, ts::TS, step::PetscInt, ptime::PetscReal, v::PetscVec, monctx::Cvoid) end

@for_petsc function TSMonitorSPEig(petsclib::$UnionPetscLib, ts::TS, step::$PetscInt, ptime::$PetscReal, v::PetscVec, monctx::Cvoid )

    @chk ccall(
               (:TSMonitorSPEig, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt, $PetscReal, CVec, Ptr{Cvoid}),
               ts, step, ptime, v, monctx,
              )


	return nothing
end 

"""
	TSRHSSplitSetIS(petsclib::PetscLibType,ts::TS, splitname::String, is::IS) 
Set the index set for the specified split

Logically Collective

Input Parameters:
- `ts`        - the `TS` context obtained from `TSCreate()`
- `splitname` - name of this split, if `NULL` the number of the split is used
- `is`        - the index set for part of the solution vector

Level: intermediate

-seealso: [](ch_ts), `TS`, `IS`, `TSRHSSplitGetIS()`

# External Links
$(_doc_external("Ts/TSRHSSplitSetIS"))
"""
function TSRHSSplitSetIS(petsclib::PetscLibType, ts::TS, splitname::String, is::IS) end

@for_petsc function TSRHSSplitSetIS(petsclib::$UnionPetscLib, ts::TS, splitname::String, is::IS )

    @chk ccall(
               (:TSRHSSplitSetIS, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{Cchar}, CIS),
               ts, splitname, is,
              )


	return nothing
end 

"""
	TSRHSSplitGetIS(petsclib::PetscLibType,ts::TS, splitname::String, is::IS) 
Retrieves the elements for a split as an `IS`

Logically Collective

Input Parameters:
- `ts`        - the `TS` context obtained from `TSCreate()`
- `splitname` - name of this split

Output Parameter:
- `is` - the index set for part of the solution vector

Level: intermediate

-seealso: [](ch_ts), `TS`, `IS`, `TSRHSSplitSetIS()`

# External Links
$(_doc_external("Ts/TSRHSSplitGetIS"))
"""
function TSRHSSplitGetIS(petsclib::PetscLibType, ts::TS, splitname::String, is::IS) end

@for_petsc function TSRHSSplitGetIS(petsclib::$UnionPetscLib, ts::TS, splitname::String, is::IS )
	is_ = Ref(is.ptr)

    @chk ccall(
               (:TSRHSSplitGetIS, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{Cchar}, Ptr{CIS}),
               ts, splitname, is_,
              )

	is.ptr = C_NULL

	return nothing
end 

"""
	TSRHSSplitSetRHSFunction(petsclib::PetscLibType,ts::TS, splitname::String, r::PetscVec, rhsfunc::TSRHSFunctionFn, ctx::Cvoid) 
Set the split right

Logically Collective

Input Parameters:
- `ts`        - the `TS` context obtained from `TSCreate()`
- `splitname` - name of this split
- `r`         - vector to hold the residual (or `NULL` to have it created internally)
- `rhsfunc`   - the RHS function evaluation routine
- `ctx`       - user-defined context for private data for the split function evaluation routine (may be `NULL`)

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSRHSFunctionFn`, `IS`, `TSRHSSplitSetIS()`

# External Links
$(_doc_external("Ts/TSRHSSplitSetRHSFunction"))
"""
function TSRHSSplitSetRHSFunction(petsclib::PetscLibType, ts::TS, splitname::String, r::PetscVec, rhsfunc::TSRHSFunctionFn, ctx::Cvoid) end

@for_petsc function TSRHSSplitSetRHSFunction(petsclib::$UnionPetscLib, ts::TS, splitname::String, r::PetscVec, rhsfunc::TSRHSFunctionFn, ctx::Cvoid )

    @chk ccall(
               (:TSRHSSplitSetRHSFunction, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{Cchar}, CVec, Ptr{TSRHSFunctionFn}, Ptr{Cvoid}),
               ts, splitname, r, rhsfunc, ctx,
              )


	return nothing
end 

"""
	TSRHSSplitSetIFunction(petsclib::PetscLibType,ts::TS, splitname::String, r::PetscVec, ifunc::TSIFunctionFn, ctx::Cvoid) 
Set the split implicit function for `TSARKIMEX`

Logically Collective

Input Parameters:
- `ts`        - the `TS` context obtained from `TSCreate()`
- `splitname` - name of this split
- `r`         - vector to hold the residual (or `NULL` to have it created internally)
- `ifunc`     - the implicit function evaluation routine
- `ctx`       - user-defined context for private data for the split function evaluation routine (may be `NULL`)

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSIFunctionFn`, `IS`, `TSRHSSplitSetIS()`, `TSARKIMEX`

# External Links
$(_doc_external("Ts/TSRHSSplitSetIFunction"))
"""
function TSRHSSplitSetIFunction(petsclib::PetscLibType, ts::TS, splitname::String, r::PetscVec, ifunc::TSIFunctionFn, ctx::Cvoid) end

@for_petsc function TSRHSSplitSetIFunction(petsclib::$UnionPetscLib, ts::TS, splitname::String, r::PetscVec, ifunc::TSIFunctionFn, ctx::Cvoid )

    @chk ccall(
               (:TSRHSSplitSetIFunction, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{Cchar}, CVec, Ptr{TSIFunctionFn}, Ptr{Cvoid}),
               ts, splitname, r, ifunc, ctx,
              )


	return nothing
end 

"""
	TSRHSSplitSetIJacobian(petsclib::PetscLibType,ts::TS, splitname::String, Amat::PetscMat, Pmat::PetscMat, ijac::TSIJacobianFn, ctx::Cvoid) 
Set the Jacobian for the split implicit function with `TSARKIMEX`

Logically Collective

Input Parameters:
- `ts`        - the `TS` context obtained from `TSCreate()`
- `splitname` - name of this split
- `Amat`      - (approximate) matrix to store Jacobian entries computed by `f`
- `Pmat`      - matrix used to compute preconditioner (usually the same as `Amat`)
- `ijac`      - the Jacobian evaluation routine
- `ctx`       - user-defined context for private data for the split function evaluation routine (may be `NULL`)

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSRHSSplitSetIFunction`, `TSIJacobianFn`, `IS`, `TSRHSSplitSetIS()`

# External Links
$(_doc_external("Ts/TSRHSSplitSetIJacobian"))
"""
function TSRHSSplitSetIJacobian(petsclib::PetscLibType, ts::TS, splitname::String, Amat::PetscMat, Pmat::PetscMat, ijac::TSIJacobianFn, ctx::Cvoid) end

@for_petsc function TSRHSSplitSetIJacobian(petsclib::$UnionPetscLib, ts::TS, splitname::String, Amat::PetscMat, Pmat::PetscMat, ijac::TSIJacobianFn, ctx::Cvoid )

    @chk ccall(
               (:TSRHSSplitSetIJacobian, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{Cchar}, CMat, CMat, Ptr{TSIJacobianFn}, Ptr{Cvoid}),
               ts, splitname, Amat, Pmat, ijac, ctx,
              )


	return nothing
end 

"""
	TSRHSSplitGetSubTS(petsclib::PetscLibType,ts::TS, splitname::String, subts::TS) 
Get the sub

Logically Collective

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Output Parameters:
- `splitname` - the number of the split
- `subts`     - the sub-`TS`

Level: advanced

-seealso: [](ch_ts), `TS`, `IS`, `TSGetRHSSplitFunction()`

# External Links
$(_doc_external("Ts/TSRHSSplitGetSubTS"))
"""
function TSRHSSplitGetSubTS(petsclib::PetscLibType, ts::TS, splitname::String, subts::TS) end

@for_petsc function TSRHSSplitGetSubTS(petsclib::$UnionPetscLib, ts::TS, splitname::String, subts::TS )
	subts_ = Ref(subts.ptr)

    @chk ccall(
               (:TSRHSSplitGetSubTS, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{Cchar}, Ptr{CTS}),
               ts, splitname, subts_,
              )

	subts.ptr = C_NULL

	return nothing
end 

"""
	n::PetscInt = TSRHSSplitGetSubTSs(petsclib::PetscLibType,ts::TS, subts::Vector{TS}) 
Get an array of all sub

Logically Collective

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Output Parameters:
- `n`     - the number of splits
- `subts` - the array of `TS` contexts

Level: advanced

-seealso: [](ch_ts), `TS`, `IS`, `TSGetRHSSplitFunction()`

# External Links
$(_doc_external("Ts/TSRHSSplitGetSubTSs"))
"""
function TSRHSSplitGetSubTSs(petsclib::PetscLibType, ts::TS, subts::Vector{TS}) end

@for_petsc function TSRHSSplitGetSubTSs(petsclib::$UnionPetscLib, ts::TS, subts::Vector{TS} )
	n_ = Ref{$PetscInt}()
	subts_ = Ref(pointer(subts))

    @chk ccall(
               (:TSRHSSplitGetSubTSs, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{$PetscInt}, Ptr{Ptr{CTS}}),
               ts, n_, subts_,
              )

	n = n_[]

	return n
end 

"""
	TSRHSSplitGetSNES(petsclib::PetscLibType,ts::TS, snes::PetscSNES) 
Returns the `SNES` (nonlinear solver) associated with
a `TS` (timestepper) context when RHS splits are used.

Not Collective, but snes is parallel if ts is parallel

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Output Parameter:
- `snes` - the nonlinear solver context

Level: intermediate

-seealso: [](ch_ts), `TS`, `SNES`, `TSCreate()`, `TSRHSSplitSetSNES()`

# External Links
$(_doc_external("Ts/TSRHSSplitGetSNES"))
"""
function TSRHSSplitGetSNES(petsclib::PetscLibType, ts::TS, snes::PetscSNES) end

@for_petsc function TSRHSSplitGetSNES(petsclib::$UnionPetscLib, ts::TS, snes::PetscSNES )
	snes_ = Ref(snes.ptr)

    @chk ccall(
               (:TSRHSSplitGetSNES, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{CSNES}),
               ts, snes_,
              )

	snes.ptr = C_NULL

	return nothing
end 

"""
	TSRHSSplitSetSNES(petsclib::PetscLibType,ts::TS, snes::PetscSNES) 
Set the `SNES` (nonlinear solver) to be used by the
timestepping context when RHS splits are used.

Collective

Input Parameters:
- `ts`   - the `TS` context obtained from `TSCreate()`
- `snes` - the nonlinear solver context

Level: intermediate

-seealso: [](ch_ts), `TS`, `SNES`, `TSCreate()`, `TSRHSSplitGetSNES()`

# External Links
$(_doc_external("Ts/TSRHSSplitSetSNES"))
"""
function TSRHSSplitSetSNES(petsclib::PetscLibType, ts::TS, snes::PetscSNES) end

@for_petsc function TSRHSSplitSetSNES(petsclib::$UnionPetscLib, ts::TS, snes::PetscSNES )

    @chk ccall(
               (:TSRHSSplitSetSNES, $petsc_library),
               PetscErrorCode,
               (CTS, CSNES),
               ts, snes,
              )


	return nothing
end 

"""
	TSSetRHSJacobianP(petsclib::PetscLibType,ts::TS, Amat::PetscMat, func::TSRHSJacobianPFn, ctx::Cvoid) 
Sets the function that computes the Jacobian of G w.r.t. the parameters p where U_t = G(U,p,t), as well as the location to store the matrix.

Logically Collective

Input Parameters:
- `ts`   - `TS` context obtained from `TSCreate()`
- `Amat` - JacobianP matrix
- `func` - function
- `ctx`  - [optional] user-defined function context

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSRHSJacobianPFn`, `TSGetRHSJacobianP()`

# External Links
$(_doc_external("Ts/TSSetRHSJacobianP"))
"""
function TSSetRHSJacobianP(petsclib::PetscLibType, ts::TS, Amat::PetscMat, func::TSRHSJacobianPFn, ctx::Cvoid) end

@for_petsc function TSSetRHSJacobianP(petsclib::$UnionPetscLib, ts::TS, Amat::PetscMat, func::TSRHSJacobianPFn, ctx::Cvoid )

    @chk ccall(
               (:TSSetRHSJacobianP, $petsc_library),
               PetscErrorCode,
               (CTS, CMat, Ptr{TSRHSJacobianPFn}, Ptr{Cvoid}),
               ts, Amat, func, ctx,
              )


	return nothing
end 

"""
	TSGetRHSJacobianP(petsclib::PetscLibType,ts::TS, Amat::PetscMat, func::TSRHSJacobianPFn, ctx::Cvoid) 
Gets the function that computes the Jacobian of G  w.r.t. the parameters p where  U_t = G(U,p,t), as well as the location to store the matrix.

Logically Collective

Input Parameter:
- `ts` - `TS` context obtained from `TSCreate()`

Output Parameters:
- `Amat` - JacobianP matrix
- `func` - function
- `ctx`  - [optional] user-defined function context

Level: intermediate

-seealso: [](ch_ts), `TSSetRHSJacobianP()`, `TS`, `TSRHSJacobianPFn`

# External Links
$(_doc_external("Ts/TSGetRHSJacobianP"))
"""
function TSGetRHSJacobianP(petsclib::PetscLibType, ts::TS, Amat::PetscMat, func::TSRHSJacobianPFn, ctx::Cvoid) end

@for_petsc function TSGetRHSJacobianP(petsclib::$UnionPetscLib, ts::TS, Amat::PetscMat, func::TSRHSJacobianPFn, ctx::Cvoid )
	Amat_ = Ref(Amat.ptr)

    @chk ccall(
               (:TSGetRHSJacobianP, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{CMat}, TSRHSJacobianPFn, Cvoid),
               ts, Amat_, func, ctx,
              )

	Amat.ptr = C_NULL

	return nothing
end 

"""
	TSComputeRHSJacobianP(petsclib::PetscLibType,ts::TS, t::PetscReal, U::PetscVec, Amat::PetscMat) 
Runs the user

Collective

Input Parameters:
- `ts` - The `TS` context obtained from `TSCreate()`
- `t`  - the time
- `U`  - the solution at which to compute the Jacobian

Output Parameter:
- `Amat` - the computed Jacobian

Level: developer

-seealso: [](ch_ts), `TSSetRHSJacobianP()`, `TS`

# External Links
$(_doc_external("Ts/TSComputeRHSJacobianP"))
"""
function TSComputeRHSJacobianP(petsclib::PetscLibType, ts::TS, t::PetscReal, U::PetscVec, Amat::PetscMat) end

@for_petsc function TSComputeRHSJacobianP(petsclib::$UnionPetscLib, ts::TS, t::$PetscReal, U::PetscVec, Amat::PetscMat )

    @chk ccall(
               (:TSComputeRHSJacobianP, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal, CVec, CMat),
               ts, t, U, Amat,
              )


	return nothing
end 

"""
	TSSetIJacobianP(petsclib::PetscLibType,ts::TS, Amat::PetscMat, func::external, ctx::Cvoid) 
Sets the function that computes the Jacobian of F w.r.t. the parameters p where F(Udot,U,p,t) = G(U,p,t), as well as the location to store the matrix.

Logically Collective

Input Parameters:
- `ts`   - `TS` context obtained from `TSCreate()`
- `Amat` - JacobianP matrix
- `func` - function
- `ctx`  - [optional] user-defined function context

Calling sequence of `func`:
- `ts`    - the `TS` context
- `t`     - current timestep
- `U`     - input vector (current ODE solution)
- `Udot`  - time derivative of state vector
- `shift` - shift to apply, see the note in `TSSetIJacobian()`
- `A`     - output matrix
- `ctx`   - [optional] user-defined function context

Level: intermediate

-seealso: [](ch_ts), `TSSetRHSJacobianP()`, `TS`

# External Links
$(_doc_external("Ts/TSSetIJacobianP"))
"""
function TSSetIJacobianP(petsclib::PetscLibType, ts::TS, Amat::PetscMat, func::external, ctx::Cvoid) end

@for_petsc function TSSetIJacobianP(petsclib::$UnionPetscLib, ts::TS, Amat::PetscMat, func::external, ctx::Cvoid )

    @chk ccall(
               (:TSSetIJacobianP, $petsc_library),
               PetscErrorCode,
               (CTS, CMat, external, Ptr{Cvoid}),
               ts, Amat, func, ctx,
              )


	return nothing
end 

"""
	TSComputeIJacobianP(petsclib::PetscLibType,ts::TS, t::PetscReal, U::PetscVec, Udot::PetscVec, shift::PetscReal, Amat::PetscMat, imex::PetscBool) 
Runs the user

Collective

Input Parameters:
- `ts`    - the `TS` context
- `t`     - current timestep
- `U`     - state vector
- `Udot`  - time derivative of state vector
- `shift` - shift to apply, see note below
- `imex`  - flag indicates if the method is IMEX so that the `RHSJacobianP` should be kept separate

Output Parameter:
- `Amat` - Jacobian matrix

Level: developer

-seealso: [](ch_ts), `TS`, `TSSetIJacobianP()`

# External Links
$(_doc_external("Ts/TSComputeIJacobianP"))
"""
function TSComputeIJacobianP(petsclib::PetscLibType, ts::TS, t::PetscReal, U::PetscVec, Udot::PetscVec, shift::PetscReal, Amat::PetscMat, imex::PetscBool) end

@for_petsc function TSComputeIJacobianP(petsclib::$UnionPetscLib, ts::TS, t::$PetscReal, U::PetscVec, Udot::PetscVec, shift::$PetscReal, Amat::PetscMat, imex::PetscBool )

    @chk ccall(
               (:TSComputeIJacobianP, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal, CVec, CVec, $PetscReal, CMat, PetscBool),
               ts, t, U, Udot, shift, Amat, imex,
              )


	return nothing
end 

"""
	TSGetCostIntegral(petsclib::PetscLibType,ts::TS, v::PetscVec) 
Returns the values of the integral term in the cost functions.
It is valid to call the routine after a backward run.

Not Collective

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Output Parameter:
- `v` - the vector containing the integrals for each cost function

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSAdjointSolve()`, `TSSetCostIntegrand()`

# External Links
$(_doc_external("Ts/TSGetCostIntegral"))
"""
function TSGetCostIntegral(petsclib::PetscLibType, ts::TS, v::PetscVec) end

@for_petsc function TSGetCostIntegral(petsclib::$UnionPetscLib, ts::TS, v::PetscVec )
	v_ = Ref(v.ptr)

    @chk ccall(
               (:TSGetCostIntegral, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{CVec}),
               ts, v_,
              )

	v.ptr = C_NULL

	return nothing
end 

"""
	TSComputeCostIntegrand(petsclib::PetscLibType,ts::TS, t::PetscReal, U::PetscVec, Q::PetscVec) 
Evaluates the integral function in the cost functions.

Input Parameters:
- `ts` - the `TS` context
- `t`  - current time
- `U`  - state vector, i.e. current solution

Output Parameter:
- `Q` - vector of size numcost to hold the outputs

Level: deprecated

-seealso: [](ch_ts), `TS`, `TSAdjointSolve()`, `TSSetCostIntegrand()`

# External Links
$(_doc_external("Ts/TSComputeCostIntegrand"))
"""
function TSComputeCostIntegrand(petsclib::PetscLibType, ts::TS, t::PetscReal, U::PetscVec, Q::PetscVec) end

@for_petsc function TSComputeCostIntegrand(petsclib::$UnionPetscLib, ts::TS, t::$PetscReal, U::PetscVec, Q::PetscVec )

    @chk ccall(
               (:TSComputeCostIntegrand, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal, CVec, CVec),
               ts, t, U, Q,
              )


	return nothing
end 

"""
	TSSetIHessianProduct(petsclib::PetscLibType,ts::TS, ihp1::PetscVec, ihessianproductfunc1::external, ihp2::PetscVec, ihessianproductfunc2::external, ihp3::PetscVec, ihessianproductfunc3::external, ihp4::PetscVec, ihessianproductfunc4::external, ctx::Cvoid) 
Sets the function that computes the vector

Logically Collective

Input Parameters:
- `ts`   - `TS` context obtained from `TSCreate()`
- `ihp1` - an array of vectors storing the result of vector-Hessian-vector product for F_UU
- `hessianproductfunc1` - vector-Hessian-vector product function for F_UU
- `ihp2` - an array of vectors storing the result of vector-Hessian-vector product for F_UP
- `hessianproductfunc2` - vector-Hessian-vector product function for F_UP
- `ihp3` - an array of vectors storing the result of vector-Hessian-vector product for F_PU
- `hessianproductfunc3` - vector-Hessian-vector product function for F_PU
- `ihp4` - an array of vectors storing the result of vector-Hessian-vector product for F_PP
- `hessianproductfunc4` - vector-Hessian-vector product function for F_PP

Calling sequence of `ihessianproductfunc1`:
- `ts`  - the `TS` context
- `t`   - current timestep
- `U`   - input vector (current ODE solution)
- `Vl`  - an array of input vectors to be left-multiplied with the Hessian
- `Vr`  - input vector to be right-multiplied with the Hessian
- `VHV` - an array of output vectors for vector-Hessian-vector product
- `ctx` - [optional] user-defined function context

Level: intermediate

-seealso: [](ch_ts), `TS`

# External Links
$(_doc_external("Ts/TSSetIHessianProduct"))
"""
function TSSetIHessianProduct(petsclib::PetscLibType, ts::TS, ihp1::PetscVec, ihessianproductfunc1::external, ihp2::PetscVec, ihessianproductfunc2::external, ihp3::PetscVec, ihessianproductfunc3::external, ihp4::PetscVec, ihessianproductfunc4::external, ctx::Cvoid) end

@for_petsc function TSSetIHessianProduct(petsclib::$UnionPetscLib, ts::TS, ihp1::PetscVec, ihessianproductfunc1::external, ihp2::PetscVec, ihessianproductfunc2::external, ihp3::PetscVec, ihessianproductfunc3::external, ihp4::PetscVec, ihessianproductfunc4::external, ctx::Cvoid )
	ihp1_ = Ref(ihp1.ptr)
	ihp2_ = Ref(ihp2.ptr)
	ihp3_ = Ref(ihp3.ptr)
	ihp4_ = Ref(ihp4.ptr)

    @chk ccall(
               (:TSSetIHessianProduct, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{CVec}, external, Ptr{CVec}, external, Ptr{CVec}, external, Ptr{CVec}, external, Ptr{Cvoid}),
               ts, ihp1_, ihessianproductfunc1, ihp2_, ihessianproductfunc2, ihp3_, ihessianproductfunc3, ihp4_, ihessianproductfunc4, ctx,
              )

	ihp1.ptr = C_NULL
	ihp2.ptr = C_NULL
	ihp3.ptr = C_NULL
	ihp4.ptr = C_NULL

	return nothing
end 

"""
	TSComputeIHessianProductFunctionUU(petsclib::PetscLibType,ts::TS, t::PetscReal, U::PetscVec, Vl::Vector{PetscVec}, Vr::PetscVec, VHV::Vector{PetscVec}) 
Runs the user

Collective

Input Parameters:
- `ts` - The `TS` context obtained from `TSCreate()`
- `t`  - the time
- `U`  - the solution at which to compute the Hessian product
- `Vl` - the array of input vectors to be multiplied with the Hessian from the left
- `Vr` - the input vector to be multiplied with the Hessian from the right

Output Parameter:
- `VHV` - the array of output vectors that store the Hessian product

Level: developer

-seealso: [](ch_ts), `TSSetIHessianProduct()`

# External Links
$(_doc_external("Ts/TSComputeIHessianProductFunctionUU"))
"""
function TSComputeIHessianProductFunctionUU(petsclib::PetscLibType, ts::TS, t::PetscReal, U::PetscVec, Vl::Vector{PetscVec}, Vr::PetscVec, VHV::Vector{PetscVec}) end

@for_petsc function TSComputeIHessianProductFunctionUU(petsclib::$UnionPetscLib, ts::TS, t::$PetscReal, U::PetscVec, Vl::Vector{PetscVec}, Vr::PetscVec, VHV::Vector{PetscVec} )

    @chk ccall(
               (:TSComputeIHessianProductFunctionUU, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal, CVec, Ptr{CVec}, CVec, Ptr{CVec}),
               ts, t, U, Vl, Vr, VHV,
              )


	return nothing
end 

"""
	TSComputeIHessianProductFunctionUP(petsclib::PetscLibType,ts::TS, t::PetscReal, U::PetscVec, Vl::Vector{PetscVec}, Vr::PetscVec, VHV::Vector{PetscVec}) 
Runs the user

Collective

Input Parameters:
- `ts` - The `TS` context obtained from `TSCreate()`
- `t`  - the time
- `U`  - the solution at which to compute the Hessian product
- `Vl` - the array of input vectors to be multiplied with the Hessian from the left
- `Vr` - the input vector to be multiplied with the Hessian from the right

Output Parameter:
- `VHV` - the array of output vectors that store the Hessian product

Level: developer

-seealso: [](ch_ts), `TSSetIHessianProduct()`

# External Links
$(_doc_external("Ts/TSComputeIHessianProductFunctionUP"))
"""
function TSComputeIHessianProductFunctionUP(petsclib::PetscLibType, ts::TS, t::PetscReal, U::PetscVec, Vl::Vector{PetscVec}, Vr::PetscVec, VHV::Vector{PetscVec}) end

@for_petsc function TSComputeIHessianProductFunctionUP(petsclib::$UnionPetscLib, ts::TS, t::$PetscReal, U::PetscVec, Vl::Vector{PetscVec}, Vr::PetscVec, VHV::Vector{PetscVec} )

    @chk ccall(
               (:TSComputeIHessianProductFunctionUP, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal, CVec, Ptr{CVec}, CVec, Ptr{CVec}),
               ts, t, U, Vl, Vr, VHV,
              )


	return nothing
end 

"""
	TSComputeIHessianProductFunctionPU(petsclib::PetscLibType,ts::TS, t::PetscReal, U::PetscVec, Vl::Vector{PetscVec}, Vr::PetscVec, VHV::Vector{PetscVec}) 
Runs the user

Collective

Input Parameters:
- `ts` - The `TS` context obtained from `TSCreate()`
- `t`  - the time
- `U`  - the solution at which to compute the Hessian product
- `Vl` - the array of input vectors to be multiplied with the Hessian from the left
- `Vr` - the input vector to be multiplied with the Hessian from the right

Output Parameter:
- `VHV` - the array of output vectors that store the Hessian product

Level: developer

-seealso: [](ch_ts), `TSSetIHessianProduct()`

# External Links
$(_doc_external("Ts/TSComputeIHessianProductFunctionPU"))
"""
function TSComputeIHessianProductFunctionPU(petsclib::PetscLibType, ts::TS, t::PetscReal, U::PetscVec, Vl::Vector{PetscVec}, Vr::PetscVec, VHV::Vector{PetscVec}) end

@for_petsc function TSComputeIHessianProductFunctionPU(petsclib::$UnionPetscLib, ts::TS, t::$PetscReal, U::PetscVec, Vl::Vector{PetscVec}, Vr::PetscVec, VHV::Vector{PetscVec} )

    @chk ccall(
               (:TSComputeIHessianProductFunctionPU, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal, CVec, Ptr{CVec}, CVec, Ptr{CVec}),
               ts, t, U, Vl, Vr, VHV,
              )


	return nothing
end 

"""
	TSComputeIHessianProductFunctionPP(petsclib::PetscLibType,ts::TS, t::PetscReal, U::PetscVec, Vl::Vector{PetscVec}, Vr::PetscVec, VHV::Vector{PetscVec}) 
Runs the user

Collective

Input Parameters:
- `ts` - The `TS` context obtained from `TSCreate()`
- `t`  - the time
- `U`  - the solution at which to compute the Hessian product
- `Vl` - the array of input vectors to be multiplied with the Hessian from the left
- `Vr` - the input vector to be multiplied with the Hessian from the right

Output Parameter:
- `VHV` - the array of output vectors that store the Hessian product

Level: developer

-seealso: [](ch_ts), `TSSetIHessianProduct()`

# External Links
$(_doc_external("Ts/TSComputeIHessianProductFunctionPP"))
"""
function TSComputeIHessianProductFunctionPP(petsclib::PetscLibType, ts::TS, t::PetscReal, U::PetscVec, Vl::Vector{PetscVec}, Vr::PetscVec, VHV::Vector{PetscVec}) end

@for_petsc function TSComputeIHessianProductFunctionPP(petsclib::$UnionPetscLib, ts::TS, t::$PetscReal, U::PetscVec, Vl::Vector{PetscVec}, Vr::PetscVec, VHV::Vector{PetscVec} )

    @chk ccall(
               (:TSComputeIHessianProductFunctionPP, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal, CVec, Ptr{CVec}, CVec, Ptr{CVec}),
               ts, t, U, Vl, Vr, VHV,
              )


	return nothing
end 

"""
	TSSetRHSHessianProduct(petsclib::PetscLibType,ts::TS, rhshp1::Vector{PetscVec}, rhshessianproductfunc1::external, rhshp2::Vector{PetscVec}, rhshessianproductfunc2::external, rhshp3::Vector{PetscVec}, rhshessianproductfunc3::external, rhshp4::Vector{PetscVec}, rhshessianproductfunc4::external, ctx::Cvoid) 
Sets the function that computes the vector
product. The Hessian is the second-order derivative of G (RHSFunction) w.r.t. the state
variable.

Logically Collective

Input Parameters:
- `ts`     - `TS` context obtained from `TSCreate()`
- `rhshp1` - an array of vectors storing the result of vector-Hessian-vector product for G_UU
- `hessianproductfunc1` - vector-Hessian-vector product function for G_UU
- `rhshp2` - an array of vectors storing the result of vector-Hessian-vector product for G_UP
- `hessianproductfunc2` - vector-Hessian-vector product function for G_UP
- `rhshp3` - an array of vectors storing the result of vector-Hessian-vector product for G_PU
- `hessianproductfunc3` - vector-Hessian-vector product function for G_PU
- `rhshp4` - an array of vectors storing the result of vector-Hessian-vector product for G_PP
- `hessianproductfunc4` - vector-Hessian-vector product function for G_PP
- `ctx`    - [optional] user-defined function context

Calling sequence of `rhshessianproductfunc1`:
- `ts`  - the `TS` context
- `t`   - current timestep
- `U`   - input vector (current ODE solution)
- `Vl`  - an array of input vectors to be left-multiplied with the Hessian
- `Vr`  - input vector to be right-multiplied with the Hessian
- `VHV` - an array of output vectors for vector-Hessian-vector product
- `ctx` - [optional] user-defined function context

Level: intermediate

-seealso: `TS`, `TSAdjoint`

# External Links
$(_doc_external("Ts/TSSetRHSHessianProduct"))
"""
function TSSetRHSHessianProduct(petsclib::PetscLibType, ts::TS, rhshp1::Vector{PetscVec}, rhshessianproductfunc1::external, rhshp2::Vector{PetscVec}, rhshessianproductfunc2::external, rhshp3::Vector{PetscVec}, rhshessianproductfunc3::external, rhshp4::Vector{PetscVec}, rhshessianproductfunc4::external, ctx::Cvoid) end

@for_petsc function TSSetRHSHessianProduct(petsclib::$UnionPetscLib, ts::TS, rhshp1::Vector{PetscVec}, rhshessianproductfunc1::external, rhshp2::Vector{PetscVec}, rhshessianproductfunc2::external, rhshp3::Vector{PetscVec}, rhshessianproductfunc3::external, rhshp4::Vector{PetscVec}, rhshessianproductfunc4::external, ctx::Cvoid )

    @chk ccall(
               (:TSSetRHSHessianProduct, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{CVec}, external, Ptr{CVec}, external, Ptr{CVec}, external, Ptr{CVec}, external, Ptr{Cvoid}),
               ts, rhshp1, rhshessianproductfunc1, rhshp2, rhshessianproductfunc2, rhshp3, rhshessianproductfunc3, rhshp4, rhshessianproductfunc4, ctx,
              )


	return nothing
end 

"""
	TSComputeRHSHessianProductFunctionUU(petsclib::PetscLibType,ts::TS, t::PetscReal, U::PetscVec, Vl::Vector{PetscVec}, Vr::PetscVec, VHV::Vector{PetscVec}) 
Runs the user

Collective

Input Parameters:
- `ts` - The `TS` context obtained from `TSCreate()`
- `t`  - the time
- `U`  - the solution at which to compute the Hessian product
- `Vl` - the array of input vectors to be multiplied with the Hessian from the left
- `Vr` - the input vector to be multiplied with the Hessian from the right

Output Parameter:
- `VHV` - the array of output vectors that store the Hessian product

Level: developer

-seealso: [](ch_ts), `TS`, `TSSetRHSHessianProduct()`

# External Links
$(_doc_external("Ts/TSComputeRHSHessianProductFunctionUU"))
"""
function TSComputeRHSHessianProductFunctionUU(petsclib::PetscLibType, ts::TS, t::PetscReal, U::PetscVec, Vl::Vector{PetscVec}, Vr::PetscVec, VHV::Vector{PetscVec}) end

@for_petsc function TSComputeRHSHessianProductFunctionUU(petsclib::$UnionPetscLib, ts::TS, t::$PetscReal, U::PetscVec, Vl::Vector{PetscVec}, Vr::PetscVec, VHV::Vector{PetscVec} )

    @chk ccall(
               (:TSComputeRHSHessianProductFunctionUU, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal, CVec, Ptr{CVec}, CVec, Ptr{CVec}),
               ts, t, U, Vl, Vr, VHV,
              )


	return nothing
end 

"""
	TSComputeRHSHessianProductFunctionUP(petsclib::PetscLibType,ts::TS, t::PetscReal, U::PetscVec, Vl::Vector{PetscVec}, Vr::PetscVec, VHV::Vector{PetscVec}) 
Runs the user

Collective

Input Parameters:
- `ts` - The `TS` context obtained from `TSCreate()`
- `t`  - the time
- `U`  - the solution at which to compute the Hessian product
- `Vl` - the array of input vectors to be multiplied with the Hessian from the left
- `Vr` - the input vector to be multiplied with the Hessian from the right

Output Parameter:
- `VHV` - the array of output vectors that store the Hessian product

Level: developer

-seealso: [](ch_ts), `TS`, `TSSetRHSHessianProduct()`

# External Links
$(_doc_external("Ts/TSComputeRHSHessianProductFunctionUP"))
"""
function TSComputeRHSHessianProductFunctionUP(petsclib::PetscLibType, ts::TS, t::PetscReal, U::PetscVec, Vl::Vector{PetscVec}, Vr::PetscVec, VHV::Vector{PetscVec}) end

@for_petsc function TSComputeRHSHessianProductFunctionUP(petsclib::$UnionPetscLib, ts::TS, t::$PetscReal, U::PetscVec, Vl::Vector{PetscVec}, Vr::PetscVec, VHV::Vector{PetscVec} )

    @chk ccall(
               (:TSComputeRHSHessianProductFunctionUP, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal, CVec, Ptr{CVec}, CVec, Ptr{CVec}),
               ts, t, U, Vl, Vr, VHV,
              )


	return nothing
end 

"""
	TSComputeRHSHessianProductFunctionPU(petsclib::PetscLibType,ts::TS, t::PetscReal, U::PetscVec, Vl::Vector{PetscVec}, Vr::PetscVec, VHV::Vector{PetscVec}) 
Runs the user

Collective

Input Parameters:
- `ts` - The `TS` context obtained from `TSCreate()`
- `t`  - the time
- `U`  - the solution at which to compute the Hessian product
- `Vl` - the array of input vectors to be multiplied with the Hessian from the left
- `Vr` - the input vector to be multiplied with the Hessian from the right

Output Parameter:
- `VHV` - the array of output vectors that store the Hessian product

Level: developer

-seealso: [](ch_ts), `TSSetRHSHessianProduct()`

# External Links
$(_doc_external("Ts/TSComputeRHSHessianProductFunctionPU"))
"""
function TSComputeRHSHessianProductFunctionPU(petsclib::PetscLibType, ts::TS, t::PetscReal, U::PetscVec, Vl::Vector{PetscVec}, Vr::PetscVec, VHV::Vector{PetscVec}) end

@for_petsc function TSComputeRHSHessianProductFunctionPU(petsclib::$UnionPetscLib, ts::TS, t::$PetscReal, U::PetscVec, Vl::Vector{PetscVec}, Vr::PetscVec, VHV::Vector{PetscVec} )

    @chk ccall(
               (:TSComputeRHSHessianProductFunctionPU, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal, CVec, Ptr{CVec}, CVec, Ptr{CVec}),
               ts, t, U, Vl, Vr, VHV,
              )


	return nothing
end 

"""
	TSComputeRHSHessianProductFunctionPP(petsclib::PetscLibType,ts::TS, t::PetscReal, U::PetscVec, Vl::Vector{PetscVec}, Vr::PetscVec, VHV::Vector{PetscVec}) 
Runs the user

Collective

Input Parameters:
- `ts` - The `TS` context obtained from `TSCreate()`
- `t`  - the time
- `U`  - the solution at which to compute the Hessian product
- `Vl` - the array of input vectors to be multiplied with the Hessian from the left
- `Vr` - the input vector to be multiplied with the Hessian from the right

Output Parameter:
- `VHV` - the array of output vectors that store the Hessian product

Level: developer

-seealso: [](ch_ts), `TSSetRHSHessianProduct()`

# External Links
$(_doc_external("Ts/TSComputeRHSHessianProductFunctionPP"))
"""
function TSComputeRHSHessianProductFunctionPP(petsclib::PetscLibType, ts::TS, t::PetscReal, U::PetscVec, Vl::Vector{PetscVec}, Vr::PetscVec, VHV::Vector{PetscVec}) end

@for_petsc function TSComputeRHSHessianProductFunctionPP(petsclib::$UnionPetscLib, ts::TS, t::$PetscReal, U::PetscVec, Vl::Vector{PetscVec}, Vr::PetscVec, VHV::Vector{PetscVec} )

    @chk ccall(
               (:TSComputeRHSHessianProductFunctionPP, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal, CVec, Ptr{CVec}, CVec, Ptr{CVec}),
               ts, t, U, Vl, Vr, VHV,
              )


	return nothing
end 

"""
	TSSetCostGradients(petsclib::PetscLibType,ts::TS, numcost::PetscInt, lambda::Vector{PetscVec}, mu::Vector{PetscVec}) 
Sets the initial value of the gradients of the cost function w.r.t. initial values and w.r.t. the problem parameters
for use by the `TS` adjoint routines.

Logically Collective

Input Parameters:
- `ts`      - the `TS` context obtained from `TSCreate()`
- `numcost` - number of gradients to be computed, this is the number of cost functions
- `lambda`  - gradients with respect to the initial condition variables, the dimension and parallel layout of these vectors is the same as the ODE solution vector
- `mu`      - gradients with respect to the parameters, the number of entries in these vectors is the same as the number of parameters

Level: beginner

-seealso: `TS`, `TSAdjointSolve()`, `TSGetCostGradients()`

# External Links
$(_doc_external("Ts/TSSetCostGradients"))
"""
function TSSetCostGradients(petsclib::PetscLibType, ts::TS, numcost::PetscInt, lambda::Vector{PetscVec}, mu::Vector{PetscVec}) end

@for_petsc function TSSetCostGradients(petsclib::$UnionPetscLib, ts::TS, numcost::$PetscInt, lambda::Vector{PetscVec}, mu::Vector{PetscVec} )

    @chk ccall(
               (:TSSetCostGradients, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt, Ptr{CVec}, Ptr{CVec}),
               ts, numcost, lambda, mu,
              )


	return nothing
end 

"""
	numcost::PetscInt = TSGetCostGradients(petsclib::PetscLibType,ts::TS, lambda::Vector{PetscVec}, mu::Vector{PetscVec}) 
Returns the gradients from the `TSAdjointSolve()`

Not Collective, but the vectors returned are parallel if `TS` is parallel

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Output Parameters:
- `numcost` - size of returned arrays
- `lambda`  - vectors containing the gradients of the cost functions with respect to the ODE/DAE solution variables
- `mu`      - vectors containing the gradients of the cost functions with respect to the problem parameters

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSAdjointSolve()`, `TSSetCostGradients()`

# External Links
$(_doc_external("Ts/TSGetCostGradients"))
"""
function TSGetCostGradients(petsclib::PetscLibType, ts::TS, lambda::Vector{PetscVec}, mu::Vector{PetscVec}) end

@for_petsc function TSGetCostGradients(petsclib::$UnionPetscLib, ts::TS, lambda::Vector{PetscVec}, mu::Vector{PetscVec} )
	numcost_ = Ref{$PetscInt}()
	lambda_ = Ref(pointer(lambda))
	mu_ = Ref(pointer(mu))

    @chk ccall(
               (:TSGetCostGradients, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{$PetscInt}, Ptr{Ptr{CVec}}, Ptr{Ptr{CVec}}),
               ts, numcost_, lambda_, mu_,
              )

	numcost = numcost_[]

	return numcost
end 

"""
	TSSetCostHessianProducts(petsclib::PetscLibType,ts::TS, numcost::PetscInt, lambda2::Vector{PetscVec}, mu2::Vector{PetscVec}, dir::PetscVec) 
Sets the initial value of the Hessian
for use by the `TS` adjoint routines.

Logically Collective

Input Parameters:
- `ts`      - the `TS` context obtained from `TSCreate()`
- `numcost` - number of cost functions
- `lambda2` - Hessian-vector product with respect to the initial condition variables, the dimension and parallel layout of these vectors is the same as the ODE solution vector
- `mu2`     - Hessian-vector product with respect to the parameters, the number of entries in these vectors is the same as the number of parameters
- `dir`     - the direction vector that are multiplied with the Hessian of the cost functions

Level: beginner

-seealso: [](ch_ts), `TS`, `TSAdjointSolve()`, `TSAdjointSetForward()`

# External Links
$(_doc_external("Ts/TSSetCostHessianProducts"))
"""
function TSSetCostHessianProducts(petsclib::PetscLibType, ts::TS, numcost::PetscInt, lambda2::Vector{PetscVec}, mu2::Vector{PetscVec}, dir::PetscVec) end

@for_petsc function TSSetCostHessianProducts(petsclib::$UnionPetscLib, ts::TS, numcost::$PetscInt, lambda2::Vector{PetscVec}, mu2::Vector{PetscVec}, dir::PetscVec )

    @chk ccall(
               (:TSSetCostHessianProducts, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt, Ptr{CVec}, Ptr{CVec}, CVec),
               ts, numcost, lambda2, mu2, dir,
              )


	return nothing
end 

"""
	numcost::PetscInt = TSGetCostHessianProducts(petsclib::PetscLibType,ts::TS, lambda2::Vector{PetscVec}, mu2::Vector{PetscVec}, dir::PetscVec) 
Returns the gradients from the `TSAdjointSolve()`

Not Collective, but vectors returned are parallel if `TS` is parallel

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Output Parameters:
- `numcost` - number of cost functions
- `lambda2` - Hessian-vector product with respect to the initial condition variables, the dimension and parallel layout of these vectors is the same as the ODE solution vector
- `mu2`     - Hessian-vector product with respect to the parameters, the number of entries in these vectors is the same as the number of parameters
- `dir`     - the direction vector that are multiplied with the Hessian of the cost functions

Level: intermediate

-seealso: [](ch_ts), `TSAdjointSolve()`, `TSSetCostHessianProducts()`

# External Links
$(_doc_external("Ts/TSGetCostHessianProducts"))
"""
function TSGetCostHessianProducts(petsclib::PetscLibType, ts::TS, lambda2::Vector{PetscVec}, mu2::Vector{PetscVec}, dir::PetscVec) end

@for_petsc function TSGetCostHessianProducts(petsclib::$UnionPetscLib, ts::TS, lambda2::Vector{PetscVec}, mu2::Vector{PetscVec}, dir::PetscVec )
	numcost_ = Ref{$PetscInt}()
	lambda2_ = Ref(pointer(lambda2))
	mu2_ = Ref(pointer(mu2))
	dir_ = Ref(dir.ptr)

    @chk ccall(
               (:TSGetCostHessianProducts, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{$PetscInt}, Ptr{Ptr{CVec}}, Ptr{Ptr{CVec}}, Ptr{CVec}),
               ts, numcost_, lambda2_, mu2_, dir_,
              )

	numcost = numcost_[]
	dir.ptr = C_NULL

	return numcost
end 

"""
	TSAdjointSetForward(petsclib::PetscLibType,ts::TS, didp::PetscMat) 
Trigger the tangent linear solver and initialize the forward sensitivities

Logically Collective

Input Parameters:
- `ts`   - the `TS` context obtained from `TSCreate()`
- `didp` - the derivative of initial values w.r.t. parameters

Level: intermediate

-seealso: [](ch_ts), `TSAdjointSolve()`, `TSSetCostHessianProducts()`, `TSAdjointResetForward()`

# External Links
$(_doc_external("Ts/TSAdjointSetForward"))
"""
function TSAdjointSetForward(petsclib::PetscLibType, ts::TS, didp::PetscMat) end

@for_petsc function TSAdjointSetForward(petsclib::$UnionPetscLib, ts::TS, didp::PetscMat )

    @chk ccall(
               (:TSAdjointSetForward, $petsc_library),
               PetscErrorCode,
               (CTS, CMat),
               ts, didp,
              )


	return nothing
end 

"""
	TSAdjointResetForward(petsclib::PetscLibType,ts::TS) 
Reset the tangent linear solver and destroy the tangent linear context

Logically Collective

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Level: intermediate

-seealso: [](ch_ts), `TSAdjointSetForward()`

# External Links
$(_doc_external("Ts/TSAdjointResetForward"))
"""
function TSAdjointResetForward(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSAdjointResetForward(petsclib::$UnionPetscLib, ts::TS )

    @chk ccall(
               (:TSAdjointResetForward, $petsc_library),
               PetscErrorCode,
               (CTS,),
               ts,
              )


	return nothing
end 

"""
	TSAdjointSetUp(petsclib::PetscLibType,ts::TS) 
Sets up the internal data structures for the later use
of an adjoint solver

Collective

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Level: advanced

-seealso: [](ch_ts), `TSCreate()`, `TSAdjointStep()`, `TSSetCostGradients()`

# External Links
$(_doc_external("Ts/TSAdjointSetUp"))
"""
function TSAdjointSetUp(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSAdjointSetUp(petsclib::$UnionPetscLib, ts::TS )

    @chk ccall(
               (:TSAdjointSetUp, $petsc_library),
               PetscErrorCode,
               (CTS,),
               ts,
              )


	return nothing
end 

"""
	TSAdjointReset(petsclib::PetscLibType,ts::TS) 
Resets a `TS` adjoint context and removes any allocated `Vec`s and `Mat`s.

Collective

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Level: beginner

-seealso: [](ch_ts), `TSCreate()`, `TSAdjointSetUp()`, `TSDestroy()`

# External Links
$(_doc_external("Ts/TSAdjointReset"))
"""
function TSAdjointReset(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSAdjointReset(petsclib::$UnionPetscLib, ts::TS )

    @chk ccall(
               (:TSAdjointReset, $petsc_library),
               PetscErrorCode,
               (CTS,),
               ts,
              )


	return nothing
end 

"""
	TSAdjointSetSteps(petsclib::PetscLibType,ts::TS, steps::PetscInt) 
Sets the number of steps the adjoint solver should take backward in time

Logically Collective

Input Parameters:
- `ts`    - the `TS` context obtained from `TSCreate()`
- `steps` - number of steps to use

Level: intermediate

-seealso: [](ch_ts), `TSAdjointSolve()`, `TS`, `TSSetExactFinalTime()`

# External Links
$(_doc_external("Ts/TSAdjointSetSteps"))
"""
function TSAdjointSetSteps(petsclib::PetscLibType, ts::TS, steps::PetscInt) end

@for_petsc function TSAdjointSetSteps(petsclib::$UnionPetscLib, ts::TS, steps::$PetscInt )

    @chk ccall(
               (:TSAdjointSetSteps, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt),
               ts, steps,
              )


	return nothing
end 

"""
	TSAdjointMonitorSetFromOptions(petsclib::PetscLibType,ts::TS, name::String, help::String, manual::String, monitor::external, monitorsetup::external) 
Sets a monitor function and viewer appropriate for the type indicated by the user

Collective

Input Parameters:
- `ts`           - `TS` object you wish to monitor
- `name`         - the monitor type one is seeking
- `help`         - message indicating what monitoring is done
- `manual`       - manual page for the monitor
- `monitor`      - the monitor function, its context must be a `PetscViewerAndFormat`
- `monitorsetup` - a function that is called once ONLY if the user selected this monitor that may set additional features of the `TS` or `PetscViewer` objects

Level: developer

-seealso: [](ch_ts), `PetscOptionsCreateViewer()`, `PetscOptionsGetReal()`, `PetscOptionsHasName()`, `PetscOptionsGetString()`,
`PetscOptionsGetIntArray()`, `PetscOptionsGetRealArray()`, `PetscOptionsBool()`
`PetscOptionsInt()`, `PetscOptionsString()`, `PetscOptionsReal()`,
`PetscOptionsName()`, `PetscOptionsBegin()`, `PetscOptionsEnd()`, `PetscOptionsHeadBegin()`,
`PetscOptionsStringArray()`, `PetscOptionsRealArray()`, `PetscOptionsScalar()`,
`PetscOptionsBoolGroupBegin()`, `PetscOptionsBoolGroup()`, `PetscOptionsBoolGroupEnd()`,
`PetscOptionsFList()`, `PetscOptionsEList()`, `PetscViewerAndFormat`

# External Links
$(_doc_external("Ts/TSAdjointMonitorSetFromOptions"))
"""
function TSAdjointMonitorSetFromOptions(petsclib::PetscLibType, ts::TS, name::String, help::String, manual::String, monitor::external, monitorsetup::external) end

@for_petsc function TSAdjointMonitorSetFromOptions(petsclib::$UnionPetscLib, ts::TS, name::String, help::String, manual::String, monitor::external, monitorsetup::external )

    @chk ccall(
               (:TSAdjointMonitorSetFromOptions, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, external, external),
               ts, name, help, manual, monitor, monitorsetup,
              )


	return nothing
end 

"""
	TSAdjointMonitorSet(petsclib::PetscLibType,ts::TS, adjointmonitor::external, adjointmctx::Cvoid, adjointmdestroy::PetscCtxDestroyFn) 
Sets an ADDITIONAL function that is to be used at every
timestep to display the iteration's  progress.

Logically Collective

Input Parameters:
- `ts`              - the `TS` context obtained from `TSCreate()`
- `adjointmonitor`  - monitoring routine
- `adjointmctx`     - [optional] user-defined context for private data for the monitor routine
(use `NULL` if no context is desired)
- `adjointmdestroy` - [optional] routine that frees monitor context (may be `NULL`), see `PetscCtxDestroyFn` for its calling sequence

Calling sequence of `adjointmonitor`:
- `ts`          - the `TS` context
- `steps`       - iteration number (after the final time step the monitor routine is called with
a step of -1, this is at the final time which may have been interpolated to)
- `time`        - current time
- `u`           - current iterate
- `numcost`     - number of cost functionos
- `lambda`      - sensitivities to initial conditions
- `mu`          - sensitivities to parameters
- `adjointmctx` - [optional] adjoint monitoring context

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSAdjointSolve()`, `TSAdjointMonitorCancel()`, `PetscCtxDestroyFn`

# External Links
$(_doc_external("Ts/TSAdjointMonitorSet"))
"""
function TSAdjointMonitorSet(petsclib::PetscLibType, ts::TS, adjointmonitor::external, adjointmctx::Cvoid, adjointmdestroy::PetscCtxDestroyFn) end

@for_petsc function TSAdjointMonitorSet(petsclib::$UnionPetscLib, ts::TS, adjointmonitor::external, adjointmctx::Cvoid, adjointmdestroy::PetscCtxDestroyFn )

    @chk ccall(
               (:TSAdjointMonitorSet, $petsc_library),
               PetscErrorCode,
               (CTS, external, Ptr{Cvoid}, Ptr{PetscCtxDestroyFn}),
               ts, adjointmonitor, adjointmctx, adjointmdestroy,
              )


	return nothing
end 

"""
	TSAdjointMonitorCancel(petsclib::PetscLibType,ts::TS) 
Clears all the adjoint monitors that have been set on a time

Logically Collective

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

-seealso: [](ch_ts), `TS`, `TSAdjointSolve()`, `TSAdjointMonitorSet()`

# External Links
$(_doc_external("Ts/TSAdjointMonitorCancel"))
"""
function TSAdjointMonitorCancel(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSAdjointMonitorCancel(petsclib::$UnionPetscLib, ts::TS )

    @chk ccall(
               (:TSAdjointMonitorCancel, $petsc_library),
               PetscErrorCode,
               (CTS,),
               ts,
              )


	return nothing
end 

"""
	TSAdjointMonitorDefault(petsclib::PetscLibType,ts::TS, step::PetscInt, time::PetscReal, v::PetscVec, numcost::PetscInt, lambda::Vector{PetscVec}, mu::Vector{PetscVec}, vf::PetscViewerAndFormat) 
the default monitor of adjoint computations

Input Parameters:
- `ts`      - the `TS` context
- `step`    - iteration number (after the final time step the monitor routine is called with a
step of -1, this is at the final time which may have been interpolated to)
- `time`    - current time
- `v`       - current iterate
- `numcost` - number of cost functionos
- `lambda`  - sensitivities to initial conditions
- `mu`      - sensitivities to parameters
- `vf`      - the viewer and format

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSAdjointSolve()`, `TSAdjointMonitorSet()`

# External Links
$(_doc_external("Ts/TSAdjointMonitorDefault"))
"""
function TSAdjointMonitorDefault(petsclib::PetscLibType, ts::TS, step::PetscInt, time::PetscReal, v::PetscVec, numcost::PetscInt, lambda::Vector{PetscVec}, mu::Vector{PetscVec}, vf::PetscViewerAndFormat) end

@for_petsc function TSAdjointMonitorDefault(petsclib::$UnionPetscLib, ts::TS, step::$PetscInt, time::$PetscReal, v::PetscVec, numcost::$PetscInt, lambda::Vector{PetscVec}, mu::Vector{PetscVec}, vf::PetscViewerAndFormat )

    @chk ccall(
               (:TSAdjointMonitorDefault, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt, $PetscReal, CVec, $PetscInt, Ptr{CVec}, Ptr{CVec}, Ptr{PetscViewerAndFormat}),
               ts, step, time, v, numcost, lambda, mu, vf,
              )


	return nothing
end 

"""
	TSAdjointMonitorDrawSensi(petsclib::PetscLibType,ts::TS, step::PetscInt, ptime::PetscReal, u::PetscVec, numcost::PetscInt, lambda::Vector{PetscVec}, mu::Vector{PetscVec}, dummy::Cvoid) 
Monitors progress of the adjoint `TS` solvers by calling
`VecView()` for the sensitivities to initial states at each timestep

Collective

Input Parameters:
- `ts`      - the `TS` context
- `step`    - current time-step
- `ptime`   - current time
- `u`       - current state
- `numcost` - number of cost functions
- `lambda`  - sensitivities to initial conditions
- `mu`      - sensitivities to parameters
- `dummy`   - either a viewer or `NULL`

Level: intermediate

-seealso: [](ch_ts), `TSAdjointSolve()`, `TSAdjointMonitorSet()`, `TSAdjointMonitorDefault()`, `VecView()`

# External Links
$(_doc_external("Ts/TSAdjointMonitorDrawSensi"))
"""
function TSAdjointMonitorDrawSensi(petsclib::PetscLibType, ts::TS, step::PetscInt, ptime::PetscReal, u::PetscVec, numcost::PetscInt, lambda::Vector{PetscVec}, mu::Vector{PetscVec}, dummy::Cvoid) end

@for_petsc function TSAdjointMonitorDrawSensi(petsclib::$UnionPetscLib, ts::TS, step::$PetscInt, ptime::$PetscReal, u::PetscVec, numcost::$PetscInt, lambda::Vector{PetscVec}, mu::Vector{PetscVec}, dummy::Cvoid )

    @chk ccall(
               (:TSAdjointMonitorDrawSensi, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt, $PetscReal, CVec, $PetscInt, Ptr{CVec}, Ptr{CVec}, Ptr{Cvoid}),
               ts, step, ptime, u, numcost, lambda, mu, dummy,
              )


	return nothing
end 

"""
	TSAdjointSetFromOptions(petsclib::PetscLibType,ts::TS, PetscOptionsObject::PetscOptionItems) 
Sets various `TS` adjoint parameters from options database.

Collective

Input Parameters:
- `ts`                 - the `TS` context
- `PetscOptionsObject` - the options context

Options Database Keys:
- `-ts_adjoint_solve <yes,no>`     - After solving the ODE/DAE solve the adjoint problem (requires `-ts_save_trajectory`)
- `-ts_adjoint_monitor`            - print information at each adjoint time step
- `-ts_adjoint_monitor_draw_sensi` - monitor the sensitivity of the first cost function wrt initial conditions (lambda[0]) graphically

Level: developer

-seealso: [](ch_ts), `TSSetSaveTrajectory()`, `TSTrajectorySetUp()`

# External Links
$(_doc_external("Ts/TSAdjointSetFromOptions"))
"""
function TSAdjointSetFromOptions(petsclib::PetscLibType, ts::TS, PetscOptionsObject::PetscOptionItems) end

@for_petsc function TSAdjointSetFromOptions(petsclib::$UnionPetscLib, ts::TS, PetscOptionsObject::PetscOptionItems )

    @chk ccall(
               (:TSAdjointSetFromOptions, $petsc_library),
               PetscErrorCode,
               (CTS, PetscOptionItems),
               ts, PetscOptionsObject,
              )


	return nothing
end 

"""
	TSAdjointStep(petsclib::PetscLibType,ts::TS) 
Steps one time step backward in the adjoint run

Collective

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Level: intermediate

-seealso: [](ch_ts), `TSAdjointSetUp()`, `TSAdjointSolve()`

# External Links
$(_doc_external("Ts/TSAdjointStep"))
"""
function TSAdjointStep(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSAdjointStep(petsclib::$UnionPetscLib, ts::TS )

    @chk ccall(
               (:TSAdjointStep, $petsc_library),
               PetscErrorCode,
               (CTS,),
               ts,
              )


	return nothing
end 

"""
	TSAdjointSolve(petsclib::PetscLibType,ts::TS) 
Solves the discrete ajoint problem for an ODE/DAE

Collective
`

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Options Database Key:
- `-ts_adjoint_view_solution <viewerinfo>` - views the first gradient with respect to the initial values

Level: intermediate

-seealso: [](ch_ts), `TSCreate()`, `TSSetCostGradients()`, `TSSetSolution()`, `TSAdjointStep()`

# External Links
$(_doc_external("Ts/TSAdjointSolve"))
"""
function TSAdjointSolve(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSAdjointSolve(petsclib::$UnionPetscLib, ts::TS )

    @chk ccall(
               (:TSAdjointSolve, $petsc_library),
               PetscErrorCode,
               (CTS,),
               ts,
              )


	return nothing
end 

"""
	TSAdjointMonitor(petsclib::PetscLibType,ts::TS, step::PetscInt, ptime::PetscReal, u::PetscVec, numcost::PetscInt, lambda::Vector{PetscVec}, mu::Vector{PetscVec}) 
Runs all user

Collective

Input Parameters:
- `ts`      - time stepping context obtained from `TSCreate()`
- `step`    - step number that has just completed
- `ptime`   - model time of the state
- `u`       - state at the current model time
- `numcost` - number of cost functions (dimension of lambda  or mu)
- `lambda`  - vectors containing the gradients of the cost functions with respect to the ODE/DAE solution variables
- `mu`      - vectors containing the gradients of the cost functions with respect to the problem parameters

Level: developer

-seealso: `TSAdjointMonitorSet()`, `TSAdjointSolve()`

# External Links
$(_doc_external("Ts/TSAdjointMonitor"))
"""
function TSAdjointMonitor(petsclib::PetscLibType, ts::TS, step::PetscInt, ptime::PetscReal, u::PetscVec, numcost::PetscInt, lambda::Vector{PetscVec}, mu::Vector{PetscVec}) end

@for_petsc function TSAdjointMonitor(petsclib::$UnionPetscLib, ts::TS, step::$PetscInt, ptime::$PetscReal, u::PetscVec, numcost::$PetscInt, lambda::Vector{PetscVec}, mu::Vector{PetscVec} )

    @chk ccall(
               (:TSAdjointMonitor, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt, $PetscReal, CVec, $PetscInt, Ptr{CVec}, Ptr{CVec}),
               ts, step, ptime, u, numcost, lambda, mu,
              )


	return nothing
end 

"""
	TSAdjointCostIntegral(petsclib::PetscLibType,ts::TS) 
Evaluate the cost integral in the adjoint run.

Collective

Input Parameter:
- `ts` - time stepping context

Level: advanced

-seealso: [](ch_ts), `TSAdjointSolve()`, `TSAdjointStep()`

# External Links
$(_doc_external("Ts/TSAdjointCostIntegral"))
"""
function TSAdjointCostIntegral(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSAdjointCostIntegral(petsclib::$UnionPetscLib, ts::TS )

    @chk ccall(
               (:TSAdjointCostIntegral, $petsc_library),
               PetscErrorCode,
               (CTS,),
               ts,
              )


	return nothing
end 

"""
	TSForwardSetUp(petsclib::PetscLibType,ts::TS) 
Sets up the internal data structures for the later use
of forward sensitivity analysis

Collective

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Level: advanced

-seealso: [](ch_ts), `TS`, `TSCreate()`, `TSDestroy()`, `TSSetUp()`

# External Links
$(_doc_external("Ts/TSForwardSetUp"))
"""
function TSForwardSetUp(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSForwardSetUp(petsclib::$UnionPetscLib, ts::TS )

    @chk ccall(
               (:TSForwardSetUp, $petsc_library),
               PetscErrorCode,
               (CTS,),
               ts,
              )


	return nothing
end 

"""
	TSForwardReset(petsclib::PetscLibType,ts::TS) 
Reset the internal data structures used by forward sensitivity analysis

Collective

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Level: advanced

-seealso: [](ch_ts), `TSCreate()`, `TSDestroy()`, `TSForwardSetUp()`

# External Links
$(_doc_external("Ts/TSForwardReset"))
"""
function TSForwardReset(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSForwardReset(petsclib::$UnionPetscLib, ts::TS )

    @chk ccall(
               (:TSForwardReset, $petsc_library),
               PetscErrorCode,
               (CTS,),
               ts,
              )


	return nothing
end 

"""
	TSForwardStep(petsclib::PetscLibType,ts::TS) 
Compute the forward sensitivity for one time step.

Collective

Input Parameter:
- `ts` - time stepping context

Level: advanced

-seealso: [](ch_ts), `TSForwardSetSensitivities()`, `TSForwardGetSensitivities()`, `TSForwardSetIntegralGradients()`, `TSForwardGetIntegralGradients()`, `TSForwardSetUp()`

# External Links
$(_doc_external("Ts/TSForwardStep"))
"""
function TSForwardStep(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSForwardStep(petsclib::$UnionPetscLib, ts::TS )

    @chk ccall(
               (:TSForwardStep, $petsc_library),
               PetscErrorCode,
               (CTS,),
               ts,
              )


	return nothing
end 

"""
	TSForwardSetSensitivities(petsclib::PetscLibType,ts::TS, nump::PetscInt, Smat::PetscMat) 
Sets the initial value of the trajectory sensitivities of solution  w.r.t. the problem parameters and initial values.

Logically Collective

Input Parameters:
- `ts`   - the `TS` context obtained from `TSCreate()`
- `nump` - number of parameters
- `Smat` - sensitivities with respect to the parameters, the number of entries in these vectors is the same as the number of parameters

Level: beginner

-seealso: [](ch_ts), `TSForwardGetSensitivities()`, `TSForwardSetIntegralGradients()`, `TSForwardGetIntegralGradients()`, `TSForwardStep()`

# External Links
$(_doc_external("Ts/TSForwardSetSensitivities"))
"""
function TSForwardSetSensitivities(petsclib::PetscLibType, ts::TS, nump::PetscInt, Smat::PetscMat) end

@for_petsc function TSForwardSetSensitivities(petsclib::$UnionPetscLib, ts::TS, nump::$PetscInt, Smat::PetscMat )

    @chk ccall(
               (:TSForwardSetSensitivities, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt, CMat),
               ts, nump, Smat,
              )


	return nothing
end 

"""
	nump::PetscInt = TSForwardGetSensitivities(petsclib::PetscLibType,ts::TS, Smat::PetscMat) 
Returns the trajectory sensitivities

Not Collective, but Smat returned is parallel if ts is parallel

Output Parameters:
- `ts`   - the `TS` context obtained from `TSCreate()`
- `nump` - number of parameters
- `Smat` - sensitivities with respect to the parameters, the number of entries in these vectors is the same as the number of parameters

Level: intermediate

-seealso: [](ch_ts), `TSForwardSetSensitivities()`, `TSForwardSetIntegralGradients()`, `TSForwardGetIntegralGradients()`, `TSForwardStep()`

# External Links
$(_doc_external("Ts/TSForwardGetSensitivities"))
"""
function TSForwardGetSensitivities(petsclib::PetscLibType, ts::TS, Smat::PetscMat) end

@for_petsc function TSForwardGetSensitivities(petsclib::$UnionPetscLib, ts::TS, Smat::PetscMat )
	nump_ = Ref{$PetscInt}()
	Smat_ = Ref(Smat.ptr)

    @chk ccall(
               (:TSForwardGetSensitivities, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{$PetscInt}, Ptr{CMat}),
               ts, nump_, Smat_,
              )

	nump = nump_[]
	Smat.ptr = C_NULL

	return nump
end 

"""
	TSForwardCostIntegral(petsclib::PetscLibType,ts::TS) 
Evaluate the cost integral in the forward run.

Collective

Input Parameter:
- `ts` - time stepping context

Level: advanced

-seealso: [](ch_ts), `TS`, `TSSolve()`, `TSAdjointCostIntegral()`

# External Links
$(_doc_external("Ts/TSForwardCostIntegral"))
"""
function TSForwardCostIntegral(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSForwardCostIntegral(petsclib::$UnionPetscLib, ts::TS )

    @chk ccall(
               (:TSForwardCostIntegral, $petsc_library),
               PetscErrorCode,
               (CTS,),
               ts,
              )


	return nothing
end 

"""
	TSForwardSetInitialSensitivities(petsclib::PetscLibType,ts::TS, didp::PetscMat) 
Set initial values for tangent linear sensitivities

Collective

Input Parameters:
- `ts`   - the `TS` context obtained from `TSCreate()`
- `didp` - parametric sensitivities of the initial condition

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSForwardSetSensitivities()`

# External Links
$(_doc_external("Ts/TSForwardSetInitialSensitivities"))
"""
function TSForwardSetInitialSensitivities(petsclib::PetscLibType, ts::TS, didp::PetscMat) end

@for_petsc function TSForwardSetInitialSensitivities(petsclib::$UnionPetscLib, ts::TS, didp::PetscMat )

    @chk ccall(
               (:TSForwardSetInitialSensitivities, $petsc_library),
               PetscErrorCode,
               (CTS, CMat),
               ts, didp,
              )


	return nothing
end 

"""
	ns::PetscInt = TSForwardGetStages(petsclib::PetscLibType,ts::TS, S::PetscMat) 
Get the number of stages and the tangent linear sensitivities at the intermediate stages

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Output Parameters:
- `ns` - number of stages
- `S`  - tangent linear sensitivities at the intermediate stages

Level: advanced

-seealso: `TS`

# External Links
$(_doc_external("Ts/TSForwardGetStages"))
"""
function TSForwardGetStages(petsclib::PetscLibType, ts::TS, S::PetscMat) end

@for_petsc function TSForwardGetStages(petsclib::$UnionPetscLib, ts::TS, S::PetscMat )
	ns_ = Ref{$PetscInt}()

    @chk ccall(
               (:TSForwardGetStages, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{$PetscInt}, CMat),
               ts, ns_, S,
              )

	ns = ns_[]

	return ns
end 

"""
	quadts::TS = TSCreateQuadratureTS(petsclib::PetscLibType,ts::TS, fwd::PetscBool) 
Create a sub

Input Parameters:
- `ts`  - the `TS` context obtained from `TSCreate()`
- `fwd` - flag indicating whether to evaluate cost integral in the forward run or the adjoint run

Output Parameter:
- `quadts` - the child `TS` context

Level: intermediate

-seealso: [](ch_ts), `TSGetQuadratureTS()`

# External Links
$(_doc_external("Ts/TSCreateQuadratureTS"))
"""
function TSCreateQuadratureTS(petsclib::PetscLibType, ts::TS, fwd::PetscBool) end

@for_petsc function TSCreateQuadratureTS(petsclib::$UnionPetscLib, ts::TS, fwd::PetscBool )
	quadts_ = Ref{CTS}()

    @chk ccall(
               (:TSCreateQuadratureTS, $petsc_library),
               PetscErrorCode,
               (CTS, PetscBool, Ptr{CTS}),
               ts, fwd, quadts_,
              )

	quadts = TS(quadts_[], petsclib)

	return quadts
end 

"""
	fwd::PetscBool = TSGetQuadratureTS(petsclib::PetscLibType,ts::TS, quadts::TS) 
Return the sub

Input Parameter:
- `ts` - the `TS` context obtained from `TSCreate()`

Output Parameters:
- `fwd`    - flag indicating whether to evaluate cost integral in the forward run or the adjoint run
- `quadts` - the child `TS` context

Level: intermediate

-seealso: [](ch_ts), `TSCreateQuadratureTS()`

# External Links
$(_doc_external("Ts/TSGetQuadratureTS"))
"""
function TSGetQuadratureTS(petsclib::PetscLibType, ts::TS, quadts::TS) end

@for_petsc function TSGetQuadratureTS(petsclib::$UnionPetscLib, ts::TS, quadts::TS )
	fwd_ = Ref{PetscBool}()
	quadts_ = Ref(quadts.ptr)

    @chk ccall(
               (:TSGetQuadratureTS, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{PetscBool}, Ptr{CTS}),
               ts, fwd_, quadts_,
              )

	fwd = fwd_[]
	quadts.ptr = C_NULL

	return fwd
end 

"""
	TSComputeSNESJacobian(petsclib::PetscLibType,ts::TS, x::PetscVec, J::PetscMat, Jpre::PetscMat) 
Compute the Jacobian needed for the `SNESSolve()` in `TS`

Collective

Input Parameters:
- `ts` - the `TS` context obtained from `TSCreate()`
- `x`  - state vector

Output Parameters:
- `J`    - Jacobian matrix
- `Jpre` - matrix used to compute the preconditioner for `J` (may be same as `J`)

Level: developer

-seealso: `SNES`, `TS`, `SNESSetJacobian()`, `TSSetRHSJacobian()`, `TSSetIJacobian()`

# External Links
$(_doc_external("Ts/TSComputeSNESJacobian"))
"""
function TSComputeSNESJacobian(petsclib::PetscLibType, ts::TS, x::PetscVec, J::PetscMat, Jpre::PetscMat) end

@for_petsc function TSComputeSNESJacobian(petsclib::$UnionPetscLib, ts::TS, x::PetscVec, J::PetscMat, Jpre::PetscMat )

    @chk ccall(
               (:TSComputeSNESJacobian, $petsc_library),
               PetscErrorCode,
               (CTS, CVec, CMat, CMat),
               ts, x, J, Jpre,
              )


	return nothing
end 

"""
	TSARKIMEXRegisterDestroy(petsclib::PetscLibType) 
Frees the list of schemes that were registered by `TSARKIMEXRegister()`.

Not Collective

Level: advanced

-seealso: [](ch_ts), `TSARKIMEX`, `TSARKIMEXRegister()`, `TSARKIMEXRegisterAll()`

# External Links
$(_doc_external("Ts/TSARKIMEXRegisterDestroy"))
"""
function TSARKIMEXRegisterDestroy(petsclib::PetscLibType) end

@for_petsc function TSARKIMEXRegisterDestroy(petsclib::$UnionPetscLib)

    @chk ccall(
               (:TSARKIMEXRegisterDestroy, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	TSARKIMEXInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the `TSARKIMEX` package. It is called
from `TSInitializePackage()`.

Level: developer

-seealso: [](ch_ts), `PetscInitialize()`, `TSARKIMEXFinalizePackage()`

# External Links
$(_doc_external("Ts/TSARKIMEXInitializePackage"))
"""
function TSARKIMEXInitializePackage(petsclib::PetscLibType) end

@for_petsc function TSARKIMEXInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:TSARKIMEXInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	TSARKIMEXFinalizePackage(petsclib::PetscLibType) 
This function destroys everything in the `TSARKIMEX` package. It is
called from `PetscFinalize()`.

Level: developer

-seealso: [](ch_ts), `PetscFinalize()`, `TSARKIMEXInitializePackage()`

# External Links
$(_doc_external("Ts/TSARKIMEXFinalizePackage"))
"""
function TSARKIMEXFinalizePackage(petsclib::PetscLibType) end

@for_petsc function TSARKIMEXFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:TSARKIMEXFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	TSARKIMEXRegister(petsclib::PetscLibType,name::TSARKIMEXType, order::PetscInt, s::PetscInt, At::Vector{PetscReal}, bt::Vector{PetscReal}, ct::Vector{PetscReal}, A::Vector{PetscReal}, b::Vector{PetscReal}, c::Vector{PetscReal}, bembedt::Vector{PetscReal}, bembed::Vector{PetscReal}, pinterp::PetscInt, binterpt::Vector{PetscReal}, binterp::Vector{PetscReal}) 
register a `TSARKIMEX` scheme by providing the entries in the Butcher tableau and optionally embedded approximations and interpolation

Logically Collective

Input Parameters:
- `name`     - identifier for method
- `order`    - approximation order of method
- `s`        - number of stages, this is the dimension of the matrices below
- `At`       - Butcher table of stage coefficients for stiff part (dimension s*s, row-major)
- `bt`       - Butcher table for completing the stiff part of the step (dimension s; NULL to use the last row of At)
- `ct`       - Abscissa of each stiff stage (dimension s, NULL to use row sums of At)
- `A`        - Non-stiff stage coefficients (dimension s*s, row-major)
- `b`        - Non-stiff step completion table (dimension s; NULL to use last row of At)
- `c`        - Non-stiff abscissa (dimension s; NULL to use row sums of A)
- `bembedt`  - Stiff part of completion table for embedded method (dimension s; NULL if not available)
- `bembed`   - Non-stiff part of completion table for embedded method (dimension s; NULL to use bembedt if provided)
- `pinterp`  - Order of the interpolation scheme, equal to the number of columns of binterpt and binterp
- `binterpt` - Coefficients of the interpolation formula for the stiff part (dimension s*pinterp)
- `binterp`  - Coefficients of the interpolation formula for the non-stiff part (dimension s*pinterp; NULL to reuse binterpt)

Level: advanced

-seealso: [](ch_ts), `TSARKIMEX`, `TSType`, `TS`

# External Links
$(_doc_external("Ts/TSARKIMEXRegister"))
"""
function TSARKIMEXRegister(petsclib::PetscLibType, name::TSARKIMEXType, order::PetscInt, s::PetscInt, At::Vector{PetscReal}, bt::Vector{PetscReal}, ct::Vector{PetscReal}, A::Vector{PetscReal}, b::Vector{PetscReal}, c::Vector{PetscReal}, bembedt::Vector{PetscReal}, bembed::Vector{PetscReal}, pinterp::PetscInt, binterpt::Vector{PetscReal}, binterp::Vector{PetscReal}) end

@for_petsc function TSARKIMEXRegister(petsclib::$UnionPetscLib, name::TSARKIMEXType, order::$PetscInt, s::$PetscInt, At::Vector{$PetscReal}, bt::Vector{$PetscReal}, ct::Vector{$PetscReal}, A::Vector{$PetscReal}, b::Vector{$PetscReal}, c::Vector{$PetscReal}, bembedt::Vector{$PetscReal}, bembed::Vector{$PetscReal}, pinterp::$PetscInt, binterpt::Vector{$PetscReal}, binterp::Vector{$PetscReal} )

    @chk ccall(
               (:TSARKIMEXRegister, $petsc_library),
               PetscErrorCode,
               (TSARKIMEXType, $PetscInt, $PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, $PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}),
               name, order, s, At, bt, ct, A, b, c, bembedt, bembed, pinterp, binterpt, binterp,
              )


	return nothing
end 

"""
	TSDIRKRegister(petsclib::PetscLibType,name::TSDIRKType, order::PetscInt, s::PetscInt, At::Vector{PetscReal}, bt::Vector{PetscReal}, ct::Vector{PetscReal}, bembedt::Vector{PetscReal}, pinterp::PetscInt, binterpt::Vector{PetscReal}) 
register a `TSDIRK` scheme by providing the entries in its Butcher tableau and, optionally, embedded approximations and interpolation

Logically Collective.

Input Parameters:
- `name`     - identifier for method
- `order`    - approximation order of method
- `s`        - number of stages, this is the dimension of the matrices below
- `At`       - Butcher table of stage coefficients (dimension `s`*`s`, row-major order)
- `bt`       - Butcher table for completing the step (dimension `s`; pass `NULL` to use the last row of `At`)
- `ct`       - Abscissa of each stage (dimension s, NULL to use row sums of At)
- `bembedt`  - Stiff part of completion table for embedded method (dimension s; `NULL` if not available)
- `pinterp`  - Order of the interpolation scheme, equal to the number of columns of `binterpt` and `binterp`
- `binterpt` - Coefficients of the interpolation formula (dimension s*pinterp)

Level: advanced

-seealso: [](ch_ts), `TSDIRK`, `TSType`, `TS`

# External Links
$(_doc_external("Ts/TSDIRKRegister"))
"""
function TSDIRKRegister(petsclib::PetscLibType, name::TSDIRKType, order::PetscInt, s::PetscInt, At::Vector{PetscReal}, bt::Vector{PetscReal}, ct::Vector{PetscReal}, bembedt::Vector{PetscReal}, pinterp::PetscInt, binterpt::Vector{PetscReal}) end

@for_petsc function TSDIRKRegister(petsclib::$UnionPetscLib, name::TSDIRKType, order::$PetscInt, s::$PetscInt, At::Vector{$PetscReal}, bt::Vector{$PetscReal}, ct::Vector{$PetscReal}, bembedt::Vector{$PetscReal}, pinterp::$PetscInt, binterpt::Vector{$PetscReal} )

    @chk ccall(
               (:TSDIRKRegister, $petsc_library),
               PetscErrorCode,
               (TSDIRKType, $PetscInt, $PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, $PetscInt, Ptr{$PetscReal}),
               name, order, s, At, bt, ct, bembedt, pinterp, binterpt,
              )


	return nothing
end 

"""
	TSARKIMEXSetType(petsclib::PetscLibType,ts::TS, arktype::TSARKIMEXType) 
Set the type of `TSARKIMEX` scheme

Logically Collective

Input Parameters:
- `ts`      - timestepping context
- `arktype` - type of `TSARKIMEX` scheme

Options Database Key:
- `-ts_arkimex_type <1bee,a2,l2,ars122,2c,2d,2e,prssp2,3,bpr3,ars443,4,5>` - set `TSARKIMEX` scheme type

Level: intermediate

-seealso: [](ch_ts), `TSARKIMEXGetType()`, `TSARKIMEX`, `TSARKIMEXType`, `TSARKIMEX1BEE`, `TSARKIMEXA2`, `TSARKIMEXL2`, `TSARKIMEXARS122`, `TSARKIMEX2C`, `TSARKIMEX2D`, `TSARKIMEX2E`, `TSARKIMEXPRSSP2`,
`TSARKIMEX3`, `TSARKIMEXBPR3`, `TSARKIMEXARS443`, `TSARKIMEX4`, `TSARKIMEX5`

# External Links
$(_doc_external("Ts/TSARKIMEXSetType"))
"""
function TSARKIMEXSetType(petsclib::PetscLibType, ts::TS, arktype::TSARKIMEXType) end

@for_petsc function TSARKIMEXSetType(petsclib::$UnionPetscLib, ts::TS, arktype::TSARKIMEXType )

    @chk ccall(
               (:TSARKIMEXSetType, $petsc_library),
               PetscErrorCode,
               (CTS, TSARKIMEXType),
               ts, arktype,
              )


	return nothing
end 

"""
	arktype::TSARKIMEXType = TSARKIMEXGetType(petsclib::PetscLibType,ts::TS) 
Get the type of `TSARKIMEX` scheme

Logically Collective

Input Parameter:
- `ts` - timestepping context

Output Parameter:
- `arktype` - type of `TSARKIMEX` scheme

Level: intermediate

-seealso: [](ch_ts), `TSARKIMEX`

# External Links
$(_doc_external("Ts/TSARKIMEXGetType"))
"""
function TSARKIMEXGetType(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSARKIMEXGetType(petsclib::$UnionPetscLib, ts::TS )
	arktype_ = Ref{TSARKIMEXType}()

    @chk ccall(
               (:TSARKIMEXGetType, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{TSARKIMEXType}),
               ts, arktype_,
              )

	arktype = unsafe_string(arktype_[])

	return arktype
end 

"""
	TSARKIMEXSetFullyImplicit(petsclib::PetscLibType,ts::TS, flg::PetscBool) 
Solve both parts of the equation implicitly, including the part that is normally solved explicitly

Logically Collective

Input Parameters:
- `ts`  - timestepping context
- `flg` - `PETSC_TRUE` for fully implicit

Level: intermediate

-seealso: [](ch_ts), `TSARKIMEX`, `TSARKIMEXGetType()`, `TSARKIMEXGetFullyImplicit()`

# External Links
$(_doc_external("Ts/TSARKIMEXSetFullyImplicit"))
"""
function TSARKIMEXSetFullyImplicit(petsclib::PetscLibType, ts::TS, flg::PetscBool) end

@for_petsc function TSARKIMEXSetFullyImplicit(petsclib::$UnionPetscLib, ts::TS, flg::PetscBool )

    @chk ccall(
               (:TSARKIMEXSetFullyImplicit, $petsc_library),
               PetscErrorCode,
               (CTS, PetscBool),
               ts, flg,
              )


	return nothing
end 

"""
	flg::PetscBool = TSARKIMEXGetFullyImplicit(petsclib::PetscLibType,ts::TS) 
Inquires if both parts of the equation are solved implicitly

Logically Collective

Input Parameter:
- `ts` - timestepping context

Output Parameter:
- `flg` - `PETSC_TRUE` for fully implicit

Level: intermediate

-seealso: [](ch_ts), `TSARKIMEXGetType()`, `TSARKIMEXSetFullyImplicit()`

# External Links
$(_doc_external("Ts/TSARKIMEXGetFullyImplicit"))
"""
function TSARKIMEXGetFullyImplicit(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSARKIMEXGetFullyImplicit(petsclib::$UnionPetscLib, ts::TS )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:TSARKIMEXGetFullyImplicit, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{PetscBool}),
               ts, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	TSDIRKSetType(petsclib::PetscLibType,ts::TS, dirktype::TSDIRKType) 
Set the type of `TSDIRK` scheme

Logically Collective

Input Parameters:
- `ts`       - timestepping context
- `dirktype` - type of `TSDIRK` scheme

Options Database Key:
- `-ts_dirkimex_type` - set `TSDIRK` scheme type

Level: intermediate

-seealso: [](ch_ts), `TSDIRKGetType()`, `TSDIRK`, `TSDIRKType`

# External Links
$(_doc_external("Ts/TSDIRKSetType"))
"""
function TSDIRKSetType(petsclib::PetscLibType, ts::TS, dirktype::TSDIRKType) end

@for_petsc function TSDIRKSetType(petsclib::$UnionPetscLib, ts::TS, dirktype::TSDIRKType )

    @chk ccall(
               (:TSDIRKSetType, $petsc_library),
               PetscErrorCode,
               (CTS, TSDIRKType),
               ts, dirktype,
              )


	return nothing
end 

"""
	dirktype::TSDIRKType = TSDIRKGetType(petsclib::PetscLibType,ts::TS) 
Get the type of `TSDIRK` scheme

Logically Collective

Input Parameter:
- `ts` - timestepping context

Output Parameter:
- `dirktype` - type of `TSDIRK` scheme

Level: intermediate

-seealso: [](ch_ts), `TSDIRKSetType()`

# External Links
$(_doc_external("Ts/TSDIRKGetType"))
"""
function TSDIRKGetType(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSDIRKGetType(petsclib::$UnionPetscLib, ts::TS )
	dirktype_ = Ref{TSDIRKType}()

    @chk ccall(
               (:TSDIRKGetType, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{TSDIRKType}),
               ts, dirktype_,
              )

	dirktype = unsafe_string(dirktype_[])

	return dirktype
end 

"""
	TSARKIMEXSetFastSlowSplit(petsclib::PetscLibType,ts::TS, fastslow::PetscBool) 
Use `TSARKIMEX` for solving a fast

Logically Collective

Input Parameters:
- `ts`       - timestepping context
- `fastslow` - `PETSC_TRUE` enables the `TSARKIMEX` solver for a fast-slow system where the RHS is split component-wise.

Options Database Key:
- `-ts_arkimex_fastslowsplit` - <true,false>

Level: intermediate

-seealso: [](ch_ts), `TSARKIMEX`, `TSARKIMEXGetFastSlowSplit()`

# External Links
$(_doc_external("Ts/TSARKIMEXSetFastSlowSplit"))
"""
function TSARKIMEXSetFastSlowSplit(petsclib::PetscLibType, ts::TS, fastslow::PetscBool) end

@for_petsc function TSARKIMEXSetFastSlowSplit(petsclib::$UnionPetscLib, ts::TS, fastslow::PetscBool )

    @chk ccall(
               (:TSARKIMEXSetFastSlowSplit, $petsc_library),
               PetscErrorCode,
               (CTS, PetscBool),
               ts, fastslow,
              )


	return nothing
end 

"""
	fastslow::PetscBool = TSARKIMEXGetFastSlowSplit(petsclib::PetscLibType,ts::TS) 
Gets whether to use `TSARKIMEX` for a fast

Not Collective

Input Parameter:
- `ts` - timestepping context

Output Parameter:
- `fastslow` - `PETSC_TRUE` if `TSARKIMEX` will be used for solving a fast-slow system, `PETSC_FALSE` otherwise

Level: intermediate

-seealso: [](ch_ts), `TSARKIMEX`, `TSARKIMEXSetFastSlowSplit()`

# External Links
$(_doc_external("Ts/TSARKIMEXGetFastSlowSplit"))
"""
function TSARKIMEXGetFastSlowSplit(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSARKIMEXGetFastSlowSplit(petsclib::$UnionPetscLib, ts::TS )
	fastslow_ = Ref{PetscBool}()

    @chk ccall(
               (:TSARKIMEXGetFastSlowSplit, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{PetscBool}),
               ts, fastslow_,
              )

	fastslow = fastslow_[]

	return fastslow
end 

"""
	TSPythonSetType(petsclib::PetscLibType,ts::TS, pyname::String) 
Initialize a `TS` object implemented in Python.

Collective

Input Parameters:
- `ts`  - the `TS` context
- `pyname`  - full dotted Python name [package].module[.{class|function}]

Options Database Key:
- `-ts_python_type <pyname>`  - python class

Level: intermediate

-seealso: [](ch_ts), `TSCreate()`, `TSSetType()`, `TSPYTHON`, `PetscPythonInitialize()`

# External Links
$(_doc_external("Ts/TSPythonSetType"))
"""
function TSPythonSetType(petsclib::PetscLibType, ts::TS, pyname::String) end

@for_petsc function TSPythonSetType(petsclib::$UnionPetscLib, ts::TS, pyname::String )

    @chk ccall(
               (:TSPythonSetType, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{Cchar}),
               ts, pyname,
              )


	return nothing
end 

"""
	pyname::String = TSPythonGetType(petsclib::PetscLibType,ts::TS) 
Get the type of a `TS` object implemented in Python.

Not Collective

Input Parameter:
- `ts`  - the `TS` context

Output Parameter:
- `pyname`  - full dotted Python name [package].module[.{class|function}]

Level: intermediate

-seealso: [](ch_ts), `TSCreate()`, `TSSetType()`, `TSPYTHON`, `PetscPythonInitialize()`, `TSPythonSetType()`

# External Links
$(_doc_external("Ts/TSPythonGetType"))
"""
function TSPythonGetType(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSPythonGetType(petsclib::$UnionPetscLib, ts::TS )
	pyname_ = Ref{Ptr{Cchar}}()

    @chk ccall(
               (:TSPythonGetType, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{Ptr{Cchar}}),
               ts, pyname_,
              )

	pyname = unsafe_wrap(Array, pyname_[], VecGetLocalSize(petsclib, x); own = false)

	return pyname
end 

"""
	TSSSPSetType(petsclib::PetscLibType,ts::TS, ssptype::TSSSPType) 
set the `TSSSP` time integration scheme to use

Logically Collective

Input Parameters:
- `ts`      - time stepping object
- `ssptype` - type of scheme to use

Options Database Keys:
- `-ts_ssp_type <rks2>`               - Type of `TSSSP` method (one of) rks2 rks3 rk104
- `-ts_ssp_nstages<rks2: 5, rks3: 9>` - Number of stages

Level: beginner

-seealso: [](ch_ts), `TSSSP`, `TSSSPGetType()`, `TSSSPSetNumStages()`, `TSSSPRKS2`, `TSSSPRKS3`, `TSSSPRK104`

# External Links
$(_doc_external("Ts/TSSSPSetType"))
"""
function TSSSPSetType(petsclib::PetscLibType, ts::TS, ssptype::TSSSPType) end

@for_petsc function TSSSPSetType(petsclib::$UnionPetscLib, ts::TS, ssptype::TSSSPType )

    @chk ccall(
               (:TSSSPSetType, $petsc_library),
               PetscErrorCode,
               (CTS, TSSSPType),
               ts, ssptype,
              )


	return nothing
end 

"""
	type::TSSSPType = TSSSPGetType(petsclib::PetscLibType,ts::TS) 
get the `TSSSP` time integration scheme

Logically Collective

Input Parameter:
- `ts` - time stepping object

Output Parameter:
- `type` - type of scheme being used

Level: beginner

-seealso: [](ch_ts), `TSSSP`, `TSSSPSetType()`, `TSSSPSetNumStages()`, `TSSSPRKS2`, `TSSSPRKS3`, `TSSSPRK104`

# External Links
$(_doc_external("Ts/TSSSPGetType"))
"""
function TSSSPGetType(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSSSPGetType(petsclib::$UnionPetscLib, ts::TS )
	type_ = Ref{TSSSPType}()

    @chk ccall(
               (:TSSSPGetType, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{TSSSPType}),
               ts, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	TSSSPSetNumStages(petsclib::PetscLibType,ts::TS, nstages::PetscInt) 
set the number of stages to use with the `TSSSP` method. Must be called after
`TSSSPSetType()`.

Logically Collective

Input Parameters:
- `ts`      - time stepping object
- `nstages` - number of stages

Options Database Keys:
- `-ts_ssp_type <rks2>`               - Type of `TSSSP` method (one of) rks2 rks3 rk104
- `-ts_ssp_nstages<rks2: 5, rks3: 9>` - Number of stages

Level: beginner

-seealso: [](ch_ts), `TSSSP`, `TSSSPGetNumStages()`, `TSSSPRKS2`, `TSSSPRKS3`, `TSSSPRK104`

# External Links
$(_doc_external("Ts/TSSSPSetNumStages"))
"""
function TSSSPSetNumStages(petsclib::PetscLibType, ts::TS, nstages::PetscInt) end

@for_petsc function TSSSPSetNumStages(petsclib::$UnionPetscLib, ts::TS, nstages::$PetscInt )

    @chk ccall(
               (:TSSSPSetNumStages, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt),
               ts, nstages,
              )


	return nothing
end 

"""
	nstages::PetscInt = TSSSPGetNumStages(petsclib::PetscLibType,ts::TS) 
get the number of stages in the `TSSSP` time integration scheme

Logically Collective

Input Parameter:
- `ts` - time stepping object

Output Parameter:
- `nstages` - number of stages

Level: beginner

-seealso: [](ch_ts), `TSSSP`, `TSSSPGetType()`, `TSSSPSetNumStages()`, `TSSSPRKS2`, `TSSSPRKS3`, `TSSSPRK104`

# External Links
$(_doc_external("Ts/TSSSPGetNumStages"))
"""
function TSSSPGetNumStages(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSSSPGetNumStages(petsclib::$UnionPetscLib, ts::TS )
	nstages_ = Ref{$PetscInt}()

    @chk ccall(
               (:TSSSPGetNumStages, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{$PetscInt}),
               ts, nstages_,
              )

	nstages = nstages_[]

	return nstages
end 

"""
	TSSSPInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the `TSSSP` package. It is called
from `TSInitializePackage()`.

Level: developer

-seealso: [](ch_ts), `PetscInitialize()`, `TSSSPFinalizePackage()`, `TSInitializePackage()`

# External Links
$(_doc_external("Ts/TSSSPInitializePackage"))
"""
function TSSSPInitializePackage(petsclib::PetscLibType) end

@for_petsc function TSSSPInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:TSSSPInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	TSSSPFinalizePackage(petsclib::PetscLibType) 
This function destroys everything in the `TSSSP` package. It is
called from `PetscFinalize()`.

Level: developer

-seealso: [](ch_ts), `PetscFinalize()`, `TSSSPInitiallizePackage()`, `TSInitializePackage()`

# External Links
$(_doc_external("Ts/TSSSPFinalizePackage"))
"""
function TSSSPFinalizePackage(petsclib::PetscLibType) end

@for_petsc function TSSSPFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:TSSSPFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	TSRKRegisterDestroy(petsclib::PetscLibType) 
Frees the list of schemes that were registered by `TSRKRegister()`.

Not Collective

Level: advanced

-seealso: [](ch_ts), `TSRK`, `TSRKRegister()`, `TSRKRegisterAll()`

# External Links
$(_doc_external("Ts/TSRKRegisterDestroy"))
"""
function TSRKRegisterDestroy(petsclib::PetscLibType) end

@for_petsc function TSRKRegisterDestroy(petsclib::$UnionPetscLib)

    @chk ccall(
               (:TSRKRegisterDestroy, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	TSRKInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the `TSRK` package. It is called
from `TSInitializePackage()`.

Level: developer

-seealso: [](ch_ts), `TSInitializePackage()`, `PetscInitialize()`, `TSRKFinalizePackage()`

# External Links
$(_doc_external("Ts/TSRKInitializePackage"))
"""
function TSRKInitializePackage(petsclib::PetscLibType) end

@for_petsc function TSRKInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:TSRKInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	TSRKFinalizePackage(petsclib::PetscLibType) 
This function destroys everything in the `TSRK` package. It is
called from `PetscFinalize()`.

Level: developer

-seealso: [](ch_ts), `PetscFinalize()`, `TSRKInitializePackage()`

# External Links
$(_doc_external("Ts/TSRKFinalizePackage"))
"""
function TSRKFinalizePackage(petsclib::PetscLibType) end

@for_petsc function TSRKFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:TSRKFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	TSRKRegister(petsclib::PetscLibType,name::TSRKType, order::PetscInt, s::PetscInt, A::Vector{PetscReal}, b::Vector{PetscReal}, c::Vector{PetscReal}, bembed::Vector{PetscReal}, p::PetscInt, binterp::Vector{PetscReal}) 
register an `TSRK` scheme by providing the entries in the Butcher tableau and optionally embedded approximations and interpolation

Not Collective, but the same schemes should be registered on all processes on which they will be used, No Fortran Support

Input Parameters:
- `name`    - identifier for method
- `order`   - approximation order of method
- `s`       - number of stages, this is the dimension of the matrices below
- `A`       - stage coefficients (dimension s*s, row-major)
- `b`       - step completion table (dimension s; NULL to use last row of A)
- `c`       - abscissa (dimension s; NULL to use row sums of A)
- `bembed`  - completion table for embedded method (dimension s; NULL if not available)
- `p`       - Order of the interpolation scheme, equal to the number of columns of binterp
- `binterp` - Coefficients of the interpolation formula (dimension s*p; NULL to reuse b with p=1)

Level: advanced

-seealso: [](ch_ts), `TSRK`

# External Links
$(_doc_external("Ts/TSRKRegister"))
"""
function TSRKRegister(petsclib::PetscLibType, name::TSRKType, order::PetscInt, s::PetscInt, A::Vector{PetscReal}, b::Vector{PetscReal}, c::Vector{PetscReal}, bembed::Vector{PetscReal}, p::PetscInt, binterp::Vector{PetscReal}) end

@for_petsc function TSRKRegister(petsclib::$UnionPetscLib, name::TSRKType, order::$PetscInt, s::$PetscInt, A::Vector{$PetscReal}, b::Vector{$PetscReal}, c::Vector{$PetscReal}, bembed::Vector{$PetscReal}, p::$PetscInt, binterp::Vector{$PetscReal} )

    @chk ccall(
               (:TSRKRegister, $petsc_library),
               PetscErrorCode,
               (TSRKType, $PetscInt, $PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, $PetscInt, Ptr{$PetscReal}),
               name, order, s, A, b, c, bembed, p, binterp,
              )


	return nothing
end 

"""
	s::PetscInt,A::PetscReal,b::PetscReal,c::PetscReal,bembed::PetscReal,p::PetscInt,binterp::PetscReal,FSAL::PetscBool = TSRKGetTableau(petsclib::PetscLibType,ts::TS) 
Get info on the `TSRK` tableau

Not Collective

Input Parameter:
- `ts` - timestepping context

Output Parameters:
- `s`       - number of stages, this is the dimension of the matrices below
- `A`       - stage coefficients (dimension s*s, row-major)
- `b`       - step completion table (dimension s)
- `c`       - abscissa (dimension s)
- `bembed`  - completion table for embedded method (dimension s; NULL if not available)
- `p`       - Order of the interpolation scheme, equal to the number of columns of binterp
- `binterp` - Coefficients of the interpolation formula (dimension s*p)
- `FSAL`    - whether or not the scheme has the First Same As Last property

Level: developer

-seealso: [](ch_ts), `TSRK`, `TSRKRegister()`, `TSRKSetType()`

# External Links
$(_doc_external("Ts/TSRKGetTableau"))
"""
function TSRKGetTableau(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSRKGetTableau(petsclib::$UnionPetscLib, ts::TS )
	s_ = Ref{$PetscInt}()
	A_ = Ref{$PetscReal}()
	b_ = Ref{$PetscReal}()
	c_ = Ref{$PetscReal}()
	bembed_ = Ref{$PetscReal}()
	p_ = Ref{$PetscInt}()
	binterp_ = Ref{$PetscReal}()
	FSAL_ = Ref{PetscBool}()

    @chk ccall(
               (:TSRKGetTableau, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{$PetscInt}, $PetscReal, $PetscReal, $PetscReal, $PetscReal, Ptr{$PetscInt}, $PetscReal, Ptr{PetscBool}),
               ts, s_, A_, b_, c_, bembed_, p_, binterp_, FSAL_,
              )

	s = s_[]
	A = A_[]
	b = b_[]
	c = c_[]
	bembed = bembed_[]
	p = p_[]
	binterp = binterp_[]
	FSAL = FSAL_[]

	return s,A,b,c,bembed,p,binterp,FSAL
end 

"""
	order::PetscInt = TSRKGetOrder(petsclib::PetscLibType,ts::TS) 
Get the order of the `TSRK` scheme

Not Collective

Input Parameter:
- `ts` - timestepping context

Output Parameter:
- `order` - order of `TSRK` scheme

Level: intermediate

-seealso: [](ch_ts), `TSRK`, `TSRKGetType()`

# External Links
$(_doc_external("Ts/TSRKGetOrder"))
"""
function TSRKGetOrder(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSRKGetOrder(petsclib::$UnionPetscLib, ts::TS )
	order_ = Ref{$PetscInt}()

    @chk ccall(
               (:TSRKGetOrder, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{$PetscInt}),
               ts, order_,
              )

	order = order_[]

	return order
end 

"""
	TSRKSetType(petsclib::PetscLibType,ts::TS, rktype::TSRKType) 
Set the type of the `TSRK` scheme

Logically Collective

Input Parameters:
- `ts`     - timestepping context
- `rktype` - type of `TSRK` scheme

Options Database Key:
- `-ts_rk_type` - <1fe,2a,3,3bs,4,5f,5dp,5bs>

Level: intermediate

-seealso: [](ch_ts), `TSRKGetType()`, `TSRK`, `TSRKType`, `TSRK1FE`, `TSRK2A`, `TSRK2B`, `TSRK3`, `TSRK3BS`, `TSRK4`, `TSRK5F`, `TSRK5DP`, `TSRK5BS`, `TSRK6VR`, `TSRK7VR`, `TSRK8VR`

# External Links
$(_doc_external("Ts/TSRKSetType"))
"""
function TSRKSetType(petsclib::PetscLibType, ts::TS, rktype::TSRKType) end

@for_petsc function TSRKSetType(petsclib::$UnionPetscLib, ts::TS, rktype::TSRKType )

    @chk ccall(
               (:TSRKSetType, $petsc_library),
               PetscErrorCode,
               (CTS, TSRKType),
               ts, rktype,
              )


	return nothing
end 

"""
	rktype::TSRKType = TSRKGetType(petsclib::PetscLibType,ts::TS) 
Get the type of `TSRK` scheme

Not Collective

Input Parameter:
- `ts` - timestepping context

Output Parameter:
- `rktype` - type of `TSRK`-scheme

Level: intermediate

-seealso: [](ch_ts), `TSRKSetType()`

# External Links
$(_doc_external("Ts/TSRKGetType"))
"""
function TSRKGetType(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSRKGetType(petsclib::$UnionPetscLib, ts::TS )
	rktype_ = Ref{TSRKType}()

    @chk ccall(
               (:TSRKGetType, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{TSRKType}),
               ts, rktype_,
              )

	rktype = unsafe_string(rktype_[])

	return rktype
end 

"""
	TSRKSetMultirate(petsclib::PetscLibType,ts::TS, use_multirate::PetscBool) 
Use the interpolation

Logically Collective

Input Parameters:
- `ts`            - timestepping context
- `use_multirate` - `PETSC_TRUE` enables the multirate `TSRK` method, sets the basic method to be RK2A and sets the ratio between slow stepsize and fast stepsize to be 2

Options Database Key:
- `-ts_rk_multirate` - <true,false>

Level: intermediate

-seealso: [](ch_ts), `TSRK`, `TSRKGetMultirate()`

# External Links
$(_doc_external("Ts/TSRKSetMultirate"))
"""
function TSRKSetMultirate(petsclib::PetscLibType, ts::TS, use_multirate::PetscBool) end

@for_petsc function TSRKSetMultirate(petsclib::$UnionPetscLib, ts::TS, use_multirate::PetscBool )

    @chk ccall(
               (:TSRKSetMultirate, $petsc_library),
               PetscErrorCode,
               (CTS, PetscBool),
               ts, use_multirate,
              )


	return nothing
end 

"""
	use_multirate::PetscBool = TSRKGetMultirate(petsclib::PetscLibType,ts::TS) 
Gets whether to use the interpolation

Not Collective

Input Parameter:
- `ts` - timestepping context

Output Parameter:
- `use_multirate` - `PETSC_TRUE` if the multirate RK method is enabled, `PETSC_FALSE` otherwise

Level: intermediate

-seealso: [](ch_ts), `TSRK`, `TSRKSetMultirate()`

# External Links
$(_doc_external("Ts/TSRKGetMultirate"))
"""
function TSRKGetMultirate(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSRKGetMultirate(petsclib::$UnionPetscLib, ts::TS )
	use_multirate_ = Ref{PetscBool}()

    @chk ccall(
               (:TSRKGetMultirate, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{PetscBool}),
               ts, use_multirate_,
              )

	use_multirate = use_multirate_[]

	return use_multirate
end 

"""
	TSBDFSetOrder(petsclib::PetscLibType,ts::TS, order::PetscInt) 
Set the order of the `TSBDF` method

Logically Collective

Input Parameters:
- `ts`    - timestepping context
- `order` - order of the method

Options Database Key:
- `-ts_bdf_order <order>` - select the order

Level: intermediate

-seealso: `TSBDFGetOrder()`, `TS`, `TSBDF`

# External Links
$(_doc_external("Ts/TSBDFSetOrder"))
"""
function TSBDFSetOrder(petsclib::PetscLibType, ts::TS, order::PetscInt) end

@for_petsc function TSBDFSetOrder(petsclib::$UnionPetscLib, ts::TS, order::$PetscInt )

    @chk ccall(
               (:TSBDFSetOrder, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt),
               ts, order,
              )


	return nothing
end 

"""
	order::PetscInt = TSBDFGetOrder(petsclib::PetscLibType,ts::TS) 
Get the order of the `TSBDF` method

Not Collective

Input Parameter:
- `ts` - timestepping context

Output Parameter:
- `order` - order of the method

Level: intermediate

-seealso: `TSBDFSetOrder()`, `TS`, `TSBDF`

# External Links
$(_doc_external("Ts/TSBDFGetOrder"))
"""
function TSBDFGetOrder(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSBDFGetOrder(petsclib::$UnionPetscLib, ts::TS )
	order_ = Ref{$PetscInt}()

    @chk ccall(
               (:TSBDFGetOrder, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{$PetscInt}),
               ts, order_,
              )

	order = order_[]

	return order
end 

"""
	TSBasicSymplecticRegisterAll(petsclib::PetscLibType) 
Registers all of the basic symplectic integration methods in `TSBASICSYMPLECTIC`

Not Collective, but should be called by all processes which will need the schemes to be registered

Level: advanced

-seealso: [](ch_ts), `TSBASICSYMPLECTIC`, `TSBasicSymplecticRegisterDestroy()`

# External Links
$(_doc_external("Ts/TSBasicSymplecticRegisterAll"))
"""
function TSBasicSymplecticRegisterAll(petsclib::PetscLibType) end

@for_petsc function TSBasicSymplecticRegisterAll(petsclib::$UnionPetscLib)

    @chk ccall(
               (:TSBasicSymplecticRegisterAll, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	TSBasicSymplecticRegisterDestroy(petsclib::PetscLibType) 
Frees the list of schemes that were registered by `TSBasicSymplecticRegister()`.

Not Collective

Level: advanced

-seealso: [](ch_ts), `TSBasicSymplecticRegister()`, `TSBasicSymplecticRegisterAll()`, `TSBASICSYMPLECTIC`

# External Links
$(_doc_external("Ts/TSBasicSymplecticRegisterDestroy"))
"""
function TSBasicSymplecticRegisterDestroy(petsclib::PetscLibType) end

@for_petsc function TSBasicSymplecticRegisterDestroy(petsclib::$UnionPetscLib)

    @chk ccall(
               (:TSBasicSymplecticRegisterDestroy, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	TSBasicSymplecticInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the `TSBASICSYMPLECTIC` package. It is called
from `TSInitializePackage()`.

Level: developer

-seealso: [](ch_ts), `PetscInitialize()`, `TSBASICSYMPLECTIC`

# External Links
$(_doc_external("Ts/TSBasicSymplecticInitializePackage"))
"""
function TSBasicSymplecticInitializePackage(petsclib::PetscLibType) end

@for_petsc function TSBasicSymplecticInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:TSBasicSymplecticInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	TSBasicSymplecticFinalizePackage(petsclib::PetscLibType) 
This function destroys everything in the `TSBASICSYMPLECTIC` package. It is
called from `PetscFinalize()`.

Level: developer

-seealso: [](ch_ts), `PetscFinalize()`, `TSBASICSYMPLECTIC`

# External Links
$(_doc_external("Ts/TSBasicSymplecticFinalizePackage"))
"""
function TSBasicSymplecticFinalizePackage(petsclib::PetscLibType) end

@for_petsc function TSBasicSymplecticFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:TSBasicSymplecticFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	TSBasicSymplecticRegister(petsclib::PetscLibType,name::TSRosWType, order::PetscInt, s::PetscInt, c::Vector{PetscReal}, d::Vector{PetscReal}) 
register a basic symplectic integration scheme by providing the coefficients.

Not Collective, but the same schemes should be registered on all processes on which they will be used

Input Parameters:
- `name`  - identifier for method
- `order` - approximation order of method
- `s`     - number of stages, this is the dimension of the matrices below
- `c`     - coefficients for updating generalized position (dimension s)
- `d`     - coefficients for updating generalized momentum (dimension s)

Level: advanced

-seealso: [](ch_ts), `TSBASICSYMPLECTIC`

# External Links
$(_doc_external("Ts/TSBasicSymplecticRegister"))
"""
function TSBasicSymplecticRegister(petsclib::PetscLibType, name::TSRosWType, order::PetscInt, s::PetscInt, c::Vector{PetscReal}, d::Vector{PetscReal}) end

@for_petsc function TSBasicSymplecticRegister(petsclib::$UnionPetscLib, name::TSRosWType, order::$PetscInt, s::$PetscInt, c::Vector{$PetscReal}, d::Vector{$PetscReal} )

    @chk ccall(
               (:TSBasicSymplecticRegister, $petsc_library),
               PetscErrorCode,
               (TSRosWType, $PetscInt, $PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}),
               name, order, s, c, d,
              )


	return nothing
end 

"""
	TSBasicSymplecticSetType(petsclib::PetscLibType,ts::TS, bsymptype::TSBasicSymplecticType) 
Set the type of the basic symplectic method

Logically Collective

Input Parameters:
- `ts`        - timestepping context
- `bsymptype` - type of the symplectic scheme

Options Database Key:
- `-ts_basicsymplectic_type <scheme>` - select the scheme

Level: intermediate

-seealso: [](ch_ts), `TSBASICSYMPLECTIC`, `TSBasicSymplecticType`

# External Links
$(_doc_external("Ts/TSBasicSymplecticSetType"))
"""
function TSBasicSymplecticSetType(petsclib::PetscLibType, ts::TS, bsymptype::TSBasicSymplecticType) end

@for_petsc function TSBasicSymplecticSetType(petsclib::$UnionPetscLib, ts::TS, bsymptype::TSBasicSymplecticType )

    @chk ccall(
               (:TSBasicSymplecticSetType, $petsc_library),
               PetscErrorCode,
               (CTS, TSBasicSymplecticType),
               ts, bsymptype,
              )


	return nothing
end 

"""
	TSBasicSymplecticGetType(petsclib::PetscLibType,ts::TS, bsymptype::TSBasicSymplecticType) 
Get the type of the basic symplectic method

Logically Collective

Input Parameters:
- `ts`        - timestepping context
- `bsymptype` - type of the basic symplectic scheme

Level: intermediate

-seealso: [](ch_ts), `TSBASICSYMPLECTIC`, `TSBasicSymplecticType`, `TSBasicSymplecticSetType()`

# External Links
$(_doc_external("Ts/TSBasicSymplecticGetType"))
"""
function TSBasicSymplecticGetType(petsclib::PetscLibType, ts::TS, bsymptype::TSBasicSymplecticType) end

@for_petsc function TSBasicSymplecticGetType(petsclib::$UnionPetscLib, ts::TS, bsymptype::TSBasicSymplecticType )

    @chk ccall(
               (:TSBasicSymplecticGetType, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{TSBasicSymplecticType}),
               ts, bsymptype,
              )


	return nothing
end 

"""
	dt::PetscReal = TSPseudoComputeTimeStep(petsclib::PetscLibType,ts::TS) 
Computes the next timestep for a currently running
pseudo-timestepping process.

Collective

Input Parameter:
- `ts` - timestep context

Output Parameter:
- `dt` - newly computed timestep

Level: developer

-seealso: [](ch_ts), `TSPSEUDO`, `TSPseudoTimeStepDefault()`, `TSPseudoSetTimeStep()`

# External Links
$(_doc_external("Ts/TSPseudoComputeTimeStep"))
"""
function TSPseudoComputeTimeStep(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSPseudoComputeTimeStep(petsclib::$UnionPetscLib, ts::TS )
	dt_ = Ref{$PetscReal}()

    @chk ccall(
               (:TSPseudoComputeTimeStep, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{$PetscReal}),
               ts, dt_,
              )

	dt = dt_[]

	return dt
end 

"""
	newdt::PetscReal,flag::PetscBool = TSPseudoVerifyTimeStepDefault(petsclib::PetscLibType,ts::TS, update::PetscVec, dtctx::Cvoid) 
Default code to verify the quality of the last timestep.

Collective, No Fortran Support

Input Parameters:
- `ts`     - the timestep context
- `dtctx`  - unused timestep context
- `update` - latest solution vector

Output Parameters:
- `newdt` - the timestep to use for the next step
- `flag`  - flag indicating whether the last time step was acceptable

Level: advanced

-seealso: [](ch_ts), `TSPSEUDO`, `TSPseudoSetVerifyTimeStep()`, `TSPseudoVerifyTimeStep()`

# External Links
$(_doc_external("Ts/TSPseudoVerifyTimeStepDefault"))
"""
function TSPseudoVerifyTimeStepDefault(petsclib::PetscLibType, ts::TS, update::PetscVec, dtctx::Cvoid) end

@for_petsc function TSPseudoVerifyTimeStepDefault(petsclib::$UnionPetscLib, ts::TS, update::PetscVec, dtctx::Cvoid )
	newdt_ = Ref{$PetscReal}()
	flag_ = Ref{PetscBool}()

    @chk ccall(
               (:TSPseudoVerifyTimeStepDefault, $petsc_library),
               PetscErrorCode,
               (CTS, CVec, Ptr{Cvoid}, Ptr{$PetscReal}, Ptr{PetscBool}),
               ts, update, dtctx, newdt_, flag_,
              )

	newdt = newdt_[]
	flag = flag_[]

	return newdt,flag
end 

"""
	dt::PetscReal,flag::PetscBool = TSPseudoVerifyTimeStep(petsclib::PetscLibType,ts::TS, update::PetscVec) 
Verifies whether the last timestep was acceptable.

Collective

Input Parameters:
- `ts`     - timestep context
- `update` - latest solution vector

Output Parameters:
- `dt`   - newly computed timestep (if it had to shrink)
- `flag` - indicates if current timestep was ok

Level: advanced

-seealso: [](ch_ts), `TSPSEUDO`, `TSPseudoSetVerifyTimeStep()`, `TSPseudoVerifyTimeStepDefault()`

# External Links
$(_doc_external("Ts/TSPseudoVerifyTimeStep"))
"""
function TSPseudoVerifyTimeStep(petsclib::PetscLibType, ts::TS, update::PetscVec) end

@for_petsc function TSPseudoVerifyTimeStep(petsclib::$UnionPetscLib, ts::TS, update::PetscVec )
	dt_ = Ref{$PetscReal}()
	flag_ = Ref{PetscBool}()

    @chk ccall(
               (:TSPseudoVerifyTimeStep, $petsc_library),
               PetscErrorCode,
               (CTS, CVec, Ptr{$PetscReal}, Ptr{PetscBool}),
               ts, update, dt_, flag_,
              )

	dt = dt_[]
	flag = flag_[]

	return dt,flag
end 

"""
	TSPseudoSetVerifyTimeStep(petsclib::PetscLibType,ts::TS, dt::external, ctx::Cvoid) 
Sets a user
last timestep.

Logically Collective

Input Parameters:
- `ts`  - timestep context
- `dt`  - user-defined function to verify timestep
- `ctx` - [optional] user-defined context for private data for the timestep verification routine (may be `NULL`)

Calling sequence of `func`:
- `ts`     - the time-step context
- `update` - latest solution vector
- `ctx`    - [optional] user-defined timestep context
- `newdt`  - the timestep to use for the next step
- `flag`   - flag indicating whether the last time step was acceptable

Level: advanced

-seealso: [](ch_ts), `TSPSEUDO`, `TSPseudoVerifyTimeStepDefault()`, `TSPseudoVerifyTimeStep()`

# External Links
$(_doc_external("Ts/TSPseudoSetVerifyTimeStep"))
"""
function TSPseudoSetVerifyTimeStep(petsclib::PetscLibType, ts::TS, dt::external, ctx::Cvoid) end

@for_petsc function TSPseudoSetVerifyTimeStep(petsclib::$UnionPetscLib, ts::TS, dt::external, ctx::Cvoid )

    @chk ccall(
               (:TSPseudoSetVerifyTimeStep, $petsc_library),
               PetscErrorCode,
               (CTS, external, Ptr{Cvoid}),
               ts, dt, ctx,
              )


	return nothing
end 

"""
	TSPseudoSetTimeStepIncrement(petsclib::PetscLibType,ts::TS, inc::PetscReal) 
Sets the scaling increment applied to
dt when using the TSPseudoTimeStepDefault() routine.

Logically Collective

Input Parameters:
- `ts`  - the timestep context
- `inc` - the scaling factor >= 1.0

Options Database Key:
- `-ts_pseudo_increment <increment>` - set pseudo increment

Level: advanced

-seealso: [](ch_ts), `TSPSEUDO`, `TSPseudoSetTimeStep()`, `TSPseudoTimeStepDefault()`

# External Links
$(_doc_external("Ts/TSPseudoSetTimeStepIncrement"))
"""
function TSPseudoSetTimeStepIncrement(petsclib::PetscLibType, ts::TS, inc::PetscReal) end

@for_petsc function TSPseudoSetTimeStepIncrement(petsclib::$UnionPetscLib, ts::TS, inc::$PetscReal )

    @chk ccall(
               (:TSPseudoSetTimeStepIncrement, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal),
               ts, inc,
              )


	return nothing
end 

"""
	TSPseudoSetMaxTimeStep(petsclib::PetscLibType,ts::TS, maxdt::PetscReal) 
Sets the maximum time step
when using the TSPseudoTimeStepDefault() routine.

Logically Collective

Input Parameters:
- `ts`    - the timestep context
- `maxdt` - the maximum time step, use a non-positive value to deactivate

Options Database Key:
- `-ts_pseudo_max_dt <increment>` - set pseudo max dt

Level: advanced

-seealso: [](ch_ts), `TSPSEUDO`, `TSPseudoSetTimeStep()`, `TSPseudoTimeStepDefault()`

# External Links
$(_doc_external("Ts/TSPseudoSetMaxTimeStep"))
"""
function TSPseudoSetMaxTimeStep(petsclib::PetscLibType, ts::TS, maxdt::PetscReal) end

@for_petsc function TSPseudoSetMaxTimeStep(petsclib::$UnionPetscLib, ts::TS, maxdt::$PetscReal )

    @chk ccall(
               (:TSPseudoSetMaxTimeStep, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal),
               ts, maxdt,
              )


	return nothing
end 

"""
	TSPseudoIncrementDtFromInitialDt(petsclib::PetscLibType,ts::TS) 
Indicates that a new timestep
is computed via the formula   dt = initial_dt*initial_fnorm/current_fnorm   rather than the default update,   dt = current_dt*previous_fnorm/current_fnorm.

Logically Collective

Input Parameter:
- `ts` - the timestep context

Options Database Key:
- `-ts_pseudo_increment_dt_from_initial_dt <true,false>` - use the initial dt to determine increment

Level: advanced

-seealso: [](ch_ts), `TSPSEUDO`, `TSPseudoSetTimeStep()`, `TSPseudoTimeStepDefault()`

# External Links
$(_doc_external("Ts/TSPseudoIncrementDtFromInitialDt"))
"""
function TSPseudoIncrementDtFromInitialDt(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSPseudoIncrementDtFromInitialDt(petsclib::$UnionPetscLib, ts::TS )

    @chk ccall(
               (:TSPseudoIncrementDtFromInitialDt, $petsc_library),
               PetscErrorCode,
               (CTS,),
               ts,
              )


	return nothing
end 

"""
	TSPseudoSetTimeStep(petsclib::PetscLibType,ts::TS, dt::external, ctx::Cvoid) 
Sets the user
called at each pseudo-timestep to update the timestep.

Logically Collective

Input Parameters:
- `ts`  - timestep context
- `dt`  - function to compute timestep
- `ctx` - [optional] user-defined context for private data required by the function (may be `NULL`)

Calling sequence of `dt`:
- `ts`    - the `TS` context
- `newdt` - the newly computed timestep
- `ctx`   - [optional] user-defined context

Level: intermediate

-seealso: [](ch_ts), `TSPSEUDO`, `TSPseudoTimeStepDefault()`, `TSPseudoComputeTimeStep()`

# External Links
$(_doc_external("Ts/TSPseudoSetTimeStep"))
"""
function TSPseudoSetTimeStep(petsclib::PetscLibType, ts::TS, dt::external, ctx::Cvoid) end

@for_petsc function TSPseudoSetTimeStep(petsclib::$UnionPetscLib, ts::TS, dt::external, ctx::Cvoid )

    @chk ccall(
               (:TSPseudoSetTimeStep, $petsc_library),
               PetscErrorCode,
               (CTS, external, Ptr{Cvoid}),
               ts, dt, ctx,
              )


	return nothing
end 

"""
	newdt::PetscReal = TSPseudoTimeStepDefault(petsclib::PetscLibType,ts::TS, dtctx::Cvoid) 
Default code to compute pseudo

Collective, No Fortran Support

Input Parameters:
- `ts`    - the timestep context
- `dtctx` - unused timestep context

Output Parameter:
- `newdt` - the timestep to use for the next step

Level: advanced

-seealso: [](ch_ts), `TSPseudoSetTimeStep()`, `TSPseudoComputeTimeStep()`, `TSPSEUDO`

# External Links
$(_doc_external("Ts/TSPseudoTimeStepDefault"))
"""
function TSPseudoTimeStepDefault(petsclib::PetscLibType, ts::TS, dtctx::Cvoid) end

@for_petsc function TSPseudoTimeStepDefault(petsclib::$UnionPetscLib, ts::TS, dtctx::Cvoid )
	newdt_ = Ref{$PetscReal}()

    @chk ccall(
               (:TSPseudoTimeStepDefault, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{$PetscReal}, Ptr{Cvoid}),
               ts, newdt_, dtctx,
              )

	newdt = newdt_[]

	return newdt
end 

"""
	TSMPRKRegisterDestroy(petsclib::PetscLibType) 
Frees the list of schemes that were registered by `TSMPRKRegister()`.

Not Collective

Level: advanced

-seealso: [](ch_ts), `TSMPRK`, `TSMPRKRegister()`, `TSMPRKRegisterAll()`

# External Links
$(_doc_external("Ts/TSMPRKRegisterDestroy"))
"""
function TSMPRKRegisterDestroy(petsclib::PetscLibType) end

@for_petsc function TSMPRKRegisterDestroy(petsclib::$UnionPetscLib)

    @chk ccall(
               (:TSMPRKRegisterDestroy, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	TSMPRKInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the `TSMPRK` package. It is called
from `PetscDLLibraryRegister()` when using dynamic libraries, and on the first call to `TSCreate_MPRK()`
when using static libraries.

Level: developer

-seealso: [](ch_ts), `TSMPRK`, `PetscInitialize()`

# External Links
$(_doc_external("Ts/TSMPRKInitializePackage"))
"""
function TSMPRKInitializePackage(petsclib::PetscLibType) end

@for_petsc function TSMPRKInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:TSMPRKInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	TSMPRKFinalizePackage(petsclib::PetscLibType) 
This function destroys everything in the `TSMPRK` package. It is
called from `PetscFinalize()`.

Level: developer

-seealso: [](ch_ts), `TSMPRK`, `PetscFinalize()`

# External Links
$(_doc_external("Ts/TSMPRKFinalizePackage"))
"""
function TSMPRKFinalizePackage(petsclib::PetscLibType) end

@for_petsc function TSMPRKFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:TSMPRKFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	TSMPRKRegister(petsclib::PetscLibType,name::TSMPRKType, order::PetscInt, sbase::PetscInt, ratio1::PetscInt, ratio2::PetscInt, Asb::Vector{PetscReal}, bsb::Vector{PetscReal}, csb::Vector{PetscReal}, rsb::Vector{PetscInt}, Amb::Vector{PetscReal}, bmb::Vector{PetscReal}, cmb::Vector{PetscReal}, rmb::Vector{PetscInt}, Af::Vector{PetscReal}, bf::Vector{PetscReal}, cf::Vector{PetscReal}) 
register a `TSMPRK` scheme by providing the entries in the Butcher tableau

Not Collective, but the same schemes should be registered on all processes on which they will be used, No Fortran Support

Input Parameters:
- `name`   - identifier for method
- `order`  - approximation order of method
- `sbase`  - number of stages in the base methods
- `ratio1` - stepsize ratio at 1st level (e.g. slow/medium)
- `ratio2` - stepsize ratio at 2nd level (e.g. medium/fast)
- `Asb`    - stage coefficients for slow components(dimension s*s, row-major)
- `bsb`    - step completion table for slow components(dimension s)
- `csb`    - abscissa for slow components(dimension s)
- `rsb`    - array of flags for repeated stages for slow components (dimension s)
- `Amb`    - stage coefficients for medium components(dimension s*s, row-major)
- `bmb`    - step completion table for medium components(dimension s)
- `cmb`    - abscissa for medium components(dimension s)
- `rmb`    - array of flags for repeated stages for medium components (dimension s)
- `Af`     - stage coefficients for fast components(dimension s*s, row-major)
- `bf`     - step completion table for fast components(dimension s)
- `cf`     - abscissa for fast components(dimension s)

Level: advanced

-seealso: [](ch_ts), `TSMPRK`

# External Links
$(_doc_external("Ts/TSMPRKRegister"))
"""
function TSMPRKRegister(petsclib::PetscLibType, name::TSMPRKType, order::PetscInt, sbase::PetscInt, ratio1::PetscInt, ratio2::PetscInt, Asb::Vector{PetscReal}, bsb::Vector{PetscReal}, csb::Vector{PetscReal}, rsb::Vector{PetscInt}, Amb::Vector{PetscReal}, bmb::Vector{PetscReal}, cmb::Vector{PetscReal}, rmb::Vector{PetscInt}, Af::Vector{PetscReal}, bf::Vector{PetscReal}, cf::Vector{PetscReal}) end

@for_petsc function TSMPRKRegister(petsclib::$UnionPetscLib, name::TSMPRKType, order::$PetscInt, sbase::$PetscInt, ratio1::$PetscInt, ratio2::$PetscInt, Asb::Vector{$PetscReal}, bsb::Vector{$PetscReal}, csb::Vector{$PetscReal}, rsb::Vector{$PetscInt}, Amb::Vector{$PetscReal}, bmb::Vector{$PetscReal}, cmb::Vector{$PetscReal}, rmb::Vector{$PetscInt}, Af::Vector{$PetscReal}, bf::Vector{$PetscReal}, cf::Vector{$PetscReal} )

    @chk ccall(
               (:TSMPRKRegister, $petsc_library),
               PetscErrorCode,
               (TSMPRKType, $PetscInt, $PetscInt, $PetscInt, $PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscInt}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscInt}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
               name, order, sbase, ratio1, ratio2, Asb, bsb, csb, rsb, Amb, bmb, cmb, rmb, Af, bf, cf,
              )


	return nothing
end 

"""
	TSMPRKSetType(petsclib::PetscLibType,ts::TS, mprktype::TSMPRKType) 
Set the type of `TSMPRK` scheme

Not Collective

Input Parameters:
- `ts`       - timestepping context
- `mprktype` - type of `TSMPRK` scheme

Options Database Key:
- `-ts_mprk_type` - <pm2,p2,p3> - select the specific scheme

Level: intermediate

-seealso: [](ch_ts), `TSMPRKGetType()`, `TSMPRK`, `TSMPRKType`

# External Links
$(_doc_external("Ts/TSMPRKSetType"))
"""
function TSMPRKSetType(petsclib::PetscLibType, ts::TS, mprktype::TSMPRKType) end

@for_petsc function TSMPRKSetType(petsclib::$UnionPetscLib, ts::TS, mprktype::TSMPRKType )

    @chk ccall(
               (:TSMPRKSetType, $petsc_library),
               PetscErrorCode,
               (CTS, TSMPRKType),
               ts, mprktype,
              )


	return nothing
end 

"""
	mprktype::TSMPRKType = TSMPRKGetType(petsclib::PetscLibType,ts::TS) 
Get the type of `TSMPRK` scheme

Not Collective

Input Parameter:
- `ts` - timestepping context

Output Parameter:
- `mprktype` - type of `TSMPRK` scheme

Level: intermediate

-seealso: [](ch_ts), `TSMPRK`

# External Links
$(_doc_external("Ts/TSMPRKGetType"))
"""
function TSMPRKGetType(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSMPRKGetType(petsclib::$UnionPetscLib, ts::TS )
	mprktype_ = Ref{TSMPRKType}()

    @chk ccall(
               (:TSMPRKGetType, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{TSMPRKType}),
               ts, mprktype_,
              )

	mprktype = unsafe_string(mprktype_[])

	return mprktype
end 

"""
	TSRosWRegisterDestroy(petsclib::PetscLibType) 
Frees the list of schemes that were registered by `TSRosWRegister()`.

Not Collective

Level: advanced

-seealso: [](ch_ts), `TSRosWRegister()`, `TSRosWRegisterAll()`

# External Links
$(_doc_external("Ts/TSRosWRegisterDestroy"))
"""
function TSRosWRegisterDestroy(petsclib::PetscLibType) end

@for_petsc function TSRosWRegisterDestroy(petsclib::$UnionPetscLib)

    @chk ccall(
               (:TSRosWRegisterDestroy, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	TSRosWInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the `TSROSW` package. It is called
from `TSInitializePackage()`.

Level: developer

-seealso: [](ch_ts), `TSROSW`, `PetscInitialize()`, `TSRosWFinalizePackage()`

# External Links
$(_doc_external("Ts/TSRosWInitializePackage"))
"""
function TSRosWInitializePackage(petsclib::PetscLibType) end

@for_petsc function TSRosWInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:TSRosWInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	TSRosWFinalizePackage(petsclib::PetscLibType) 
This function destroys everything in the `TSROSW` package. It is
called from `PetscFinalize()`.

Level: developer

-seealso: [](ch_ts), `TSROSW`, `PetscFinalize()`, `TSRosWInitializePackage()`

# External Links
$(_doc_external("Ts/TSRosWFinalizePackage"))
"""
function TSRosWFinalizePackage(petsclib::PetscLibType) end

@for_petsc function TSRosWFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:TSRosWFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	TSRosWRegister(petsclib::PetscLibType,name::TSRosWType, order::PetscInt, s::PetscInt, A::Vector{PetscReal}, Gamma::Vector{PetscReal}, b::Vector{PetscReal}, bembed::Vector{PetscReal}, pinterp::PetscInt, binterpt::Vector{PetscReal}) 
register a `TSROSW`, Rosenbrock W scheme by providing the entries in the Butcher tableau and optionally embedded approximations and interpolation

Not Collective, but the same schemes should be registered on all processes on which they will be used

Input Parameters:
- `name`     - identifier for method
- `order`    - approximation order of method
- `s`        - number of stages, this is the dimension of the matrices below
- `A`        - Table of propagated stage coefficients (dimension s*s, row-major), strictly lower triangular
- `Gamma`    - Table of coefficients in implicit stage equations (dimension s*s, row-major), lower triangular with nonzero diagonal
- `b`        - Step completion table (dimension s)
- `bembed`   - Step completion table for a scheme of order one less (dimension s, NULL if no embedded scheme is available)
- `pinterp`  - Order of the interpolation scheme, equal to the number of columns of binterpt
- `binterpt` - Coefficients of the interpolation formula (dimension s*pinterp)

Level: advanced

-seealso: [](ch_ts), `TSROSW`

# External Links
$(_doc_external("Ts/TSRosWRegister"))
"""
function TSRosWRegister(petsclib::PetscLibType, name::TSRosWType, order::PetscInt, s::PetscInt, A::Vector{PetscReal}, Gamma::Vector{PetscReal}, b::Vector{PetscReal}, bembed::Vector{PetscReal}, pinterp::PetscInt, binterpt::Vector{PetscReal}) end

@for_petsc function TSRosWRegister(petsclib::$UnionPetscLib, name::TSRosWType, order::$PetscInt, s::$PetscInt, A::Vector{$PetscReal}, Gamma::Vector{$PetscReal}, b::Vector{$PetscReal}, bembed::Vector{$PetscReal}, pinterp::$PetscInt, binterpt::Vector{$PetscReal} )

    @chk ccall(
               (:TSRosWRegister, $petsc_library),
               PetscErrorCode,
               (TSRosWType, $PetscInt, $PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, $PetscInt, Ptr{$PetscReal}),
               name, order, s, A, Gamma, b, bembed, pinterp, binterpt,
              )


	return nothing
end 

"""
	TSRosWRegisterRos4(petsclib::PetscLibType,name::TSRosWType, gamma::PetscReal, a2::PetscReal, a3::PetscReal, b3::PetscReal, e4::PetscReal) 
register a fourth order Rosenbrock scheme by providing parameter choices

Not Collective, but the same schemes should be registered on all processes on which they will be used

Input Parameters:
- `name`  - identifier for method
- `gamma` - leading coefficient (diagonal entry)
- `a2`    - design parameter, see Table 7.2 of {cite}`wanner1996solving`
- `a3`    - design parameter or `PETSC_DETERMINE` to satisfy one of the order five conditions (Eq 7.22)
- `b3`    - design parameter, see Table 7.2 of {cite}`wanner1996solving`
- `e4`    - design parameter for embedded method, see coefficient E4 in ros4.f code from Hairer

Level: developer

-seealso: [](ch_ts), `TSROSW`, `TSRosWRegister()`

# External Links
$(_doc_external("Ts/TSRosWRegisterRos4"))
"""
function TSRosWRegisterRos4(petsclib::PetscLibType, name::TSRosWType, gamma::PetscReal, a2::PetscReal, a3::PetscReal, b3::PetscReal, e4::PetscReal) end

@for_petsc function TSRosWRegisterRos4(petsclib::$UnionPetscLib, name::TSRosWType, gamma::$PetscReal, a2::$PetscReal, a3::$PetscReal, b3::$PetscReal, e4::$PetscReal )

    @chk ccall(
               (:TSRosWRegisterRos4, $petsc_library),
               PetscErrorCode,
               (TSRosWType, $PetscReal, $PetscReal, $PetscReal, $PetscReal, $PetscReal),
               name, gamma, a2, a3, b3, e4,
              )


	return nothing
end 

"""
	TSRosWSetType(petsclib::PetscLibType,ts::TS, roswtype::TSRosWType) 
Set the type of Rosenbrock

Logically Collective

Input Parameters:
- `ts`       - timestepping context
- `roswtype` - type of Rosenbrock-W scheme

Level: beginner

-seealso: [](ch_ts), `TSRosWGetType()`, `TSROSW`, `TSROSW2M`, `TSROSW2P`, `TSROSWRA3PW`, `TSROSWRA34PW2`, `TSROSWRODAS3`, `TSROSWSANDU3`, `TSROSWASSP3P3S1C`, `TSROSWLASSP3P4S2C`, `TSROSWLLSSP3P4S2C`, `TSROSWARK3`

# External Links
$(_doc_external("Ts/TSRosWSetType"))
"""
function TSRosWSetType(petsclib::PetscLibType, ts::TS, roswtype::TSRosWType) end

@for_petsc function TSRosWSetType(petsclib::$UnionPetscLib, ts::TS, roswtype::TSRosWType )

    @chk ccall(
               (:TSRosWSetType, $petsc_library),
               PetscErrorCode,
               (CTS, TSRosWType),
               ts, roswtype,
              )


	return nothing
end 

"""
	rostype::TSRosWType = TSRosWGetType(petsclib::PetscLibType,ts::TS) 
Get the type of Rosenbrock

Logically Collective

Input Parameter:
- `ts` - timestepping context

Output Parameter:
- `rostype` - type of Rosenbrock-W scheme

Level: intermediate

-seealso: [](ch_ts), `TSRosWType`, `TSRosWSetType()`

# External Links
$(_doc_external("Ts/TSRosWGetType"))
"""
function TSRosWGetType(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSRosWGetType(petsclib::$UnionPetscLib, ts::TS )
	rostype_ = Ref{TSRosWType}()

    @chk ccall(
               (:TSRosWGetType, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{TSRosWType}),
               ts, rostype_,
              )

	rostype = unsafe_string(rostype_[])

	return rostype
end 

"""
	TSRosWSetRecomputeJacobian(petsclib::PetscLibType,ts::TS, flg::PetscBool) 
Set whether to recompute the Jacobian at each stage. The default is to update the Jacobian once per step.

Logically Collective

Input Parameters:
- `ts`  - timestepping context
- `flg` - `PETSC_TRUE` to recompute the Jacobian at each stage

Level: intermediate

-seealso: [](ch_ts), `TSRosWType`, `TSRosWGetType()`

# External Links
$(_doc_external("Ts/TSRosWSetRecomputeJacobian"))
"""
function TSRosWSetRecomputeJacobian(petsclib::PetscLibType, ts::TS, flg::PetscBool) end

@for_petsc function TSRosWSetRecomputeJacobian(petsclib::$UnionPetscLib, ts::TS, flg::PetscBool )

    @chk ccall(
               (:TSRosWSetRecomputeJacobian, $petsc_library),
               PetscErrorCode,
               (CTS, PetscBool),
               ts, flg,
              )


	return nothing
end 

"""
	TSIRKTableauCreate(petsclib::PetscLibType,ts::TS, nstages::PetscInt, A::PetscReal, b::PetscReal, c::PetscReal, binterp::PetscReal, A_inv::PetscScalar, A_inv_rowsum::PetscScalar, I_s::PetscScalar) 
create the tableau for `TSIRK` and provide the entries

Not Collective

Input Parameters:
- `ts`           - timestepping context
- `nstages`      - number of stages, this is the dimension of the matrices below
- `A`            - stage coefficients (dimension nstages*nstages, row-major)
- `b`            - step completion table (dimension nstages)
- `c`            - abscissa (dimension nstages)
- `binterp`      - coefficients of the interpolation formula (dimension nstages)
- `A_inv`        - inverse of A (dimension nstages*nstages, row-major)
- `A_inv_rowsum` - row sum of the inverse of A (dimension nstages)
- `I_s`          - identity matrix (dimension nstages*nstages)

Level: advanced

-seealso: [](ch_ts), `TSIRK`, `TSIRKRegister()`

# External Links
$(_doc_external("Ts/TSIRKTableauCreate"))
"""
function TSIRKTableauCreate(petsclib::PetscLibType, ts::TS, nstages::PetscInt, A::PetscReal, b::PetscReal, c::PetscReal, binterp::PetscReal, A_inv::PetscScalar, A_inv_rowsum::PetscScalar, I_s::PetscScalar) end

@for_petsc function TSIRKTableauCreate(petsclib::$UnionPetscLib, ts::TS, nstages::$PetscInt, A::$PetscReal, b::$PetscReal, c::$PetscReal, binterp::$PetscReal, A_inv::$PetscScalar, A_inv_rowsum::$PetscScalar, I_s::$PetscScalar )

    @chk ccall(
               (:TSIRKTableauCreate, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscScalar}, Ptr{$PetscScalar}, Ptr{$PetscScalar}),
               ts, nstages, A, b, c, binterp, A_inv, A_inv_rowsum, I_s,
              )


	return nothing
end 

"""
	TSIRKRegister(petsclib::PetscLibType,sname::String, fnc::external) 
adds a `TSIRK` implementation

Not Collective, No Fortran Support

Input Parameters:
- `sname`    - name of user-defined IRK scheme
- `function` - function to create method context

Level: advanced

-seealso: [](ch_ts), `TSIRK`, `TSIRKRegisterAll()`

# External Links
$(_doc_external("Ts/TSIRKRegister"))
"""
function TSIRKRegister(petsclib::PetscLibType, sname::String, fnc::external) end

@for_petsc function TSIRKRegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:TSIRKRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	TSIRKRegisterDestroy(petsclib::PetscLibType) 
Frees the list of schemes that were registered by `TSIRKRegister()`.

Not Collective

Level: advanced

-seealso: [](ch_ts), `TSIRK`, `TSIRKRegister()`, `TSIRKRegisterAll()`

# External Links
$(_doc_external("Ts/TSIRKRegisterDestroy"))
"""
function TSIRKRegisterDestroy(petsclib::PetscLibType) end

@for_petsc function TSIRKRegisterDestroy(petsclib::$UnionPetscLib)

    @chk ccall(
               (:TSIRKRegisterDestroy, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	TSIRKInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the `TSIRK` package. It is called
from `TSInitializePackage()`.

Level: developer

-seealso: [](ch_ts), `TSIRK`, `PetscInitialize()`, `TSIRKFinalizePackage()`, `TSInitializePackage()`

# External Links
$(_doc_external("Ts/TSIRKInitializePackage"))
"""
function TSIRKInitializePackage(petsclib::PetscLibType) end

@for_petsc function TSIRKInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:TSIRKInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	TSIRKFinalizePackage(petsclib::PetscLibType) 
This function destroys everything in the `TSIRK` package. It is
called from `PetscFinalize()`.

Level: developer

-seealso: [](ch_ts), `TSIRK`, `PetscFinalize()`, `TSInitializePackage()`

# External Links
$(_doc_external("Ts/TSIRKFinalizePackage"))
"""
function TSIRKFinalizePackage(petsclib::PetscLibType) end

@for_petsc function TSIRKFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:TSIRKFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	TSIRKSetType(petsclib::PetscLibType,ts::TS, irktype::TSIRKType) 
Set the type of `TSIRK` scheme to use

Logically Collective

Input Parameters:
- `ts`      - timestepping context
- `irktype` - type of `TSIRK` scheme

Options Database Key:
- `-ts_irk_type <gauss>` - set irk type

Level: intermediate

-seealso: [](ch_ts), `TSIRKGetType()`, `TSIRK`, `TSIRKType`, `TSIRKGAUSS`

# External Links
$(_doc_external("Ts/TSIRKSetType"))
"""
function TSIRKSetType(petsclib::PetscLibType, ts::TS, irktype::TSIRKType) end

@for_petsc function TSIRKSetType(petsclib::$UnionPetscLib, ts::TS, irktype::TSIRKType )

    @chk ccall(
               (:TSIRKSetType, $petsc_library),
               PetscErrorCode,
               (CTS, TSIRKType),
               ts, irktype,
              )


	return nothing
end 

"""
	irktype::TSIRKType = TSIRKGetType(petsclib::PetscLibType,ts::TS) 
Get the type of `TSIRK` IMEX scheme being used

Logically Collective

Input Parameter:
- `ts` - timestepping context

Output Parameter:
- `irktype` - type of `TSIRK` IMEX scheme

Level: intermediate

-seealso: [](ch_ts), `TSIRK`, `TSIRKType`, `TSIRKGAUSS`

# External Links
$(_doc_external("Ts/TSIRKGetType"))
"""
function TSIRKGetType(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSIRKGetType(petsclib::$UnionPetscLib, ts::TS )
	irktype_ = Ref{TSIRKType}()

    @chk ccall(
               (:TSIRKGetType, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{TSIRKType}),
               ts, irktype_,
              )

	irktype = unsafe_string(irktype_[])

	return irktype
end 

"""
	TSIRKSetNumStages(petsclib::PetscLibType,ts::TS, nstages::PetscInt) 
Set the number of stages of `TSIRK` scheme to use

Logically Collective

Input Parameters:
- `ts`      - timestepping context
- `nstages` - number of stages of `TSIRK` scheme

Options Database Key:
- `-ts_irk_nstages <int>` - set number of stages

Level: intermediate

-seealso: [](ch_ts), `TSIRKGetNumStages()`, `TSIRK`

# External Links
$(_doc_external("Ts/TSIRKSetNumStages"))
"""
function TSIRKSetNumStages(petsclib::PetscLibType, ts::TS, nstages::PetscInt) end

@for_petsc function TSIRKSetNumStages(petsclib::$UnionPetscLib, ts::TS, nstages::$PetscInt )

    @chk ccall(
               (:TSIRKSetNumStages, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt),
               ts, nstages,
              )


	return nothing
end 

"""
	TSIRKGetNumStages(petsclib::PetscLibType,ts::TS, nstages::PetscInt) 
Get the number of stages of `TSIRK` scheme

Logically Collective

Input Parameters:
- `ts`      - timestepping context
- `nstages` - number of stages of `TSIRK` scheme

Level: intermediate

-seealso: [](ch_ts), `TSIRKSetNumStages()`, `TSIRK`

# External Links
$(_doc_external("Ts/TSIRKGetNumStages"))
"""
function TSIRKGetNumStages(petsclib::PetscLibType, ts::TS, nstages::PetscInt) end

@for_petsc function TSIRKGetNumStages(petsclib::$UnionPetscLib, ts::TS, nstages::$PetscInt )

    @chk ccall(
               (:TSIRKGetNumStages, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{$PetscInt}),
               ts, nstages,
              )


	return nothing
end 

"""
	TSDiscGradSetFormulation(petsclib::PetscLibType,ts::TS, Sfunc::external, Ffunc::external, Gfunc::external, ctx::Cvoid) 
Set the construction method for S, F, and grad F from the
formulation u_t = S(u) \nabla F(u) for `TSDISCGRAD`

Not Collective

Input Parameters:
- `ts`    - timestepping context
- `Sfunc` - constructor for the S matrix from the formulation
- `Ffunc` - functional F from the formulation
- `Gfunc` - constructor for the gradient of F from the formulation
- `ctx`   - optional context for the functions

Calling sequence of `Sfunc`:
- `ts`   - the integrator
- `time` - the current time
- `u`    - the solution
- `S`    - the S-matrix from the formulation
- `ctx`  - the user context

Calling sequence of `Ffunc`:
- `ts`   - the integrator
- `time` - the current time
- `u`    - the solution
- `F`    - the computed function from the formulation
- `ctx`  - the user context

Calling sequence of `Gfunc`:
- `ts`   - the integrator
- `time` - the current time
- `u`    - the solution
- `G`    - the gradient of the computed function from the formulation
- `ctx`  - the user context

Level: intermediate

-seealso: [](ch_ts), `TSDISCGRAD`, `TSDiscGradGetFormulation()`

# External Links
$(_doc_external("Ts/TSDiscGradSetFormulation"))
"""
function TSDiscGradSetFormulation(petsclib::PetscLibType, ts::TS, Sfunc::external, Ffunc::external, Gfunc::external, ctx::Cvoid) end

@for_petsc function TSDiscGradSetFormulation(petsclib::$UnionPetscLib, ts::TS, Sfunc::external, Ffunc::external, Gfunc::external, ctx::Cvoid )

    @chk ccall(
               (:TSDiscGradSetFormulation, $petsc_library),
               PetscErrorCode,
               (CTS, external, external, external, Ptr{Cvoid}),
               ts, Sfunc, Ffunc, Gfunc, ctx,
              )


	return nothing
end 

"""
	dgtype::TSDGType = TSDiscGradGetType(petsclib::PetscLibType,ts::TS) 
Checks for which discrete gradient to use in formulation for `TSDISCGRAD`

Not Collective

Input Parameter:
- `ts` - timestepping context

Output Parameter:
- `dgtype` - Discrete gradient type <none, gonzalez, average>

Level: advanced

-seealso: [](ch_ts), `TSDISCGRAD`, `TSDiscGradSetType()`

# External Links
$(_doc_external("Ts/TSDiscGradGetType"))
"""
function TSDiscGradGetType(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSDiscGradGetType(petsclib::$UnionPetscLib, ts::TS )
	dgtype_ = Ref{TSDGType}()

    @chk ccall(
               (:TSDiscGradGetType, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{TSDGType}),
               ts, dgtype_,
              )

	dgtype = unsafe_string(dgtype_[])

	return dgtype
end 

"""
	TSDiscGradSetType(petsclib::PetscLibType,ts::TS, dgtype::TSDGType) 
Sets discrete gradient formulation.

Not Collective

Input Parameters:
- `ts`     - timestepping context
- `dgtype` - Discrete gradient type <none, gonzalez, average>

Options Database Key:
- `-ts_discgrad_type <type>` - flag to choose discrete gradient type

Level: intermediate

-seealso: [](ch_ts), `TSDISCGRAD`

# External Links
$(_doc_external("Ts/TSDiscGradSetType"))
"""
function TSDiscGradSetType(petsclib::PetscLibType, ts::TS, dgtype::TSDGType) end

@for_petsc function TSDiscGradSetType(petsclib::$UnionPetscLib, ts::TS, dgtype::TSDGType )

    @chk ccall(
               (:TSDiscGradSetType, $petsc_library),
               PetscErrorCode,
               (CTS, TSDGType),
               ts, dgtype,
              )


	return nothing
end 

"""
	TSSundialsGetIterations(petsclib::PetscLibType,ts::TS, nonlin::Cint, lin::Cint) 
Gets the number of nonlinear and linear iterations used so far by `TSSUNDIALS`.

Not Collective

Input Parameter:
- `ts` - the time-step context

Output Parameters:
- `nonlin` - number of nonlinear iterations
- `lin`    - number of linear iterations

Level: advanced

-seealso: [](ch_ts), `TSSundialsSetType()`, `TSSundialsSetMaxl()`,
`TSSundialsSetLinearTolerance()`, `TSSundialsSetGramSchmidtType()`, `TSSundialsSetTolerance()`,
`TSSundialsGetPC()`, `TSSetExactFinalTime()`

# External Links
$(_doc_external("Ts/TSSundialsGetIterations"))
"""
function TSSundialsGetIterations(petsclib::PetscLibType, ts::TS, nonlin::Cint, lin::Cint) end

@for_petsc function TSSundialsGetIterations(petsclib::$UnionPetscLib, ts::TS, nonlin::Cint, lin::Cint )

    @chk ccall(
               (:TSSundialsGetIterations, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{Cint}, Ptr{Cint}),
               ts, nonlin, lin,
              )


	return nothing
end 

"""
	TSSundialsSetType(petsclib::PetscLibType,ts::TS, type::TSSundialsLmmType) 
Sets the method that `TSSUNDIALS` will use for integration.

Logically Collective

Input Parameters:
- `ts`   - the time-step context
- `type` - one of  `SUNDIALS_ADAMS` or `SUNDIALS_BDF`

Level: intermediate

-seealso: [](ch_ts), `TSSundialsGetIterations()`, `TSSundialsSetMaxl()`,
`TSSundialsSetLinearTolerance()`, `TSSundialsSetGramSchmidtType()`,
`TSSundialsSetTolerance()`, `TSSundialsGetPC()`, `TSSetExactFinalTime()`

# External Links
$(_doc_external("Ts/TSSundialsSetType"))
"""
function TSSundialsSetType(petsclib::PetscLibType, ts::TS, type::TSSundialsLmmType) end

@for_petsc function TSSundialsSetType(petsclib::$UnionPetscLib, ts::TS, type::TSSundialsLmmType )

    @chk ccall(
               (:TSSundialsSetType, $petsc_library),
               PetscErrorCode,
               (CTS, TSSundialsLmmType),
               ts, type,
              )


	return nothing
end 

"""
	TSSundialsSetMaxord(petsclib::PetscLibType,ts::TS, maxord::PetscInt) 
Sets the maximum order for BDF/Adams method used by `TSSUNDIALS`.

Logically Collective

Input Parameters:
- `ts`     - the time-step context
- `maxord` - maximum order of BDF / Adams method

Level: advanced

-seealso: [](ch_ts), `TSSundialsGetIterations()`, `TSSundialsSetType()`,
`TSSundialsSetLinearTolerance()`, `TSSundialsSetGramSchmidtType()`,
`TSSundialsSetTolerance()`, `TSSundialsGetPC()`, `TSSetExactFinalTime()`

# External Links
$(_doc_external("Ts/TSSundialsSetMaxord"))
"""
function TSSundialsSetMaxord(petsclib::PetscLibType, ts::TS, maxord::PetscInt) end

@for_petsc function TSSundialsSetMaxord(petsclib::$UnionPetscLib, ts::TS, maxord::$PetscInt )

    @chk ccall(
               (:TSSundialsSetMaxord, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt),
               ts, maxord,
              )


	return nothing
end 

"""
	TSSundialsSetMaxl(petsclib::PetscLibType,ts::TS, maxl::PetscInt) 
Sets the dimension of the Krylov space used by
GMRES in the linear solver in `TSSUNDIALS`. `TSSUNDIALS` DOES NOT use restarted GMRES so
this is the maximum number of GMRES steps that will be used.

Logically Collective

Input Parameters:
- `ts`   - the time-step context
- `maxl` - number of direction vectors (the dimension of Krylov subspace).

Level: advanced

-seealso: [](ch_ts), `TSSundialsGetIterations()`, `TSSundialsSetType()`,
`TSSundialsSetLinearTolerance()`, `TSSundialsSetGramSchmidtType()`, `TSSundialsSetTolerance()`,
`TSSundialsGetPC()`, `TSSetExactFinalTime()`

# External Links
$(_doc_external("Ts/TSSundialsSetMaxl"))
"""
function TSSundialsSetMaxl(petsclib::PetscLibType, ts::TS, maxl::PetscInt) end

@for_petsc function TSSundialsSetMaxl(petsclib::$UnionPetscLib, ts::TS, maxl::$PetscInt )

    @chk ccall(
               (:TSSundialsSetMaxl, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt),
               ts, maxl,
              )


	return nothing
end 

"""
	TSSundialsSetLinearTolerance(petsclib::PetscLibType,ts::TS, tol::PetscReal) 
Sets the tolerance used to solve the linear
system by `TSSUNDIALS`.

Logically Collective

Input Parameters:
- `ts`  - the time-step context
- `tol` - the factor by which the tolerance on the nonlinear solver is
multiplied to get the tolerance on the linear solver, .05 by default.

Level: advanced

-seealso: [](ch_ts), `TSSundialsGetIterations()`, `TSSundialsSetType()`, `TSSundialsSetMaxl()`,
`TSSundialsSetGramSchmidtType()`, `TSSundialsSetTolerance()`,
`TSSundialsGetPC()`,
`TSSetExactFinalTime()`

# External Links
$(_doc_external("Ts/TSSundialsSetLinearTolerance"))
"""
function TSSundialsSetLinearTolerance(petsclib::PetscLibType, ts::TS, tol::PetscReal) end

@for_petsc function TSSundialsSetLinearTolerance(petsclib::$UnionPetscLib, ts::TS, tol::$PetscReal )

    @chk ccall(
               (:TSSundialsSetLinearTolerance, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal),
               ts, tol,
              )


	return nothing
end 

"""
	TSSundialsSetGramSchmidtType(petsclib::PetscLibType,ts::TS, type::TSSundialsGramSchmidtType) 
Sets type of orthogonalization used
in GMRES method by `TSSUNDIALS` linear solver.

Logically Collective

Input Parameters:
- `ts`   - the time-step context
- `type` - either `SUNDIALS_MODIFIED_GS` or `SUNDIALS_CLASSICAL_GS`

Level: advanced

-seealso: [](ch_ts), `TSSundialsGetIterations()`, `TSSundialsSetType()`, `TSSundialsSetMaxl()`,
`TSSundialsSetLinearTolerance()`, `TSSundialsSetTolerance()`,
`TSSundialsGetPC()`,
`TSSetExactFinalTime()`

# External Links
$(_doc_external("Ts/TSSundialsSetGramSchmidtType"))
"""
function TSSundialsSetGramSchmidtType(petsclib::PetscLibType, ts::TS, type::TSSundialsGramSchmidtType) end

@for_petsc function TSSundialsSetGramSchmidtType(petsclib::$UnionPetscLib, ts::TS, type::TSSundialsGramSchmidtType )

    @chk ccall(
               (:TSSundialsSetGramSchmidtType, $petsc_library),
               PetscErrorCode,
               (CTS, TSSundialsGramSchmidtType),
               ts, type,
              )


	return nothing
end 

"""
	TSSundialsSetTolerance(petsclib::PetscLibType,ts::TS, aabs::PetscReal, rel::PetscReal) 
Sets the absolute and relative tolerance used by
`TSSUNDIALS` for error control.

Logically Collective

Input Parameters:
- `ts`   - the time-step context
- `aabs` - the absolute tolerance
- `rel`  - the relative tolerance

See the CVODE/SUNDIALS users manual for exact details on these parameters. Essentially
these regulate the size of the error for a SINGLE timestep.

Level: intermediate

-seealso: [](ch_ts), `TSSundialsGetIterations()`, `TSSundialsSetType()`, `TSSundialsSetGMRESMaxl()`,
`TSSundialsSetLinearTolerance()`, `TSSundialsSetGramSchmidtType()`,
`TSSundialsGetPC()`,
`TSSetExactFinalTime()`

# External Links
$(_doc_external("Ts/TSSundialsSetTolerance"))
"""
function TSSundialsSetTolerance(petsclib::PetscLibType, ts::TS, aabs::PetscReal, rel::PetscReal) end

@for_petsc function TSSundialsSetTolerance(petsclib::$UnionPetscLib, ts::TS, aabs::$PetscReal, rel::$PetscReal )

    @chk ccall(
               (:TSSundialsSetTolerance, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal, $PetscReal),
               ts, aabs, rel,
              )


	return nothing
end 

"""
	TSSundialsGetPC(petsclib::PetscLibType,ts::TS, pc::PC) 
Extract the PC context from a time

Input Parameter:
- `ts` - the time-step context

Output Parameter:
- `pc` - the preconditioner context

Level: advanced

-seealso: [](ch_ts), `TSSundialsGetIterations()`, `TSSundialsSetType()`, `TSSundialsSetMaxl()`,
`TSSundialsSetLinearTolerance()`, `TSSundialsSetGramSchmidtType()`, `TSSundialsSetTolerance()`

# External Links
$(_doc_external("Ts/TSSundialsGetPC"))
"""
function TSSundialsGetPC(petsclib::PetscLibType, ts::TS, pc::PC) end

@for_petsc function TSSundialsGetPC(petsclib::$UnionPetscLib, ts::TS, pc::PC )

    @chk ccall(
               (:TSSundialsGetPC, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{PC}),
               ts, pc,
              )


	return nothing
end 

"""
	TSSundialsSetMinTimeStep(petsclib::PetscLibType,ts::TS, mindt::PetscReal) 
Smallest time step to be chosen by the adaptive controller.

Input Parameters:
- `ts`    - the time-step context
- `mindt` - lowest time step if positive, negative to deactivate

-seealso: [](ch_ts), `TSSundialsSetType()`, `TSSundialsSetTolerance()`,

# External Links
$(_doc_external("Ts/TSSundialsSetMinTimeStep"))
"""
function TSSundialsSetMinTimeStep(petsclib::PetscLibType, ts::TS, mindt::PetscReal) end

@for_petsc function TSSundialsSetMinTimeStep(petsclib::$UnionPetscLib, ts::TS, mindt::$PetscReal )

    @chk ccall(
               (:TSSundialsSetMinTimeStep, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal),
               ts, mindt,
              )


	return nothing
end 

"""
	TSSundialsSetMaxTimeStep(petsclib::PetscLibType,ts::TS, maxdt::PetscReal) 
Largest time step to be chosen by the adaptive controller.

Input Parameters:
- `ts`    - the time-step context
- `maxdt` - lowest time step if positive, negative to deactivate

Level: beginner

-seealso: [](ch_ts), `TSSundialsSetType()`, `TSSundialsSetTolerance()`,

# External Links
$(_doc_external("Ts/TSSundialsSetMaxTimeStep"))
"""
function TSSundialsSetMaxTimeStep(petsclib::PetscLibType, ts::TS, maxdt::PetscReal) end

@for_petsc function TSSundialsSetMaxTimeStep(petsclib::$UnionPetscLib, ts::TS, maxdt::$PetscReal )

    @chk ccall(
               (:TSSundialsSetMaxTimeStep, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal),
               ts, maxdt,
              )


	return nothing
end 

"""
	TSSundialsMonitorInternalSteps(petsclib::PetscLibType,ts::TS, ft::PetscBool) 
Monitor `TSSUNDIALS` internal steps (Defaults to false).

Input Parameters:
- `ts` - the time-step context
- `ft` - `PETSC_TRUE` if monitor, else `PETSC_FALSE`

Level: beginner

-seealso: [](ch_ts), `TSSundialsGetIterations()`, `TSSundialsSetType()`, `TSSundialsSetMaxl()`,
`TSSundialsSetLinearTolerance()`, `TSSundialsSetGramSchmidtType()`, `TSSundialsSetTolerance()`,
`TSSundialsGetPC()`

# External Links
$(_doc_external("Ts/TSSundialsMonitorInternalSteps"))
"""
function TSSundialsMonitorInternalSteps(petsclib::PetscLibType, ts::TS, ft::PetscBool) end

@for_petsc function TSSundialsMonitorInternalSteps(petsclib::$UnionPetscLib, ts::TS, ft::PetscBool )

    @chk ccall(
               (:TSSundialsMonitorInternalSteps, $petsc_library),
               PetscErrorCode,
               (CTS, PetscBool),
               ts, ft,
              )


	return nothing
end 

"""
	TSSundialsSetUseDense(petsclib::PetscLibType,ts::TS, use_dense::PetscBool) 
Set a flag to use a dense linear solver in `TSSUNDIALS` (serial only)

Logically Collective

Input Parameters:
- `ts`        - the time-step context
- `use_dense` - `PETSC_TRUE` to use the dense solver

Level: advanced

-seealso: [](ch_ts), `TSSUNDIALS`

# External Links
$(_doc_external("Ts/TSSundialsSetUseDense"))
"""
function TSSundialsSetUseDense(petsclib::PetscLibType, ts::TS, use_dense::PetscBool) end

@for_petsc function TSSundialsSetUseDense(petsclib::$UnionPetscLib, ts::TS, use_dense::PetscBool )

    @chk ccall(
               (:TSSundialsSetUseDense, $petsc_library),
               PetscErrorCode,
               (CTS, PetscBool),
               ts, use_dense,
              )


	return nothing
end 

"""
	theta::PetscReal = TSThetaGetTheta(petsclib::PetscLibType,ts::TS) 
Get the abscissa of the stage in (0,1] for `TSTHETA`

Not Collective

Input Parameter:
- `ts` - timestepping context

Output Parameter:
- `theta` - stage abscissa

Level: advanced

-seealso: [](ch_ts), `TSThetaSetTheta()`, `TSTHETA`

# External Links
$(_doc_external("Ts/TSThetaGetTheta"))
"""
function TSThetaGetTheta(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSThetaGetTheta(petsclib::$UnionPetscLib, ts::TS )
	theta_ = Ref{$PetscReal}()

    @chk ccall(
               (:TSThetaGetTheta, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{$PetscReal}),
               ts, theta_,
              )

	theta = theta_[]

	return theta
end 

"""
	TSThetaSetTheta(petsclib::PetscLibType,ts::TS, theta::PetscReal) 
Set the abscissa of the stage in (0,1]  for `TSTHETA`

Not Collective

Input Parameters:
- `ts`    - timestepping context
- `theta` - stage abscissa

Options Database Key:
- `-ts_theta_theta <theta>` - set theta

Level: intermediate

-seealso: [](ch_ts), `TSThetaGetTheta()`, `TSTHETA`, `TSCN`

# External Links
$(_doc_external("Ts/TSThetaSetTheta"))
"""
function TSThetaSetTheta(petsclib::PetscLibType, ts::TS, theta::PetscReal) end

@for_petsc function TSThetaSetTheta(petsclib::$UnionPetscLib, ts::TS, theta::$PetscReal )

    @chk ccall(
               (:TSThetaSetTheta, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal),
               ts, theta,
              )


	return nothing
end 

"""
	endpoint::PetscBool = TSThetaGetEndpoint(petsclib::PetscLibType,ts::TS) 
Gets whether to use the endpoint variant of the method (e.g. trapezoid/Crank

Not Collective

Input Parameter:
- `ts` - timestepping context

Output Parameter:
- `endpoint` - `PETSC_TRUE` when using the endpoint variant

Level: advanced

-seealso: [](ch_ts), `TSThetaSetEndpoint()`, `TSTHETA`, `TSCN`

# External Links
$(_doc_external("Ts/TSThetaGetEndpoint"))
"""
function TSThetaGetEndpoint(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSThetaGetEndpoint(petsclib::$UnionPetscLib, ts::TS )
	endpoint_ = Ref{PetscBool}()

    @chk ccall(
               (:TSThetaGetEndpoint, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{PetscBool}),
               ts, endpoint_,
              )

	endpoint = endpoint_[]

	return endpoint
end 

"""
	TSThetaSetEndpoint(petsclib::PetscLibType,ts::TS, flg::PetscBool) 
Sets whether to use the endpoint variant of the method (e.g. trapezoid/Crank

Not Collective

Input Parameters:
- `ts`  - timestepping context
- `flg` - `PETSC_TRUE` to use the endpoint variant

Options Database Key:
- `-ts_theta_endpoint <flg>` - use the endpoint variant

Level: intermediate

-seealso: [](ch_ts), `TSTHETA`, `TSCN`

# External Links
$(_doc_external("Ts/TSThetaSetEndpoint"))
"""
function TSThetaSetEndpoint(petsclib::PetscLibType, ts::TS, flg::PetscBool) end

@for_petsc function TSThetaSetEndpoint(petsclib::$UnionPetscLib, ts::TS, flg::PetscBool )

    @chk ccall(
               (:TSThetaSetEndpoint, $petsc_library),
               PetscErrorCode,
               (CTS, PetscBool),
               ts, flg,
              )


	return nothing
end 

"""
	TSGLLESetType(petsclib::PetscLibType,ts::TS, type::TSGLLEType) 
sets the class of general linear method, `TSGLLE` to use for time

Collective

Input Parameters:
- `ts`   - the `TS` context
- `type` - a method

Options Database Key:
- `-ts_gl_type <type>` - sets the method, use -help for a list of available method (e.g. irks)

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSGLLEType`, `TSGLLE`

# External Links
$(_doc_external("Ts/TSGLLESetType"))
"""
function TSGLLESetType(petsclib::PetscLibType, ts::TS, type::TSGLLEType) end

@for_petsc function TSGLLESetType(petsclib::$UnionPetscLib, ts::TS, type::TSGLLEType )

    @chk ccall(
               (:TSGLLESetType, $petsc_library),
               PetscErrorCode,
               (CTS, TSGLLEType),
               ts, type,
              )


	return nothing
end 

"""
	TSGLLESetAcceptType(petsclib::PetscLibType,ts::TS, type::TSGLLEAcceptType) 
sets the acceptance test for `TSGLLE`

Logically Collective

Input Parameters:
- `ts`   - the `TS` context
- `type` - the type

Options Database Key:
- `-ts_gl_accept_type <type>` - sets the method used to determine whether to accept or reject a step

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSGLLE`, `TSGLLEAcceptRegister()`, `TSGLLEAdapt`

# External Links
$(_doc_external("Ts/TSGLLESetAcceptType"))
"""
function TSGLLESetAcceptType(petsclib::PetscLibType, ts::TS, type::TSGLLEAcceptType) end

@for_petsc function TSGLLESetAcceptType(petsclib::$UnionPetscLib, ts::TS, type::TSGLLEAcceptType )

    @chk ccall(
               (:TSGLLESetAcceptType, $petsc_library),
               PetscErrorCode,
               (CTS, TSGLLEAcceptType),
               ts, type,
              )


	return nothing
end 

"""
	TSGLLEGetAdapt(petsclib::PetscLibType,ts::TS, adapt::TSGLLEAdapt) 
gets the `TSGLLEAdapt` object from the `TS`

Not Collective

Input Parameter:
- `ts` - the `TS` context

Output Parameter:
- `adapt` - the `TSGLLEAdapt` context

Level: advanced

-seealso: [](ch_ts), `TS`, `TSGLLE`, `TSGLLEAdapt`, `TSGLLEAdaptRegister()`

# External Links
$(_doc_external("Ts/TSGLLEGetAdapt"))
"""
function TSGLLEGetAdapt(petsclib::PetscLibType, ts::TS, adapt::TSGLLEAdapt) end

@for_petsc function TSGLLEGetAdapt(petsclib::$UnionPetscLib, ts::TS, adapt::TSGLLEAdapt )

    @chk ccall(
               (:TSGLLEGetAdapt, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{TSGLLEAdapt}),
               ts, adapt,
              )


	return nothing
end 

"""
	TSGLLERegister(petsclib::PetscLibType,sname::String, fnc::external) 
adds a `TSGLLE` implementation

Not Collective, No Fortran Support

Input Parameters:
- `sname`    - name of user-defined general linear scheme
- `function` - routine to create method context

Level: advanced

-seealso: [](ch_ts), `TSGLLE`, `TSGLLEType`, `TSGLLERegisterAll()`

# External Links
$(_doc_external("Ts/TSGLLERegister"))
"""
function TSGLLERegister(petsclib::PetscLibType, sname::String, fnc::external) end

@for_petsc function TSGLLERegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:TSGLLERegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	TSGLLEAcceptRegister(petsclib::PetscLibType,sname::String, fnc::TSGLLEAcceptFn) 
adds a `TSGLLE` acceptance scheme

Not Collective

Input Parameters:
- `sname`    - name of user-defined acceptance scheme
- `function` - routine to create method context, see `TSGLLEAcceptFn` for the calling sequence

Level: advanced

-seealso: [](ch_ts), `TSGLLE`, `TSGLLEType`, `TSGLLERegisterAll()`, `TSGLLEAcceptFn`

# External Links
$(_doc_external("Ts/TSGLLEAcceptRegister"))
"""
function TSGLLEAcceptRegister(petsclib::PetscLibType, sname::String, fnc::TSGLLEAcceptFn) end

@for_petsc function TSGLLEAcceptRegister(petsclib::$UnionPetscLib, sname::String, fnc::TSGLLEAcceptFn )

    @chk ccall(
               (:TSGLLEAcceptRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Ptr{TSGLLEAcceptFn}),
               sname, fnc,
              )


	return nothing
end 

"""
	TSGLLEInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the `TSGLLE` package. It is called
from `TSInitializePackage()`.

Level: developer

-seealso: [](ch_ts), `PetscInitialize()`, `TSInitializePackage()`, `TSGLLEFinalizePackage()`

# External Links
$(_doc_external("Ts/TSGLLEInitializePackage"))
"""
function TSGLLEInitializePackage(petsclib::PetscLibType) end

@for_petsc function TSGLLEInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:TSGLLEInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	TSGLLEFinalizePackage(petsclib::PetscLibType) 
This function destroys everything in the `TSGLLE` package. It is
called from `PetscFinalize()`.

Level: developer

-seealso: [](ch_ts), `PetscFinalize()`, `TSGLLEInitializePackage()`, `TSInitializePackage()`

# External Links
$(_doc_external("Ts/TSGLLEFinalizePackage"))
"""
function TSGLLEFinalizePackage(petsclib::PetscLibType) end

@for_petsc function TSGLLEFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:TSGLLEFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	TSAlpha2SetPredictor(petsclib::PetscLibType,ts::TS, predictor::TSAlpha2PredictorFn, ctx::Cvoid) 
sets the callback for computing a predictor (i.e., initial guess
for the nonlinear solver).

Input Parameters:
- `ts`        - timestepping context
- `predictor` - callback to set the predictor in each step
- `ctx`       - the application context, which may be set to `NULL` if not used

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSALPHA2`, `TSAlpha2PredictorFn`

# External Links
$(_doc_external("Ts/TSAlpha2SetPredictor"))
"""
function TSAlpha2SetPredictor(petsclib::PetscLibType, ts::TS, predictor::TSAlpha2PredictorFn, ctx::Cvoid) end

@for_petsc function TSAlpha2SetPredictor(petsclib::$UnionPetscLib, ts::TS, predictor::TSAlpha2PredictorFn, ctx::Cvoid )

    @chk ccall(
               (:TSAlpha2SetPredictor, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{TSAlpha2PredictorFn}, Ptr{Cvoid}),
               ts, predictor, ctx,
              )


	return nothing
end 

"""
	TSAlpha2SetRadius(petsclib::PetscLibType,ts::TS, radius::PetscReal) 
sets the desired spectral radius of the method for `TSALPHA2`
(i.e. high-frequency numerical damping)

Logically Collective

Input Parameters:
- `ts`     - timestepping context
- `radius` - the desired spectral radius

Options Database Key:
- `-ts_alpha_radius <radius>` - set the desired spectral radius

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSALPHA2`, `TSAlpha2SetParams()`, `TSAlpha2GetParams()`

# External Links
$(_doc_external("Ts/TSAlpha2SetRadius"))
"""
function TSAlpha2SetRadius(petsclib::PetscLibType, ts::TS, radius::PetscReal) end

@for_petsc function TSAlpha2SetRadius(petsclib::$UnionPetscLib, ts::TS, radius::$PetscReal )

    @chk ccall(
               (:TSAlpha2SetRadius, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal),
               ts, radius,
              )


	return nothing
end 

"""
	TSAlpha2SetParams(petsclib::PetscLibType,ts::TS, alpha_m::PetscReal, alpha_f::PetscReal, gamma::PetscReal, beta::PetscReal) 
sets the algorithmic parameters for `TSALPHA2`

Logically Collective

Input Parameters:
- `ts`      - timestepping context
- `alpha_m` - algorithmic parameter
- `alpha_f` - algorithmic parameter
- `gamma`   - algorithmic parameter
- `beta`    - algorithmic parameter

Options Database Keys:
- `-ts_alpha_alpha_m <alpha_m>` - set alpha_m
- `-ts_alpha_alpha_f <alpha_f>` - set alpha_f
- `-ts_alpha_gamma   <gamma>`   - set gamma
- `-ts_alpha_beta    <beta>`    - set beta

Level: advanced

-seealso: [](ch_ts), `TS`, `TSALPHA2`, `TSAlpha2SetRadius()`, `TSAlpha2GetParams()`

# External Links
$(_doc_external("Ts/TSAlpha2SetParams"))
"""
function TSAlpha2SetParams(petsclib::PetscLibType, ts::TS, alpha_m::PetscReal, alpha_f::PetscReal, gamma::PetscReal, beta::PetscReal) end

@for_petsc function TSAlpha2SetParams(petsclib::$UnionPetscLib, ts::TS, alpha_m::$PetscReal, alpha_f::$PetscReal, gamma::$PetscReal, beta::$PetscReal )

    @chk ccall(
               (:TSAlpha2SetParams, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal, $PetscReal, $PetscReal, $PetscReal),
               ts, alpha_m, alpha_f, gamma, beta,
              )


	return nothing
end 

"""
	alpha_m::PetscReal,alpha_f::PetscReal,gamma::PetscReal,beta::PetscReal = TSAlpha2GetParams(petsclib::PetscLibType,ts::TS) 
gets the algorithmic parameters for `TSALPHA2`

Not Collective

Input Parameter:
- `ts` - timestepping context

Output Parameters:
- `alpha_m` - algorithmic parameter
- `alpha_f` - algorithmic parameter
- `gamma`   - algorithmic parameter
- `beta`    - algorithmic parameter

Level: advanced

-seealso: [](ch_ts), `TS`, `TSALPHA2`, `TSAlpha2SetRadius()`, `TSAlpha2SetParams()`

# External Links
$(_doc_external("Ts/TSAlpha2GetParams"))
"""
function TSAlpha2GetParams(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSAlpha2GetParams(petsclib::$UnionPetscLib, ts::TS )
	alpha_m_ = Ref{$PetscReal}()
	alpha_f_ = Ref{$PetscReal}()
	gamma_ = Ref{$PetscReal}()
	beta_ = Ref{$PetscReal}()

    @chk ccall(
               (:TSAlpha2GetParams, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
               ts, alpha_m_, alpha_f_, gamma_, beta_,
              )

	alpha_m = alpha_m_[]
	alpha_f = alpha_f_[]
	gamma = gamma_[]
	beta = beta_[]

	return alpha_m,alpha_f,gamma,beta
end 

"""
	TSAlphaSetRadius(petsclib::PetscLibType,ts::TS, radius::PetscReal) 
sets the desired spectral radius of the method for `TSALPHA`
(i.e. high-frequency numerical damping)

Logically Collective

Input Parameters:
- `ts`     - timestepping context
- `radius` - the desired spectral radius

Options Database Key:
- `-ts_alpha_radius <radius>` - set alpha radius

Level: intermediate

-seealso: [](ch_ts), `TS`, `TSALPHA`, `TSAlphaSetParams()`, `TSAlphaGetParams()`

# External Links
$(_doc_external("Ts/TSAlphaSetRadius"))
"""
function TSAlphaSetRadius(petsclib::PetscLibType, ts::TS, radius::PetscReal) end

@for_petsc function TSAlphaSetRadius(petsclib::$UnionPetscLib, ts::TS, radius::$PetscReal )

    @chk ccall(
               (:TSAlphaSetRadius, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal),
               ts, radius,
              )


	return nothing
end 

"""
	TSAlphaSetParams(petsclib::PetscLibType,ts::TS, alpha_m::PetscReal, alpha_f::PetscReal, gamma::PetscReal) 
sets the algorithmic parameters for `TSALPHA`

Logically Collective

Input Parameters:
- `ts`      - timestepping context
- `alpha_m` - algorithmic parameter
- `alpha_f` - algorithmic parameter
- `gamma`   - algorithmic parameter

Options Database Keys:
- `-ts_alpha_alpha_m <alpha_m>` - set alpha_m
- `-ts_alpha_alpha_f <alpha_f>` - set alpha_f
- `-ts_alpha_gamma   <gamma>`   - set gamma

Level: advanced

-seealso: [](ch_ts), `TS`, `TSALPHA`, `TSAlphaSetRadius()`, `TSAlphaGetParams()`

# External Links
$(_doc_external("Ts/TSAlphaSetParams"))
"""
function TSAlphaSetParams(petsclib::PetscLibType, ts::TS, alpha_m::PetscReal, alpha_f::PetscReal, gamma::PetscReal) end

@for_petsc function TSAlphaSetParams(petsclib::$UnionPetscLib, ts::TS, alpha_m::$PetscReal, alpha_f::$PetscReal, gamma::$PetscReal )

    @chk ccall(
               (:TSAlphaSetParams, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal, $PetscReal, $PetscReal),
               ts, alpha_m, alpha_f, gamma,
              )


	return nothing
end 

"""
	alpha_m::PetscReal,alpha_f::PetscReal,gamma::PetscReal = TSAlphaGetParams(petsclib::PetscLibType,ts::TS) 
gets the algorithmic parameters for `TSALPHA`

Not Collective

Input Parameter:
- `ts` - timestepping context

Output Parameters:
- `alpha_m` - algorithmic parameter
- `alpha_f` - algorithmic parameter
- `gamma`   - algorithmic parameter

Level: advanced

-seealso: [](ch_ts), `TS`, `TSALPHA`, `TSAlphaSetRadius()`, `TSAlphaSetParams()`

# External Links
$(_doc_external("Ts/TSAlphaGetParams"))
"""
function TSAlphaGetParams(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSAlphaGetParams(petsclib::$UnionPetscLib, ts::TS )
	alpha_m_ = Ref{$PetscReal}()
	alpha_f_ = Ref{$PetscReal}()
	gamma_ = Ref{$PetscReal}()

    @chk ccall(
               (:TSAlphaGetParams, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
               ts, alpha_m_, alpha_f_, gamma_,
              )

	alpha_m = alpha_m_[]
	alpha_f = alpha_f_[]
	gamma = gamma_[]

	return alpha_m,alpha_f,gamma
end 

"""
	TSGLEERegisterAll(petsclib::PetscLibType) 
Registers all of the General Linear with Error Estimation methods in `TSGLEE`

Not Collective, but should be called by all processes which will need the schemes to be registered

Level: advanced

-seealso: [](ch_ts), `TSGLEERegisterDestroy()`

# External Links
$(_doc_external("Ts/TSGLEERegisterAll"))
"""
function TSGLEERegisterAll(petsclib::PetscLibType) end

@for_petsc function TSGLEERegisterAll(petsclib::$UnionPetscLib)

    @chk ccall(
               (:TSGLEERegisterAll, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	TSGLEERegisterDestroy(petsclib::PetscLibType) 
Frees the list of schemes that were registered by `TSGLEERegister()`.

Not Collective

Level: advanced

-seealso: [](ch_ts), `TSGLEERegister()`, `TSGLEERegisterAll()`

# External Links
$(_doc_external("Ts/TSGLEERegisterDestroy"))
"""
function TSGLEERegisterDestroy(petsclib::PetscLibType) end

@for_petsc function TSGLEERegisterDestroy(petsclib::$UnionPetscLib)

    @chk ccall(
               (:TSGLEERegisterDestroy, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	TSGLEEInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the `TSGLEE` package. It is called
from `TSInitializePackage()`.

Level: developer

-seealso: [](ch_ts), `PetscInitialize()`

# External Links
$(_doc_external("Ts/TSGLEEInitializePackage"))
"""
function TSGLEEInitializePackage(petsclib::PetscLibType) end

@for_petsc function TSGLEEInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:TSGLEEInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	TSGLEEFinalizePackage(petsclib::PetscLibType) 
This function destroys everything in the `TSGLEE` package. It is
called from `PetscFinalize()`.

Level: developer

-seealso: [](ch_ts), `PetscFinalize()`

# External Links
$(_doc_external("Ts/TSGLEEFinalizePackage"))
"""
function TSGLEEFinalizePackage(petsclib::PetscLibType) end

@for_petsc function TSGLEEFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:TSGLEEFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	TSGLEERegister(petsclib::PetscLibType,name::TSGLEEType, order::PetscInt, s::PetscInt, r::PetscInt, gamma::PetscReal, A::Vector{PetscReal}, B::Vector{PetscReal}, U::Vector{PetscReal}, V::Vector{PetscReal}, S::Vector{PetscReal}, F::Vector{PetscReal}, c::Vector{PetscReal}, Fembed::Vector{PetscReal}, Ferror::Vector{PetscReal}, Serror::Vector{PetscReal}, pinterp::PetscInt, binterp::Vector{PetscReal}) 
register a new `TSGLEE` scheme by providing the entries in the Butcher tableau

Not Collective, but the same schemes should be registered on all processes on which they will be used, No Fortran Support

Input Parameters:
- `name`    - identifier for method
- `order`   - order of method
- `s`       - number of stages
- `r`       - number of steps
- `gamma`   - LTE ratio
- `A`       - stage coefficients (dimension s*s, row-major)
- `B`       - step completion coefficients (dimension r*s, row-major)
- `U`       - method coefficients (dimension s*r, row-major)
- `V`       - method coefficients (dimension r*r, row-major)
- `S`       - starting coefficients
- `F`       - finishing coefficients
- `c`       - abscissa (dimension s; NULL to use row sums of A)
- `Fembed`  - step completion coefficients for embedded method
- `Ferror`  - error computation coefficients
- `Serror`  - error initialization coefficients
- `pinterp` - order of interpolation (0 if unavailable)
- `binterp` - array of interpolation coefficients (NULL if unavailable)

Level: advanced

-seealso: [](ch_ts), `TSGLEE`

# External Links
$(_doc_external("Ts/TSGLEERegister"))
"""
function TSGLEERegister(petsclib::PetscLibType, name::TSGLEEType, order::PetscInt, s::PetscInt, r::PetscInt, gamma::PetscReal, A::Vector{PetscReal}, B::Vector{PetscReal}, U::Vector{PetscReal}, V::Vector{PetscReal}, S::Vector{PetscReal}, F::Vector{PetscReal}, c::Vector{PetscReal}, Fembed::Vector{PetscReal}, Ferror::Vector{PetscReal}, Serror::Vector{PetscReal}, pinterp::PetscInt, binterp::Vector{PetscReal}) end

@for_petsc function TSGLEERegister(petsclib::$UnionPetscLib, name::TSGLEEType, order::$PetscInt, s::$PetscInt, r::$PetscInt, gamma::$PetscReal, A::Vector{$PetscReal}, B::Vector{$PetscReal}, U::Vector{$PetscReal}, V::Vector{$PetscReal}, S::Vector{$PetscReal}, F::Vector{$PetscReal}, c::Vector{$PetscReal}, Fembed::Vector{$PetscReal}, Ferror::Vector{$PetscReal}, Serror::Vector{$PetscReal}, pinterp::$PetscInt, binterp::Vector{$PetscReal} )

    @chk ccall(
               (:TSGLEERegister, $petsc_library),
               PetscErrorCode,
               (TSGLEEType, $PetscInt, $PetscInt, $PetscInt, $PetscReal, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, $PetscInt, Ptr{$PetscReal}),
               name, order, s, r, gamma, A, B, U, V, S, F, c, Fembed, Ferror, Serror, pinterp, binterp,
              )


	return nothing
end 

"""
	TSGLEESetType(petsclib::PetscLibType,ts::TS, gleetype::TSGLEEType) 
Set the type of `TSGLEE` scheme

Logically Collective

Input Parameters:
- `ts`       - timestepping context
- `gleetype` - type of `TSGLEE` scheme

Level: intermediate

-seealso: [](ch_ts), `TSGLEEGetType()`, `TSGLEE`

# External Links
$(_doc_external("Ts/TSGLEESetType"))
"""
function TSGLEESetType(petsclib::PetscLibType, ts::TS, gleetype::TSGLEEType) end

@for_petsc function TSGLEESetType(petsclib::$UnionPetscLib, ts::TS, gleetype::TSGLEEType )

    @chk ccall(
               (:TSGLEESetType, $petsc_library),
               PetscErrorCode,
               (CTS, TSGLEEType),
               ts, gleetype,
              )


	return nothing
end 

"""
	gleetype::TSGLEEType = TSGLEEGetType(petsclib::PetscLibType,ts::TS) 
Get the type of `TSGLEE` scheme

Logically Collective

Input Parameter:
- `ts` - timestepping context

Output Parameter:
- `gleetype` - type of `TSGLEE` scheme

Level: intermediate

-seealso: [](ch_ts), `TSGLEE`, `TSGLEESetType()`

# External Links
$(_doc_external("Ts/TSGLEEGetType"))
"""
function TSGLEEGetType(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSGLEEGetType(petsclib::$UnionPetscLib, ts::TS )
	gleetype_ = Ref{TSGLEEType}()

    @chk ccall(
               (:TSGLEEGetType, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{TSGLEEType}),
               ts, gleetype_,
              )

	gleetype = unsafe_string(gleetype_[])

	return gleetype
end 

"""
	TSEIMEXSetMaxRows(petsclib::PetscLibType,ts::TS, nrows::PetscInt) 
Set the maximum number of rows for `TSEIMEX` schemes

Logically Collective

Input Parameters:
- `ts`    - timestepping context
- `nrows` - maximum number of rows

Level: intermediate

-seealso: [](ch_ts), `TSEIMEXSetRowCol()`, `TSEIMEXSetOrdAdapt()`, `TSEIMEX`

# External Links
$(_doc_external("Ts/TSEIMEXSetMaxRows"))
"""
function TSEIMEXSetMaxRows(petsclib::PetscLibType, ts::TS, nrows::PetscInt) end

@for_petsc function TSEIMEXSetMaxRows(petsclib::$UnionPetscLib, ts::TS, nrows::$PetscInt )

    @chk ccall(
               (:TSEIMEXSetMaxRows, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt),
               ts, nrows,
              )


	return nothing
end 

"""
	TSEIMEXSetRowCol(petsclib::PetscLibType,ts::TS, row::PetscInt, col::PetscInt) 
Set the number of rows and the number of columns for the tableau that represents the T solution in the `TSEIMEX` scheme

Logically Collective

Input Parameters:
- `ts`  - timestepping context
- `row` - the row
- `col` - the column

Level: intermediate

-seealso: [](ch_ts), `TSEIMEXSetMaxRows()`, `TSEIMEXSetOrdAdapt()`, `TSEIMEX`

# External Links
$(_doc_external("Ts/TSEIMEXSetRowCol"))
"""
function TSEIMEXSetRowCol(petsclib::PetscLibType, ts::TS, row::PetscInt, col::PetscInt) end

@for_petsc function TSEIMEXSetRowCol(petsclib::$UnionPetscLib, ts::TS, row::$PetscInt, col::$PetscInt )

    @chk ccall(
               (:TSEIMEXSetRowCol, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt, $PetscInt),
               ts, row, col,
              )


	return nothing
end 

"""
	TSEIMEXSetOrdAdapt(petsclib::PetscLibType,ts::TS, flg::PetscBool) 
Set the order adaptativity for the `TSEIMEX` schemes

Logically Collective

Input Parameters:
- `ts`  - timestepping context
- `flg` - index in the T table

Level: intermediate

-seealso: [](ch_ts), `TSEIMEXSetRowCol()`, `TSEIMEX`

# External Links
$(_doc_external("Ts/TSEIMEXSetOrdAdapt"))
"""
function TSEIMEXSetOrdAdapt(petsclib::PetscLibType, ts::TS, flg::PetscBool) end

@for_petsc function TSEIMEXSetOrdAdapt(petsclib::$UnionPetscLib, ts::TS, flg::PetscBool )

    @chk ccall(
               (:TSEIMEXSetOrdAdapt, $petsc_library),
               PetscErrorCode,
               (CTS, PetscBool),
               ts, flg,
              )


	return nothing
end 

"""
	TSMonitorDMDARayDestroy(petsclib::PetscLibType,mctx::Cvoid) 

# External Links
$(_doc_external("Ts/TSMonitorDMDARayDestroy"))
"""
function TSMonitorDMDARayDestroy(petsclib::PetscLibType, mctx::Cvoid) end

@for_petsc function TSMonitorDMDARayDestroy(petsclib::$UnionPetscLib, mctx::Cvoid )

    @chk ccall(
               (:TSMonitorDMDARayDestroy, $petsc_library),
               PetscErrorCode,
               (Cvoid,),
               mctx,
              )


	return nothing
end 

"""
	TSMonitorDMDARay(petsclib::PetscLibType,ts::TS, steps::PetscInt, time::PetscReal, u::PetscVec, mctx::Cvoid) 

# External Links
$(_doc_external("Ts/TSMonitorDMDARay"))
"""
function TSMonitorDMDARay(petsclib::PetscLibType, ts::TS, steps::PetscInt, time::PetscReal, u::PetscVec, mctx::Cvoid) end

@for_petsc function TSMonitorDMDARay(petsclib::$UnionPetscLib, ts::TS, steps::$PetscInt, time::$PetscReal, u::PetscVec, mctx::Cvoid )

    @chk ccall(
               (:TSMonitorDMDARay, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt, $PetscReal, CVec, Ptr{Cvoid}),
               ts, steps, time, u, mctx,
              )


	return nothing
end 

"""
	TSMonitorLGDMDARay(petsclib::PetscLibType,ts::TS, step::PetscInt, ptime::PetscReal, u::PetscVec, ctx::Cvoid) 

# External Links
$(_doc_external("Ts/TSMonitorLGDMDARay"))
"""
function TSMonitorLGDMDARay(petsclib::PetscLibType, ts::TS, step::PetscInt, ptime::PetscReal, u::PetscVec, ctx::Cvoid) end

@for_petsc function TSMonitorLGDMDARay(petsclib::$UnionPetscLib, ts::TS, step::$PetscInt, ptime::$PetscReal, u::PetscVec, ctx::Cvoid )

    @chk ccall(
               (:TSMonitorLGDMDARay, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt, $PetscReal, CVec, Ptr{Cvoid}),
               ts, step, ptime, u, ctx,
              )


	return nothing
end 

"""
	TSSetPostEventStep(petsclib::PetscLibType,ts::TS, dt1::PetscReal) 
Set the first time step to use after the event

Logically Collective

Input Parameters:
- `ts`  - time integration context
- `dt1` - first post event step

Options Database Key:
- `-ts_event_post_event_step <dt1>` - first time step after the event

Level: advanced

-seealso: [](ch_ts), `TS`, `TSEvent`, `TSSetEventHandler()`, `TSSetPostEventSecondStep()`

# External Links
$(_doc_external("Ts/TSSetPostEventStep"))
"""
function TSSetPostEventStep(petsclib::PetscLibType, ts::TS, dt1::PetscReal) end

@for_petsc function TSSetPostEventStep(petsclib::$UnionPetscLib, ts::TS, dt1::$PetscReal )

    @chk ccall(
               (:TSSetPostEventStep, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal),
               ts, dt1,
              )


	return nothing
end 

"""
	TSSetPostEventSecondStep(petsclib::PetscLibType,ts::TS, dt2::PetscReal) 
Set the second time step to use after the event

Logically Collective

Input Parameters:
- `ts`  - time integration context
- `dt2` - second post event step

Options Database Key:
- `-ts_event_post_event_second_step <dt2>` - second time step after the event

Level: advanced

-seealso: [](ch_ts), `TS`, `TSEvent`, `TSSetEventHandler()`, `TSSetPostEventStep()`

# External Links
$(_doc_external("Ts/TSSetPostEventSecondStep"))
"""
function TSSetPostEventSecondStep(petsclib::PetscLibType, ts::TS, dt2::PetscReal) end

@for_petsc function TSSetPostEventSecondStep(petsclib::$UnionPetscLib, ts::TS, dt2::$PetscReal )

    @chk ccall(
               (:TSSetPostEventSecondStep, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal),
               ts, dt2,
              )


	return nothing
end 

"""
	TSSetEventTolerances(petsclib::PetscLibType,ts::TS, tol::PetscReal, vtol::Vector{PetscReal}) 
Set tolerances for event (indicator function) zero crossings

Logically Collective

Input Parameters:
- `ts`   - time integration context
- `tol`  - tolerance, `PETSC_CURRENT` to leave the current value
- `vtol` - array of tolerances or `NULL`, used in preference to `tol` if present

Options Database Key:
- `-ts_event_tol <tol>` - tolerance for event (indicator function) zero crossing

Level: beginner

-seealso: [](ch_ts), `TS`, `TSEvent`, `TSSetEventHandler()`

# External Links
$(_doc_external("Ts/TSSetEventTolerances"))
"""
function TSSetEventTolerances(petsclib::PetscLibType, ts::TS, tol::PetscReal, vtol::Vector{PetscReal}) end

@for_petsc function TSSetEventTolerances(petsclib::$UnionPetscLib, ts::TS, tol::$PetscReal, vtol::Vector{$PetscReal} )

    @chk ccall(
               (:TSSetEventTolerances, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscReal, Ptr{$PetscReal}),
               ts, tol, vtol,
              )


	return nothing
end 

"""
	TSSetEventHandler(petsclib::PetscLibType,ts::TS, nevents::PetscInt, direction::Vector{PetscInt}, terminate::Vector{PetscBool}, indicator::external, postevent::external, ctx::Cvoid) 
Sets functions and parameters used for indicating events and handling them

Logically Collective

Input Parameters:
- `ts`        - the `TS` context obtained from `TSCreate()`
- `nevents`   - number of local events (i.e. managed by the given MPI process)
- `direction` - direction of zero crossing to be detected (one for each local event).
`-1` => zero crossing in negative direction,
`+1` => zero crossing in positive direction, `0` => both ways
- `terminate` - flag to indicate whether time stepping should be terminated after
an event is detected (one for each local event)
- `indicator` - callback defininig the user indicator functions whose sign changes (see `direction`) mark presence of the events
- `postevent` - [optional] user post-event callback; it can change the solution, ODE etc at the time of the event
- `ctx`       - [optional] user-defined context for private data for the
`indicator()` and `postevent()` routines (use `NULL` if no
context is desired)

Calling sequence of `indicator`:
- `ts`     - the `TS` context
- `t`      - current time
- `U`      - current solution
- `fvalue` - output array with values of local indicator functions (length == `nevents`) for time t and state-vector U
- `ctx`    - the context passed as the final argument to `TSSetEventHandler()`

Calling sequence of `postevent`:
- `ts`           - the `TS` context
- `nevents_zero` - number of triggered local events (whose indicator function is marked as crossing zero, and direction is appropriate)
- `events_zero`  - indices of the triggered local events
- `t`            - current time
- `U`            - current solution
- `forwardsolve` - flag to indicate whether `TS` is doing a forward solve (`PETSC_TRUE`) or adjoint solve (`PETSC_FALSE`)
- `ctx`          - the context passed as the final argument to `TSSetEventHandler()`

Options Database Keys:
- `-ts_event_tol <tol>`                       - tolerance for zero crossing check of indicator functions
- `-ts_event_monitor`                         - print choices made by event handler
- `-ts_event_recorder_initial_size <recsize>` - initial size of event recorder
- `-ts_event_post_event_step <dt1>`           - first time step after event
- `-ts_event_post_event_second_step <dt2>`    - second time step after event
- `-ts_event_dt_min <dt>`                     - minimum time step considered for TSEvent

Level: intermediate

-seealso: [](ch_ts), `TSEvent`, `TSCreate()`, `TSSetTimeStep()`, `TSSetConvergedReason()`

# External Links
$(_doc_external("Ts/TSSetEventHandler"))
"""
function TSSetEventHandler(petsclib::PetscLibType, ts::TS, nevents::PetscInt, direction::Vector{PetscInt}, terminate::Vector{PetscBool}, indicator::external, postevent::external, ctx::Cvoid) end

@for_petsc function TSSetEventHandler(petsclib::$UnionPetscLib, ts::TS, nevents::$PetscInt, direction::Vector{$PetscInt}, terminate::Vector{PetscBool}, indicator::external, postevent::external, ctx::Cvoid )

    @chk ccall(
               (:TSSetEventHandler, $petsc_library),
               PetscErrorCode,
               (CTS, $PetscInt, Ptr{$PetscInt}, Ptr{PetscBool}, external, external, Ptr{Cvoid}),
               ts, nevents, direction, terminate, indicator, postevent, ctx,
              )


	return nothing
end 

"""
	nevents::PetscInt = TSGetNumEvents(petsclib::PetscLibType,ts::TS) 
Get the number of events defined on the given MPI process

Logically Collective

Input Parameter:
- `ts` - the `TS` context

Output Parameter:
- `nevents` - the number of local events on each MPI process

Level: intermediate

-seealso: [](ch_ts), `TSEvent`, `TSSetEventHandler()`

# External Links
$(_doc_external("Ts/TSGetNumEvents"))
"""
function TSGetNumEvents(petsclib::PetscLibType, ts::TS) end

@for_petsc function TSGetNumEvents(petsclib::$UnionPetscLib, ts::TS )
	nevents_ = Ref{$PetscInt}()

    @chk ccall(
               (:TSGetNumEvents, $petsc_library),
               PetscErrorCode,
               (CTS, Ptr{$PetscInt}),
               ts, nevents_,
              )

	nevents = nevents_[]

	return nevents
end 

