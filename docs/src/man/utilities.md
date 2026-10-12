# Utilities

This page documents utility functions for initialization, options handling, system information, and code auditing.

## Initialization

The initialization functions manage PETSc's lifecycle, including MPI initialization and cleanup.

```@autodocs
Modules = [PETSc]
Pages   = ["init.jl"]
```

## Destroying objects together

On more than one process nothing frees a PETSc object but `destroy!`, called on every rank in the same order. `with_destroy` does that for every object registered with `track!` when its block ends.

```@autodocs
Modules = [PETSc]
Pages   = ["with_destroy.jl"]
```

## Options

PETSc uses an options database to configure solvers and other objects at runtime. These functions provide access to the options system.

```@autodocs
Modules = [PETSc]
Pages   = ["options.jl"]
```

## System Utilities

General system-level utilities for working with PETSc objects.

```@autodocs
Modules = [PETSc]
Pages   = ["sys.jl"]
```

## Code Auditing

The audit utilities help identify potential memory leaks by tracking PETSc object creation and destruction.

```@autodocs
Modules = [PETSc]
Pages   = ["audit.jl", "audit_names.jl"]
```
