# Contributing

Contributions are highly welcome, in particular since only part of the PETSc functionality is currently being tested and has high-level interfaces. 

You can thus help in many ways:
1) Add more examples
2) Add new tests
3) Update documentation
4) Report bugs
5) Add a high-level interface for unsupported features
6) Keep the routines up to date with future PETSc versions
7) Keep the precompiled binaries in [PETSc_jll](https://github.com/JuliaBinaryWrappers/PETSc_jll.jl) up to date.


#### Autowrappers
The low-level `LibPETSc` wrappers are generated, once per PETSc release, by the maintainer-only generator in `wrapping/generator/`; see `wrapping/WRAPPING.md` for how it works and what to watch for. Never edit files in `src/autowrapped/` by hand.

Note, however, that a range of additional changes were necessary and we thus had manually fix a number of things. It is therefore *not* recommended to rerun these autowrappers for newer versions of PETSc. 
Since there are usually only a limited number of new or updated functions between PETSc releases, it is recommended to run the a wrapper only for these new functions and replace those affected accordingly.   

Make sure that the tests work!

#### Running the tests
The tests live in `test/`, one folder per group: `core` (initialization, handles, lifetimes, the API register), `vecmat`, `solvers`, `dm`, `lowlevel` (the generated `LibPETSc` wrappers), `regression`, `examples` (the scripts in `examples/` and the manual's examples) and `mpi` (files run on 4 ranks). A new test file goes in its group's folder and is listed in that group's block in `test/runtests.jl`.

The whole suite takes about 15 minutes. To run only some groups:

```julia
using Pkg
Pkg.test("PETSc"; test_args = ["dm", "mpi"])
```

#### Adding new functionality
Please open a pull request to add any of the above contributions.

New high-level functions must follow the [naming conventions](naming.md). The low-level `LibPETSc` layer is exempt: it keeps the PETSc C names verbatim.