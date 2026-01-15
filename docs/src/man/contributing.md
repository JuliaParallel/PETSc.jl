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
We originally used the function `/wrapping/generatejuliabindings.jl` to wrap the whole PETSc library, which borrows python routines by Barry Smith. The PETSc version for which we generated these original wrappers was 3.23.6.

Note, however, that a range of additional changes were necessary and we thus had manually fix a number of things. It is therefore *not* recommended to rerun these autowrappers for newer versions of PETSc. 
Since there are usually only a limited number of new or updated functions between PETSc releases, it is recommended to run the a wrapper only for these new functions and replace those affected accordingly.   

Make sure that the tests work!

#### Adding new functionality
Please open a pull request to add any of the above contributions.