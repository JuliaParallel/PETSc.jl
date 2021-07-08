const CDM = Ptr{Cvoid}

abstract type AbstractDM{PetscLib} end
mutable struct DM{PetscLib} <: AbstractDM{PetscLib}
    ptr::CDM
    opts::Options{PetscLib}
    DM{PetscLib}(ptr, opts = Options(PetscLib)) where {PetscLib} =
        new{PetscLib}(ptr, opts)
end

"""
    DMSetUp!(da::DM)

# External Links
$(_doc_external("DM/DMSetUp"))
"""
function DMSetUp! end

@for_petsc function DMSetUp!(da::DM{$PetscLib})
    with(da.opts) do
        @chk ccall(
            (:DMSetFromOptions, $petsc_library),
            PetscErrorCode,
            (CDM,),
            da,
        )

        @chk ccall((:DMSetUp, $petsc_library), PetscErrorCode, (CDM,), da)
    end
end

@for_petsc begin
    function destroy(da::DM{$PetscLib})
        finalized($PetscScalar) || @chk ccall(
            (:DMDestroy, $petsc_library),
            PetscErrorCode,
            (Ptr{CDM},),
            da,
        )
        da.ptr = C_NULL
        return nothing
    end
end

"""
    DMGetDimension(dm::DM)

Return the topological dimension of the `dm`

# External Links
$(_doc_external("DM/DMGetDimension"))
"""
function DMGetDimension end

@for_petsc function DMGetDimension(dm::DM{$PetscLib})
    dim = Ref{$PetscInt}()

    @chk ccall(
        (:DMGetDimension, $petsc_library),
        PetscErrorCode,
        (CDM, Ptr{$PetscInt}),
        dm,
        dim,
    )

    return dim[]
end

"""
    gettype(dm::DM)

Gets type name of the `dm`

# External Links
$(_doc_external("DM/DMGetType"))
"""
function gettype(::DM) end

@for_petsc function gettype(dm::DM{$PetscLib})
    t_r = Ref{Cstring}()
    @chk ccall(
        (:DMGetType, $petsc_library),
        PetscErrorCode,
        (CDM, Ptr{Cstring}),
        dm,
        t_r,
    )
    return unsafe_string(t_r[])
end

"""
    view(dm::DM, viewer::Viewer=ViewerStdout(petsclib, getcomm(dm)))

view a `dm` with `viewer`

# External Links
$(_doc_external("DM/DMView"))
"""
function view(::DM) end

@for_petsc function view(
    dm::DM{$PetscLib},
    viewer::AbstractViewer{$PetscLib} = ViewerStdout($petsclib, getcomm(dm)),
)
    @chk ccall(
        (:DMView, $petsc_library),
        PetscErrorCode,
        (CDM, CPetscViewer),
        dm,
        viewer,
    )
    return nothing
end
