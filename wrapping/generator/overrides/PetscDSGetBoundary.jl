# override for PetscDSGetBoundary; C signature: PetscDSGetBoundary(PetscDS ds, PetscInt bd, PetscWeakForm* wf, DMBoundaryConditionType* type, char* name[], DMLabel* label, PetscInt* Nv, PetscInt* values[], PetscInt* field, PetscInt* Nc, PetscInt* comps[], PetscVoidFn** func, PetscVoidFn** func_t, void** ctx)
"""
	Nv::PetscInt,values::Vector{PetscInt},field::PetscInt,Nc::PetscInt,comps::Vector{PetscInt} = PetscDSGetBoundary(petsclib::PetscLibType,ds::PetscDS, bd::PetscInt, wf::PetscWeakForm, type::DMBoundaryConditionType, name::String, label::DMLabel, func::PetscVoidFn, func_t::PetscVoidFn, ctx::Cvoid) 
Gets a boundary condition from the model

Input Parameters:
- `ds` - The `PetscDS` object
- `bd` - The boundary condition number

Output Parameters:
- `wf`     - The `PetscWeakForm` holding the pointwise functions
- `type`   - The type of condition, e.g. `DM_BC_ESSENTIAL`/`DM_BC_ESSENTIAL_FIELD` (Dirichlet), or `DM_BC_NATURAL` (Neumann)
- `name`   - The boundary condition name
- `label`  - The label defining constrained points
- `Nv`     - The number of `DMLabel` ids for constrained points
- `values` - An array of ids for constrained points
- `field`  - The field to constrain
- `Nc`     - The number of constrained field components
- `comps`  - An array of constrained component numbers
- `func`   - A pointwise function giving boundary values
- `func_t` - A pointwise function giving the time derivative of the boundary values
- `ctx`    - An optional user context for `bcFunc`

Options Database Keys:
- `-bc_<boundary name> <num>`      - Overrides the boundary ids
- `-bc_<boundary name>_comp <num>` - Overrides the boundary components

Level: developer

-seealso: `PetscDS`, `PetscWeakForm`, `DMBoundaryConditionType`, `PetscDSAddBoundary()`, `DMLabel`

# External Links
$(_doc_external("Dm/PetscDSGetBoundary"))
"""
function PetscDSGetBoundary(petsclib::PetscLibType, ds::PetscDS, bd::PetscInt) end

@for_petsc function PetscDSGetBoundary(petsclib::$UnionPetscLib, ds::PetscDS, bd::$PetscInt)
	wf_ref    = Ref{PetscWeakForm}(C_NULL)
	type_ref  = Ref{DMBoundaryConditionType}()
	name_ref  = Ref{Ptr{Cchar}}(C_NULL)
	label_ref = Ref{DMLabel}(C_NULL)
	Nv_ref    = Ref{$PetscInt}(0)
	values_ref = Ref{Ptr{$PetscInt}}(C_NULL)
	field_ref = Ref{$PetscInt}(0)
	Nc_ref    = Ref{$PetscInt}(0)
	comps_ref = Ref{Ptr{$PetscInt}}(C_NULL)
	func_ref  = Ref{Ptr{Cvoid}}(C_NULL)
	func_t_ref = Ref{Ptr{Cvoid}}(C_NULL)
	ctx_ref   = Ref{Ptr{Cvoid}}(C_NULL)

    @chk ccall(
               (:PetscDSGetBoundary, $petsc_library),
               PetscErrorCode,
               (PetscDS, $PetscInt,
                Ptr{PetscWeakForm}, Ptr{DMBoundaryConditionType}, Ptr{Ptr{Cchar}}, Ptr{DMLabel},
                Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}}, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{$PetscInt}},
                Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}, Ptr{Ptr{Cvoid}}),
               ds, bd,
               wf_ref, type_ref, name_ref, label_ref,
               Nv_ref, values_ref, field_ref, Nc_ref, comps_ref,
               func_ref, func_t_ref, ctx_ref,
              )

	Nv = Nv_ref[]
	values = Nv > 0 ? unsafe_wrap(Array, values_ref[], Nv; own = false) : $PetscInt[]
	Nc = Nc_ref[]
	comps = Nc > 0 ? unsafe_wrap(Array, comps_ref[], Nc; own = false) : $PetscInt[]

	return (wf = wf_ref[], type = type_ref[], name = name_ref[] == C_NULL ? "" : unsafe_string(name_ref[]),
	        label = label_ref[], Nv = Nv, values = values, field = field_ref[],
	        Nc = Nc, comps = comps, func = func_ref[], func_t = func_t_ref[], ctx = ctx_ref[])
end

