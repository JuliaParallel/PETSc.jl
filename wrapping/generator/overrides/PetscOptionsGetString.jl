# override for PetscOptionsGetString; C signature: PetscOptionsGetString(PetscOptions options, char pre[], char name[], char string[], size_t len, PetscBool* set)
"""
	string::Union{Bool,String} = PetscOptionsGetString(petsclib::PetscLibType,options::AbstractPetscOptions, pre::String, name::String, string::String, len::Csize_t) 
Gets the string value for a particular option in
the database.

Not Collective

Input Parameters:
- `options` - options database, use `NULL` for default global database
- `pre`     - string to prepend to name or `NULL`
- `name`    - the option one is seeking
- `len`     - maximum length of the string including null termination

Output Parameters:
- `string` - returns the value of the parameter ifn set, otherwise `false`

Level: beginner

-seealso: `PetscOptionsGetInt()`, `PetscOptionsGetReal()`,
`PetscOptionsHasName()`, `PetscOptionsGetIntArray()`, `PetscOptionsGetRealArray()`, `PetscOptionsBool()`,
`PetscOptionsName()`, `PetscOptionsBegin()`, `PetscOptionsEnd()`, `PetscOptionsHeadBegin()`,
`PetscOptionsStringArray()`, `PetscOptionsRealArray()`, `PetscOptionsScalar()`,
`PetscOptionsBoolGroupBegin()`, `PetscOptionsBoolGroup()`, `PetscOptionsBoolGroupEnd()`,
`PetscOptionsFList()`, `PetscOptionsEList()`

# External Links
$(_doc_external("Sys/PetscOptionsGetString"))
"""
function PetscOptionsGetString(petsclib::PetscLibType, options::AbstractPetscOptions, pre::Union{Ptr, String}, name::String) end

@for_petsc function PetscOptionsGetString(petsclib::$UnionPetscLib, options::AbstractPetscOptions, pre::Union{Ptr,String}, name::String)
	set_ = Ref{PetscBool}()
    val = Vector{UInt8}(undef, 256)

    @chk ccall(
               (:PetscOptionsGetString, $petsc_library),
               PetscErrorCode,
               (COptions, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}, Csize_t, Ptr{PetscBool}),
               options, pre, name, val, sizeof(val), set_,
              )

	set = set_[]
    if set
        val = GC.@preserve val unsafe_string(pointer(val))
    else
        val = false
    end
  

	return val
end

