mutable struct _p_PetscSectionSym end
const PetscSectionSym = Ptr{_p_PetscSectionSym}


"""
	PetscSectionCreate(petsclib::PetscLibType,comm::MPI_Comm, s::PetscSection) 
Allocates a `PetscSection` and sets the map contents to the default.

Collective

Input Parameters:
- `comm` - the MPI communicator
- `s`    - pointer to the section

Level: beginner

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionSetChart()`, `PetscSectionDestroy()`, `PetscSectionCreateGlobalSection()`

# External Links
$(_doc_external("Vec/PetscSectionCreate"))
"""
function PetscSectionCreate(petsclib::PetscLibType, comm::MPI_Comm, s::PetscSection) end

# Convenience method: accept a Ref so callers can pass `Ref{PetscSection}()` directly
@for_petsc function PetscSectionCreate(petsclib::$UnionPetscLib, comm::MPI_Comm, s::Ref{PetscSection} )

    @chk ccall(
               (:PetscSectionCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PetscSection}),
               comm, s,
              )


	return nothing
end

@for_petsc function PetscSectionCreate(petsclib::$UnionPetscLib, comm::MPI_Comm, s::PetscSection )

    @chk ccall(
               (:PetscSectionCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PetscSection}),
               comm, s,
              )


	return nothing
end 

"""
	PetscSectionCopy(petsclib::PetscLibType,section::PetscSection, newSection::PetscSection) 
Creates a shallow (if possible) copy of the `PetscSection`

Collective

Input Parameter:
- `section` - the `PetscSection`

Output Parameter:
- `newSection` - the copy

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionCreate()`, `PetscSectionDestroy()`

# External Links
$(_doc_external("Vec/PetscSectionCopy"))
"""
function PetscSectionCopy(petsclib::PetscLibType, section::PetscSection, newSection::PetscSection) end

@for_petsc function PetscSectionCopy(petsclib::$UnionPetscLib, section::PetscSection, newSection::PetscSection )

    @chk ccall(
               (:PetscSectionCopy, $petsc_library),
               PetscErrorCode,
               (PetscSection, PetscSection),
               section, newSection,
              )


	return nothing
end 

"""
	PetscSectionClone(petsclib::PetscLibType,section::PetscSection, newSection::PetscSection) 
Creates a shallow (if possible) copy of the `PetscSection`

Collective

Input Parameter:
- `section` - the `PetscSection`

Output Parameter:
- `newSection` - the copy

Level: beginner

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionCreate()`, `PetscSectionDestroy()`, `PetscSectionCopy()`

# External Links
$(_doc_external("Vec/PetscSectionClone"))
"""
function PetscSectionClone(petsclib::PetscLibType, section::PetscSection, newSection::PetscSection) end

@for_petsc function PetscSectionClone(petsclib::$UnionPetscLib, section::PetscSection, newSection::PetscSection )

    @chk ccall(
               (:PetscSectionClone, $petsc_library),
               PetscErrorCode,
               (PetscSection, Ptr{PetscSection}),
               section, newSection,
              )


	return nothing
end 

"""
	PetscSectionSetFromOptions(petsclib::PetscLibType,s::PetscSection) 
sets parameters in a `PetscSection` from the options database

Collective

Input Parameter:
- `s` - the `PetscSection`

Options Database Key:
- `-petscsection_point_major` - `PETSC_TRUE` for point-major order

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionCreate()`, `PetscSectionDestroy()`

# External Links
$(_doc_external("Vec/PetscSectionSetFromOptions"))
"""
function PetscSectionSetFromOptions(petsclib::PetscLibType, s::PetscSection) end

@for_petsc function PetscSectionSetFromOptions(petsclib::$UnionPetscLib, s::PetscSection )

    @chk ccall(
               (:PetscSectionSetFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscSection,),
               s,
              )


	return nothing
end 

"""
	congruent::PetscBool = PetscSectionCompare(petsclib::PetscLibType,s1::PetscSection, s2::PetscSection) 
Compares two sections

Collective

Input Parameters:
- `s1` - the first `PetscSection`
- `s2` - the second `PetscSection`

Output Parameter:
- `congruent` - `PETSC_TRUE` if the two sections are congruent, `PETSC_FALSE` otherwise

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionCreate()`, `PetscSectionCopy()`, `PetscSectionClone()`

# External Links
$(_doc_external("Vec/PetscSectionCompare"))
"""
function PetscSectionCompare(petsclib::PetscLibType, s1::PetscSection, s2::PetscSection) end

@for_petsc function PetscSectionCompare(petsclib::$UnionPetscLib, s1::PetscSection, s2::PetscSection )
	congruent_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscSectionCompare, $petsc_library),
               PetscErrorCode,
               (PetscSection, PetscSection, Ptr{PetscBool}),
               s1, s2, congruent_,
              )

	congruent = congruent_[]

	return congruent
end 

"""
	numFields::PetscInt = PetscSectionGetNumFields(petsclib::PetscLibType,s::PetscSection) 
Returns the number of fields in a `PetscSection`, or 0 if no fields were defined.

Not Collective

Input Parameter:
- `s` - the `PetscSection`

Output Parameter:
- `numFields` - the number of fields defined, or 0 if none were defined

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionSetNumFields()`

# External Links
$(_doc_external("Vec/PetscSectionGetNumFields"))
"""
function PetscSectionGetNumFields(petsclib::PetscLibType, s::PetscSection) end

@for_petsc function PetscSectionGetNumFields(petsclib::$UnionPetscLib, s::PetscSection )
	numFields_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscSectionGetNumFields, $petsc_library),
               PetscErrorCode,
               (PetscSection, Ptr{$PetscInt}),
               s, numFields_,
              )

	numFields = numFields_[]

	return numFields
end 

"""
	PetscSectionSetNumFields(petsclib::PetscLibType,s::PetscSection, numFields::PetscInt) 
Sets the number of fields in a `PetscSection`

Not Collective

Input Parameters:
- `s`         - the `PetscSection`
- `numFields` - the number of fields

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionGetNumFields()`, `PetscSectionSetChart()`, `PetscSectionReset()`

# External Links
$(_doc_external("Vec/PetscSectionSetNumFields"))
"""
function PetscSectionSetNumFields(petsclib::PetscLibType, s::PetscSection, numFields::PetscInt) end

@for_petsc function PetscSectionSetNumFields(petsclib::$UnionPetscLib, s::PetscSection, numFields::$PetscInt )

    @chk ccall(
               (:PetscSectionSetNumFields, $petsc_library),
               PetscErrorCode,
               (PetscSection, $PetscInt),
               s, numFields,
              )


	return nothing
end 

"""
	PetscSectionGetFieldName(petsclib::PetscLibType,s::PetscSection, field::PetscInt, fieldName::String) 
Returns the name of a field in the `PetscSection`

Not Collective

Input Parameters:
- `s`     - the `PetscSection`
- `field` - the field number

Output Parameter:
- `fieldName` - the field name

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionSetFieldName()`, `PetscSectionSetNumFields()`, `PetscSectionGetNumFields()`

# External Links
$(_doc_external("Vec/PetscSectionGetFieldName"))
"""
function PetscSectionGetFieldName(petsclib::PetscLibType, s::PetscSection, field::PetscInt, fieldName::String) end

@for_petsc function PetscSectionGetFieldName(petsclib::$UnionPetscLib, s::PetscSection, field::$PetscInt, fieldName::String )
	fieldName_ = Ref(pointer(fieldName))

    @chk ccall(
               (:PetscSectionGetFieldName, $petsc_library),
               PetscErrorCode,
               (PetscSection, $PetscInt, Ptr{Ptr{Cchar}}),
               s, field, fieldName_,
              )


	return nothing
end 

"""
	PetscSectionSetFieldName(petsclib::PetscLibType,s::PetscSection, field::PetscInt, fieldName::String) 
Sets the name of a field in the `PetscSection`

Not Collective

Input Parameters:
- `s`         - the `PetscSection`
- `field`     - the field number
- `fieldName` - the field name

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSectionGetFieldName()`, `PetscSectionSetNumFields()`, `PetscSectionGetNumFields()`

# External Links
$(_doc_external("Vec/PetscSectionSetFieldName"))
"""
function PetscSectionSetFieldName(petsclib::PetscLibType, s::PetscSection, field::PetscInt, fieldName::String) end

@for_petsc function PetscSectionSetFieldName(petsclib::$UnionPetscLib, s::PetscSection, field::$PetscInt, fieldName::String )

    @chk ccall(
               (:PetscSectionSetFieldName, $petsc_library),
               PetscErrorCode,
               (PetscSection, $PetscInt, Ptr{Cchar}),
               s, field, fieldName,
              )


	return nothing
end 

"""
	PetscSectionGetComponentName(petsclib::PetscLibType,s::PetscSection, field::PetscInt, comp::PetscInt, compName::String) 
Gets the name of a field component in the `PetscSection`

Not Collective

Input Parameters:
- `s`     - the `PetscSection`
- `field` - the field number
- `comp`  - the component number

Output Parameter:
- `compName` - the component name

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionGetFieldName()`, `PetscSectionSetNumFields()`, `PetscSectionGetNumFields()`,
`PetscSectionSetComponentName()`, `PetscSectionSetFieldName()`, `PetscSectionGetFieldComponents()`, `PetscSectionSetFieldComponents()`

# External Links
$(_doc_external("Vec/PetscSectionGetComponentName"))
"""
function PetscSectionGetComponentName(petsclib::PetscLibType, s::PetscSection, field::PetscInt, comp::PetscInt, compName::String) end

@for_petsc function PetscSectionGetComponentName(petsclib::$UnionPetscLib, s::PetscSection, field::$PetscInt, comp::$PetscInt, compName::String )
	compName_ = Ref(pointer(compName))

    @chk ccall(
               (:PetscSectionGetComponentName, $petsc_library),
               PetscErrorCode,
               (PetscSection, $PetscInt, $PetscInt, Ptr{Ptr{Cchar}}),
               s, field, comp, compName_,
              )


	return nothing
end 

"""
	PetscSectionSetComponentName(petsclib::PetscLibType,s::PetscSection, field::PetscInt, comp::PetscInt, compName::String) 
Sets the name of a field component in the `PetscSection`

Not Collective

Input Parameters:
- `s`        - the `PetscSection`
- `field`    - the field number
- `comp`     - the component number
- `compName` - the component name

Level: advanced

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionGetComponentName()`, `PetscSectionSetNumFields()`, `PetscSectionGetNumFields()`,
`PetscSectionSetFieldName()`, `PetscSectionGetFieldComponents()`, `PetscSectionSetFieldComponents()`

# External Links
$(_doc_external("Vec/PetscSectionSetComponentName"))
"""
function PetscSectionSetComponentName(petsclib::PetscLibType, s::PetscSection, field::PetscInt, comp::PetscInt, compName::String) end

@for_petsc function PetscSectionSetComponentName(petsclib::$UnionPetscLib, s::PetscSection, field::$PetscInt, comp::$PetscInt, compName::String )

    @chk ccall(
               (:PetscSectionSetComponentName, $petsc_library),
               PetscErrorCode,
               (PetscSection, $PetscInt, $PetscInt, Ptr{Cchar}),
               s, field, comp, compName,
              )


	return nothing
end 

"""
	numComp::PetscInt = PetscSectionGetFieldComponents(petsclib::PetscLibType,s::PetscSection, field::PetscInt) 
Returns the number of field components for the given field.

Not Collective

Input Parameters:
- `s`     - the `PetscSection`
- `field` - the field number

Output Parameter:
- `numComp` - the number of field components

Level: advanced

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionSetFieldComponents()`, `PetscSectionGetNumFields()`,
`PetscSectionSetComponentName()`, `PetscSectionGetComponentName()`

# External Links
$(_doc_external("Vec/PetscSectionGetFieldComponents"))
"""
function PetscSectionGetFieldComponents(petsclib::PetscLibType, s::PetscSection, field::PetscInt) end

@for_petsc function PetscSectionGetFieldComponents(petsclib::$UnionPetscLib, s::PetscSection, field::$PetscInt )
	numComp_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscSectionGetFieldComponents, $petsc_library),
               PetscErrorCode,
               (PetscSection, $PetscInt, Ptr{$PetscInt}),
               s, field, numComp_,
              )

	numComp = numComp_[]

	return numComp
end 

"""
	PetscSectionSetFieldComponents(petsclib::PetscLibType,s::PetscSection, field::PetscInt, numComp::PetscInt) 
Sets the number of field components for the given field.

Not Collective

Input Parameters:
- `s`       - the `PetscSection`
- `field`   - the field number
- `numComp` - the number of field components

Level: advanced

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionGetFieldComponents()`, `PetscSectionSetComponentName()`,
`PetscSectionGetComponentName()`, `PetscSectionGetNumFields()`

# External Links
$(_doc_external("Vec/PetscSectionSetFieldComponents"))
"""
function PetscSectionSetFieldComponents(petsclib::PetscLibType, s::PetscSection, field::PetscInt, numComp::PetscInt) end

@for_petsc function PetscSectionSetFieldComponents(petsclib::$UnionPetscLib, s::PetscSection, field::$PetscInt, numComp::$PetscInt )

    @chk ccall(
               (:PetscSectionSetFieldComponents, $petsc_library),
               PetscErrorCode,
               (PetscSection, $PetscInt, $PetscInt),
               s, field, numComp,
              )


	return nothing
end 

"""
	pStart::PetscInt,pEnd::PetscInt = PetscSectionGetChart(petsclib::PetscLibType,s::PetscSection) 
Returns the range [`pStart`, `pEnd`) in which points (indices) lie for this `PetscSection` on this MPI process

Not Collective

Input Parameter:
- `s` - the `PetscSection`

Output Parameters:
- `pStart` - the first point
- `pEnd`   - one past the last point

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionSetChart()`, `PetscSectionCreate()`

# External Links
$(_doc_external("Vec/PetscSectionGetChart"))
"""
function PetscSectionGetChart(petsclib::PetscLibType, s::PetscSection) end

@for_petsc function PetscSectionGetChart(petsclib::$UnionPetscLib, s::PetscSection )
	pStart_ = Ref{$PetscInt}()
	pEnd_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscSectionGetChart, $petsc_library),
               PetscErrorCode,
               (PetscSection, Ptr{$PetscInt}, Ptr{$PetscInt}),
               s, pStart_, pEnd_,
              )

	pStart = pStart_[]
	pEnd = pEnd_[]

	return pStart,pEnd
end 

"""
	PetscSectionSetChart(petsclib::PetscLibType,s::PetscSection, pStart::PetscInt, pEnd::PetscInt) 
Sets the range [`pStart`, `pEnd`) in which points (indices) lie for this `PetscSection` on this MPI process

Not Collective

Input Parameters:
- `s`      - the `PetscSection`
- `pStart` - the first point
- `pEnd`   - one past the last point, `pStart`  ≤  `pEnd`

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionGetChart()`, `PetscSectionCreate()`, `PetscSectionSetNumFields()`

# External Links
$(_doc_external("Vec/PetscSectionSetChart"))
"""
function PetscSectionSetChart(petsclib::PetscLibType, s::PetscSection, pStart::PetscInt, pEnd::PetscInt) end

@for_petsc function PetscSectionSetChart(petsclib::$UnionPetscLib, s::PetscSection, pStart::$PetscInt, pEnd::$PetscInt )

    @chk ccall(
               (:PetscSectionSetChart, $petsc_library),
               PetscErrorCode,
               (PetscSection, $PetscInt, $PetscInt),
               s, pStart, pEnd,
              )


	return nothing
end 

"""
	PetscSectionGetPermutation(petsclib::PetscLibType,s::PetscSection, perm::IS) 
Returns the permutation of [0, `pEnd`

Not Collective

Input Parameter:
- `s` - the `PetscSection`

Output Parameter:
- `perm` - The permutation as an `IS`

Level: intermediate

-seealso: [](sec_scatter), `IS`, `PetscSection`, `PetscSectionSetPermutation()`, `PetscSectionCreate()`

# External Links
$(_doc_external("Vec/PetscSectionGetPermutation"))
"""
function PetscSectionGetPermutation(petsclib::PetscLibType, s::PetscSection, perm::IS) end

@for_petsc function PetscSectionGetPermutation(petsclib::$UnionPetscLib, s::PetscSection, perm::IS )
	perm_ = Ref(perm.ptr)

    @chk ccall(
               (:PetscSectionGetPermutation, $petsc_library),
               PetscErrorCode,
               (PetscSection, Ptr{CIS}),
               s, perm_,
              )

	perm.ptr = C_NULL

	return nothing
end 

"""
	PetscSectionSetPermutation(petsclib::PetscLibType,s::PetscSection, perm::IS) 
Sets a permutation of the chart for this section, [0, `pEnd`

Not Collective

Input Parameters:
- `s`    - the `PetscSection`
- `perm` - the permutation of points

Level: intermediate

-seealso: [](sec_scatter), `IS`, `PetscSection`, `PetscSectionSetUp()`, `PetscSectionGetPermutation()`, `PetscSectionPermute()`, `PetscSectionCreate()`

# External Links
$(_doc_external("Vec/PetscSectionSetPermutation"))
"""
function PetscSectionSetPermutation(petsclib::PetscLibType, s::PetscSection, perm::IS) end

@for_petsc function PetscSectionSetPermutation(petsclib::$UnionPetscLib, s::PetscSection, perm::IS )

    @chk ccall(
               (:PetscSectionSetPermutation, $petsc_library),
               PetscErrorCode,
               (PetscSection, CIS),
               s, perm,
              )


	return nothing
end 

"""
	PetscSectionGetBlockStarts(petsclib::PetscLibType,s::PetscSection, blockStarts::PetscBT) 
Returns a table indicating which points start new blocks

Not Collective, No Fortran Support

Input Parameter:
- `s` - the `PetscSection`

Output Parameter:
- `blockStarts` - The `PetscBT` with a 1 for each point that begins a block

-seealso: [](sec_scatter), `IS`, `PetscSection`, `PetscSectionSetBlockStarts()`, `PetscSectionCreate()`, `DMCreateMatrix()`, `MatSetVariableBlockSizes()`

# External Links
$(_doc_external("Vec/PetscSectionGetBlockStarts"))
"""
function PetscSectionGetBlockStarts(petsclib::PetscLibType, s::PetscSection, blockStarts::PetscBT) end

@for_petsc function PetscSectionGetBlockStarts(petsclib::$UnionPetscLib, s::PetscSection, blockStarts::PetscBT )

    @chk ccall(
               (:PetscSectionGetBlockStarts, $petsc_library),
               PetscErrorCode,
               (PetscSection, Ptr{PetscBT}),
               s, blockStarts,
              )


	return nothing
end 

"""
	PetscSectionSetBlockStarts(petsclib::PetscLibType,s::PetscSection, blockStarts::PetscBT) 
Sets a table indicating which points start new blocks

Not Collective, No Fortran Support

Input Parameters:
- `s`           - the `PetscSection`
- `blockStarts` - The `PetscBT` with a 1 for each point that begins a block

Level: intermediate

-seealso: [](sec_scatter), `IS`, `PetscSection`, `PetscSectionGetBlockStarts()`, `PetscSectionCreate()`, `DMCreateMatrix()`, `MatSetVariableBlockSizes()`

# External Links
$(_doc_external("Vec/PetscSectionSetBlockStarts"))
"""
function PetscSectionSetBlockStarts(petsclib::PetscLibType, s::PetscSection, blockStarts::PetscBT) end

@for_petsc function PetscSectionSetBlockStarts(petsclib::$UnionPetscLib, s::PetscSection, blockStarts::PetscBT )

    @chk ccall(
               (:PetscSectionSetBlockStarts, $petsc_library),
               PetscErrorCode,
               (PetscSection, PetscBT),
               s, blockStarts,
              )


	return nothing
end 

"""
	pm::PetscBool = PetscSectionGetPointMajor(petsclib::PetscLibType,s::PetscSection) 
Returns the flag for dof ordering, `PETSC_TRUE` if it is point major, `PETSC_FALSE` if it is field major

Not Collective

Input Parameter:
- `s` - the `PetscSection`

Output Parameter:
- `pm` - the flag for point major ordering

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionSetPointMajor()`

# External Links
$(_doc_external("Vec/PetscSectionGetPointMajor"))
"""
function PetscSectionGetPointMajor(petsclib::PetscLibType, s::PetscSection) end

@for_petsc function PetscSectionGetPointMajor(petsclib::$UnionPetscLib, s::PetscSection )
	pm_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscSectionGetPointMajor, $petsc_library),
               PetscErrorCode,
               (PetscSection, Ptr{PetscBool}),
               s, pm_,
              )

	pm = pm_[]

	return pm
end 

"""
	PetscSectionSetPointMajor(petsclib::PetscLibType,s::PetscSection, pm::PetscBool) 
Sets the flag for dof ordering, `PETSC_TRUE` for point major, otherwise it will be field major

Not Collective

Input Parameters:
- `s`  - the `PetscSection`
- `pm` - the flag for point major ordering

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionGetPointMajor()`, `PetscSectionSetPermutation()`

# External Links
$(_doc_external("Vec/PetscSectionSetPointMajor"))
"""
function PetscSectionSetPointMajor(petsclib::PetscLibType, s::PetscSection, pm::PetscBool) end

@for_petsc function PetscSectionSetPointMajor(petsclib::$UnionPetscLib, s::PetscSection, pm::PetscBool )

    @chk ccall(
               (:PetscSectionSetPointMajor, $petsc_library),
               PetscErrorCode,
               (PetscSection, PetscBool),
               s, pm,
              )


	return nothing
end 

"""
	includesConstraints::PetscBool = PetscSectionGetIncludesConstraints(petsclib::PetscLibType,s::PetscSection) 
Returns the flag indicating if constrained dofs were included when computing offsets in the `PetscSection`.
The value is set with `PetscSectionSetIncludesConstraints()`

Not Collective

Input Parameter:
- `s` - the `PetscSection`

Output Parameter:
- `includesConstraints` - the flag indicating if constrained dofs were included when computing offsets

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionSetIncludesConstraints()`

# External Links
$(_doc_external("Vec/PetscSectionGetIncludesConstraints"))
"""
function PetscSectionGetIncludesConstraints(petsclib::PetscLibType, s::PetscSection) end

@for_petsc function PetscSectionGetIncludesConstraints(petsclib::$UnionPetscLib, s::PetscSection )
	includesConstraints_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscSectionGetIncludesConstraints, $petsc_library),
               PetscErrorCode,
               (PetscSection, Ptr{PetscBool}),
               s, includesConstraints_,
              )

	includesConstraints = includesConstraints_[]

	return includesConstraints
end 

"""
	PetscSectionSetIncludesConstraints(petsclib::PetscLibType,s::PetscSection, includesConstraints::PetscBool) 
Sets the flag indicating if constrained dofs are to be included when computing offsets

Not Collective

Input Parameters:
- `s`                   - the `PetscSection`
- `includesConstraints` - the flag indicating if constrained dofs are to be included when computing offsets

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionGetIncludesConstraints()`

# External Links
$(_doc_external("Vec/PetscSectionSetIncludesConstraints"))
"""
function PetscSectionSetIncludesConstraints(petsclib::PetscLibType, s::PetscSection, includesConstraints::PetscBool) end

@for_petsc function PetscSectionSetIncludesConstraints(petsclib::$UnionPetscLib, s::PetscSection, includesConstraints::PetscBool )

    @chk ccall(
               (:PetscSectionSetIncludesConstraints, $petsc_library),
               PetscErrorCode,
               (PetscSection, PetscBool),
               s, includesConstraints,
              )


	return nothing
end 

"""
	numDof::PetscInt = PetscSectionGetDof(petsclib::PetscLibType,s::PetscSection, point::PetscInt) 
Return the total number of degrees of freedom associated with a given point.

Not Collective

Input Parameters:
- `s`     - the `PetscSection`
- `point` - the point

Output Parameter:
- `numDof` - the number of dof

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionSetDof()`, `PetscSectionCreate()`

# External Links
$(_doc_external("Vec/PetscSectionGetDof"))
"""
function PetscSectionGetDof(petsclib::PetscLibType, s::PetscSection, point::PetscInt) end

@for_petsc function PetscSectionGetDof(petsclib::$UnionPetscLib, s::PetscSection, point::$PetscInt )
	numDof_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscSectionGetDof, $petsc_library),
               PetscErrorCode,
               (PetscSection, $PetscInt, Ptr{$PetscInt}),
               s, point, numDof_,
              )

	numDof = numDof_[]

	return numDof
end 

"""
	PetscSectionSetDof(petsclib::PetscLibType,s::PetscSection, point::PetscInt, numDof::PetscInt) 
Sets the total number of degrees of freedom associated with a given point.

Not Collective

Input Parameters:
- `s`      - the `PetscSection`
- `point`  - the point
- `numDof` - the number of dof, these values may be negative -(dof+1) to indicate they are off process

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionGetDof()`, `PetscSectionAddDof()`, `PetscSectionCreate()`

# External Links
$(_doc_external("Vec/PetscSectionSetDof"))
"""
function PetscSectionSetDof(petsclib::PetscLibType, s::PetscSection, point::PetscInt, numDof::PetscInt) end

@for_petsc function PetscSectionSetDof(petsclib::$UnionPetscLib, s::PetscSection, point::$PetscInt, numDof::$PetscInt )

    @chk ccall(
               (:PetscSectionSetDof, $petsc_library),
               PetscErrorCode,
               (PetscSection, $PetscInt, $PetscInt),
               s, point, numDof,
              )


	return nothing
end 

"""
	PetscSectionAddDof(petsclib::PetscLibType,s::PetscSection, point::PetscInt, numDof::PetscInt) 
Adds to the total number of degrees of freedom associated with a given point.

Not Collective

Input Parameters:
- `s`      - the `PetscSection`
- `point`  - the point
- `numDof` - the number of additional dof

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionGetDof()`, `PetscSectionSetDof()`, `PetscSectionCreate()`

# External Links
$(_doc_external("Vec/PetscSectionAddDof"))
"""
function PetscSectionAddDof(petsclib::PetscLibType, s::PetscSection, point::PetscInt, numDof::PetscInt) end

@for_petsc function PetscSectionAddDof(petsclib::$UnionPetscLib, s::PetscSection, point::$PetscInt, numDof::$PetscInt )

    @chk ccall(
               (:PetscSectionAddDof, $petsc_library),
               PetscErrorCode,
               (PetscSection, $PetscInt, $PetscInt),
               s, point, numDof,
              )


	return nothing
end 

"""
	numDof::PetscInt = PetscSectionGetFieldDof(petsclib::PetscLibType,s::PetscSection, point::PetscInt, field::PetscInt) 
Return the number of degrees of freedom associated with a field on a given point.

Not Collective

Input Parameters:
- `s`     - the `PetscSection`
- `point` - the point
- `field` - the field

Output Parameter:
- `numDof` - the number of dof

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionSetFieldDof()`, `PetscSectionCreate()`

# External Links
$(_doc_external("Vec/PetscSectionGetFieldDof"))
"""
function PetscSectionGetFieldDof(petsclib::PetscLibType, s::PetscSection, point::PetscInt, field::PetscInt) end

@for_petsc function PetscSectionGetFieldDof(petsclib::$UnionPetscLib, s::PetscSection, point::$PetscInt, field::$PetscInt )
	numDof_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscSectionGetFieldDof, $petsc_library),
               PetscErrorCode,
               (PetscSection, $PetscInt, $PetscInt, Ptr{$PetscInt}),
               s, point, field, numDof_,
              )

	numDof = numDof_[]

	return numDof
end 

"""
	PetscSectionSetFieldDof(petsclib::PetscLibType,s::PetscSection, point::PetscInt, field::PetscInt, numDof::PetscInt) 
Sets the number of degrees of freedom associated with a field on a given point.

Not Collective

Input Parameters:
- `s`      - the `PetscSection`
- `point`  - the point
- `field`  - the field
- `numDof` - the number of dof, these values may be negative -(dof+1) to indicate they are off process

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionGetFieldDof()`, `PetscSectionCreate()`, `PetscSectionAddDof()`, `PetscSectionSetDof()`

# External Links
$(_doc_external("Vec/PetscSectionSetFieldDof"))
"""
function PetscSectionSetFieldDof(petsclib::PetscLibType, s::PetscSection, point::PetscInt, field::PetscInt, numDof::PetscInt) end

@for_petsc function PetscSectionSetFieldDof(petsclib::$UnionPetscLib, s::PetscSection, point::$PetscInt, field::$PetscInt, numDof::$PetscInt )

    @chk ccall(
               (:PetscSectionSetFieldDof, $petsc_library),
               PetscErrorCode,
               (PetscSection, $PetscInt, $PetscInt, $PetscInt),
               s, point, field, numDof,
              )


	return nothing
end 

"""
	PetscSectionAddFieldDof(petsclib::PetscLibType,s::PetscSection, point::PetscInt, field::PetscInt, numDof::PetscInt) 
Adds a number of degrees of freedom associated with a field on a given point.

Not Collective

Input Parameters:
- `s`      - the `PetscSection`
- `point`  - the point
- `field`  - the field
- `numDof` - the number of dof

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionSetFieldDof()`, `PetscSectionGetFieldDof()`, `PetscSectionCreate()`

# External Links
$(_doc_external("Vec/PetscSectionAddFieldDof"))
"""
function PetscSectionAddFieldDof(petsclib::PetscLibType, s::PetscSection, point::PetscInt, field::PetscInt, numDof::PetscInt) end

@for_petsc function PetscSectionAddFieldDof(petsclib::$UnionPetscLib, s::PetscSection, point::$PetscInt, field::$PetscInt, numDof::$PetscInt )

    @chk ccall(
               (:PetscSectionAddFieldDof, $petsc_library),
               PetscErrorCode,
               (PetscSection, $PetscInt, $PetscInt, $PetscInt),
               s, point, field, numDof,
              )


	return nothing
end 

"""
	numDof::PetscInt = PetscSectionGetConstraintDof(petsclib::PetscLibType,s::PetscSection, point::PetscInt) 
Return the number of constrained degrees of freedom associated with a given point.

Not Collective

Input Parameters:
- `s`     - the `PetscSection`
- `point` - the point

Output Parameter:
- `numDof` - the number of dof which are fixed by constraints

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionGetDof()`, `PetscSectionSetConstraintDof()`, `PetscSectionCreate()`

# External Links
$(_doc_external("Vec/PetscSectionGetConstraintDof"))
"""
function PetscSectionGetConstraintDof(petsclib::PetscLibType, s::PetscSection, point::PetscInt) end

@for_petsc function PetscSectionGetConstraintDof(petsclib::$UnionPetscLib, s::PetscSection, point::$PetscInt )
	numDof_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscSectionGetConstraintDof, $petsc_library),
               PetscErrorCode,
               (PetscSection, $PetscInt, Ptr{$PetscInt}),
               s, point, numDof_,
              )

	numDof = numDof_[]

	return numDof
end 

"""
	PetscSectionSetConstraintDof(petsclib::PetscLibType,s::PetscSection, point::PetscInt, numDof::PetscInt) 
Set the number of constrained degrees of freedom associated with a given point.

Not Collective

Input Parameters:
- `s`      - the `PetscSection`
- `point`  - the point
- `numDof` - the number of dof which are fixed by constraints

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionSetDof()`, `PetscSectionGetConstraintDof()`, `PetscSectionCreate()`

# External Links
$(_doc_external("Vec/PetscSectionSetConstraintDof"))
"""
function PetscSectionSetConstraintDof(petsclib::PetscLibType, s::PetscSection, point::PetscInt, numDof::PetscInt) end

@for_petsc function PetscSectionSetConstraintDof(petsclib::$UnionPetscLib, s::PetscSection, point::$PetscInt, numDof::$PetscInt )

    @chk ccall(
               (:PetscSectionSetConstraintDof, $petsc_library),
               PetscErrorCode,
               (PetscSection, $PetscInt, $PetscInt),
               s, point, numDof,
              )


	return nothing
end 

"""
	PetscSectionAddConstraintDof(petsclib::PetscLibType,s::PetscSection, point::PetscInt, numDof::PetscInt) 
Increment the number of constrained degrees of freedom associated with a given point.

Not Collective

Input Parameters:
- `s`      - the `PetscSection`
- `point`  - the point
- `numDof` - the number of additional dof which are fixed by constraints

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionAddDof()`, `PetscSectionGetConstraintDof()`, `PetscSectionCreate()`

# External Links
$(_doc_external("Vec/PetscSectionAddConstraintDof"))
"""
function PetscSectionAddConstraintDof(petsclib::PetscLibType, s::PetscSection, point::PetscInt, numDof::PetscInt) end

@for_petsc function PetscSectionAddConstraintDof(petsclib::$UnionPetscLib, s::PetscSection, point::$PetscInt, numDof::$PetscInt )

    @chk ccall(
               (:PetscSectionAddConstraintDof, $petsc_library),
               PetscErrorCode,
               (PetscSection, $PetscInt, $PetscInt),
               s, point, numDof,
              )


	return nothing
end 

"""
	numDof::PetscInt = PetscSectionGetFieldConstraintDof(petsclib::PetscLibType,s::PetscSection, point::PetscInt, field::PetscInt) 
Return the number of constrained degrees of freedom associated with a given field on a point.

Not Collective

Input Parameters:
- `s`     - the `PetscSection`
- `point` - the point
- `field` - the field

Output Parameter:
- `numDof` - the number of dof which are fixed by constraints

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionGetDof()`, `PetscSectionSetFieldConstraintDof()`, `PetscSectionCreate()`

# External Links
$(_doc_external("Vec/PetscSectionGetFieldConstraintDof"))
"""
function PetscSectionGetFieldConstraintDof(petsclib::PetscLibType, s::PetscSection, point::PetscInt, field::PetscInt) end

@for_petsc function PetscSectionGetFieldConstraintDof(petsclib::$UnionPetscLib, s::PetscSection, point::$PetscInt, field::$PetscInt )
	numDof_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscSectionGetFieldConstraintDof, $petsc_library),
               PetscErrorCode,
               (PetscSection, $PetscInt, $PetscInt, Ptr{$PetscInt}),
               s, point, field, numDof_,
              )

	numDof = numDof_[]

	return numDof
end 

"""
	PetscSectionSetFieldConstraintDof(petsclib::PetscLibType,s::PetscSection, point::PetscInt, field::PetscInt, numDof::PetscInt) 
Set the number of constrained degrees of freedom associated with a given field on a point.

Not Collective

Input Parameters:
- `s`      - the `PetscSection`
- `point`  - the point
- `field`  - the field
- `numDof` - the number of dof which are fixed by constraints

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionSetDof()`, `PetscSectionGetFieldConstraintDof()`, `PetscSectionCreate()`

# External Links
$(_doc_external("Vec/PetscSectionSetFieldConstraintDof"))
"""
function PetscSectionSetFieldConstraintDof(petsclib::PetscLibType, s::PetscSection, point::PetscInt, field::PetscInt, numDof::PetscInt) end

@for_petsc function PetscSectionSetFieldConstraintDof(petsclib::$UnionPetscLib, s::PetscSection, point::$PetscInt, field::$PetscInt, numDof::$PetscInt )

    @chk ccall(
               (:PetscSectionSetFieldConstraintDof, $petsc_library),
               PetscErrorCode,
               (PetscSection, $PetscInt, $PetscInt, $PetscInt),
               s, point, field, numDof,
              )


	return nothing
end 

"""
	PetscSectionAddFieldConstraintDof(petsclib::PetscLibType,s::PetscSection, point::PetscInt, field::PetscInt, numDof::PetscInt) 
Increment the number of constrained degrees of freedom associated with a given field on a point.

Not Collective

Input Parameters:
- `s`      - the `PetscSection`
- `point`  - the point
- `field`  - the field
- `numDof` - the number of additional dof which are fixed by constraints

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionAddDof()`, `PetscSectionGetFieldConstraintDof()`, `PetscSectionCreate()`

# External Links
$(_doc_external("Vec/PetscSectionAddFieldConstraintDof"))
"""
function PetscSectionAddFieldConstraintDof(petsclib::PetscLibType, s::PetscSection, point::PetscInt, field::PetscInt, numDof::PetscInt) end

@for_petsc function PetscSectionAddFieldConstraintDof(petsclib::$UnionPetscLib, s::PetscSection, point::$PetscInt, field::$PetscInt, numDof::$PetscInt )

    @chk ccall(
               (:PetscSectionAddFieldConstraintDof, $petsc_library),
               PetscErrorCode,
               (PetscSection, $PetscInt, $PetscInt, $PetscInt),
               s, point, field, numDof,
              )


	return nothing
end 

"""
	PetscSectionSetUpBC(petsclib::PetscLibType,s::PetscSection) 
Setup the subsections describing boundary conditions.

Not Collective

Input Parameter:
- `s` - the `PetscSection`

Level: advanced

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionSetUp()`, `PetscSectionCreate()`

# External Links
$(_doc_external("Vec/PetscSectionSetUpBC"))
"""
function PetscSectionSetUpBC(petsclib::PetscLibType, s::PetscSection) end

@for_petsc function PetscSectionSetUpBC(petsclib::$UnionPetscLib, s::PetscSection )

    @chk ccall(
               (:PetscSectionSetUpBC, $petsc_library),
               PetscErrorCode,
               (PetscSection,),
               s,
              )


	return nothing
end 

"""
	PetscSectionSetUp(petsclib::PetscLibType,s::PetscSection) 
Calculate offsets based upon the number of degrees of freedom for each point in preparation for use of the `PetscSection`

Not Collective

Input Parameter:
- `s` - the `PetscSection`

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionCreate()`, `PetscSectionSetPermutation()`

# External Links
$(_doc_external("Vec/PetscSectionSetUp"))
"""
function PetscSectionSetUp(petsclib::PetscLibType, s::PetscSection) end

@for_petsc function PetscSectionSetUp(petsclib::$UnionPetscLib, s::PetscSection )

    @chk ccall(
               (:PetscSectionSetUp, $petsc_library),
               PetscErrorCode,
               (PetscSection,),
               s,
              )


	return nothing
end 

"""
	maxDof::PetscInt = PetscSectionGetMaxDof(petsclib::PetscLibType,s::PetscSection) 
Return the maximum number of degrees of freedom on any point in the `PetscSection`

Not Collective

Input Parameter:
- `s` - the `PetscSection`

Output Parameter:
- `maxDof` - the maximum dof

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionGetDof()`, `PetscSectionSetDof()`, `PetscSectionAddDof()`, `PetscSectionCreate()`

# External Links
$(_doc_external("Vec/PetscSectionGetMaxDof"))
"""
function PetscSectionGetMaxDof(petsclib::PetscLibType, s::PetscSection) end

@for_petsc function PetscSectionGetMaxDof(petsclib::$UnionPetscLib, s::PetscSection )
	maxDof_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscSectionGetMaxDof, $petsc_library),
               PetscErrorCode,
               (PetscSection, Ptr{$PetscInt}),
               s, maxDof_,
              )

	maxDof = maxDof_[]

	return maxDof
end 

"""
	size::PetscInt = PetscSectionGetStorageSize(petsclib::PetscLibType,s::PetscSection) 
Return the size of an array or local `Vec` capable of holding all the degrees of freedom defined in a `PetscSection`

Not Collective

Input Parameter:
- `s` - the `PetscSection`

Output Parameter:
- `size` - the size of an array which can hold all the dofs

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionGetOffset()`, `PetscSectionGetConstrainedStorageSize()`, `PetscSectionCreate()`

# External Links
$(_doc_external("Vec/PetscSectionGetStorageSize"))
"""
function PetscSectionGetStorageSize(petsclib::PetscLibType, s::PetscSection) end

@for_petsc function PetscSectionGetStorageSize(petsclib::$UnionPetscLib, s::PetscSection )
	size_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscSectionGetStorageSize, $petsc_library),
               PetscErrorCode,
               (PetscSection, Ptr{$PetscInt}),
               s, size_,
              )

	size = size_[]

	return size
end 

"""
	size::PetscInt = PetscSectionGetConstrainedStorageSize(petsclib::PetscLibType,s::PetscSection) 
Return the size of an array or local `Vec` capable of holding all unconstrained degrees of freedom in a `PetscSection`

Not Collective

Input Parameter:
- `s` - the `PetscSection`

Output Parameter:
- `size` - the size of an array which can hold all unconstrained dofs

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionGetStorageSize()`, `PetscSectionGetOffset()`, `PetscSectionCreate()`

# External Links
$(_doc_external("Vec/PetscSectionGetConstrainedStorageSize"))
"""
function PetscSectionGetConstrainedStorageSize(petsclib::PetscLibType, s::PetscSection) end

@for_petsc function PetscSectionGetConstrainedStorageSize(petsclib::$UnionPetscLib, s::PetscSection )
	size_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscSectionGetConstrainedStorageSize, $petsc_library),
               PetscErrorCode,
               (PetscSection, Ptr{$PetscInt}),
               s, size_,
              )

	size = size_[]

	return size
end 

"""
	gsection::PetscSection = PetscSectionCreateGlobalSection(petsclib::PetscLibType,s::PetscSection, sf::PetscSF, usePermutation::PetscBool, includeConstraints::PetscBool, locOffsets::PetscBool) 
Create a parallel section describing the global layout using
a local (sequential) `PetscSection` on each MPI process and a `PetscSF` describing the section point overlap.

Input Parameters:
- `s`                  - The `PetscSection` for the local field layout
- `sf`                 - The `PetscSF` describing parallel layout of the section points (leaves are unowned local points)
- `usePermutation`     - By default this is `PETSC_TRUE`, meaning any permutation of the local section is transferred to the global section
- `includeConstraints` - By default this is `PETSC_FALSE`, meaning that the global field vector will not possess constrained dofs
- `localOffsets`       - If `PETSC_TRUE`, use local rather than global offsets for the points

Output Parameter:
- `gsection` - The `PetscSection` for the global field layout

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionCreate()`, `PetscSectionCreateGlobalSectionCensored()`

# External Links
$(_doc_external("Vec/PetscSectionCreateGlobalSection"))
"""
function PetscSectionCreateGlobalSection(petsclib::PetscLibType, s::PetscSection, sf::PetscSF, usePermutation::PetscBool, includeConstraints::PetscBool, locOffsets::PetscBool) end

@for_petsc function PetscSectionCreateGlobalSection(petsclib::$UnionPetscLib, s::PetscSection, sf::PetscSF, usePermutation::PetscBool, includeConstraints::PetscBool, locOffsets::PetscBool )
	gsection_ = Ref{PetscSection}()

    @chk ccall(
               (:PetscSectionCreateGlobalSection, $petsc_library),
               PetscErrorCode,
               (PetscSection, PetscSF, PetscBool, PetscBool, PetscBool, Ptr{PetscSection}),
               s, sf, usePermutation, includeConstraints, locOffsets, gsection_,
              )

	gsection = gsection_[]

	return gsection
end 

"""
	gsection::PetscSection = PetscSectionCreateGlobalSectionCensored(petsclib::PetscLibType,s::PetscSection, sf::PetscSF, includeConstraints::PetscBool, numExcludes::PetscInt, excludes::Vector{PetscInt}) 
Create a `PetscSection` describing the globallayout using
a local (sequential) `PetscSection` on each MPI process and an `PetscSF` describing the section point overlap.

Input Parameters:
- `s`                  - The `PetscSection` for the local field layout
- `sf`                 - The `PetscSF` describing parallel layout of the section points
- `includeConstraints` - By default this is `PETSC_FALSE`, meaning that the global vector will not possess constrained dofs
- `numExcludes`        - The number of exclusion ranges, this must have the same value on all MPI processes
- `excludes`           - An array [start_0, end_0, start_1, end_1, ...] where there are `numExcludes` pairs and must have the same values on all MPI processes

Output Parameter:
- `gsection` - The `PetscSection` for the global field layout

Level: advanced

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionCreate()`

# External Links
$(_doc_external("Vec/PetscSectionCreateGlobalSectionCensored"))
"""
function PetscSectionCreateGlobalSectionCensored(petsclib::PetscLibType, s::PetscSection, sf::PetscSF, includeConstraints::PetscBool, numExcludes::PetscInt, excludes::Vector{PetscInt}) end

@for_petsc function PetscSectionCreateGlobalSectionCensored(petsclib::$UnionPetscLib, s::PetscSection, sf::PetscSF, includeConstraints::PetscBool, numExcludes::$PetscInt, excludes::Vector{$PetscInt} )
	gsection_ = Ref{PetscSection}()

    @chk ccall(
               (:PetscSectionCreateGlobalSectionCensored, $petsc_library),
               PetscErrorCode,
               (PetscSection, PetscSF, PetscBool, $PetscInt, Ptr{$PetscInt}, Ptr{PetscSection}),
               s, sf, includeConstraints, numExcludes, excludes, gsection_,
              )

	gsection = gsection_[]

	return gsection
end 

"""
	PetscSectionGetPointLayout(petsclib::PetscLibType,comm::MPI_Comm, s::PetscSection, layout::PetscLayout) 
Get a `PetscLayout` for the points with nonzero dof counts of the unnamed default field within this `PetscSection`s local chart

Collective

Input Parameters:
- `comm` - The `MPI_Comm`
- `s`    - The `PetscSection`

Output Parameter:
- `layout` - The point layout for the data that defines the section

Level: advanced

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionGetValueLayout()`, `PetscSectionCreate()`

# External Links
$(_doc_external("Vec/PetscSectionGetPointLayout"))
"""
function PetscSectionGetPointLayout(petsclib::PetscLibType, comm::MPI_Comm, s::PetscSection, layout::PetscLayout) end

@for_petsc function PetscSectionGetPointLayout(petsclib::$UnionPetscLib, comm::MPI_Comm, s::PetscSection, layout::PetscLayout )

    @chk ccall(
               (:PetscSectionGetPointLayout, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, PetscSection, Ptr{PetscLayout}),
               comm, s, layout,
              )


	return nothing
end 

"""
	PetscSectionGetValueLayout(petsclib::PetscLibType,comm::MPI_Comm, s::PetscSection, layout::PetscLayout) 
Get the `PetscLayout` associated with the section dofs of a `PetscSection`

Collective

Input Parameters:
- `comm` - The `MPI_Comm`
- `s`    - The `PetscSection`

Output Parameter:
- `layout` - The dof layout for the section

Level: advanced

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionGetPointLayout()`, `PetscSectionCreate()`

# External Links
$(_doc_external("Vec/PetscSectionGetValueLayout"))
"""
function PetscSectionGetValueLayout(petsclib::PetscLibType, comm::MPI_Comm, s::PetscSection, layout::PetscLayout) end

@for_petsc function PetscSectionGetValueLayout(petsclib::$UnionPetscLib, comm::MPI_Comm, s::PetscSection, layout::PetscLayout )

    @chk ccall(
               (:PetscSectionGetValueLayout, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, PetscSection, Ptr{PetscLayout}),
               comm, s, layout,
              )


	return nothing
end 

"""
	offset::PetscInt = PetscSectionGetOffset(petsclib::PetscLibType,s::PetscSection, point::PetscInt) 
Return the offset into an array or `Vec` for the dof associated with the given point.

Not Collective

Input Parameters:
- `s`     - the `PetscSection`
- `point` - the point

Output Parameter:
- `offset` - the offset

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionGetFieldOffset()`, `PetscSectionCreate()`, `PetscSectionSetPointMajor()`

# External Links
$(_doc_external("Vec/PetscSectionGetOffset"))
"""
function PetscSectionGetOffset(petsclib::PetscLibType, s::PetscSection, point::PetscInt) end

@for_petsc function PetscSectionGetOffset(petsclib::$UnionPetscLib, s::PetscSection, point::$PetscInt )
	offset_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscSectionGetOffset, $petsc_library),
               PetscErrorCode,
               (PetscSection, $PetscInt, Ptr{$PetscInt}),
               s, point, offset_,
              )

	offset = offset_[]

	return offset
end 

"""
	PetscSectionSetOffset(petsclib::PetscLibType,s::PetscSection, point::PetscInt, offset::PetscInt) 
Set the offset into an array or `Vec` for the dof associated with the given point.

Not Collective

Input Parameters:
- `s`      - the `PetscSection`
- `point`  - the point
- `offset` - the offset, these values may be negative indicating the values are off process

Level: developer

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionGetFieldOffset()`, `PetscSectionCreate()`, `PetscSectionSetUp()`

# External Links
$(_doc_external("Vec/PetscSectionSetOffset"))
"""
function PetscSectionSetOffset(petsclib::PetscLibType, s::PetscSection, point::PetscInt, offset::PetscInt) end

@for_petsc function PetscSectionSetOffset(petsclib::$UnionPetscLib, s::PetscSection, point::$PetscInt, offset::$PetscInt )

    @chk ccall(
               (:PetscSectionSetOffset, $petsc_library),
               PetscErrorCode,
               (PetscSection, $PetscInt, $PetscInt),
               s, point, offset,
              )


	return nothing
end 

"""
	offset::PetscInt = PetscSectionGetFieldOffset(petsclib::PetscLibType,s::PetscSection, point::PetscInt, field::PetscInt) 
Return the offset into an array or `Vec` for the field dof associated with the given point.

Not Collective

Input Parameters:
- `s`     - the `PetscSection`
- `point` - the point
- `field` - the field

Output Parameter:
- `offset` - the offset

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionGetOffset()`, `PetscSectionCreate()`, `PetscSectionGetFieldPointOffset()`

# External Links
$(_doc_external("Vec/PetscSectionGetFieldOffset"))
"""
function PetscSectionGetFieldOffset(petsclib::PetscLibType, s::PetscSection, point::PetscInt, field::PetscInt) end

@for_petsc function PetscSectionGetFieldOffset(petsclib::$UnionPetscLib, s::PetscSection, point::$PetscInt, field::$PetscInt )
	offset_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscSectionGetFieldOffset, $petsc_library),
               PetscErrorCode,
               (PetscSection, $PetscInt, $PetscInt, Ptr{$PetscInt}),
               s, point, field, offset_,
              )

	offset = offset_[]

	return offset
end 

"""
	PetscSectionSetFieldOffset(petsclib::PetscLibType,s::PetscSection, point::PetscInt, field::PetscInt, offset::PetscInt) 
Set the offset into an array or `Vec` for the dof associated with the given field at a point.

Not Collective

Input Parameters:
- `s`      - the `PetscSection`
- `point`  - the point
- `field`  - the field
- `offset` - the offset, these values may be negative indicating the values are off process

Level: developer

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionGetFieldOffset()`, `PetscSectionSetOffset()`, `PetscSectionCreate()`, `PetscSectionSetUp()`

# External Links
$(_doc_external("Vec/PetscSectionSetFieldOffset"))
"""
function PetscSectionSetFieldOffset(petsclib::PetscLibType, s::PetscSection, point::PetscInt, field::PetscInt, offset::PetscInt) end

@for_petsc function PetscSectionSetFieldOffset(petsclib::$UnionPetscLib, s::PetscSection, point::$PetscInt, field::$PetscInt, offset::$PetscInt )

    @chk ccall(
               (:PetscSectionSetFieldOffset, $petsc_library),
               PetscErrorCode,
               (PetscSection, $PetscInt, $PetscInt, $PetscInt),
               s, point, field, offset,
              )


	return nothing
end 

"""
	offset::PetscInt = PetscSectionGetFieldPointOffset(petsclib::PetscLibType,s::PetscSection, point::PetscInt, field::PetscInt) 
Return the offset for the first field dof associated with the given point relative to the offset for that point for the
unnamed default field's first dof

Not Collective

Input Parameters:
- `s`     - the `PetscSection`
- `point` - the point
- `field` - the field

Output Parameter:
- `offset` - the offset

Level: advanced

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionGetOffset()`, `PetscSectionCreate()`, `PetscSectionGetFieldOffset()`

# External Links
$(_doc_external("Vec/PetscSectionGetFieldPointOffset"))
"""
function PetscSectionGetFieldPointOffset(petsclib::PetscLibType, s::PetscSection, point::PetscInt, field::PetscInt) end

@for_petsc function PetscSectionGetFieldPointOffset(petsclib::$UnionPetscLib, s::PetscSection, point::$PetscInt, field::$PetscInt )
	offset_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscSectionGetFieldPointOffset, $petsc_library),
               PetscErrorCode,
               (PetscSection, $PetscInt, $PetscInt, Ptr{$PetscInt}),
               s, point, field, offset_,
              )

	offset = offset_[]

	return offset
end 

"""
	start::PetscInt,end_::PetscInt = PetscSectionGetOffsetRange(petsclib::PetscLibType,s::PetscSection) 
Return the full range of offsets [`start`, `end`) for a `PetscSection`

Not Collective

Input Parameter:
- `s` - the `PetscSection`

Output Parameters:
- `start` - the minimum offset
- `end`   - one more than the maximum offset

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionGetOffset()`, `PetscSectionCreate()`

# External Links
$(_doc_external("Vec/PetscSectionGetOffsetRange"))
"""
function PetscSectionGetOffsetRange(petsclib::PetscLibType, s::PetscSection) end

@for_petsc function PetscSectionGetOffsetRange(petsclib::$UnionPetscLib, s::PetscSection )
	start_ = Ref{$PetscInt}()
	end__ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscSectionGetOffsetRange, $petsc_library),
               PetscErrorCode,
               (PetscSection, Ptr{$PetscInt}, Ptr{$PetscInt}),
               s, start_, end__,
              )

	start = start_[]
	end_ = end__[]

	return start,end_
end 

"""
	subs::PetscSection = PetscSectionCreateSubsection(petsclib::PetscLibType,s::PetscSection, len::PetscInt, fields::Vector{PetscInt}) 
Create a new, smaller `PetscSection` composed of only selected fields

Collective

Input Parameters:
- `s`      - the `PetscSection`
- `len`    - the number of subfields
- `fields` - the subfield numbers

Output Parameter:
- `subs` - the subsection

Level: advanced

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionCreateSupersection()`, `PetscSectionCreate()`

# External Links
$(_doc_external("Vec/PetscSectionCreateSubsection"))
"""
function PetscSectionCreateSubsection(petsclib::PetscLibType, s::PetscSection, len::PetscInt, fields::Vector{PetscInt}) end

@for_petsc function PetscSectionCreateSubsection(petsclib::$UnionPetscLib, s::PetscSection, len::$PetscInt, fields::Vector{$PetscInt} )
	subs_ = Ref{PetscSection}()

    @chk ccall(
               (:PetscSectionCreateSubsection, $petsc_library),
               PetscErrorCode,
               (PetscSection, $PetscInt, Ptr{$PetscInt}, Ptr{PetscSection}),
               s, len, fields, subs_,
              )

	subs = subs_[]

	return subs
end 

"""
	subs::PetscSection = PetscSectionCreateComponentSubsection(petsclib::PetscLibType,s::PetscSection, len::PetscInt, comps::Vector{PetscInt}) 
Create a new, smaller `PetscSection` composed of only selected components

Collective

Input Parameters:
- `s`     - the `PetscSection`
- `len`   - the number of components
- `comps` - the component numbers

Output Parameter:
- `subs` - the subsection

Level: advanced

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionCreateSupersection()`, `PetscSectionCreate()`

# External Links
$(_doc_external("Vec/PetscSectionCreateComponentSubsection"))
"""
function PetscSectionCreateComponentSubsection(petsclib::PetscLibType, s::PetscSection, len::PetscInt, comps::Vector{PetscInt}) end

@for_petsc function PetscSectionCreateComponentSubsection(petsclib::$UnionPetscLib, s::PetscSection, len::$PetscInt, comps::Vector{$PetscInt} )
	subs_ = Ref{PetscSection}()

    @chk ccall(
               (:PetscSectionCreateComponentSubsection, $petsc_library),
               PetscErrorCode,
               (PetscSection, $PetscInt, Ptr{$PetscInt}, Ptr{PetscSection}),
               s, len, comps, subs_,
              )

	subs = subs_[]

	return subs
end 

"""
	supers::PetscSection = PetscSectionCreateSupersection(petsclib::PetscLibType,s::Vector{PetscSection}, len::PetscInt) 
Create a new, larger section composed of multiple `PetscSection`s

Collective

Input Parameters:
- `s`   - the input sections
- `len` - the number of input sections

Output Parameter:
- `supers` - the supersection

Level: advanced

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionCreateSubsection()`, `PetscSectionCreate()`

# External Links
$(_doc_external("Vec/PetscSectionCreateSupersection"))
"""
function PetscSectionCreateSupersection(petsclib::PetscLibType, s::Vector{PetscSection}, len::PetscInt) end

@for_petsc function PetscSectionCreateSupersection(petsclib::$UnionPetscLib, s::Vector{PetscSection}, len::$PetscInt )
	supers_ = Ref{PetscSection}()

    @chk ccall(
               (:PetscSectionCreateSupersection, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscSection}, $PetscInt, Ptr{PetscSection}),
               s, len, supers_,
              )

	supers = supers_[]

	return supers
end 

"""
	subs::PetscSection = PetscSectionCreateSubmeshSection(petsclib::PetscLibType,s::PetscSection, subpointIS::IS) 
Create a new, smaller section with support on the submesh

Collective

Input Parameters:
- `s`          - the `PetscSection`
- `subpointIS` - a sorted list of points in the original mesh which are in the submesh

Output Parameter:
- `subs` - the subsection

Level: advanced

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionCreateSubdomainSection()`, `PetscSectionCreateSubsection()`, `DMPlexGetSubpointMap()`, `PetscSectionCreate()`

# External Links
$(_doc_external("Vec/PetscSectionCreateSubmeshSection"))
"""
function PetscSectionCreateSubmeshSection(petsclib::PetscLibType, s::PetscSection, subpointIS::IS) end

@for_petsc function PetscSectionCreateSubmeshSection(petsclib::$UnionPetscLib, s::PetscSection, subpointIS::IS )
	subs_ = Ref{PetscSection}()

    @chk ccall(
               (:PetscSectionCreateSubmeshSection, $petsc_library),
               PetscErrorCode,
               (PetscSection, CIS, Ptr{PetscSection}),
               s, subpointIS, subs_,
              )

	subs = subs_[]

	return subs
end 

"""
	subs::PetscSection = PetscSectionCreateSubdomainSection(petsclib::PetscLibType,s::PetscSection, subpointMap::IS) 
Create a new, smaller section with support on a subdomain of the mesh

Collective

Input Parameters:
- `s`           - the `PetscSection`
- `subpointMap` - a sorted list of points in the original mesh which are in the subdomain

Output Parameter:
- `subs` - the subsection

Level: advanced

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionCreateSubmeshSection()`, `PetscSectionCreateSubsection()`, `DMPlexGetSubpointMap()`, `PetscSectionCreate()`

# External Links
$(_doc_external("Vec/PetscSectionCreateSubdomainSection"))
"""
function PetscSectionCreateSubdomainSection(petsclib::PetscLibType, s::PetscSection, subpointMap::IS) end

@for_petsc function PetscSectionCreateSubdomainSection(petsclib::$UnionPetscLib, s::PetscSection, subpointMap::IS )
	subs_ = Ref{PetscSection}()

    @chk ccall(
               (:PetscSectionCreateSubdomainSection, $petsc_library),
               PetscErrorCode,
               (PetscSection, CIS, Ptr{PetscSection}),
               s, subpointMap, subs_,
              )

	subs = subs_[]

	return subs
end 

"""
	PetscSectionViewFromOptions(petsclib::PetscLibType,A::PetscSection, obj::PetscObject, name::String) 
View the `PetscSection` based on values in the options database

Collective

Input Parameters:
- `A`    - the `PetscSection` object to view
- `obj`  - Optional object that provides the options prefix used for the options
- `name` - command line option

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionView`, `PetscObjectViewFromOptions()`, `PetscSectionCreate()`, `PetscSectionView()`

# External Links
$(_doc_external("Vec/PetscSectionViewFromOptions"))
"""
function PetscSectionViewFromOptions(petsclib::PetscLibType, A::PetscSection, obj::PetscObject, name::String) end

@for_petsc function PetscSectionViewFromOptions(petsclib::$UnionPetscLib, A::PetscSection, obj::PetscObject, name::String )

    @chk ccall(
               (:PetscSectionViewFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscSection, PetscObject, Ptr{Cchar}),
               A, obj, name,
              )


	return nothing
end 

"""
	PetscSectionView(petsclib::PetscLibType,s::PetscSection, viewer::PetscViewer) 
Views a `PetscSection`

Collective

Input Parameters:
- `s`      - the `PetscSection` object to view
- `viewer` - the viewer

Level: beginner

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionCreate()`, `PetscSectionDestroy()`, `PetscSectionLoad()`, `PetscViewer`

# External Links
$(_doc_external("Vec/PetscSectionView"))
"""
function PetscSectionView(petsclib::PetscLibType, s::PetscSection, viewer::PetscViewer) end

@for_petsc function PetscSectionView(petsclib::$UnionPetscLib, s::PetscSection, viewer::PetscViewer )

    @chk ccall(
               (:PetscSectionView, $petsc_library),
               PetscErrorCode,
               (PetscSection, PetscViewer),
               s, viewer,
              )


	return nothing
end 

"""
	PetscSectionLoad(petsclib::PetscLibType,s::PetscSection, viewer::PetscViewer) 
Loads a `PetscSection`

Collective

Input Parameters:
- `s`      - the `PetscSection` object to load
- `viewer` - the viewer

Level: beginner

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionCreate()`, `PetscSectionDestroy()`, `PetscSectionView()`

# External Links
$(_doc_external("Vec/PetscSectionLoad"))
"""
function PetscSectionLoad(petsclib::PetscLibType, s::PetscSection, viewer::PetscViewer) end

@for_petsc function PetscSectionLoad(petsclib::$UnionPetscLib, s::PetscSection, viewer::PetscViewer )

    @chk ccall(
               (:PetscSectionLoad, $petsc_library),
               PetscErrorCode,
               (PetscSection, PetscViewer),
               s, viewer,
              )


	return nothing
end 

"""
	PetscSectionArrayView(petsclib::PetscLibType,s::PetscSection, array::Cvoid, data_type::PetscDataType, viewer::PetscViewer) 
View an array, using the section to structure the values

Collective

Input Parameters:
- `s`         - the organizing `PetscSection`
- `array`     - the array of values
- `data_type` - the `PetscDataType` of the array
- `viewer`    - the `PetscViewer`

Level: developer

-seealso: `PetscSection`, `PetscViewer`, `PetscSectionCreate()`, `VecSetValuesSection()`, `PetscSectionVecView()`

# External Links
$(_doc_external("Vec/PetscSectionArrayView"))
"""
function PetscSectionArrayView(petsclib::PetscLibType, s::PetscSection, array::Cvoid, data_type::PetscDataType, viewer::PetscViewer) end

@for_petsc function PetscSectionArrayView(petsclib::$UnionPetscLib, s::PetscSection, array::Cvoid, data_type::PetscDataType, viewer::PetscViewer )

    @chk ccall(
               (:PetscSectionArrayView, $petsc_library),
               PetscErrorCode,
               (PetscSection, Ptr{Cvoid}, PetscDataType, PetscViewer),
               s, array, data_type, viewer,
              )


	return nothing
end 

"""
	PetscSectionResetClosurePermutation(petsclib::PetscLibType,section::PetscSection) 
Remove any existing closure permutation

Input Parameter:
- `section` - The `PetscSection`

Level: intermediate

-seealso: `PetscSectionSetClosurePermutation()`, `PetscSectionSetClosureIndex()`, `PetscSectionReset()`

# External Links
$(_doc_external("Vec/PetscSectionResetClosurePermutation"))
"""
function PetscSectionResetClosurePermutation(petsclib::PetscLibType, section::PetscSection) end

@for_petsc function PetscSectionResetClosurePermutation(petsclib::$UnionPetscLib, section::PetscSection )

    @chk ccall(
               (:PetscSectionResetClosurePermutation, $petsc_library),
               PetscErrorCode,
               (PetscSection,),
               section,
              )


	return nothing
end 

"""
	PetscSectionReset(petsclib::PetscLibType,s::PetscSection) 
Frees all section data, the section is then as if `PetscSectionCreate()` had just been called.

Not Collective

Input Parameter:
- `s` - the `PetscSection`

Level: beginner

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionCreate()`

# External Links
$(_doc_external("Vec/PetscSectionReset"))
"""
function PetscSectionReset(petsclib::PetscLibType, s::PetscSection) end

@for_petsc function PetscSectionReset(petsclib::$UnionPetscLib, s::PetscSection )

    @chk ccall(
               (:PetscSectionReset, $petsc_library),
               PetscErrorCode,
               (PetscSection,),
               s,
              )


	return nothing
end 

"""
	PetscSectionDestroy(petsclib::PetscLibType,s::PetscSection) 
Frees a `PetscSection`

Not Collective

Input Parameter:
- `s` - the `PetscSection`

Level: beginner

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionCreate()`, `PetscSectionReset()`

# External Links
$(_doc_external("Vec/PetscSectionDestroy"))
"""
function PetscSectionDestroy(petsclib::PetscLibType, s::PetscSection) end

# Accept a Ref{PetscSection} (pointer-to-pointer) so callers can pass a Ref directly
@for_petsc function PetscSectionDestroy(petsclib::$UnionPetscLib, s::Ref{PetscSection} )

    @chk ccall(
               (:PetscSectionDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscSection},),
               s,
              )


	return nothing
end

@for_petsc function PetscSectionDestroy(petsclib::$UnionPetscLib, s::PetscSection )

    @chk ccall(
               (:PetscSectionDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscSection},),
               s,
              )


	return nothing
end 

"""
	hasConstraints::PetscBool = PetscSectionHasConstraints(petsclib::PetscLibType,s::PetscSection) 
Determine whether a `PetscSection` has constrained dofs

Not Collective

Input Parameter:
- `s` - The `PetscSection`

Output Parameter:
- `hasConstraints` - flag indicating that the section has constrained dofs

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSectionSetConstraintIndices()`, `PetscSectionGetConstraintDof()`, `PetscSection`

# External Links
$(_doc_external("Vec/PetscSectionHasConstraints"))
"""
function PetscSectionHasConstraints(petsclib::PetscLibType, s::PetscSection) end

@for_petsc function PetscSectionHasConstraints(petsclib::$UnionPetscLib, s::PetscSection )
	hasConstraints_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscSectionHasConstraints, $petsc_library),
               PetscErrorCode,
               (PetscSection, Ptr{PetscBool}),
               s, hasConstraints_,
              )

	hasConstraints = hasConstraints_[]

	return hasConstraints
end 

"""
	indices::Vector{PetscInt} = PetscSectionGetConstraintIndices(petsclib::PetscLibType,s::PetscSection, point::PetscInt) 
Get the point dof numbers, in [0, dof), which are constrained for a given point

Not Collective

Input Parameters:
- `s`     - The `PetscSection`
- `point` - The point

Output Parameter:
- `indices` - The constrained dofs

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSectionSetConstraintIndices()`, `PetscSectionGetConstraintDof()`, `PetscSection`

# External Links
$(_doc_external("Vec/PetscSectionGetConstraintIndices"))
"""
function PetscSectionGetConstraintIndices(petsclib::PetscLibType, s::PetscSection, point::PetscInt) end

@for_petsc function PetscSectionGetConstraintIndices(petsclib::$UnionPetscLib, s::PetscSection, point::$PetscInt )
	indices_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:PetscSectionGetConstraintIndices, $petsc_library),
               PetscErrorCode,
               (PetscSection, $PetscInt, Ptr{Ptr{$PetscInt}}),
               s, point, indices_,
              )

	indices = unsafe_wrap(Array, indices_[], VecGetLocalSize(petsclib, x); own = false)

	return indices
end 

"""
	PetscSectionSetConstraintIndices(petsclib::PetscLibType,s::PetscSection, point::PetscInt, indices::Vector{PetscInt}) 
Set the point dof numbers, in [0, dof), which are constrained

Not Collective

Input Parameters:
- `s`       - The `PetscSection`
- `point`   - The point
- `indices` - The constrained dofs

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSectionGetConstraintIndices()`, `PetscSectionGetConstraintDof()`, `PetscSection`

# External Links
$(_doc_external("Vec/PetscSectionSetConstraintIndices"))
"""
function PetscSectionSetConstraintIndices(petsclib::PetscLibType, s::PetscSection, point::PetscInt, indices::Vector{PetscInt}) end

@for_petsc function PetscSectionSetConstraintIndices(petsclib::$UnionPetscLib, s::PetscSection, point::$PetscInt, indices::Vector{$PetscInt} )

    @chk ccall(
               (:PetscSectionSetConstraintIndices, $petsc_library),
               PetscErrorCode,
               (PetscSection, $PetscInt, Ptr{$PetscInt}),
               s, point, indices,
              )


	return nothing
end 

"""
	indices::Vector{PetscInt} = PetscSectionGetFieldConstraintIndices(petsclib::PetscLibType,s::PetscSection, point::PetscInt, field::PetscInt) 
Get the field dof numbers, in [0, fdof), which are constrained

Not Collective

Input Parameters:
- `s`     - The `PetscSection`
- `field` - The field number
- `point` - The point

Output Parameter:
- `indices` - The constrained dofs sorted in ascending order, the length is returned by `PetscSectionGetConstraintDof()`.

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSectionSetFieldConstraintIndices()`, `PetscSectionGetConstraintIndices()`, `PetscSectionGetConstraintDof()`, `PetscSection`

# External Links
$(_doc_external("Vec/PetscSectionGetFieldConstraintIndices"))
"""
function PetscSectionGetFieldConstraintIndices(petsclib::PetscLibType, s::PetscSection, point::PetscInt, field::PetscInt) end

@for_petsc function PetscSectionGetFieldConstraintIndices(petsclib::$UnionPetscLib, s::PetscSection, point::$PetscInt, field::$PetscInt )
	indices_ = Ref{Ptr{$PetscInt}}()

    @chk ccall(
               (:PetscSectionGetFieldConstraintIndices, $petsc_library),
               PetscErrorCode,
               (PetscSection, $PetscInt, $PetscInt, Ptr{Ptr{$PetscInt}}),
               s, point, field, indices_,
              )

	indices = unsafe_wrap(Array, indices_[], VecGetLocalSize(petsclib, x); own = false)

	return indices
end 

"""
	PetscSectionSetFieldConstraintIndices(petsclib::PetscLibType,s::PetscSection, point::PetscInt, field::PetscInt, indices::Vector{PetscInt}) 
Set the field dof numbers, in [0, fdof), which are constrained

Not Collective

Input Parameters:
- `s`       - The `PetscSection`
- `point`   - The point
- `field`   - The field number
- `indices` - The constrained dofs

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSectionSetConstraintIndices()`, `PetscSectionGetFieldConstraintIndices()`, `PetscSectionGetConstraintDof()`, `PetscSection`

# External Links
$(_doc_external("Vec/PetscSectionSetFieldConstraintIndices"))
"""
function PetscSectionSetFieldConstraintIndices(petsclib::PetscLibType, s::PetscSection, point::PetscInt, field::PetscInt, indices::Vector{PetscInt}) end

@for_petsc function PetscSectionSetFieldConstraintIndices(petsclib::$UnionPetscLib, s::PetscSection, point::$PetscInt, field::$PetscInt, indices::Vector{$PetscInt} )

    @chk ccall(
               (:PetscSectionSetFieldConstraintIndices, $petsc_library),
               PetscErrorCode,
               (PetscSection, $PetscInt, $PetscInt, Ptr{$PetscInt}),
               s, point, field, indices,
              )


	return nothing
end 

"""
	PetscSectionPermute(petsclib::PetscLibType,section::PetscSection, permutation::IS, sectionNew::PetscSection) 
Reorder the section according to the input point permutation

Collective

Input Parameters:
- `section`     - The `PetscSection` object
- `permutation` - The point permutation, old point p becomes new point perm[p]

Output Parameter:
- `sectionNew` - The permuted `PetscSection`

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `IS`, `PetscSection`, `MatPermute()`, `PetscSectionSetPermutation()`

# External Links
$(_doc_external("Vec/PetscSectionPermute"))
"""
function PetscSectionPermute(petsclib::PetscLibType, section::PetscSection, permutation::IS, sectionNew::PetscSection) end

@for_petsc function PetscSectionPermute(petsclib::$UnionPetscLib, section::PetscSection, permutation::IS, sectionNew::PetscSection )

    @chk ccall(
               (:PetscSectionPermute, $petsc_library),
               PetscErrorCode,
               (PetscSection, CIS, Ptr{PetscSection}),
               section, permutation, sectionNew,
              )


	return nothing
end 

"""
	PetscSectionSetClosureIndex(petsclib::PetscLibType,section::PetscSection, obj::PetscObject, clSection::PetscSection, clPoints::IS) 
Create an internal data structure to speed up closure queries.

Collective

Input Parameters:
- `section`   - The `PetscSection`
- `obj`       - A `PetscObject` which serves as the key for this index
- `clSection` - `PetscSection` giving the size of the closure of each point
- `clPoints`  - `IS` giving the points in each closure

Level: advanced

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionGetClosureIndex()`, `DMPlexCreateClosureIndex()`

# External Links
$(_doc_external("Vec/PetscSectionSetClosureIndex"))
"""
function PetscSectionSetClosureIndex(petsclib::PetscLibType, section::PetscSection, obj::PetscObject, clSection::PetscSection, clPoints::IS) end

@for_petsc function PetscSectionSetClosureIndex(petsclib::$UnionPetscLib, section::PetscSection, obj::PetscObject, clSection::PetscSection, clPoints::IS )

    @chk ccall(
               (:PetscSectionSetClosureIndex, $petsc_library),
               PetscErrorCode,
               (PetscSection, PetscObject, PetscSection, CIS),
               section, obj, clSection, clPoints,
              )


	return nothing
end 

"""
	PetscSectionGetClosureIndex(petsclib::PetscLibType,section::PetscSection, obj::PetscObject, clSection::PetscSection, clPoints::IS) 
Get the cache of points in the closure of each point in the section set with `PetscSectionSetClosureIndex()`

Collective

Input Parameters:
- `section` - The `PetscSection`
- `obj`     - A `PetscObject` which serves as the key for this index

Output Parameters:
- `clSection` - `PetscSection` giving the size of the closure of each point
- `clPoints`  - `IS` giving the points in each closure

Level: advanced

-seealso: [PetscSection](ch_petscsection), `PetscSectionSetClosureIndex()`, `DMPlexCreateClosureIndex()`

# External Links
$(_doc_external("Vec/PetscSectionGetClosureIndex"))
"""
function PetscSectionGetClosureIndex(petsclib::PetscLibType, section::PetscSection, obj::PetscObject, clSection::PetscSection, clPoints::IS) end

@for_petsc function PetscSectionGetClosureIndex(petsclib::$UnionPetscLib, section::PetscSection, obj::PetscObject, clSection::PetscSection, clPoints::IS )
	clPoints_ = Ref(clPoints.ptr)

    @chk ccall(
               (:PetscSectionGetClosureIndex, $petsc_library),
               PetscErrorCode,
               (PetscSection, PetscObject, Ptr{PetscSection}, Ptr{CIS}),
               section, obj, clSection, clPoints_,
              )

	clPoints.ptr = C_NULL

	return nothing
end 

"""
	PetscSectionSetClosurePermutation(petsclib::PetscLibType,section::PetscSection, obj::PetscObject, depth::PetscInt, perm::IS) 
Set the dof permutation for the closure of each cell in the section, meaning clPerm[newIndex] = oldIndex.

Not Collective

Input Parameters:
- `section` - The `PetscSection`
- `obj`     - A `PetscObject` which serves as the key for this index (usually a `DM`)
- `depth`   - Depth of points on which to apply the given permutation
- `perm`    - Permutation of the cell dof closure

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `IS`, `PetscSectionGetClosurePermutation()`, `PetscSectionGetClosureIndex()`, `DMPlexCreateClosureIndex()`, `PetscCopyMode`

# External Links
$(_doc_external("Vec/PetscSectionSetClosurePermutation"))
"""
function PetscSectionSetClosurePermutation(petsclib::PetscLibType, section::PetscSection, obj::PetscObject, depth::PetscInt, perm::IS) end

@for_petsc function PetscSectionSetClosurePermutation(petsclib::$UnionPetscLib, section::PetscSection, obj::PetscObject, depth::$PetscInt, perm::IS )

    @chk ccall(
               (:PetscSectionSetClosurePermutation, $petsc_library),
               PetscErrorCode,
               (PetscSection, PetscObject, $PetscInt, CIS),
               section, obj, depth, perm,
              )


	return nothing
end 

"""
	PetscSectionGetClosurePermutation(petsclib::PetscLibType,section::PetscSection, obj::PetscObject, depth::PetscInt, clSize::PetscInt, perm::IS) 
Get the dof permutation for the closure of each cell in the section, meaning clPerm[newIndex] = oldIndex.

Not Collective

Input Parameters:
- `section` - The `PetscSection`
- `obj`     - A `PetscObject` which serves as the key for this index (usually a DM)
- `depth`   - Depth stratum on which to obtain closure permutation
- `clSize`  - Closure size to be permuted (e.g., may vary with element topology and degree)

Output Parameter:
- `perm` - The dof closure permutation

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `IS`, `PetscSectionSetClosurePermutation()`, `PetscSectionGetClosureInversePermutation()`, `PetscSectionGetClosureIndex()`, `PetscSectionSetClosureIndex()`, `DMPlexCreateClosureIndex()`

# External Links
$(_doc_external("Vec/PetscSectionGetClosurePermutation"))
"""
function PetscSectionGetClosurePermutation(petsclib::PetscLibType, section::PetscSection, obj::PetscObject, depth::PetscInt, clSize::PetscInt, perm::IS) end

@for_petsc function PetscSectionGetClosurePermutation(petsclib::$UnionPetscLib, section::PetscSection, obj::PetscObject, depth::$PetscInt, clSize::$PetscInt, perm::IS )
	perm_ = Ref(perm.ptr)

    @chk ccall(
               (:PetscSectionGetClosurePermutation, $petsc_library),
               PetscErrorCode,
               (PetscSection, PetscObject, $PetscInt, $PetscInt, Ptr{CIS}),
               section, obj, depth, clSize, perm_,
              )

	perm.ptr = C_NULL

	return nothing
end 

"""
	PetscSectionGetClosureInversePermutation(petsclib::PetscLibType,section::PetscSection, obj::PetscObject, depth::PetscInt, clSize::PetscInt, perm::IS) 
Get the inverse dof permutation for the closure of each cell in the section, meaning clPerm[oldIndex] = newIndex.

Not Collective

Input Parameters:
- `section` - The `PetscSection`
- `obj`     - A `PetscObject` which serves as the key for this index (usually a `DM`)
- `depth`   - Depth stratum on which to obtain closure permutation
- `clSize`  - Closure size to be permuted (e.g., may vary with element topology and degree)

Output Parameter:
- `perm` - The dof closure permutation

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `IS`, `PetscSectionSetClosurePermutation()`, `PetscSectionGetClosureIndex()`, `PetscSectionSetClosureIndex()`, `DMPlexCreateClosureIndex()`

# External Links
$(_doc_external("Vec/PetscSectionGetClosureInversePermutation"))
"""
function PetscSectionGetClosureInversePermutation(petsclib::PetscLibType, section::PetscSection, obj::PetscObject, depth::PetscInt, clSize::PetscInt, perm::IS) end

@for_petsc function PetscSectionGetClosureInversePermutation(petsclib::$UnionPetscLib, section::PetscSection, obj::PetscObject, depth::$PetscInt, clSize::$PetscInt, perm::IS )
	perm_ = Ref(perm.ptr)

    @chk ccall(
               (:PetscSectionGetClosureInversePermutation, $petsc_library),
               PetscErrorCode,
               (PetscSection, PetscObject, $PetscInt, $PetscInt, Ptr{CIS}),
               section, obj, depth, clSize, perm_,
              )

	perm.ptr = C_NULL

	return nothing
end 

"""
	PetscSectionGetField(petsclib::PetscLibType,s::PetscSection, field::PetscInt, subs::PetscSection) 
Get the `PetscSection` associated with a single field

Input Parameters:
- `s`     - The `PetscSection`
- `field` - The field number

Output Parameter:
- `subs` - The `PetscSection` for the given field, note the chart of `subs` is not set

Level: intermediate

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `IS`, `PetscSectionSetNumFields()`

# External Links
$(_doc_external("Vec/PetscSectionGetField"))
"""
function PetscSectionGetField(petsclib::PetscLibType, s::PetscSection, field::PetscInt, subs::PetscSection) end

@for_petsc function PetscSectionGetField(petsclib::$UnionPetscLib, s::PetscSection, field::$PetscInt, subs::PetscSection )

    @chk ccall(
               (:PetscSectionGetField, $petsc_library),
               PetscErrorCode,
               (PetscSection, $PetscInt, Ptr{PetscSection}),
               s, field, subs,
              )


	return nothing
end 

"""
	PetscSectionSetSym(petsclib::PetscLibType,section::PetscSection, sym::PetscSectionSym) 
Set the symmetries for the data referred to by the section

Collective

Input Parameters:
- `section` - the section describing data layout
- `sym`     - the symmetry describing the affect of orientation on the access of the data

Level: developer

-seealso: [PetscSection](ch_petscsection), `PetscSectionSym`, `PetscSectionGetSym()`, `PetscSectionSymCreate()`

# External Links
$(_doc_external("Vec/PetscSectionSetSym"))
"""
function PetscSectionSetSym(petsclib::PetscLibType, section::PetscSection, sym::PetscSectionSym) end

@for_petsc function PetscSectionSetSym(petsclib::$UnionPetscLib, section::PetscSection, sym::PetscSectionSym )

    @chk ccall(
               (:PetscSectionSetSym, $petsc_library),
               PetscErrorCode,
               (PetscSection, PetscSectionSym),
               section, sym,
              )


	return nothing
end 

"""
	PetscSectionGetSym(petsclib::PetscLibType,section::PetscSection, sym::PetscSectionSym) 
Get the symmetries for the data referred to by the section

Not Collective

Input Parameter:
- `section` - the section describing data layout

Output Parameter:
- `sym` - the symmetry describing the affect of orientation on the access of the data, provided previously by `PetscSectionSetSym()`

Level: developer

-seealso: [PetscSection](ch_petscsection), `PetscSectionSym`, `PetscSectionSetSym()`, `PetscSectionSymCreate()`

# External Links
$(_doc_external("Vec/PetscSectionGetSym"))
"""
function PetscSectionGetSym(petsclib::PetscLibType, section::PetscSection, sym::PetscSectionSym) end

@for_petsc function PetscSectionGetSym(petsclib::$UnionPetscLib, section::PetscSection, sym::PetscSectionSym )

    @chk ccall(
               (:PetscSectionGetSym, $petsc_library),
               PetscErrorCode,
               (PetscSection, Ptr{PetscSectionSym}),
               section, sym,
              )


	return nothing
end 

"""
	PetscSectionSetFieldSym(petsclib::PetscLibType,section::PetscSection, field::PetscInt, sym::PetscSectionSym) 
Set the symmetries for the data referred to by a field of the section

Collective

Input Parameters:
- `section` - the section describing data layout
- `field`   - the field number
- `sym`     - the symmetry describing the affect of orientation on the access of the data

Level: developer

-seealso: [PetscSection](ch_petscsection), `PetscSectionSym`, `PetscSectionGetFieldSym()`, `PetscSectionSymCreate()`

# External Links
$(_doc_external("Vec/PetscSectionSetFieldSym"))
"""
function PetscSectionSetFieldSym(petsclib::PetscLibType, section::PetscSection, field::PetscInt, sym::PetscSectionSym) end

@for_petsc function PetscSectionSetFieldSym(petsclib::$UnionPetscLib, section::PetscSection, field::$PetscInt, sym::PetscSectionSym )

    @chk ccall(
               (:PetscSectionSetFieldSym, $petsc_library),
               PetscErrorCode,
               (PetscSection, $PetscInt, PetscSectionSym),
               section, field, sym,
              )


	return nothing
end 

"""
	PetscSectionGetFieldSym(petsclib::PetscLibType,section::PetscSection, field::PetscInt, sym::PetscSectionSym) 
Get the symmetries for the data referred to by a field of the section

Collective

Input Parameters:
- `section` - the section describing data layout
- `field`   - the field number

Output Parameter:
- `sym` - the symmetry describing the affect of orientation on the access of the data

Level: developer

-seealso: [PetscSection](ch_petscsection), `PetscSectionSym`, `PetscSectionSetFieldSym()`, `PetscSectionSymCreate()`

# External Links
$(_doc_external("Vec/PetscSectionGetFieldSym"))
"""
function PetscSectionGetFieldSym(petsclib::PetscLibType, section::PetscSection, field::PetscInt, sym::PetscSectionSym) end

@for_petsc function PetscSectionGetFieldSym(petsclib::$UnionPetscLib, section::PetscSection, field::$PetscInt, sym::PetscSectionSym )

    @chk ccall(
               (:PetscSectionGetFieldSym, $petsc_library),
               PetscErrorCode,
               (PetscSection, $PetscInt, Ptr{PetscSectionSym}),
               section, field, sym,
              )


	return nothing
end 

"""
	perms::PetscInt,rots::PetscScalar = PetscSectionGetPointSyms(petsclib::PetscLibType,section::PetscSection, numPoints::PetscInt, points::PetscInt) 
Get the symmetries for a set of points in a `PetscSection` under specific orientations.

Not Collective

Input Parameters:
- `section`   - the section
- `numPoints` - the number of points
- `points`    - an array of size 2 * `numPoints`, containing a list of (point, orientation) pairs. (An orientation is an
arbitrary integer: its interpretation is up to sym.  Orientations are used by `DM`: for their interpretation in that
context, see `DMPlexGetConeOrientation()`).

Output Parameters:
- `perms` - The permutations for the given orientations (or `NULL` if there is no symmetry or the permutation is the identity).
- `rots`  - The field rotations symmetries for the given orientations (or `NULL` if there is no symmetry or the rotations are all
identity).

Example of usage, gathering dofs into a local array (lArray) from a section array (sArray):
-seealso: [PetscSection](ch_petscsection), `PetscSectionSym`, `PetscSectionRestorePointSyms()`, `PetscSectionSymCreate()`, `PetscSectionSetSym()`, `PetscSectionGetSym()`

# External Links
$(_doc_external("Vec/PetscSectionGetPointSyms"))
"""
function PetscSectionGetPointSyms(petsclib::PetscLibType, section::PetscSection, numPoints::PetscInt, points::PetscInt) end

@for_petsc function PetscSectionGetPointSyms(petsclib::$UnionPetscLib, section::PetscSection, numPoints::$PetscInt, points::$PetscInt )
	perms_ = Ref{$PetscInt}()
	rots_ = Ref{$PetscScalar}()

    @chk ccall(
               (:PetscSectionGetPointSyms, $petsc_library),
               PetscErrorCode,
               (PetscSection, $PetscInt, Ptr{$PetscInt}, $PetscInt, $PetscScalar),
               section, numPoints, points, perms_, rots_,
              )

	perms = perms_[]
	rots = rots_[]

	return perms,rots
end 

"""
	PetscSectionRestorePointSyms(petsclib::PetscLibType,section::PetscSection, numPoints::PetscInt, points::PetscInt, perms::PetscInt, rots::PetscScalar) 
Restore the symmetries returned by `PetscSectionGetPointSyms()`

Not Collective

Input Parameters:
- `section`   - the section
- `numPoints` - the number of points
- `points`    - an array of size 2 * `numPoints`, containing a list of (point, orientation) pairs. (An orientation is an
arbitrary integer: its interpretation is up to sym.  Orientations are used by `DM`: for their interpretation in that
context, see `DMPlexGetConeOrientation()`).
- `perms`     - The permutations for the given orientations: set to `NULL` at conclusion
- `rots`      - The field rotations symmetries for the given orientations: set to `NULL` at conclusion

Level: developer

-seealso: [PetscSection](ch_petscsection), `PetscSectionSym`, `PetscSectionGetPointSyms()`, `PetscSectionSymCreate()`, `PetscSectionSetSym()`, `PetscSectionGetSym()`

# External Links
$(_doc_external("Vec/PetscSectionRestorePointSyms"))
"""
function PetscSectionRestorePointSyms(petsclib::PetscLibType, section::PetscSection, numPoints::PetscInt, points::PetscInt, perms::PetscInt, rots::PetscScalar) end

@for_petsc function PetscSectionRestorePointSyms(petsclib::$UnionPetscLib, section::PetscSection, numPoints::$PetscInt, points::$PetscInt, perms::$PetscInt, rots::$PetscScalar )

    @chk ccall(
               (:PetscSectionRestorePointSyms, $petsc_library),
               PetscErrorCode,
               (PetscSection, $PetscInt, Ptr{$PetscInt}, $PetscInt, $PetscScalar),
               section, numPoints, points, perms, rots,
              )


	return nothing
end 

"""
	perms::PetscInt,rots::PetscScalar = PetscSectionGetFieldPointSyms(petsclib::PetscLibType,section::PetscSection, field::PetscInt, numPoints::PetscInt, points::PetscInt) 
Get the symmetries for a set of points in a field of a `PetscSection` under specific orientations.

Not Collective

Input Parameters:
- `section`   - the section
- `field`     - the field of the section
- `numPoints` - the number of points
- `points`    - an array of size 2 * `numPoints`, containing a list of (point, orientation) pairs. (An orientation is an
arbitrary integer: its interpretation is up to sym.  Orientations are used by `DM`: for their interpretation in that
context, see `DMPlexGetConeOrientation()`).

Output Parameters:
- `perms` - The permutations for the given orientations (or `NULL` if there is no symmetry or the permutation is the identity).
- `rots`  - The field rotations symmetries for the given orientations (or `NULL` if there is no symmetry or the rotations are all
identity).

Level: developer

-seealso: [PetscSection](ch_petscsection), `PetscSectionSym`, `PetscSectionGetPointSyms()`, `PetscSectionRestoreFieldPointSyms()`

# External Links
$(_doc_external("Vec/PetscSectionGetFieldPointSyms"))
"""
function PetscSectionGetFieldPointSyms(petsclib::PetscLibType, section::PetscSection, field::PetscInt, numPoints::PetscInt, points::PetscInt) end

@for_petsc function PetscSectionGetFieldPointSyms(petsclib::$UnionPetscLib, section::PetscSection, field::$PetscInt, numPoints::$PetscInt, points::$PetscInt )
	perms_ = Ref{$PetscInt}()
	rots_ = Ref{$PetscScalar}()

    @chk ccall(
               (:PetscSectionGetFieldPointSyms, $petsc_library),
               PetscErrorCode,
               (PetscSection, $PetscInt, $PetscInt, Ptr{$PetscInt}, $PetscInt, $PetscScalar),
               section, field, numPoints, points, perms_, rots_,
              )

	perms = perms_[]
	rots = rots_[]

	return perms,rots
end 

"""
	PetscSectionRestoreFieldPointSyms(petsclib::PetscLibType,section::PetscSection, field::PetscInt, numPoints::PetscInt, points::PetscInt, perms::PetscInt, rots::PetscScalar) 
Restore the symmetries returned by `PetscSectionGetFieldPointSyms()`

Not Collective

Input Parameters:
- `section`   - the section
- `field`     - the field number
- `numPoints` - the number of points
- `points`    - an array of size 2 * `numPoints`, containing a list of (point, orientation) pairs. (An orientation is an
arbitrary integer: its interpretation is up to sym.  Orientations are used by `DM`: for their interpretation in that
context, see `DMPlexGetConeOrientation()`).
- `perms`     - The permutations for the given orientations: set to NULL at conclusion
- `rots`      - The field rotations symmetries for the given orientations: set to NULL at conclusion

Level: developer

-seealso: [PetscSection](ch_petscsection), `PetscSectionSym`, `PetscSectionRestorePointSyms()`, `petscSectionGetFieldPointSyms()`, `PetscSectionSymCreate()`, `PetscSectionSetSym()`, `PetscSectionGetSym()`

# External Links
$(_doc_external("Vec/PetscSectionRestoreFieldPointSyms"))
"""
function PetscSectionRestoreFieldPointSyms(petsclib::PetscLibType, section::PetscSection, field::PetscInt, numPoints::PetscInt, points::PetscInt, perms::PetscInt, rots::PetscScalar) end

@for_petsc function PetscSectionRestoreFieldPointSyms(petsclib::$UnionPetscLib, section::PetscSection, field::$PetscInt, numPoints::$PetscInt, points::$PetscInt, perms::$PetscInt, rots::$PetscScalar )

    @chk ccall(
               (:PetscSectionRestoreFieldPointSyms, $petsc_library),
               PetscErrorCode,
               (PetscSection, $PetscInt, $PetscInt, Ptr{$PetscInt}, $PetscInt, $PetscScalar),
               section, field, numPoints, points, perms, rots,
              )


	return nothing
end 

"""
	flg::PetscBool = PetscSectionGetUseFieldOffsets(petsclib::PetscLibType,s::PetscSection) 
Get the flag indicating if field offsets are used directly in a global section, rather than just the point offset

Not Collective

Input Parameter:
- `s` - the global `PetscSection`

Output Parameter:
- `flg` - the flag

Level: developer

-seealso: [PetscSection](ch_petscsection), `PetscSectionSym`, `PetscSectionSetChart()`, `PetscSectionCreate()`

# External Links
$(_doc_external("Vec/PetscSectionGetUseFieldOffsets"))
"""
function PetscSectionGetUseFieldOffsets(petsclib::PetscLibType, s::PetscSection) end

@for_petsc function PetscSectionGetUseFieldOffsets(petsclib::$UnionPetscLib, s::PetscSection )
	flg_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscSectionGetUseFieldOffsets, $petsc_library),
               PetscErrorCode,
               (PetscSection, Ptr{PetscBool}),
               s, flg_,
              )

	flg = flg_[]

	return flg
end 

"""
	PetscSectionSetUseFieldOffsets(petsclib::PetscLibType,s::PetscSection, flg::PetscBool) 
Set the flag to use field offsets directly in a global section, rather than just the point offset

Not Collective

Input Parameters:
- `s`   - the global `PetscSection`
- `flg` - the flag

Level: developer

-seealso: [PetscSection](ch_petscsection), `PetscSectionSym`, `PetscSectionGetUseFieldOffsets()`, `PetscSectionSetChart()`, `PetscSectionCreate()`

# External Links
$(_doc_external("Vec/PetscSectionSetUseFieldOffsets"))
"""
function PetscSectionSetUseFieldOffsets(petsclib::PetscLibType, s::PetscSection, flg::PetscBool) end

@for_petsc function PetscSectionSetUseFieldOffsets(petsclib::$UnionPetscLib, s::PetscSection, flg::PetscBool )

    @chk ccall(
               (:PetscSectionSetUseFieldOffsets, $petsc_library),
               PetscErrorCode,
               (PetscSection, PetscBool),
               s, flg,
              )


	return nothing
end 

"""
	PetscSectionVecView(petsclib::PetscLibType,s::PetscSection, v::PetscVec, viewer::PetscViewer) 
View a vector, using the section to structure the values

Collective

Input Parameters:
- `s`      - the organizing `PetscSection`
- `v`      - the `Vec`
- `viewer` - the `PetscViewer`

Level: developer

-seealso: `PetscSection`, `PetscViewer`, `PetscSectionCreate()`, `VecSetValuesSection()`, `PetscSectionArrayView()`

# External Links
$(_doc_external("Vec/PetscSectionVecView"))
"""
function PetscSectionVecView(petsclib::PetscLibType, s::PetscSection, v::PetscVec, viewer::PetscViewer) end

@for_petsc function PetscSectionVecView(petsclib::$UnionPetscLib, s::PetscSection, v::PetscVec, viewer::PetscViewer )

    @chk ccall(
               (:PetscSectionVecView, $petsc_library),
               PetscErrorCode,
               (PetscSection, CVec, PetscViewer),
               s, v, viewer,
              )


	return nothing
end 

"""
	val::Vector{PetscReal} = PetscSectionVecNorm(petsclib::PetscLibType,s::PetscSection, gs::PetscSection, x::PetscVec, type::NormType) 
Computes the vector norm of each field

Input Parameters:
- `s`    - the local Section
- `gs`   - the global section
- `x`    - the vector
- `type` - one of `NORM_1`, `NORM_2`, `NORM_INFINITY`.

Output Parameter:
- `val` - the array of norms

Level: intermediate

-seealso: `VecNorm()`, `PetscSectionCreate()`

# External Links
$(_doc_external("Vec/PetscSectionVecNorm"))
"""
function PetscSectionVecNorm(petsclib::PetscLibType, s::PetscSection, gs::PetscSection, x::PetscVec, type::NormType) end

@for_petsc function PetscSectionVecNorm(petsclib::$UnionPetscLib, s::PetscSection, gs::PetscSection, x::PetscVec, type::NormType )
	val = Vector{$PetscReal}(undef, ni);  # CHECK SIZE!!

    @chk ccall(
               (:PetscSectionVecNorm, $petsc_library),
               PetscErrorCode,
               (PetscSection, PetscSection, CVec, NormType, Ptr{$PetscReal}),
               s, gs, x, type, val,
              )


	return val
end 

"""
	gsection::PetscSection = PetscSectionCreateGlobalSectionLabel(petsclib::PetscLibType,s::PetscSection, sf::PetscSF, includeConstraints::PetscBool, label::DMLabel, labelValue::PetscInt) 
Create a section describing the global field layout using
the local section and an `PetscSF` describing the section point overlap.

Collective

Input Parameters:
- `s`                  - The `PetscSection` for the local field layout
- `sf`                 - The `PetscSF` describing parallel layout of the section points
- `includeConstraints` - By default this is `PETSC_FALSE`, meaning that the global field vector will not possess constrained dofs
- `label`              - The label specifying the points
- `labelValue`         - The label stratum specifying the points

Output Parameter:
- `gsection` - The `PetscSection` for the global field layout

Level: developer

-seealso: `DMLabel`, `DM`, `PetscSectionCreate()`

# External Links
$(_doc_external("Dm/PetscSectionCreateGlobalSectionLabel"))
"""
function PetscSectionCreateGlobalSectionLabel(petsclib::PetscLibType, s::PetscSection, sf::PetscSF, includeConstraints::PetscBool, label::DMLabel, labelValue::PetscInt) end

@for_petsc function PetscSectionCreateGlobalSectionLabel(petsclib::$UnionPetscLib, s::PetscSection, sf::PetscSF, includeConstraints::PetscBool, label::DMLabel, labelValue::$PetscInt )
	gsection_ = Ref{PetscSection}()

    @chk ccall(
               (:PetscSectionCreateGlobalSectionLabel, $petsc_library),
               PetscErrorCode,
               (PetscSection, PetscSF, PetscBool, DMLabel, $PetscInt, Ptr{PetscSection}),
               s, sf, includeConstraints, label, labelValue, gsection_,
              )

	gsection = gsection_[]

	return gsection
end 

"""
	sym::PetscSectionSym = PetscSectionSymCreate(petsclib::PetscLibType,comm::MPI_Comm) 
Creates an empty `PetscSectionSym` object.

Collective

Input Parameter:
- `comm` - the MPI communicator

Output Parameter:
- `sym` - pointer to the new set of symmetries

Level: developer

-seealso: [PetscSection](ch_petscsection), `PetscSection`, `PetscSectionSym`, `PetscSectionSymDestroy()`

# External Links
$(_doc_external("Vec/PetscSectionSymCreate"))
"""
function PetscSectionSymCreate(petsclib::PetscLibType, comm::MPI_Comm) end

@for_petsc function PetscSectionSymCreate(petsclib::$UnionPetscLib, comm::MPI_Comm )
	sym_ = Ref{PetscSectionSym}()

    @chk ccall(
               (:PetscSectionSymCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PetscSectionSym}),
               comm, sym_,
              )

	sym = sym_[]

	return sym
end 

"""
	PetscSectionSymSetType(petsclib::PetscLibType,sym::PetscSectionSym, method::PetscSectionSymType) 
Builds a `PetscSectionSym`, for a particular implementation.

Collective

Input Parameters:
- `sym`    - The section symmetry object
- `method` - The name of the section symmetry type

Level: developer

-seealso: [PetscSection](ch_petscsection), `PetscSectionSym`, `PetscSectionSymType`, `PetscSectionSymGetType()`, `PetscSectionSymCreate()`

# External Links
$(_doc_external("Vec/PetscSectionSymSetType"))
"""
function PetscSectionSymSetType(petsclib::PetscLibType, sym::PetscSectionSym, method::PetscSectionSymType) end

@for_petsc function PetscSectionSymSetType(petsclib::$UnionPetscLib, sym::PetscSectionSym, method::PetscSectionSymType )

    @chk ccall(
               (:PetscSectionSymSetType, $petsc_library),
               PetscErrorCode,
               (PetscSectionSym, PetscSectionSymType),
               sym, method,
              )


	return nothing
end 

"""
	type::PetscSectionSymType = PetscSectionSymGetType(petsclib::PetscLibType,sym::PetscSectionSym) 
Gets the section symmetry type name (as a string) from the `PetscSectionSym`.

Not Collective

Input Parameter:
- `sym` - The section symmetry

Output Parameter:
- `type` - The index set type name

Level: developer

-seealso: [PetscSection](ch_petscsection), `PetscSectionSym`, `PetscSectionSymType`, `PetscSectionSymSetType()`, `PetscSectionSymCreate()`

# External Links
$(_doc_external("Vec/PetscSectionSymGetType"))
"""
function PetscSectionSymGetType(petsclib::PetscLibType, sym::PetscSectionSym) end

@for_petsc function PetscSectionSymGetType(petsclib::$UnionPetscLib, sym::PetscSectionSym )
	type_ = Ref{PetscSectionSymType}()

    @chk ccall(
               (:PetscSectionSymGetType, $petsc_library),
               PetscErrorCode,
               (PetscSectionSym, Ptr{PetscSectionSymType}),
               sym, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	PetscSectionSymRegister(petsclib::PetscLibType,sname::String, fnc::external) 
Registers a new section symmetry implementation

Not Collective, No Fortran Support

Input Parameters:
- `sname`    - The name of a new user-defined creation routine
- `function` - The creation routine itself

Level: developer

-seealso: [PetscSection](ch_petscsection), `PetscSectionSym`, `PetscSectionSymType`, `PetscSectionSymCreate()`, `PetscSectionSymSetType()`

# External Links
$(_doc_external("Vec/PetscSectionSymRegister"))
"""
function PetscSectionSymRegister(petsclib::PetscLibType, sname::String, fnc::external) end

@for_petsc function PetscSectionSymRegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:PetscSectionSymRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	PetscSectionSymDestroy(petsclib::PetscLibType,sym::PetscSectionSym) 
Destroys a section symmetry.

Collective

Input Parameter:
- `sym` - the section symmetry

Level: developer

-seealso: [PetscSection](ch_petscsection), `PetscSectionSym`, `PetscSectionSymCreate()`

# External Links
$(_doc_external("Vec/PetscSectionSymDestroy"))
"""
function PetscSectionSymDestroy(petsclib::PetscLibType, sym::PetscSectionSym) end

@for_petsc function PetscSectionSymDestroy(petsclib::$UnionPetscLib, sym::PetscSectionSym )

    @chk ccall(
               (:PetscSectionSymDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscSectionSym},),
               sym,
              )


	return nothing
end 

"""
	PetscSectionSymView(petsclib::PetscLibType,sym::PetscSectionSym, viewer::PetscViewer) 
Displays a section symmetry

Collective

Input Parameters:
- `sym`    - the index set
- `viewer` - viewer used to display the set, for example `PETSC_VIEWER_STDOUT_SELF`.

Level: developer

-seealso: `PetscSectionSym`, `PetscViewer`, `PetscViewerASCIIOpen()`

# External Links
$(_doc_external("Vec/PetscSectionSymView"))
"""
function PetscSectionSymView(petsclib::PetscLibType, sym::PetscSectionSym, viewer::PetscViewer) end

@for_petsc function PetscSectionSymView(petsclib::$UnionPetscLib, sym::PetscSectionSym, viewer::PetscViewer )

    @chk ccall(
               (:PetscSectionSymView, $petsc_library),
               PetscErrorCode,
               (PetscSectionSym, PetscViewer),
               sym, viewer,
              )


	return nothing
end 

"""
	PetscSectionSymCopy(petsclib::PetscLibType,sym::PetscSectionSym, nsym::PetscSectionSym) 
Copy the symmetries, assuming that the point structure is compatible

Not Collective

Input Parameter:
- `sym` - the `PetscSectionSym`

Output Parameter:
- `nsym` - the equivalent symmetries

Level: developer

-seealso: [PetscSection](ch_petscsection), `PetscSectionSym`, `PetscSectionSymCreate()`, `PetscSectionSetSym()`, `PetscSectionGetSym()`, `PetscSectionSymLabelSetStratum()`, `PetscSectionGetPointSyms()`

# External Links
$(_doc_external("Vec/PetscSectionSymCopy"))
"""
function PetscSectionSymCopy(petsclib::PetscLibType, sym::PetscSectionSym, nsym::PetscSectionSym) end

@for_petsc function PetscSectionSymCopy(petsclib::$UnionPetscLib, sym::PetscSectionSym, nsym::PetscSectionSym )

    @chk ccall(
               (:PetscSectionSymCopy, $petsc_library),
               PetscErrorCode,
               (PetscSectionSym, PetscSectionSym),
               sym, nsym,
              )


	return nothing
end 

"""
	PetscSectionSymDistribute(petsclib::PetscLibType,sym::PetscSectionSym, migrationSF::PetscSF, dsym::PetscSectionSym) 
Distribute the symmetries in accordance with the input `PetscSF`

Collective

Input Parameters:
- `sym`         - the `PetscSectionSym`
- `migrationSF` - the distribution map from roots to leaves

Output Parameter:
- `dsym` - the redistributed symmetries

Level: developer

-seealso: [PetscSection](ch_petscsection), `PetscSectionSym`, `PetscSectionSymCreate()`, `PetscSectionSetSym()`, `PetscSectionGetSym()`, `PetscSectionSymLabelSetStratum()`, `PetscSectionGetPointSyms()`

# External Links
$(_doc_external("Vec/PetscSectionSymDistribute"))
"""
function PetscSectionSymDistribute(petsclib::PetscLibType, sym::PetscSectionSym, migrationSF::PetscSF, dsym::PetscSectionSym) end

@for_petsc function PetscSectionSymDistribute(petsclib::$UnionPetscLib, sym::PetscSectionSym, migrationSF::PetscSF, dsym::PetscSectionSym )

    @chk ccall(
               (:PetscSectionSymDistribute, $petsc_library),
               PetscErrorCode,
               (PetscSectionSym, PetscSF, Ptr{PetscSectionSym}),
               sym, migrationSF, dsym,
              )


	return nothing
end 

"""
	PetscSectionSymLabelSetLabel(petsclib::PetscLibType,sym::PetscSectionSym, label::DMLabel) 
set the label whose strata will define the points that receive symmetries

Logically

Input Parameters:
- `sym`   - the section symmetries
- `label` - the `DMLabel` describing the types of points

Level: developer:

-seealso: `DMLabel`, `DM`, `PetscSectionSymLabelSetStratum()`, `PetscSectionSymCreateLabel()`, `PetscSectionGetPointSyms()`

# External Links
$(_doc_external("Dm/PetscSectionSymLabelSetLabel"))
"""
function PetscSectionSymLabelSetLabel(petsclib::PetscLibType, sym::PetscSectionSym, label::DMLabel) end

@for_petsc function PetscSectionSymLabelSetLabel(petsclib::$UnionPetscLib, sym::PetscSectionSym, label::DMLabel )

    @chk ccall(
               (:PetscSectionSymLabelSetLabel, $petsc_library),
               PetscErrorCode,
               (PetscSectionSym, DMLabel),
               sym, label,
              )


	return nothing
end 

"""
	size::PetscInt,minOrient::PetscInt,maxOrient::PetscInt,perms::PetscInt,rots::PetscScalar = PetscSectionSymLabelGetStratum(petsclib::PetscLibType,sym::PetscSectionSym, stratum::PetscInt) 
get the symmetries for the orientations of a stratum

Logically Collective

Input Parameters:
- `sym`     - the section symmetries
- `stratum` - the stratum value in the label that we are assigning symmetries for

Output Parameters:
- `size`      - the number of dofs for points in the `stratum` of the label
- `minOrient` - the smallest orientation for a point in this `stratum`
- `maxOrient` - one greater than the largest orientation for a ppoint in this `stratum` (i.e., orientations are in the range [`minOrient`, `maxOrient`))
- `perms`     - `NULL` if there are no permutations, or (`maxOrient` - `minOrient`) permutations, one for each orientation.  A `NULL` permutation is the identity
- `rots`      - `NULL` if there are no rotations, or (`maxOrient` - `minOrient`) sets of rotations, one for each orientation.  A `NULL` set of orientations is the identity

Level: developer

-seealso: `DMLabel`, `DM`, `PetscSectionSymLabelSetStratum()`, `PetscSectionSymCreate()`, `PetscSectionSetSym()`, `PetscSectionGetPointSyms()`, `PetscSectionSymCreateLabel()`

# External Links
$(_doc_external("Dm/PetscSectionSymLabelGetStratum"))
"""
function PetscSectionSymLabelGetStratum(petsclib::PetscLibType, sym::PetscSectionSym, stratum::PetscInt) end

@for_petsc function PetscSectionSymLabelGetStratum(petsclib::$UnionPetscLib, sym::PetscSectionSym, stratum::$PetscInt )
	size_ = Ref{$PetscInt}()
	minOrient_ = Ref{$PetscInt}()
	maxOrient_ = Ref{$PetscInt}()
	perms_ = Ref{$PetscInt}()
	rots_ = Ref{$PetscScalar}()

    @chk ccall(
               (:PetscSectionSymLabelGetStratum, $petsc_library),
               PetscErrorCode,
               (PetscSectionSym, $PetscInt, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{$PetscInt}, $PetscInt, $PetscScalar),
               sym, stratum, size_, minOrient_, maxOrient_, perms_, rots_,
              )

	size = size_[]
	minOrient = minOrient_[]
	maxOrient = maxOrient_[]
	perms = perms_[]
	rots = rots_[]

	return size,minOrient,maxOrient,perms,rots
end 

"""
	PetscSectionSymLabelSetStratum(petsclib::PetscLibType,sym::PetscSectionSym, stratum::PetscInt, size::PetscInt, minOrient::PetscInt, maxOrient::PetscInt, mode::PetscCopyMode, perms::PetscInt, rots::PetscScalar) 
set the symmetries for the orientations of a stratum

Logically

Input Parameters:
- `sym`       - the section symmetries
- `stratum`   - the stratum value in the label that we are assigning symmetries for
- `size`      - the number of dofs for points in the `stratum` of the label
- `minOrient` - the smallest orientation for a point in this `stratum`
- `maxOrient` - one greater than the largest orientation for a point in this `stratum` (i.e., orientations are in the range [`minOrient`, `maxOrient`))
- `mode`      - how `sym` should copy the `perms` and `rots` arrays
- `perms`     - `NULL` if there are no permutations, or (`maxOrient` - `minOrient`) permutations, one for each orientation.  A `NULL` permutation is the identity
- `rots`      - `NULL` if there are no rotations, or (`maxOrient` - `minOrient`) sets of rotations, one for each orientation.  A `NULL` set of orientations is the identity

Level: developer

-seealso: `DMLabel`, `DM`, `PetscSectionSymLabelGetStratum()`, `PetscSectionSymCreate()`, `PetscSectionSetSym()`, `PetscSectionGetPointSyms()`, `PetscSectionSymCreateLabel()`

# External Links
$(_doc_external("Dm/PetscSectionSymLabelSetStratum"))
"""
function PetscSectionSymLabelSetStratum(petsclib::PetscLibType, sym::PetscSectionSym, stratum::PetscInt, size::PetscInt, minOrient::PetscInt, maxOrient::PetscInt, mode::PetscCopyMode, perms::PetscInt, rots::PetscScalar) end

@for_petsc function PetscSectionSymLabelSetStratum(petsclib::$UnionPetscLib, sym::PetscSectionSym, stratum::$PetscInt, size::$PetscInt, minOrient::$PetscInt, maxOrient::$PetscInt, mode::PetscCopyMode, perms::$PetscInt, rots::$PetscScalar )

    @chk ccall(
               (:PetscSectionSymLabelSetStratum, $petsc_library),
               PetscErrorCode,
               (PetscSectionSym, $PetscInt, $PetscInt, $PetscInt, $PetscInt, PetscCopyMode, $PetscInt, $PetscScalar),
               sym, stratum, size, minOrient, maxOrient, mode, perms, rots,
              )


	return nothing
end 

"""
	sym::PetscSectionSym = PetscSectionSymCreateLabel(petsclib::PetscLibType,comm::MPI_Comm, label::DMLabel) 
Create a section symmetry that assigns one symmetry to each stratum of a label

Collective

Input Parameters:
- `comm`  - the MPI communicator for the new symmetry
- `label` - the label defining the strata

Output Parameter:
- `sym` - the section symmetries

Level: developer

-seealso: `DMLabel`, `DM`, `PetscSectionSymCreate()`, `PetscSectionSetSym()`, `PetscSectionGetSym()`, `PetscSectionSymLabelSetStratum()`, `PetscSectionGetPointSyms()`

# External Links
$(_doc_external("Dm/PetscSectionSymCreateLabel"))
"""
function PetscSectionSymCreateLabel(petsclib::PetscLibType, comm::MPI_Comm, label::DMLabel) end

@for_petsc function PetscSectionSymCreateLabel(petsclib::$UnionPetscLib, comm::MPI_Comm, label::DMLabel )
	sym_ = Ref{PetscSectionSym}()

    @chk ccall(
               (:PetscSectionSymCreateLabel, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, DMLabel, Ptr{PetscSectionSym}),
               comm, label, sym_,
              )

	sym = sym_[]

	return sym
end 

