# -------------------------------------------------------
# autodefined type arguments for class ------
mutable struct _n_PetscDrawAxis end
const PetscDrawAxis = Ptr{_n_PetscDrawAxis}

mutable struct _n_PetscDrawLG end
const PetscDrawLG = Ptr{_n_PetscDrawLG}

mutable struct _n_PetscDrawSP end
const PetscDrawSP = Ptr{_n_PetscDrawSP}

mutable struct _n_PetscDrawHG end
const PetscDrawHG = Ptr{_n_PetscDrawHG}

mutable struct _n_PetscDrawBar end
const PetscDrawBar = Ptr{_n_PetscDrawBar}

# -------------------------------------------------------
"""
	PetscDrawSetSave(petsclib::PetscLibType,draw::PetscDraw, filename::String) 
Saves images produced in a `PetscDraw` into a file

Collective

Input Parameters:
- `draw`     - the graphics context
- `filename` - name of the file, if .ext then uses name of draw object plus .ext using .ext to determine the image type

Options Database Keys:
- `-draw_save <filename>`                      - filename could be name.ext or .ext (where .ext determines the type of graphics file to save, for example .png)
- `-draw_save_final_image [optional filename]` - saves the final image displayed in a window
- `-draw_save_single_file`                     - saves each new image in the same file, normally each new image is saved in a new file with filename/filename_%d.ext

Level: intermediate

-seealso: `PetscDraw`, `PetscDrawOpenX()`, `PetscDrawOpenImage()`, `PetscDrawSetFromOptions()`, `PetscDrawCreate()`, `PetscDrawDestroy()`, `PetscDrawSetSaveFinalImage()`

# External Links
$(_doc_external("Sys/PetscDrawSetSave"))
"""
function PetscDrawSetSave(petsclib::PetscLibType, draw::PetscDraw, filename::String) end

@for_petsc function PetscDrawSetSave(petsclib::$UnionPetscLib, draw::PetscDraw, filename::String )

    @chk ccall(
               (:PetscDrawSetSave, $petsc_library),
               PetscErrorCode,
               (PetscDraw, Ptr{Cchar}),
               draw, filename,
              )


	return nothing
end 

"""
	PetscDrawSetSaveMovie(petsclib::PetscLibType,draw::PetscDraw, movieext::String) 
Saves a movie produced from a `PetscDraw` into a file

Collective

Input Parameters:
- `draw`     - the graphics context
- `movieext` - optional extension defining the movie format

Options Database Key:
- `-draw_save_movie <.ext>` - saves a movie with extension .ext

Level: intermediate

-seealso: `PetscDraw`, `PetscDrawSetSave()`, `PetscDrawSetFromOptions()`, `PetscDrawCreate()`, `PetscDrawDestroy()`

# External Links
$(_doc_external("Sys/PetscDrawSetSaveMovie"))
"""
function PetscDrawSetSaveMovie(petsclib::PetscLibType, draw::PetscDraw, movieext::String) end

@for_petsc function PetscDrawSetSaveMovie(petsclib::$UnionPetscLib, draw::PetscDraw, movieext::String )

    @chk ccall(
               (:PetscDrawSetSaveMovie, $petsc_library),
               PetscErrorCode,
               (PetscDraw, Ptr{Cchar}),
               draw, movieext,
              )


	return nothing
end 

"""
	PetscDrawSetSaveFinalImage(petsclib::PetscLibType,draw::PetscDraw, filename::String) 
Saves the final image produced in a `PetscDraw` into a file

Collective

Input Parameters:
- `draw`     - the graphics context
- `filename` - name of the file, if NULL or empty uses name set with `PetscDrawSetSave()` or the name of the draw object

Options Database Key:
- `-draw_save_final_image  <filename>` - filename could be name.ext or .ext (where .ext determines the type of graphics file to save, for example .png)

Level: intermediate

-seealso: `PetscDraw`, `PetscDrawSetSave()`, `PetscDrawSetFromOptions()`, `PetscDrawCreate()`, `PetscDrawDestroy()`

# External Links
$(_doc_external("Sys/PetscDrawSetSaveFinalImage"))
"""
function PetscDrawSetSaveFinalImage(petsclib::PetscLibType, draw::PetscDraw, filename::String) end

@for_petsc function PetscDrawSetSaveFinalImage(petsclib::$UnionPetscLib, draw::PetscDraw, filename::String )

    @chk ccall(
               (:PetscDrawSetSaveFinalImage, $petsc_library),
               PetscErrorCode,
               (PetscDraw, Ptr{Cchar}),
               draw, filename,
              )


	return nothing
end 

"""
	PetscDrawSave(petsclib::PetscLibType,draw::PetscDraw) 
Saves a drawn image

Collective

Input Parameter:
- `draw` - the drawing context

Level: advanced

-seealso: `PetscDraw`, `PetscDrawSetSave()`

# External Links
$(_doc_external("Sys/PetscDrawSave"))
"""
function PetscDrawSave(petsclib::PetscLibType, draw::PetscDraw) end

@for_petsc function PetscDrawSave(petsclib::$UnionPetscLib, draw::PetscDraw )

    @chk ccall(
               (:PetscDrawSave, $petsc_library),
               PetscErrorCode,
               (PetscDraw,),
               draw,
              )


	return nothing
end 

"""
	PetscDrawSaveMovie(petsclib::PetscLibType,draw::PetscDraw) 
Saves a movie from previously saved images

Collective

Input Parameter:
- `draw` - the drawing context

Level: advanced

-seealso: `PetscDraw`, `PetscDrawSetSave()`, `PetscDrawSetSaveMovie()`

# External Links
$(_doc_external("Sys/PetscDrawSaveMovie"))
"""
function PetscDrawSaveMovie(petsclib::PetscLibType, draw::PetscDraw) end

@for_petsc function PetscDrawSaveMovie(petsclib::$UnionPetscLib, draw::PetscDraw )

    @chk ccall(
               (:PetscDrawSaveMovie, $petsc_library),
               PetscErrorCode,
               (PetscDraw,),
               draw,
              )


	return nothing
end 

"""
	PetscDrawEllipse(petsclib::PetscLibType,draw::PetscDraw, x::PetscReal, y::PetscReal, a::PetscReal, b::PetscReal, c::Cint) 
Draws an ellipse onto a drawable.

Not Collective

Input Parameters:
- `draw` - The drawing context
- `x`    - The x coordinate of the center
- `y`    - The y coordinate of the center
- `a`    - The major axes length
- `b`    - The minor axes length
- `c`    - The color

Level: beginner

-seealso: `PetscDraw`, `PetscDrawRectangle()`, `PetscDrawTriangle()`, `PetscDrawMarker()`, `PetscDrawPoint()`, `PetscDrawString()`, `PetscDrawArrow()`

# External Links
$(_doc_external("Sys/PetscDrawEllipse"))
"""
function PetscDrawEllipse(petsclib::PetscLibType, draw::PetscDraw, x::PetscReal, y::PetscReal, a::PetscReal, b::PetscReal, c::Cint) end

@for_petsc function PetscDrawEllipse(petsclib::$UnionPetscLib, draw::PetscDraw, x::$PetscReal, y::$PetscReal, a::$PetscReal, b::$PetscReal, c::Cint )

    @chk ccall(
               (:PetscDrawEllipse, $petsc_library),
               PetscErrorCode,
               (PetscDraw, $PetscReal, $PetscReal, $PetscReal, $PetscReal, Cint),
               draw, x, y, a, b, c,
              )


	return nothing
end 

"""
	PetscDrawFlush(petsclib::PetscLibType,draw::PetscDraw) 
Flushes graphical output.

Collective

Input Parameter:
- `draw` - the drawing context

Level: beginner

-seealso: `PetscDraw`, `PetscDrawClear()`

# External Links
$(_doc_external("Sys/PetscDrawFlush"))
"""
function PetscDrawFlush(petsclib::PetscLibType, draw::PetscDraw) end

@for_petsc function PetscDrawFlush(petsclib::$UnionPetscLib, draw::PetscDraw )

    @chk ccall(
               (:PetscDrawFlush, $petsc_library),
               PetscErrorCode,
               (PetscDraw,),
               draw,
              )


	return nothing
end 

"""
	x_user::PetscReal,y_user::PetscReal,x_phys::PetscReal,y_phys::PetscReal = PetscDrawGetMouseButton(petsclib::PetscLibType,draw::PetscDraw, button::PetscDrawButton) 
Returns location of mouse and which button was
pressed. Waits for button to be pressed.

Collective

Input Parameter:
- `draw` - the window to be used

Output Parameters:
- `button` - one of `PETSC_BUTTON_LEFT`, `PETSC_BUTTON_CENTER`, `PETSC_BUTTON_RIGHT`, `PETSC_BUTTON_WHEEL_UP`, `PETSC_BUTTON_WHEEL_DOWN`
- `x_user` - horizontal user coordinate of location (user may pass in NULL).
- `y_user` - vertical user coordinate of location (user may pass in NULL).
- `x_phys` - horizontal window coordinate (user may pass in NULL).
- `y_phys` - vertical window coordinate (user may pass in NULL).

-seealso: `PetscDraw`, `PetscDrawButton`

# External Links
$(_doc_external("Sys/PetscDrawGetMouseButton"))
"""
function PetscDrawGetMouseButton(petsclib::PetscLibType, draw::PetscDraw, button::PetscDrawButton) end

@for_petsc function PetscDrawGetMouseButton(petsclib::$UnionPetscLib, draw::PetscDraw, button::PetscDrawButton )
	x_user_ = Ref{$PetscReal}()
	y_user_ = Ref{$PetscReal}()
	x_phys_ = Ref{$PetscReal}()
	y_phys_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscDrawGetMouseButton, $petsc_library),
               PetscErrorCode,
               (PetscDraw, Ptr{PetscDrawButton}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
               draw, button, x_user_, y_user_, x_phys_, y_phys_,
              )

	x_user = x_user_[]
	y_user = y_user_[]
	x_phys = x_phys_[]
	y_phys = y_phys_[]

	return x_user,y_user,x_phys,y_phys
end 

"""
	PetscDrawView(petsclib::PetscLibType,indraw::PetscDraw, viewer::PetscViewer) 
Prints the `PetscDraw` data structure.

Collective

Input Parameters:
- `indraw` - the `PetscDraw` context
- `viewer` - visualization context

See PetscDrawSetFromOptions() for options database keys

-seealso: `PetscDraw`, `PetscViewerASCIIOpen()`, `PetscViewer`

# External Links
$(_doc_external("Sys/PetscDrawView"))
"""
function PetscDrawView(petsclib::PetscLibType, indraw::PetscDraw, viewer::PetscViewer) end

@for_petsc function PetscDrawView(petsclib::$UnionPetscLib, indraw::PetscDraw, viewer::PetscViewer )

    @chk ccall(
               (:PetscDrawView, $petsc_library),
               PetscErrorCode,
               (PetscDraw, PetscViewer),
               indraw, viewer,
              )


	return nothing
end 

"""
	PetscDrawViewFromOptions(petsclib::PetscLibType,A::PetscDraw, obj::PetscObject, name::String) 
View a `PetscDraw` from the option database

Collective

Input Parameters:
- `A`    - the `PetscDraw` context
- `obj`  - Optional object
- `name` - command line option

Level: intermediate

-seealso: `PetscDraw`, `PetscDrawView`, `PetscObjectViewFromOptions()`, `PetscDrawCreate()`

# External Links
$(_doc_external("Sys/PetscDrawViewFromOptions"))
"""
function PetscDrawViewFromOptions(petsclib::PetscLibType, A::PetscDraw, obj::PetscObject, name::String) end

@for_petsc function PetscDrawViewFromOptions(petsclib::$UnionPetscLib, A::PetscDraw, obj::PetscObject, name::String )

    @chk ccall(
               (:PetscDrawViewFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscDraw, PetscObject, Ptr{Cchar}),
               A, obj, name,
              )


	return nothing
end 

"""
	indraw::PetscDraw = PetscDrawCreate(petsclib::PetscLibType,comm::MPI_Comm, display::String, title::String, x::Cint, y::Cint, w::Cint, h::Cint) 
Creates a graphics context.

Collective

Input Parameters:
- `comm`    - MPI communicator
- `display` - X display when using X Windows
- `title`   - optional title added to top of window
- `x`       - horizonatl coordinate of lower left corner of window or `PETSC_DECIDE`
- `y`       - vertical coordinate of lower left corner of window or `PETSC_DECIDE`
- `w`       - width of window, `PETSC_DECIDE`, `PETSC_DRAW_HALF_SIZE`, `PETSC_DRAW_FULL_SIZE`, `PETSC_DRAW_THIRD_SIZE` or `PETSC_DRAW_QUARTER_SIZE`
- `h`       - height of window, `PETSC_DECIDE`, `PETSC_DRAW_HALF_SIZE`, `PETSC_DRAW_FULL_SIZE`, `PETSC_DRAW_THIRD_SIZE` or `PETSC_DRAW_QUARTER_SIZE`

Output Parameter:
- `indraw` - location to put the `PetscDraw` context

Level: beginner

-seealso: `PetscDrawSetType()`, `PetscDrawSetFromOptions()`, `PetscDrawDestroy()`, `PetscDrawLGCreate()`, `PetscDrawSPCreate()`,
`PetscDrawViewPortsCreate()`, `PetscDrawViewPortsSet()`, `PetscDrawAxisCreate()`, `PetscDrawHGCreate()`, `PetscDrawBarCreate()`,
`PetscViewerDrawGetDraw()`, `PetscDrawSetSave()`, `PetscDrawSetSaveMovie()`, `PetscDrawSetSaveFinalImage()`,
`PetscDrawOpenX()`, `PetscDrawOpenImage()`, `PetscDrawIsNull()`, `PetscDrawGetPopup()`, `PetscDrawCheckResizedWindow()`, `PetscDrawResizeWindow()`,
`PetscDrawGetWindowSize()`, `PetscDrawLine()`, `PetscDrawArrow()`, `PetscDrawLineSetWidth()`, `PetscDrawLineGetWidth()`, `PetscDrawMarker()`,
`PetscDrawPoint()`, `PetscDrawRectangle()`, `PetscDrawTriangle()`, `PetscDrawEllipse()`, `PetscDrawString()`, `PetscDrawStringCentered()`,
`PetscDrawStringBoxed()`, `PetscDrawStringVertical()`, `PetscDrawSetViewPort()`, `PetscDrawGetViewPort()`,
`PetscDrawSplitViewPort()`, `PetscDrawSetTitle()`, `PetscDrawAppendTitle()`, `PetscDrawGetTitle()`, `PetscDrawSetPause()`, `PetscDrawGetPause()`,
`PetscDrawPause()`, `PetscDrawSetDoubleBuffer()`, `PetscDrawClear()`, `PetscDrawFlush()`, `PetscDrawGetSingleton()`, `PetscDrawGetMouseButton()`,
`PetscDrawZoom()`, `PetscDrawGetBoundingBox()`

# External Links
$(_doc_external("Sys/PetscDrawCreate"))
"""
function PetscDrawCreate(petsclib::PetscLibType, comm::MPI_Comm, display::String, title::String, x::Cint, y::Cint, w::Cint, h::Cint) end

@for_petsc function PetscDrawCreate(petsclib::$UnionPetscLib, comm::MPI_Comm, display::String, title::String, x::Cint, y::Cint, w::Cint, h::Cint )
	indraw_ = Ref{PetscDraw}()

    @chk ccall(
               (:PetscDrawCreate, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}, Cint, Cint, Cint, Cint, Ptr{PetscDraw}),
               comm, display, title, x, y, w, h, indraw_,
              )

	indraw = indraw_[]

	return indraw
end 

"""
	PetscDrawSetType(petsclib::PetscLibType,draw::PetscDraw, type::PetscDrawType) 
Builds graphics object for a particular implementation

Collective

Input Parameters:
- `draw` - the graphics context
- `type` - for example, `PETSC_DRAW_X`

Options Database Key:
- `-draw_type  <type>` - Sets the type; use -help for a list of available methods (for instance, x)

Level: intermediate

-seealso: `PetscDraw`, `PETSC_DRAW_X`, `PETSC_DRAW_TIKZ`, `PETSC_DRAW_IMAGE`, `PetscDrawSetFromOptions()`, `PetscDrawCreate()`, `PetscDrawDestroy()`, `PetscDrawType`

# External Links
$(_doc_external("Sys/PetscDrawSetType"))
"""
function PetscDrawSetType(petsclib::PetscLibType, draw::PetscDraw, type::PetscDrawType) end

@for_petsc function PetscDrawSetType(petsclib::$UnionPetscLib, draw::PetscDraw, type::PetscDrawType )

    @chk ccall(
               (:PetscDrawSetType, $petsc_library),
               PetscErrorCode,
               (PetscDraw, PetscDrawType),
               draw, type,
              )


	return nothing
end 

"""
	type::PetscDrawType = PetscDrawGetType(petsclib::PetscLibType,draw::PetscDraw) 
Gets the `PetscDraw` type as a string from the `PetscDraw` object.

Not Collective

Input Parameter:
- `draw` - Krylov context

Output Parameter:
- `type` - name of PetscDraw method

Level: advanced

-seealso: `PetscDraw`, `PetscDrawType`, `PetscDrawSetType()`, `PetscDrawCreate()`

# External Links
$(_doc_external("Sys/PetscDrawGetType"))
"""
function PetscDrawGetType(petsclib::PetscLibType, draw::PetscDraw) end

@for_petsc function PetscDrawGetType(petsclib::$UnionPetscLib, draw::PetscDraw )
	type_ = Ref{PetscDrawType}()

    @chk ccall(
               (:PetscDrawGetType, $petsc_library),
               PetscErrorCode,
               (PetscDraw, Ptr{PetscDrawType}),
               draw, type_,
              )

	type = unsafe_string(type_[])

	return type
end 

"""
	PetscDrawRegister(petsclib::PetscLibType,sname::String, fnc::external) 
Adds a method to the graphics package.

Not Collective, No Fortran Support

Input Parameters:
- `sname`    - name of a new user-defined graphics class
- `function` - routine to create method context

Level: developer

-seealso: `PetscDraw`, `PetscDrawRegisterAll()`, `PetscDrawRegisterDestroy()`, `PetscDrawType`, `PetscDrawSetType()`

# External Links
$(_doc_external("Sys/PetscDrawRegister"))
"""
function PetscDrawRegister(petsclib::PetscLibType, sname::String, fnc::external) end

@for_petsc function PetscDrawRegister(petsclib::$UnionPetscLib, sname::String, fnc::external )

    @chk ccall(
               (:PetscDrawRegister, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, external),
               sname, fnc,
              )


	return nothing
end 

"""
	PetscDrawSetOptionsPrefix(petsclib::PetscLibType,draw::PetscDraw, prefix::String) 
Sets the prefix used for searching for all
`PetscDraw` options in the database.

Logically Collective

Input Parameters:
- `draw`   - the draw context
- `prefix` - the prefix to prepend to all option names

Level: advanced

-seealso: `PetscDraw`, `PetscDrawSetFromOptions()`, `PetscDrawCreate()`

# External Links
$(_doc_external("Sys/PetscDrawSetOptionsPrefix"))
"""
function PetscDrawSetOptionsPrefix(petsclib::PetscLibType, draw::PetscDraw, prefix::String) end

@for_petsc function PetscDrawSetOptionsPrefix(petsclib::$UnionPetscLib, draw::PetscDraw, prefix::String )

    @chk ccall(
               (:PetscDrawSetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (PetscDraw, Ptr{Cchar}),
               draw, prefix,
              )


	return nothing
end 

"""
	PetscDrawSetFromOptions(petsclib::PetscLibType,draw::PetscDraw) 
Sets the graphics type from the options database.
Defaults to a PETSc X Windows graphics.

Collective

Input Parameter:
- `draw` - the graphics context

Options Database Keys:
- `-nox`                                       - do not use X graphics (ignore graphics calls, but run program correctly)
- `-nox_warning`                               - when X Windows support is not installed this prevents the warning message from being printed
- `-draw_pause <pause amount>`                 - - -1 indicates wait for mouse input, -2 indicates pause when window is to be destroyed
- `-draw_marker_type`                          - <x,point>
- `-draw_save [optional filename]`             - (X Windows only) saves each image before it is cleared to a file
- `-draw_save_final_image [optional filename]` - (X Windows only) saves the final image displayed in a window
- `-draw_save_movie`                           - converts image files to a movie  at the end of the run. See PetscDrawSetSave()
- `-draw_save_single_file`                     - saves each new image in the same file, normally each new image is saved in a new file with 'filename/filename_%d.ext'
- `-draw_save_on_clear`                        - saves an image on each clear, mainly for debugging
- `-draw_save_on_flush`                        - saves an image on each flush, mainly for debugging

Level: intermediate

-seealso: `PetscDraw`, `PetscDrawCreate()`, `PetscDrawSetType()`, `PetscDrawSetSave()`, `PetscDrawSetSaveFinalImage()`, `PetscDrawPause()`, `PetscDrawSetPause()`

# External Links
$(_doc_external("Sys/PetscDrawSetFromOptions"))
"""
function PetscDrawSetFromOptions(petsclib::PetscLibType, draw::PetscDraw) end

@for_petsc function PetscDrawSetFromOptions(petsclib::$UnionPetscLib, draw::PetscDraw )

    @chk ccall(
               (:PetscDrawSetFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscDraw,),
               draw,
              )


	return nothing
end 

"""
	PetscDrawSetViewPort(petsclib::PetscLibType,draw::PetscDraw, xl::PetscReal, yl::PetscReal, xr::PetscReal, yr::PetscReal) 
Sets the portion of the window (page) to which draw
routines will write.

Collective

Input Parameters:
- `xl`   - the horizontal coordinate of the lower left corner of the subwindow.
- `yl`   - the vertical coordinate of the lower left corner of the subwindow.
- `xr`   - the horizontal coordinate of the upper right corner of the subwindow.
- `yr`   - the vertical coordinate of the upper right corner of the subwindow.
- `draw` - the drawing context

Level: advanced

-seealso: `PetscDrawGetViewPort()`, `PetscDraw`, `PetscDrawSplitViewPort()`, `PetscDrawViewPortsCreate()`

# External Links
$(_doc_external("Sys/PetscDrawSetViewPort"))
"""
function PetscDrawSetViewPort(petsclib::PetscLibType, draw::PetscDraw, xl::PetscReal, yl::PetscReal, xr::PetscReal, yr::PetscReal) end

@for_petsc function PetscDrawSetViewPort(petsclib::$UnionPetscLib, draw::PetscDraw, xl::$PetscReal, yl::$PetscReal, xr::$PetscReal, yr::$PetscReal )

    @chk ccall(
               (:PetscDrawSetViewPort, $petsc_library),
               PetscErrorCode,
               (PetscDraw, $PetscReal, $PetscReal, $PetscReal, $PetscReal),
               draw, xl, yl, xr, yr,
              )


	return nothing
end 

"""
	xl::PetscReal,yl::PetscReal,xr::PetscReal,yr::PetscReal = PetscDrawGetViewPort(petsclib::PetscLibType,draw::PetscDraw) 
Gets the portion of the window (page) to which draw
routines will write.

Collective

Input Parameter:
- `draw` - the drawing context

Output Parameters:
- `xl` - the horizontal coordinate of the lower left corner of the subwindow.
- `yl` - the vertical coordinate of the lower left corner of the subwindow.
- `xr` - the horizontal coordinate of the upper right corner of the subwindow.
- `yr` - the vertical coordinate of the upper right corner of the subwindow.

Level: advanced

-seealso: `PetscDraw`, `PetscDrawSplitViewPort()`, `PetscDrawSetViewPort()`

# External Links
$(_doc_external("Sys/PetscDrawGetViewPort"))
"""
function PetscDrawGetViewPort(petsclib::PetscLibType, draw::PetscDraw) end

@for_petsc function PetscDrawGetViewPort(petsclib::$UnionPetscLib, draw::PetscDraw )
	xl_ = Ref{$PetscReal}()
	yl_ = Ref{$PetscReal}()
	xr_ = Ref{$PetscReal}()
	yr_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscDrawGetViewPort, $petsc_library),
               PetscErrorCode,
               (PetscDraw, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
               draw, xl_, yl_, xr_, yr_,
              )

	xl = xl_[]
	yl = yl_[]
	xr = xr_[]
	yr = yr_[]

	return xl,yl,xr,yr
end 

"""
	newports::Vector{PetscDrawViewPorts} = PetscDrawViewPortsCreate(petsclib::PetscLibType,draw::PetscDraw, nports::PetscInt) 
Splits a window into smaller view ports. Each processor shares all the viewports.

Collective

Input Parameters:
- `draw`   - the drawing context
- `nports` - the number of ports

Output Parameter:
- `newports` - a `PetscDrawViewPorts` context (C structure)

Options Database Key:
- `-draw_ports` - display multiple fields in the same window with PetscDrawPorts() instead of in separate windows

Level: advanced

-seealso: `PetscDrawSplitViewPort()`, `PetscDrawSetViewPort()`, `PetscDrawViewPortsSet()`, `PetscDrawViewPortsDestroy()`

# External Links
$(_doc_external("Sys/PetscDrawViewPortsCreate"))
"""
function PetscDrawViewPortsCreate(petsclib::PetscLibType, draw::PetscDraw, nports::PetscInt) end

@for_petsc function PetscDrawViewPortsCreate(petsclib::$UnionPetscLib, draw::PetscDraw, nports::$PetscInt )
	newports_ = Ref{Ptr{PetscDrawViewPorts}}()

    @chk ccall(
               (:PetscDrawViewPortsCreate, $petsc_library),
               PetscErrorCode,
               (PetscDraw, $PetscInt, Ptr{Ptr{PetscDrawViewPorts}}),
               draw, nports, newports_,
              )

	newports = unsafe_wrap(Array, newports_[], VecGetLocalSize(petsclib, x); own = false)

	return newports
end 

"""
	newports::Vector{PetscDrawViewPorts} = PetscDrawViewPortsCreateRect(petsclib::PetscLibType,draw::PetscDraw, nx::PetscInt, ny::PetscInt) 
Splits a window into smaller
view ports. Each processor shares all the viewports. The number
of views in the x- and y-directions is specified.

Collective

Input Parameters:
- `draw` - the drawing context
- `nx`   - the number of x divisions
- `ny`   - the number of y divisions

Output Parameter:
- `newports` - a `PetscDrawViewPorts` context (C structure)

Level: advanced

-seealso: `PetscDrawSplitViewPort()`, `PetscDrawSetViewPort()`, `PetscDrawViewPortsSet()`, `PetscDrawViewPortsDestroy()`, `PetscDrawViewPorts`

# External Links
$(_doc_external("Sys/PetscDrawViewPortsCreateRect"))
"""
function PetscDrawViewPortsCreateRect(petsclib::PetscLibType, draw::PetscDraw, nx::PetscInt, ny::PetscInt) end

@for_petsc function PetscDrawViewPortsCreateRect(petsclib::$UnionPetscLib, draw::PetscDraw, nx::$PetscInt, ny::$PetscInt )
	newports_ = Ref{Ptr{PetscDrawViewPorts}}()

    @chk ccall(
               (:PetscDrawViewPortsCreateRect, $petsc_library),
               PetscErrorCode,
               (PetscDraw, $PetscInt, $PetscInt, Ptr{Ptr{PetscDrawViewPorts}}),
               draw, nx, ny, newports_,
              )

	newports = unsafe_wrap(Array, newports_[], VecGetLocalSize(petsclib, x); own = false)

	return newports
end 

"""
	PetscDrawViewPortsDestroy(petsclib::PetscLibType,ports::PetscDrawViewPorts) 
frees a `PetscDrawViewPorts` object

Collective on the `PetscDraw` inside `ports`

Input Parameter:
- `ports` - the `PetscDrawViewPorts` object

Level: advanced

-seealso: `PetscDrawViewPorts`, `PetscDrawSplitViewPort()`, `PetscDrawSetViewPort()`, `PetscDrawViewPortsSet()`, `PetscDrawViewPortsCreate()`

# External Links
$(_doc_external("Sys/PetscDrawViewPortsDestroy"))
"""
function PetscDrawViewPortsDestroy(petsclib::PetscLibType, ports::PetscDrawViewPorts) end

@for_petsc function PetscDrawViewPortsDestroy(petsclib::$UnionPetscLib, ports::PetscDrawViewPorts )

    @chk ccall(
               (:PetscDrawViewPortsDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscDrawViewPorts},),
               ports,
              )


	return nothing
end 

"""
	PetscDrawViewPortsSet(petsclib::PetscLibType,ports::PetscDrawViewPorts, port::PetscInt) 
sets a draw object to use a particular subport

Logically Collective on the `PetscDraw` inside `ports`

Input Parameters:
- `ports` - the `PetscDrawViewPorts` object
- `port`  - the port number, from 0 to nports-1

Level: advanced

-seealso: `PetscDrawViewPorts`, `PetscDrawSplitViewPort()`, `PetscDrawSetViewPort()`, `PetscDrawViewPortsDestroy()`, `PetscDrawViewPortsCreate()`

# External Links
$(_doc_external("Sys/PetscDrawViewPortsSet"))
"""
function PetscDrawViewPortsSet(petsclib::PetscLibType, ports::PetscDrawViewPorts, port::PetscInt) end

@for_petsc function PetscDrawViewPortsSet(petsclib::$UnionPetscLib, ports::PetscDrawViewPorts, port::$PetscInt )

    @chk ccall(
               (:PetscDrawViewPortsSet, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscDrawViewPorts}, $PetscInt),
               ports, port,
              )


	return nothing
end 

"""
	PetscDrawString(petsclib::PetscLibType,draw::PetscDraw, xl::PetscReal, yl::PetscReal, cl::Cint, text::String) 
draws text onto a drawable.

Not Collective

Input Parameters:
- `draw` - the drawing context
- `xl`   - coordinate of lower left corner of text
- `yl`   - coordinate of lower left corner of text
- `cl`   - the color of the text
- `text` - the text to draw

Level: beginner

-seealso: `PetscDraw`, `PetscDrawStringVertical()`, `PetscDrawStringCentered()`, `PetscDrawStringBoxed()`, `PetscDrawStringSetSize()`,
`PetscDrawStringGetSize()`, `PetscDrawLine()`, `PetscDrawRectangle()`, `PetscDrawTriangle()`, `PetscDrawEllipse()`,
`PetscDrawMarker()`, `PetscDrawPoint()`

# External Links
$(_doc_external("Sys/PetscDrawString"))
"""
function PetscDrawString(petsclib::PetscLibType, draw::PetscDraw, xl::PetscReal, yl::PetscReal, cl::Cint, text::String) end

@for_petsc function PetscDrawString(petsclib::$UnionPetscLib, draw::PetscDraw, xl::$PetscReal, yl::$PetscReal, cl::Cint, text::String )

    @chk ccall(
               (:PetscDrawString, $petsc_library),
               PetscErrorCode,
               (PetscDraw, $PetscReal, $PetscReal, Cint, Ptr{Cchar}),
               draw, xl, yl, cl, text,
              )


	return nothing
end 

"""
	PetscDrawStringVertical(petsclib::PetscLibType,draw::PetscDraw, xl::PetscReal, yl::PetscReal, cl::Cint, text::String) 
draws text onto a drawable.

Not Collective

Input Parameters:
- `draw` - the drawing context
- `xl`   - coordinate of upper left corner of text
- `yl`   - coordinate of upper left corner of text
- `cl`   - the color of the text
- `text` - the text to draw

Level: beginner

-seealso: `PetscDraw`, `PetscDrawString()`, `PetscDrawStringCentered()`, `PetscDrawStringBoxed()`, `PetscDrawStringSetSize()`,
`PetscDrawStringGetSize()`

# External Links
$(_doc_external("Sys/PetscDrawStringVertical"))
"""
function PetscDrawStringVertical(petsclib::PetscLibType, draw::PetscDraw, xl::PetscReal, yl::PetscReal, cl::Cint, text::String) end

@for_petsc function PetscDrawStringVertical(petsclib::$UnionPetscLib, draw::PetscDraw, xl::$PetscReal, yl::$PetscReal, cl::Cint, text::String )

    @chk ccall(
               (:PetscDrawStringVertical, $petsc_library),
               PetscErrorCode,
               (PetscDraw, $PetscReal, $PetscReal, Cint, Ptr{Cchar}),
               draw, xl, yl, cl, text,
              )


	return nothing
end 

"""
	PetscDrawStringCentered(petsclib::PetscLibType,draw::PetscDraw, xc::PetscReal, yl::PetscReal, cl::Cint, text::String) 
draws text onto a drawable centered at a point

Not Collective

Input Parameters:
- `draw` - the drawing context
- `xc`   - the coordinates of right-left center of text
- `yl`   - the coordinates of lower edge of text
- `cl`   - the color of the text
- `text` - the text to draw

Level: beginner

-seealso: `PetscDraw`, `PetscDrawStringVertical()`, `PetscDrawString()`, `PetscDrawStringBoxed()`, `PetscDrawStringSetSize()`,
`PetscDrawStringGetSize()`

# External Links
$(_doc_external("Sys/PetscDrawStringCentered"))
"""
function PetscDrawStringCentered(petsclib::PetscLibType, draw::PetscDraw, xc::PetscReal, yl::PetscReal, cl::Cint, text::String) end

@for_petsc function PetscDrawStringCentered(petsclib::$UnionPetscLib, draw::PetscDraw, xc::$PetscReal, yl::$PetscReal, cl::Cint, text::String )

    @chk ccall(
               (:PetscDrawStringCentered, $petsc_library),
               PetscErrorCode,
               (PetscDraw, $PetscReal, $PetscReal, Cint, Ptr{Cchar}),
               draw, xc, yl, cl, text,
              )


	return nothing
end 

"""
	w::PetscReal,h::PetscReal = PetscDrawStringBoxed(petsclib::PetscLibType,draw::PetscDraw, sxl::PetscReal, syl::PetscReal, sc::Cint, bc::Cint, text::String) 
Draws a string with a box around it

Not Collective

Input Parameters:
- `draw` - the drawing context
- `sxl`  - the coordinates of center of the box
- `syl`  - the coordinates of top line of box
- `sc`   - the color of the text
- `bc`   - the color of the bounding box
- `text` - the text to draw

Output Parameters:
- `w` - the width of the resulting box (optional)
- `h` - the height of resulting box (optional)

Level: beginner

-seealso: `PetscDraw`, `PetscDrawStringVertical()`, `PetscDrawString()`, `PetscDrawStringCentered()`, `PetscDrawStringSetSize()`,
`PetscDrawStringGetSize()`

# External Links
$(_doc_external("Sys/PetscDrawStringBoxed"))
"""
function PetscDrawStringBoxed(petsclib::PetscLibType, draw::PetscDraw, sxl::PetscReal, syl::PetscReal, sc::Cint, bc::Cint, text::String) end

@for_petsc function PetscDrawStringBoxed(petsclib::$UnionPetscLib, draw::PetscDraw, sxl::$PetscReal, syl::$PetscReal, sc::Cint, bc::Cint, text::String )
	w_ = Ref{$PetscReal}()
	h_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscDrawStringBoxed, $petsc_library),
               PetscErrorCode,
               (PetscDraw, $PetscReal, $PetscReal, Cint, Cint, Ptr{Cchar}, Ptr{$PetscReal}, Ptr{$PetscReal}),
               draw, sxl, syl, sc, bc, text, w_, h_,
              )

	w = w_[]
	h = h_[]

	return w,h
end 

"""
	PetscDrawStringSetSize(petsclib::PetscLibType,draw::PetscDraw, width::PetscReal, height::PetscReal) 
Sets the size for character text.

Not Collective

Input Parameters:
- `draw`   - the drawing context
- `width`  - the width in user coordinates
- `height` - the character height in user coordinates

Level: advanced

-seealso: `PetscDraw`, `PetscDrawStringVertical()`, `PetscDrawString()`, `PetscDrawStringCentered()`, `PetscDrawStringBoxed()`,
`PetscDrawStringGetSize()`

# External Links
$(_doc_external("Sys/PetscDrawStringSetSize"))
"""
function PetscDrawStringSetSize(petsclib::PetscLibType, draw::PetscDraw, width::PetscReal, height::PetscReal) end

@for_petsc function PetscDrawStringSetSize(petsclib::$UnionPetscLib, draw::PetscDraw, width::$PetscReal, height::$PetscReal )

    @chk ccall(
               (:PetscDrawStringSetSize, $petsc_library),
               PetscErrorCode,
               (PetscDraw, $PetscReal, $PetscReal),
               draw, width, height,
              )


	return nothing
end 

"""
	PetscDrawStringGetSize(petsclib::PetscLibType,draw::PetscDraw, width::PetscReal, height::PetscReal) 
Gets the size for character text.  The width is
relative to the user coordinates of the window.

Not Collective

Input Parameters:
- `draw`   - the drawing context
- `width`  - the width in user coordinates
- `height` - the character height

Level: advanced

-seealso: `PetscDraw`, `PetscDrawStringVertical()`, `PetscDrawString()`, `PetscDrawStringCentered()`, `PetscDrawStringBoxed()`,
`PetscDrawStringSetSize()`

# External Links
$(_doc_external("Sys/PetscDrawStringGetSize"))
"""
function PetscDrawStringGetSize(petsclib::PetscLibType, draw::PetscDraw, width::PetscReal, height::PetscReal) end

@for_petsc function PetscDrawStringGetSize(petsclib::$UnionPetscLib, draw::PetscDraw, width::$PetscReal, height::$PetscReal )

    @chk ccall(
               (:PetscDrawStringGetSize, $petsc_library),
               PetscErrorCode,
               (PetscDraw, Ptr{$PetscReal}, Ptr{$PetscReal}),
               draw, width, height,
              )


	return nothing
end 

"""
	PetscDrawFinalizePackage(petsclib::PetscLibType) 
This function destroys everything in the PETSc interface to the `PetscDraw` package. It is
called from `PetscFinalize()`.

Level: developer

-seealso: `PetscDraw`, `PetscFinalize()`

# External Links
$(_doc_external("Sys/PetscDrawFinalizePackage"))
"""
function PetscDrawFinalizePackage(petsclib::PetscLibType) end

@for_petsc function PetscDrawFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscDrawFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscDrawInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the `PetscDraw` package. It is called
from PetscDLLibraryRegister_petsc() when using dynamic libraries, and on the call to `PetscInitialize()`
when using shared or static libraries.

Level: developer

-seealso: `PetscDraw`, `PetscInitialize()`

# External Links
$(_doc_external("Sys/PetscDrawInitializePackage"))
"""
function PetscDrawInitializePackage(petsclib::PetscLibType) end

@for_petsc function PetscDrawInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscDrawInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscDrawResizeWindow(petsclib::PetscLibType,draw::PetscDraw, w::Cint, h::Cint) 
Allows one to resize a window from a program.

Collective

Input Parameters:
- `draw` - the window
- `w`    - the new width of the window
- `h`    - the new height of the window

Level: intermediate

-seealso: `PetscDraw`, `PetscDrawCheckResizedWindow()`

# External Links
$(_doc_external("Sys/PetscDrawResizeWindow"))
"""
function PetscDrawResizeWindow(petsclib::PetscLibType, draw::PetscDraw, w::Cint, h::Cint) end

@for_petsc function PetscDrawResizeWindow(petsclib::$UnionPetscLib, draw::PetscDraw, w::Cint, h::Cint )

    @chk ccall(
               (:PetscDrawResizeWindow, $petsc_library),
               PetscErrorCode,
               (PetscDraw, Cint, Cint),
               draw, w, h,
              )


	return nothing
end 

"""
	PetscDrawGetWindowSize(petsclib::PetscLibType,draw::PetscDraw, w::Cint, h::Cint) 
Gets the size of the window.

Not Collective

Input Parameter:
- `draw` - the window

Output Parameters:
- `w` - the window width
- `h` - the window height

Level: intermediate

-seealso: `PetscDraw`, `PetscDrawResizeWindow()`, `PetscDrawCheckResizedWindow()`

# External Links
$(_doc_external("Sys/PetscDrawGetWindowSize"))
"""
function PetscDrawGetWindowSize(petsclib::PetscLibType, draw::PetscDraw, w::Cint, h::Cint) end

@for_petsc function PetscDrawGetWindowSize(petsclib::$UnionPetscLib, draw::PetscDraw, w::Cint, h::Cint )

    @chk ccall(
               (:PetscDrawGetWindowSize, $petsc_library),
               PetscErrorCode,
               (PetscDraw, Ptr{Cint}, Ptr{Cint}),
               draw, w, h,
              )


	return nothing
end 

"""
	PetscDrawCheckResizedWindow(petsclib::PetscLibType,draw::PetscDraw) 
Checks if the user has resized the window.

Collective

Input Parameter:
- `draw` - the window

Level: advanced

-seealso: `PetscDraw`, `PetscDrawResizeWindow()`

# External Links
$(_doc_external("Sys/PetscDrawCheckResizedWindow"))
"""
function PetscDrawCheckResizedWindow(petsclib::PetscLibType, draw::PetscDraw) end

@for_petsc function PetscDrawCheckResizedWindow(petsclib::$UnionPetscLib, draw::PetscDraw )

    @chk ccall(
               (:PetscDrawCheckResizedWindow, $petsc_library),
               PetscErrorCode,
               (PetscDraw,),
               draw,
              )


	return nothing
end 

"""
	PetscDrawGetTitle(petsclib::PetscLibType,draw::PetscDraw, title::String) 
Gets pointer to title of a `PetscDraw` context.

Not Collective

Input Parameter:
- `draw` - the graphics context

Output Parameter:
- `title` - the title

Level: intermediate

-seealso: `PetscDraw`, `PetscDrawSetTitle()`

# External Links
$(_doc_external("Sys/PetscDrawGetTitle"))
"""
function PetscDrawGetTitle(petsclib::PetscLibType, draw::PetscDraw, title::String) end

@for_petsc function PetscDrawGetTitle(petsclib::$UnionPetscLib, draw::PetscDraw, title::String )
	title_ = Ref(pointer(title))

    @chk ccall(
               (:PetscDrawGetTitle, $petsc_library),
               PetscErrorCode,
               (PetscDraw, Ptr{Ptr{Cchar}}),
               draw, title_,
              )


	return nothing
end 

"""
	PetscDrawSetTitle(petsclib::PetscLibType,draw::PetscDraw, title::String) 
Sets the title of a `PetscDraw` context.

Collective

Input Parameters:
- `draw`  - the graphics context
- `title` - the title

Level: intermediate

-seealso: `PetscDraw`, `PetscDrawGetTitle()`, `PetscDrawAppendTitle()`

# External Links
$(_doc_external("Sys/PetscDrawSetTitle"))
"""
function PetscDrawSetTitle(petsclib::PetscLibType, draw::PetscDraw, title::String) end

@for_petsc function PetscDrawSetTitle(petsclib::$UnionPetscLib, draw::PetscDraw, title::String )

    @chk ccall(
               (:PetscDrawSetTitle, $petsc_library),
               PetscErrorCode,
               (PetscDraw, Ptr{Cchar}),
               draw, title,
              )


	return nothing
end 

"""
	PetscDrawAppendTitle(petsclib::PetscLibType,draw::PetscDraw, title::String) 
Appends to the title of a `PetscDraw` context.

Collective

Input Parameters:
- `draw`  - the graphics context
- `title` - the title

Level: advanced

-seealso: `PetscDraw`, `PetscDrawSetTitle()`, `PetscDrawGetTitle()`

# External Links
$(_doc_external("Sys/PetscDrawAppendTitle"))
"""
function PetscDrawAppendTitle(petsclib::PetscLibType, draw::PetscDraw, title::String) end

@for_petsc function PetscDrawAppendTitle(petsclib::$UnionPetscLib, draw::PetscDraw, title::String )

    @chk ccall(
               (:PetscDrawAppendTitle, $petsc_library),
               PetscErrorCode,
               (PetscDraw, Ptr{Cchar}),
               draw, title,
              )


	return nothing
end 

"""
	PetscDrawDestroy(petsclib::PetscLibType,draw::PetscDraw) 
Deletes a draw context.

Collective

Input Parameter:
- `draw` - the drawing context

Level: beginner

-seealso: `PetscDraw`, `PetscDrawCreate()`

# External Links
$(_doc_external("Sys/PetscDrawDestroy"))
"""
function PetscDrawDestroy(petsclib::PetscLibType, draw::PetscDraw) end

@for_petsc function PetscDrawDestroy(petsclib::$UnionPetscLib, draw::PetscDraw )

    @chk ccall(
               (:PetscDrawDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscDraw},),
               draw,
              )


	return nothing
end 

"""
	PetscDrawGetPopup(petsclib::PetscLibType,draw::PetscDraw, popup::PetscDraw) 
Creates a popup window associated with a `PetscDraw` window.

Collective

Input Parameter:
- `draw` - the original window

Output Parameter:
- `popup` - the new popup window

Level: advanced

-seealso: `PetscDraw`, `PetscDrawScalePopup()`, `PetscDrawCreate()`

# External Links
$(_doc_external("Sys/PetscDrawGetPopup"))
"""
function PetscDrawGetPopup(petsclib::PetscLibType, draw::PetscDraw, popup::PetscDraw) end

@for_petsc function PetscDrawGetPopup(petsclib::$UnionPetscLib, draw::PetscDraw, popup::PetscDraw )

    @chk ccall(
               (:PetscDrawGetPopup, $petsc_library),
               PetscErrorCode,
               (PetscDraw, Ptr{PetscDraw}),
               draw, popup,
              )


	return nothing
end 

"""
	PetscDrawSetDisplay(petsclib::PetscLibType,draw::PetscDraw, display::String) 
Sets the display where a `PetscDraw` object will be displayed

Input Parameters:
- `draw`    - the drawing context
- `display` - the X windows display

Level: advanced

-seealso: `PetscDraw`, `PetscDrawOpenX()`, `PetscDrawCreate()`

# External Links
$(_doc_external("Sys/PetscDrawSetDisplay"))
"""
function PetscDrawSetDisplay(petsclib::PetscLibType, draw::PetscDraw, display::String) end

@for_petsc function PetscDrawSetDisplay(petsclib::$UnionPetscLib, draw::PetscDraw, display::String )

    @chk ccall(
               (:PetscDrawSetDisplay, $petsc_library),
               PetscErrorCode,
               (PetscDraw, Ptr{Cchar}),
               draw, display,
              )


	return nothing
end 

"""
	PetscDrawSetDoubleBuffer(petsclib::PetscLibType,draw::PetscDraw) 
Sets a window to be double buffered.

Logically Collective

Input Parameter:
- `draw` - the drawing context

Level: intermediate

-seealso: `PetscDraw`, `PetscDrawOpenX()`, `PetscDrawCreate()`

# External Links
$(_doc_external("Sys/PetscDrawSetDoubleBuffer"))
"""
function PetscDrawSetDoubleBuffer(petsclib::PetscLibType, draw::PetscDraw) end

@for_petsc function PetscDrawSetDoubleBuffer(petsclib::$UnionPetscLib, draw::PetscDraw )

    @chk ccall(
               (:PetscDrawSetDoubleBuffer, $petsc_library),
               PetscErrorCode,
               (PetscDraw,),
               draw,
              )


	return nothing
end 

"""
	PetscDrawGetSingleton(petsclib::PetscLibType,draw::PetscDraw, sdraw::PetscDraw) 
Gain access to a `PetscDraw` object as if it were owned
by the one process.

Collective

Input Parameter:
- `draw` - the original window

Output Parameter:
- `sdraw` - the singleton window

Level: advanced

-seealso: `PetscDraw`, `PetscDrawRestoreSingleton()`, `PetscViewerGetSingleton()`, `PetscViewerRestoreSingleton()`

# External Links
$(_doc_external("Sys/PetscDrawGetSingleton"))
"""
function PetscDrawGetSingleton(petsclib::PetscLibType, draw::PetscDraw, sdraw::PetscDraw) end

@for_petsc function PetscDrawGetSingleton(petsclib::$UnionPetscLib, draw::PetscDraw, sdraw::PetscDraw )

    @chk ccall(
               (:PetscDrawGetSingleton, $petsc_library),
               PetscErrorCode,
               (PetscDraw, Ptr{PetscDraw}),
               draw, sdraw,
              )


	return nothing
end 

"""
	PetscDrawRestoreSingleton(petsclib::PetscLibType,draw::PetscDraw, sdraw::PetscDraw) 
Remove access to a `PetscDraw` object obtained with `PetscDrawGetSingleton()`
by the one process.

Collective

Input Parameters:
- `draw`  - the original window
- `sdraw` - the singleton window

Level: advanced

-seealso: `PetscDraw`, `PetscDrawGetSingleton()`, `PetscViewerGetSingleton()`, `PetscViewerRestoreSingleton()`

# External Links
$(_doc_external("Sys/PetscDrawRestoreSingleton"))
"""
function PetscDrawRestoreSingleton(petsclib::PetscLibType, draw::PetscDraw, sdraw::PetscDraw) end

@for_petsc function PetscDrawRestoreSingleton(petsclib::$UnionPetscLib, draw::PetscDraw, sdraw::PetscDraw )

    @chk ccall(
               (:PetscDrawRestoreSingleton, $petsc_library),
               PetscErrorCode,
               (PetscDraw, Ptr{PetscDraw}),
               draw, sdraw,
              )


	return nothing
end 

"""
	PetscDrawSetVisible(petsclib::PetscLibType,draw::PetscDraw, visible::PetscBool) 
Sets if the drawing surface (the 'window') is visible on its display.

Input Parameters:
- `draw`    - the drawing window
- `visible` - if the surface should be visible

Level: intermediate

-seealso: `PetscDraw`

# External Links
$(_doc_external("Sys/PetscDrawSetVisible"))
"""
function PetscDrawSetVisible(petsclib::PetscLibType, draw::PetscDraw, visible::PetscBool) end

@for_petsc function PetscDrawSetVisible(petsclib::$UnionPetscLib, draw::PetscDraw, visible::PetscBool )

    @chk ccall(
               (:PetscDrawSetVisible, $petsc_library),
               PetscErrorCode,
               (PetscDraw, PetscBool),
               draw, visible,
              )


	return nothing
end 

"""
	PetscDrawSetCoordinates(petsclib::PetscLibType,draw::PetscDraw, xl::PetscReal, yl::PetscReal, xr::PetscReal, yr::PetscReal) 
Sets the application coordinates of the corners of
the window (or page).

Not Collective

Input Parameters:
- `draw` - the drawing object
- `xl`   - the lower left x coordinate
- `yl`   - the lower left y coordinate
- `xr`   - the upper right x coordinate
- `yr`   - the upper right y coordinate

Level: advanced

-seealso: `PetscDraw`, `PetscDrawGetCoordinates()`

# External Links
$(_doc_external("Sys/PetscDrawSetCoordinates"))
"""
function PetscDrawSetCoordinates(petsclib::PetscLibType, draw::PetscDraw, xl::PetscReal, yl::PetscReal, xr::PetscReal, yr::PetscReal) end

@for_petsc function PetscDrawSetCoordinates(petsclib::$UnionPetscLib, draw::PetscDraw, xl::$PetscReal, yl::$PetscReal, xr::$PetscReal, yr::$PetscReal )

    @chk ccall(
               (:PetscDrawSetCoordinates, $petsc_library),
               PetscErrorCode,
               (PetscDraw, $PetscReal, $PetscReal, $PetscReal, $PetscReal),
               draw, xl, yl, xr, yr,
              )


	return nothing
end 

"""
	xl::PetscReal,yl::PetscReal,xr::PetscReal,yr::PetscReal = PetscDrawGetCoordinates(petsclib::PetscLibType,draw::PetscDraw) 
Gets the application coordinates of the corners of
the window (or page).

Not Collective

Input Parameter:
- `draw` - the drawing object

Output Parameters:
- `xl` - the horizontal coordinate of the lower left corner of the drawing region.
- `yl` - the vertical coordinate of the lower left corner of the drawing region.
- `xr` - the horizontal coordinate of the upper right corner of the drawing region.
- `yr` - the vertical coordinate of the upper right corner of the drawing region.

Level: advanced

-seealso: `PetscDraw`, `PetscDrawSetCoordinates()`

# External Links
$(_doc_external("Sys/PetscDrawGetCoordinates"))
"""
function PetscDrawGetCoordinates(petsclib::PetscLibType, draw::PetscDraw) end

@for_petsc function PetscDrawGetCoordinates(petsclib::$UnionPetscLib, draw::PetscDraw )
	xl_ = Ref{$PetscReal}()
	yl_ = Ref{$PetscReal}()
	xr_ = Ref{$PetscReal}()
	yr_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscDrawGetCoordinates, $petsc_library),
               PetscErrorCode,
               (PetscDraw, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
               draw, xl_, yl_, xr_, yr_,
              )

	xl = xl_[]
	yl = yl_[]
	xr = xr_[]
	yr = yr_[]

	return xl,yl,xr,yr
end 

"""
	xl::PetscReal,yl::PetscReal,xr::PetscReal,yr::PetscReal = PetscDrawGetBoundingBox(petsclib::PetscLibType,draw::PetscDraw) 
Gets the bounding box of all `PetscDrawStringBoxed()` commands

Not Collective

Input Parameter:
- `draw` - the drawing context

Output Parameters:
- `xl` - horizontal coordinate of lower left corner of bounding box
- `yl` - vertical coordinate of lower left corner of bounding box
- `xr` - horizontal coordinate of upper right corner of bounding box
- `yr` - vertical coordinate of upper right corner of bounding box

Level: intermediate

-seealso: `PetscDraw`, `PetscDrawPushCurrentPoint()`, `PetscDrawPopCurrentPoint()`, `PetscDrawSetCurrentPoint()`

# External Links
$(_doc_external("Sys/PetscDrawGetBoundingBox"))
"""
function PetscDrawGetBoundingBox(petsclib::PetscLibType, draw::PetscDraw) end

@for_petsc function PetscDrawGetBoundingBox(petsclib::$UnionPetscLib, draw::PetscDraw )
	xl_ = Ref{$PetscReal}()
	yl_ = Ref{$PetscReal}()
	xr_ = Ref{$PetscReal}()
	yr_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscDrawGetBoundingBox, $petsc_library),
               PetscErrorCode,
               (PetscDraw, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
               draw, xl_, yl_, xr_, yr_,
              )

	xl = xl_[]
	yl = yl_[]
	xr = xr_[]
	yr = yr_[]

	return xl,yl,xr,yr
end 

"""
	x::PetscReal,y::PetscReal = PetscDrawGetCurrentPoint(petsclib::PetscLibType,draw::PetscDraw) 
Gets the current draw point, some codes use this point to determine where to draw next

Not Collective

Input Parameter:
- `draw` - the drawing context

Output Parameters:
- `x` - horizontal coordinate of the current point
- `y` - vertical coordinate of the current point

Level: intermediate

-seealso: `PetscDraw`, `PetscDrawPushCurrentPoint()`, `PetscDrawPopCurrentPoint()`, `PetscDrawSetCurrentPoint()`

# External Links
$(_doc_external("Sys/PetscDrawGetCurrentPoint"))
"""
function PetscDrawGetCurrentPoint(petsclib::PetscLibType, draw::PetscDraw) end

@for_petsc function PetscDrawGetCurrentPoint(petsclib::$UnionPetscLib, draw::PetscDraw )
	x_ = Ref{$PetscReal}()
	y_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscDrawGetCurrentPoint, $petsc_library),
               PetscErrorCode,
               (PetscDraw, Ptr{$PetscReal}, Ptr{$PetscReal}),
               draw, x_, y_,
              )

	x = x_[]
	y = y_[]

	return x,y
end 

"""
	PetscDrawSetCurrentPoint(petsclib::PetscLibType,draw::PetscDraw, x::PetscReal, y::PetscReal) 
Sets the current draw point, some codes use this point to determine where to draw next

Not Collective

Input Parameters:
- `draw` - the drawing context
- `x`    - horizontal coordinate of the current point
- `y`    - vertical coordinate of the current point

Level: intermediate

-seealso: `PetscDraw`, `PetscDrawPushCurrentPoint()`, `PetscDrawPopCurrentPoint()`, `PetscDrawGetCurrentPoint()`

# External Links
$(_doc_external("Sys/PetscDrawSetCurrentPoint"))
"""
function PetscDrawSetCurrentPoint(petsclib::PetscLibType, draw::PetscDraw, x::PetscReal, y::PetscReal) end

@for_petsc function PetscDrawSetCurrentPoint(petsclib::$UnionPetscLib, draw::PetscDraw, x::$PetscReal, y::$PetscReal )

    @chk ccall(
               (:PetscDrawSetCurrentPoint, $petsc_library),
               PetscErrorCode,
               (PetscDraw, $PetscReal, $PetscReal),
               draw, x, y,
              )


	return nothing
end 

"""
	PetscDrawPushCurrentPoint(petsclib::PetscLibType,draw::PetscDraw, x::PetscReal, y::PetscReal) 
Pushes a new current draw point, retaining the old one, some codes use this point to determine where to draw next

Not Collective

Input Parameters:
- `draw` - the drawing context
- `x`    - horizontal coordinate of the current point
- `y`    - vertical coordinate of the current point

Level: intermediate

-seealso: `PetscDraw`, `PetscDrawPopCurrentPoint()`, `PetscDrawGetCurrentPoint()`

# External Links
$(_doc_external("Sys/PetscDrawPushCurrentPoint"))
"""
function PetscDrawPushCurrentPoint(petsclib::PetscLibType, draw::PetscDraw, x::PetscReal, y::PetscReal) end

@for_petsc function PetscDrawPushCurrentPoint(petsclib::$UnionPetscLib, draw::PetscDraw, x::$PetscReal, y::$PetscReal )

    @chk ccall(
               (:PetscDrawPushCurrentPoint, $petsc_library),
               PetscErrorCode,
               (PetscDraw, $PetscReal, $PetscReal),
               draw, x, y,
              )


	return nothing
end 

"""
	PetscDrawPopCurrentPoint(petsclib::PetscLibType,draw::PetscDraw) 
Pops a current draw point (discarding it)

Not Collective

Input Parameter:
- `draw` - the drawing context

Level: intermediate

-seealso: `PetscDraw`, `PetscDrawPushCurrentPoint()`, `PetscDrawSetCurrentPoint()`, `PetscDrawGetCurrentPoint()`

# External Links
$(_doc_external("Sys/PetscDrawPopCurrentPoint"))
"""
function PetscDrawPopCurrentPoint(petsclib::PetscLibType, draw::PetscDraw) end

@for_petsc function PetscDrawPopCurrentPoint(petsclib::$UnionPetscLib, draw::PetscDraw )

    @chk ccall(
               (:PetscDrawPopCurrentPoint, $petsc_library),
               PetscErrorCode,
               (PetscDraw,),
               draw,
              )


	return nothing
end 

"""
	PetscDrawLine(petsclib::PetscLibType,draw::PetscDraw, xl::PetscReal, yl::PetscReal, xr::PetscReal, yr::PetscReal, cl::Cint) 
draws a line onto a drawable.

Not Collective

Input Parameters:
- `draw` - the drawing context
- `xl`   - horizontal coordinate of first end point
- `yl`   - vertical coordinate of first end point
- `xr`   - horizontal coordinate of second end point
- `yr`   - vertical coordinate of second end point
- `cl`   - the colors of the endpoints

Level: beginner

-seealso: `PetscDraw`, `PetscDrawArrow()`, `PetscDrawLineSetWidth()`, `PetscDrawLineGetWidth()`, `PetscDrawRectangle()`, `PetscDrawTriangle()`, `PetscDrawEllipse()`,
`PetscDrawMarker()`, `PetscDrawPoint()`

# External Links
$(_doc_external("Sys/PetscDrawLine"))
"""
function PetscDrawLine(petsclib::PetscLibType, draw::PetscDraw, xl::PetscReal, yl::PetscReal, xr::PetscReal, yr::PetscReal, cl::Cint) end

@for_petsc function PetscDrawLine(petsclib::$UnionPetscLib, draw::PetscDraw, xl::$PetscReal, yl::$PetscReal, xr::$PetscReal, yr::$PetscReal, cl::Cint )

    @chk ccall(
               (:PetscDrawLine, $petsc_library),
               PetscErrorCode,
               (PetscDraw, $PetscReal, $PetscReal, $PetscReal, $PetscReal, Cint),
               draw, xl, yl, xr, yr, cl,
              )


	return nothing
end 

"""
	PetscDrawArrow(petsclib::PetscLibType,draw::PetscDraw, xl::PetscReal, yl::PetscReal, xr::PetscReal, yr::PetscReal, cl::Cint) 
draws a line with arrow head at end if the line is long enough

Not Collective

Input Parameters:
- `draw` - the drawing context
- `xl`   - horizontal coordinate of first end point
- `yl`   - vertical coordinate of first end point
- `xr`   - horizontal coordinate of second end point
- `yr`   - vertical coordinate of second end point
- `cl`   - the colors of the endpoints

Level: beginner

-seealso: `PetscDraw`, `PetscDrawLine()`, `PetscDrawLineSetWidth()`, `PetscDrawLineGetWidth()`, `PetscDrawRectangle()`, `PetscDrawTriangle()`, `PetscDrawEllipse()`,
`PetscDrawMarker()`, `PetscDrawPoint()`

# External Links
$(_doc_external("Sys/PetscDrawArrow"))
"""
function PetscDrawArrow(petsclib::PetscLibType, draw::PetscDraw, xl::PetscReal, yl::PetscReal, xr::PetscReal, yr::PetscReal, cl::Cint) end

@for_petsc function PetscDrawArrow(petsclib::$UnionPetscLib, draw::PetscDraw, xl::$PetscReal, yl::$PetscReal, xr::$PetscReal, yr::$PetscReal, cl::Cint )

    @chk ccall(
               (:PetscDrawArrow, $petsc_library),
               PetscErrorCode,
               (PetscDraw, $PetscReal, $PetscReal, $PetscReal, $PetscReal, Cint),
               draw, xl, yl, xr, yr, cl,
              )


	return nothing
end 

"""
	PetscDrawLineSetWidth(petsclib::PetscLibType,draw::PetscDraw, width::PetscReal) 
Sets the line width for future draws.  The width is
relative to the user coordinates of the window; 0.0 denotes the natural
width; 1.0 denotes the entire viewport.

Not Collective

Input Parameters:
- `draw`  - the drawing context
- `width` - the width in user coordinates

Level: advanced

-seealso: `PetscDraw`, `PetscDrawLineGetWidth()`, `PetscDrawLine()`, `PetscDrawArrow()`

# External Links
$(_doc_external("Sys/PetscDrawLineSetWidth"))
"""
function PetscDrawLineSetWidth(petsclib::PetscLibType, draw::PetscDraw, width::PetscReal) end

@for_petsc function PetscDrawLineSetWidth(petsclib::$UnionPetscLib, draw::PetscDraw, width::$PetscReal )

    @chk ccall(
               (:PetscDrawLineSetWidth, $petsc_library),
               PetscErrorCode,
               (PetscDraw, $PetscReal),
               draw, width,
              )


	return nothing
end 

"""
	width::PetscReal = PetscDrawLineGetWidth(petsclib::PetscLibType,draw::PetscDraw) 
Gets the line width for future draws.  The width is
relative to the user coordinates of the window; 0.0 denotes the natural
width; 1.0 denotes the interior viewport.

Not Collective

Input Parameter:
- `draw` - the drawing context

Output Parameter:
- `width` - the width in user coordinates

Level: advanced

-seealso: `PetscDraw`, `PetscDrawLineSetWidth()`, `PetscDrawLine()`, `PetscDrawArrow()`

# External Links
$(_doc_external("Sys/PetscDrawLineGetWidth"))
"""
function PetscDrawLineGetWidth(petsclib::PetscLibType, draw::PetscDraw) end

@for_petsc function PetscDrawLineGetWidth(petsclib::$UnionPetscLib, draw::PetscDraw )
	width_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscDrawLineGetWidth, $petsc_library),
               PetscErrorCode,
               (PetscDraw, Ptr{$PetscReal}),
               draw, width_,
              )

	width = width_[]

	return width
end 

"""
	PetscDrawTriangle(petsclib::PetscLibType,draw::PetscDraw, x1::PetscReal, y_1::PetscReal, x2::PetscReal, y2::PetscReal, x3::PetscReal, y3::PetscReal, c1::Cint, c2::Cint, c3::Cint) 
draws a triangle  onto a drawable.

Not Collective

Input Parameters:
- `draw` - the drawing context
- `x1`   - coordinate of the first vertex
- `y_1`  - coordinate of the first vertex
- `x2`   - coordinate of the second vertex
- `y2`   - coordinate of the second vertex
- `x3`   - coordinate of the third vertex
- `y3`   - coordinate of the third vertex
- `c1`   - color of the first vertex
- `c2`   - color of the second vertex
- `c3`   - color of the third vertext

Level: beginner

-seealso: `PetscDraw`, `PetscDrawLine()`, `PetscDrawRectangle()`, `PetscDrawEllipse()`, `PetscDrawMarker()`, `PetscDrawPoint()`, `PetscDrawArrow()`

# External Links
$(_doc_external("Sys/PetscDrawTriangle"))
"""
function PetscDrawTriangle(petsclib::PetscLibType, draw::PetscDraw, x1::PetscReal, y_1::PetscReal, x2::PetscReal, y2::PetscReal, x3::PetscReal, y3::PetscReal, c1::Cint, c2::Cint, c3::Cint) end

@for_petsc function PetscDrawTriangle(petsclib::$UnionPetscLib, draw::PetscDraw, x1::$PetscReal, y_1::$PetscReal, x2::$PetscReal, y2::$PetscReal, x3::$PetscReal, y3::$PetscReal, c1::Cint, c2::Cint, c3::Cint )

    @chk ccall(
               (:PetscDrawTriangle, $petsc_library),
               PetscErrorCode,
               (PetscDraw, $PetscReal, $PetscReal, $PetscReal, $PetscReal, $PetscReal, $PetscReal, Cint, Cint, Cint),
               draw, x1, y_1, x2, y2, x3, y3, c1, c2, c3,
              )


	return nothing
end 

"""
	PetscDrawScalePopup(petsclib::PetscLibType,popup::PetscDraw, min::PetscReal, max::PetscReal) 
draws a contour scale window.

Collective

Input Parameters:
- `popup` - the window (often a window obtained via `PetscDrawGetPopup()`
- `min`   - minimum value being plotted
- `max`   - maximum value being plotted

Level: intermediate

-seealso: `PetscDraw`, `PetscDrawGetPopup()`, `PetscDrawTensorContour()`

# External Links
$(_doc_external("Sys/PetscDrawScalePopup"))
"""
function PetscDrawScalePopup(petsclib::PetscLibType, popup::PetscDraw, min::PetscReal, max::PetscReal) end

@for_petsc function PetscDrawScalePopup(petsclib::$UnionPetscLib, popup::PetscDraw, min::$PetscReal, max::$PetscReal )

    @chk ccall(
               (:PetscDrawScalePopup, $petsc_library),
               PetscErrorCode,
               (PetscDraw, $PetscReal, $PetscReal),
               popup, min, max,
              )


	return nothing
end 

"""
	PetscDrawTensorContour(petsclib::PetscLibType,draw::PetscDraw, m::Cint, n::Cint, xi::Vector{PetscReal}, yi::Vector{PetscReal}, v::Vector{PetscReal}) 
draws a contour plot for a two

Collective, but `draw` must be sequential

Input Parameters:
- `draw` - the draw context
- `m`    - the number of local mesh points in the x direction
- `n`    - the number of local mesh points in the y direction
- `xi`   - the locations of the global mesh points in the horizontal direction (optional, use `NULL` to indicate uniform spacing on [0,1])
- `yi`   - the locations of the global mesh points in the vertical direction (optional, use `NULL` to indicate uniform spacing on [0,1])
- `v`    - the values

Options Database Keys:
- `-draw_x_shared_colormap` - Indicates use of private colormap
- `-draw_contour_grid`      - draws grid contour

Level: intermediate

-seealso: `PetscDraw`, `PetscDrawTensorContourPatch()`, `PetscDrawScalePopup()`

# External Links
$(_doc_external("Sys/PetscDrawTensorContour"))
"""
function PetscDrawTensorContour(petsclib::PetscLibType, draw::PetscDraw, m::Cint, n::Cint, xi::Vector{PetscReal}, yi::Vector{PetscReal}, v::Vector{PetscReal}) end

@for_petsc function PetscDrawTensorContour(petsclib::$UnionPetscLib, draw::PetscDraw, m::Cint, n::Cint, xi::Vector{$PetscReal}, yi::Vector{$PetscReal}, v::Vector{$PetscReal} )

    @chk ccall(
               (:PetscDrawTensorContour, $petsc_library),
               PetscErrorCode,
               (PetscDraw, Cint, Cint, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
               draw, m, n, xi, yi, v,
              )


	return nothing
end 

"""
	PetscDrawTensorContourPatch(petsclib::PetscLibType,draw::PetscDraw, m::Cint, n::Cint, x::PetscReal, y::PetscReal, min::PetscReal, max::PetscReal, v::PetscReal) 
draws a rectangular patch of a contour plot
for a two-dimensional array.

Not Collective

Input Parameters:
- `draw` - the draw context
- `m`    - the number of local mesh points in the x direction
- `n`    - the number of local mesh points in the y direction
- `x`    - the horizontal locations of the local mesh points
- `y`    - the vertical locations of the local mesh points
- `min`  - the minimum value in the entire contour
- `max`  - the maximum value in the entire contour
- `v`    - the data

Options Database Key:
- `-draw_x_shared_colormap` - Activates private colormap

Level: advanced

-seealso: `PetscDraw`, `PetscDrawTensorContour()`

# External Links
$(_doc_external("Sys/PetscDrawTensorContourPatch"))
"""
function PetscDrawTensorContourPatch(petsclib::PetscLibType, draw::PetscDraw, m::Cint, n::Cint, x::PetscReal, y::PetscReal, min::PetscReal, max::PetscReal, v::PetscReal) end

@for_petsc function PetscDrawTensorContourPatch(petsclib::$UnionPetscLib, draw::PetscDraw, m::Cint, n::Cint, x::$PetscReal, y::$PetscReal, min::$PetscReal, max::$PetscReal, v::$PetscReal )

    @chk ccall(
               (:PetscDrawTensorContourPatch, $petsc_library),
               PetscErrorCode,
               (PetscDraw, Cint, Cint, Ptr{$PetscReal}, Ptr{$PetscReal}, $PetscReal, $PetscReal, Ptr{$PetscReal}),
               draw, m, n, x, y, min, max, v,
              )


	return nothing
end 

"""
	PetscDrawPoint(petsclib::PetscLibType,draw::PetscDraw, xl::PetscReal, yl::PetscReal, cl::Cint) 
draws a point onto a drawable.

Not Collective

Input Parameters:
- `draw` - the drawing context
- `xl`   - horizatonal coordinate of the point
- `yl`   - vertical coordinate of the point
- `cl`   - the color of the point

Level: beginner

-seealso: `PetscDraw`, `PetscDrawPointPixel()`, `PetscDrawPointSetSize()`, `PetscDrawLine()`, `PetscDrawRectangle()`, `PetscDrawTriangle()`, `PetscDrawEllipse()`,
`PetscDrawMarker()`, `PetscDrawString()`, `PetscDrawArrow()`

# External Links
$(_doc_external("Sys/PetscDrawPoint"))
"""
function PetscDrawPoint(petsclib::PetscLibType, draw::PetscDraw, xl::PetscReal, yl::PetscReal, cl::Cint) end

@for_petsc function PetscDrawPoint(petsclib::$UnionPetscLib, draw::PetscDraw, xl::$PetscReal, yl::$PetscReal, cl::Cint )

    @chk ccall(
               (:PetscDrawPoint, $petsc_library),
               PetscErrorCode,
               (PetscDraw, $PetscReal, $PetscReal, Cint),
               draw, xl, yl, cl,
              )


	return nothing
end 

"""
	PetscDrawPointPixel(petsclib::PetscLibType,draw::PetscDraw, x::Cint, y::Cint, c::Cint) 
draws a point onto a drawable, in pixel coordinates

Not Collective

Input Parameters:
- `draw` - the drawing context
- `x`    - horizontal pixel coordinates of the point
- `y`    - vertical pixel coordinates of the point
- `c`    - the color of the point

Level: beginner

-seealso: `PetscDraw`, `PetscDrawPoint()`, `PetscDrawPointSetSize()`

# External Links
$(_doc_external("Sys/PetscDrawPointPixel"))
"""
function PetscDrawPointPixel(petsclib::PetscLibType, draw::PetscDraw, x::Cint, y::Cint, c::Cint) end

@for_petsc function PetscDrawPointPixel(petsclib::$UnionPetscLib, draw::PetscDraw, x::Cint, y::Cint, c::Cint )

    @chk ccall(
               (:PetscDrawPointPixel, $petsc_library),
               PetscErrorCode,
               (PetscDraw, Cint, Cint, Cint),
               draw, x, y, c,
              )


	return nothing
end 

"""
	PetscDrawPointSetSize(petsclib::PetscLibType,draw::PetscDraw, width::PetscReal) 
Sets the point size for future draws.  The size is
relative to the user coordinates of the window; 0.0 denotes the natural
width, 1.0 denotes the entire viewport.

Not Collective

Input Parameters:
- `draw`  - the drawing context
- `width` - the width in user coordinates

Level: advanced

-seealso: `PetscDraw`, `PetscDrawPoint()`, `PetscDrawMarker()`

# External Links
$(_doc_external("Sys/PetscDrawPointSetSize"))
"""
function PetscDrawPointSetSize(petsclib::PetscLibType, draw::PetscDraw, width::PetscReal) end

@for_petsc function PetscDrawPointSetSize(petsclib::$UnionPetscLib, draw::PetscDraw, width::$PetscReal )

    @chk ccall(
               (:PetscDrawPointSetSize, $petsc_library),
               PetscErrorCode,
               (PetscDraw, $PetscReal),
               draw, width,
              )


	return nothing
end 

"""
	PetscDrawClear(petsclib::PetscLibType,draw::PetscDraw) 
Clears graphical output. All processors must call this routine.
Does not return until the draw in context is clear.

Collective

Input Parameter:
- `draw` - the drawing context

Level: intermediate

-seealso: `PetscDrawBOP()`, `PetscDrawEOP()`

# External Links
$(_doc_external("Sys/PetscDrawClear"))
"""
function PetscDrawClear(petsclib::PetscLibType, draw::PetscDraw) end

@for_petsc function PetscDrawClear(petsclib::$UnionPetscLib, draw::PetscDraw )

    @chk ccall(
               (:PetscDrawClear, $petsc_library),
               PetscErrorCode,
               (PetscDraw,),
               draw,
              )


	return nothing
end 

"""
	PetscDrawBOP(petsclib::PetscLibType,draw::PetscDraw) 
Begins a new page or frame on the selected graphical device.

Logically Collective

Input Parameter:
- `draw` - the drawing context

Level: advanced

-seealso: `PetscDrawEOP()`, `PetscDrawClear()`

# External Links
$(_doc_external("Sys/PetscDrawBOP"))
"""
function PetscDrawBOP(petsclib::PetscLibType, draw::PetscDraw) end

@for_petsc function PetscDrawBOP(petsclib::$UnionPetscLib, draw::PetscDraw )

    @chk ccall(
               (:PetscDrawBOP, $petsc_library),
               PetscErrorCode,
               (PetscDraw,),
               draw,
              )


	return nothing
end 

"""
	PetscDrawEOP(petsclib::PetscLibType,draw::PetscDraw) 
Ends a page or frame on the selected graphical device.

Logically Collective

Input Parameter:
- `draw` - the drawing context

Level: advanced

-seealso: `PetscDrawBOP()`, `PetscDrawClear()`

# External Links
$(_doc_external("Sys/PetscDrawEOP"))
"""
function PetscDrawEOP(petsclib::PetscLibType, draw::PetscDraw) end

@for_petsc function PetscDrawEOP(petsclib::$UnionPetscLib, draw::PetscDraw )

    @chk ccall(
               (:PetscDrawEOP, $petsc_library),
               PetscErrorCode,
               (PetscDraw,),
               draw,
              )


	return nothing
end 

"""
	PetscDrawPause(petsclib::PetscLibType,draw::PetscDraw) 
Waits n seconds or until user input, depending on input
to `PetscDrawSetPause()`.

Collective

Input Parameter:
- `draw` - the drawing context

Level: beginner

-seealso: `PetscDraw`, `PetscDrawSetPause()`, `PetscDrawGetPause()`

# External Links
$(_doc_external("Sys/PetscDrawPause"))
"""
function PetscDrawPause(petsclib::PetscLibType, draw::PetscDraw) end

@for_petsc function PetscDrawPause(petsclib::$UnionPetscLib, draw::PetscDraw )

    @chk ccall(
               (:PetscDrawPause, $petsc_library),
               PetscErrorCode,
               (PetscDraw,),
               draw,
              )


	return nothing
end 

"""
	PetscDrawSetPause(petsclib::PetscLibType,draw::PetscDraw, lpause::PetscReal) 
Sets the amount of time that program pauses after
a `PetscDrawPause()` is called.

Logically Collective

Input Parameters:
- `draw`   - the drawing object
- `lpause` - number of seconds to pause, -1 implies until user input, -2 pauses only on the `PetscDrawDestroy()`

Options Database Key:
- `-draw_pause value` - set the time to pause

Level: intermediate

-seealso: `PetscDraw`, `PetscDrawGetPause()`, `PetscDrawPause()`

# External Links
$(_doc_external("Sys/PetscDrawSetPause"))
"""
function PetscDrawSetPause(petsclib::PetscLibType, draw::PetscDraw, lpause::PetscReal) end

@for_petsc function PetscDrawSetPause(petsclib::$UnionPetscLib, draw::PetscDraw, lpause::$PetscReal )

    @chk ccall(
               (:PetscDrawSetPause, $petsc_library),
               PetscErrorCode,
               (PetscDraw, $PetscReal),
               draw, lpause,
              )


	return nothing
end 

"""
	PetscDrawGetPause(petsclib::PetscLibType,draw::PetscDraw, lpause::PetscReal) 
Gets the amount of time that program pauses after
a `PetscDrawPause()` is called.

Not Collective

Input Parameters:
- `draw`   - the drawing object
- `lpause` - number of seconds to pause, -1 implies until user input

Level: intermediate

-seealso: `PetscDraw`, `PetscDrawSetPause()`, `PetscDrawPause()`

# External Links
$(_doc_external("Sys/PetscDrawGetPause"))
"""
function PetscDrawGetPause(petsclib::PetscLibType, draw::PetscDraw, lpause::PetscReal) end

@for_petsc function PetscDrawGetPause(petsclib::$UnionPetscLib, draw::PetscDraw, lpause::$PetscReal )

    @chk ccall(
               (:PetscDrawGetPause, $petsc_library),
               PetscErrorCode,
               (PetscDraw, Ptr{$PetscReal}),
               draw, lpause,
              )


	return nothing
end 

"""
	PetscDrawMarker(petsclib::PetscLibType,draw::PetscDraw, xl::PetscReal, yl::PetscReal, cl::Cint) 
draws a marker onto a drawable.

Not Collective

Input Parameters:
- `draw` - the drawing context
- `xl`   - horizontal coordinate of the marker
- `yl`   - vertical coordinate of the marker
- `cl`   - the color of the marker

Level: beginner

-seealso: `PetscDraw`, `PetscDrawPoint()`, `PetscDrawString()`, `PetscDrawSetMarkerType()`, `PetscDrawGetMarkerType()`

# External Links
$(_doc_external("Sys/PetscDrawMarker"))
"""
function PetscDrawMarker(petsclib::PetscLibType, draw::PetscDraw, xl::PetscReal, yl::PetscReal, cl::Cint) end

@for_petsc function PetscDrawMarker(petsclib::$UnionPetscLib, draw::PetscDraw, xl::$PetscReal, yl::$PetscReal, cl::Cint )

    @chk ccall(
               (:PetscDrawMarker, $petsc_library),
               PetscErrorCode,
               (PetscDraw, $PetscReal, $PetscReal, Cint),
               draw, xl, yl, cl,
              )


	return nothing
end 

"""
	PetscDrawSetMarkerType(petsclib::PetscLibType,draw::PetscDraw, mtype::PetscDrawMarkerType) 
sets the type of marker to display with `PetscDrawMarker()`

Not Collective

Input Parameters:
- `draw`  - the drawing context
- `mtype` - either `PETSC_DRAW_MARKER_CROSS` (default) or `PETSC_DRAW_MARKER_POINT`

Options Database Key:
- `-draw_marker_type` - x or point

Level: beginner

-seealso: `PetscDraw`, `PetscDrawPoint()`, `PetscDrawMarker()`, `PetscDrawGetMarkerType()`, `PetscDrawMarkerType`

# External Links
$(_doc_external("Sys/PetscDrawSetMarkerType"))
"""
function PetscDrawSetMarkerType(petsclib::PetscLibType, draw::PetscDraw, mtype::PetscDrawMarkerType) end

@for_petsc function PetscDrawSetMarkerType(petsclib::$UnionPetscLib, draw::PetscDraw, mtype::PetscDrawMarkerType )

    @chk ccall(
               (:PetscDrawSetMarkerType, $petsc_library),
               PetscErrorCode,
               (PetscDraw, PetscDrawMarkerType),
               draw, mtype,
              )


	return nothing
end 

"""
	PetscDrawGetMarkerType(petsclib::PetscLibType,draw::PetscDraw, mtype::PetscDrawMarkerType) 
gets the type of marker to display with `PetscDrawMarker()`

Not Collective

Input Parameters:
- `draw`  - the drawing context
- `mtype` - either `PETSC_DRAW_MARKER_CROSS` (default) or `PETSC_DRAW_MARKER_POINT`

Level: beginner

-seealso: `PetscDraw`, `PetscDrawPoint()`, `PetscDrawMarker()`, `PetscDrawSetMarkerType()`, `PetscDrawMarkerType`

# External Links
$(_doc_external("Sys/PetscDrawGetMarkerType"))
"""
function PetscDrawGetMarkerType(petsclib::PetscLibType, draw::PetscDraw, mtype::PetscDrawMarkerType) end

@for_petsc function PetscDrawGetMarkerType(petsclib::$UnionPetscLib, draw::PetscDraw, mtype::PetscDrawMarkerType )

    @chk ccall(
               (:PetscDrawGetMarkerType, $petsc_library),
               PetscErrorCode,
               (PetscDraw, Ptr{PetscDrawMarkerType}),
               draw, mtype,
              )


	return nothing
end 

"""
	PetscDrawIndicatorFunction(petsclib::PetscLibType,draw::PetscDraw, xmin::PetscReal, xmax::PetscReal, ymin::PetscReal, ymax::PetscReal, c::Cint, indicator::external, ctx::Cvoid) 
Draws an indicator function (where a relationship is true) on a `PetscDraw`

Not Collective

Input Parameters:
- `draw`      - a `PetscDraw`
- `xmin`      - region to draw indicator function
- `xmax`      - region to draw indicator function
- `ymin`      - region to draw indicator function
- `ymax`      - region to draw indicator function
- `c`         - the color of the region
- `indicator` - the indicator function
- `ctx`       - the context to pass to the indicator function

Level: developer

-seealso: `PetscDraw`

# External Links
$(_doc_external("Sys/PetscDrawIndicatorFunction"))
"""
function PetscDrawIndicatorFunction(petsclib::PetscLibType, draw::PetscDraw, xmin::PetscReal, xmax::PetscReal, ymin::PetscReal, ymax::PetscReal, c::Cint, indicator::external, ctx::Cvoid) end

@for_petsc function PetscDrawIndicatorFunction(petsclib::$UnionPetscLib, draw::PetscDraw, xmin::$PetscReal, xmax::$PetscReal, ymin::$PetscReal, ymax::$PetscReal, c::Cint, indicator::external, ctx::Cvoid )

    @chk ccall(
               (:PetscDrawIndicatorFunction, $petsc_library),
               PetscErrorCode,
               (PetscDraw, $PetscReal, $PetscReal, $PetscReal, $PetscReal, Cint, external, Ptr{Cvoid}),
               draw, xmin, xmax, ymin, ymax, c, indicator, ctx,
              )


	return nothing
end 

"""
	PetscDrawCoordinateToPixel(petsclib::PetscLibType,draw::PetscDraw, x::PetscReal, y::PetscReal, i::Cint, j::Cint) 
given a coordinate in a `PetscDraw` returns the pixel location

Not Collective

Input Parameters:
- `draw` - the draw where the coordinates are defined
- `x`    - the horizontal coordinate
- `y`    - the vertical coordinate

Output Parameters:
- `i` - the horizontal pixel location
- `j` - the vertical pixel location

Level: developer

-seealso: `PetscDraw`

# External Links
$(_doc_external("Sys/PetscDrawCoordinateToPixel"))
"""
function PetscDrawCoordinateToPixel(petsclib::PetscLibType, draw::PetscDraw, x::PetscReal, y::PetscReal, i::Cint, j::Cint) end

@for_petsc function PetscDrawCoordinateToPixel(petsclib::$UnionPetscLib, draw::PetscDraw, x::$PetscReal, y::$PetscReal, i::Cint, j::Cint )

    @chk ccall(
               (:PetscDrawCoordinateToPixel, $petsc_library),
               PetscErrorCode,
               (PetscDraw, $PetscReal, $PetscReal, Ptr{Cint}, Ptr{Cint}),
               draw, x, y, i, j,
              )


	return nothing
end 

"""
	x::PetscReal,y::PetscReal = PetscDrawPixelToCoordinate(petsclib::PetscLibType,draw::PetscDraw, i::Cint, j::Cint) 
given a pixel in a `PetscDraw` returns the coordinate

Not Collective

Input Parameters:
- `draw` - the draw where the coordinates are defined
- `i`    - the horizontal pixel location
- `j`    - the vertical pixel location

Output Parameters:
- `x` - the horizontal coordinate
- `y` - the vertical coordinate

Level: developer

-seealso: `PetscDraw`

# External Links
$(_doc_external("Sys/PetscDrawPixelToCoordinate"))
"""
function PetscDrawPixelToCoordinate(petsclib::PetscLibType, draw::PetscDraw, i::Cint, j::Cint) end

@for_petsc function PetscDrawPixelToCoordinate(petsclib::$UnionPetscLib, draw::PetscDraw, i::Cint, j::Cint )
	x_ = Ref{$PetscReal}()
	y_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscDrawPixelToCoordinate, $petsc_library),
               PetscErrorCode,
               (PetscDraw, Cint, Cint, Ptr{$PetscReal}, Ptr{$PetscReal}),
               draw, i, j, x_, y_,
              )

	x = x_[]
	y = y_[]

	return x,y
end 

"""
	PetscDrawRectangle(petsclib::PetscLibType,draw::PetscDraw, xl::PetscReal, yl::PetscReal, xr::PetscReal, yr::PetscReal, c1::Cint, c2::Cint, c3::Cint, c4::Cint) 
draws a rectangle onto a `PetscDraw` object

Not Collective

Input Parameters:
- `draw` - the drawing context
- `xl`   - coordinates of the lower left corner
- `yl`   - coordinates of the lower left corner
- `xr`   - coordinate of the upper right corner
- `yr`   - coordinate of the upper right corner
- `c1`   - the color of the first corner
- `c2`   - the color of the second corner
- `c3`   - the color of the third corner
- `c4`   - the color of the fourth corner

Level: beginner

-seealso: `PetscDraw`, `PetscDrawLine()`, `PetscDrawTriangle()`, `PetscDrawEllipse()`,
`PetscDrawMarker()`, `PetscDrawPoint()`, `PetscDrawString()`, `PetscDrawArrow()`

# External Links
$(_doc_external("Sys/PetscDrawRectangle"))
"""
function PetscDrawRectangle(petsclib::PetscLibType, draw::PetscDraw, xl::PetscReal, yl::PetscReal, xr::PetscReal, yr::PetscReal, c1::Cint, c2::Cint, c3::Cint, c4::Cint) end

@for_petsc function PetscDrawRectangle(petsclib::$UnionPetscLib, draw::PetscDraw, xl::$PetscReal, yl::$PetscReal, xr::$PetscReal, yr::$PetscReal, c1::Cint, c2::Cint, c3::Cint, c4::Cint )

    @chk ccall(
               (:PetscDrawRectangle, $petsc_library),
               PetscErrorCode,
               (PetscDraw, $PetscReal, $PetscReal, $PetscReal, $PetscReal, Cint, Cint, Cint, Cint),
               draw, xl, yl, xr, yr, c1, c2, c3, c4,
              )


	return nothing
end 

"""
	PetscDrawOpenImage(petsclib::PetscLibType,comm::MPI_Comm, filename::String, w::Cint, h::Cint, draw::PetscDraw) 
Opens an image for use with the `PetscDraw` routines.

Collective

Input Parameters:
- `comm`     - the communicator that will share image
- `filename` - optional name of the file where the image will be stored
- `w`        - the image width in pixels
- `h`        - the image height in pixels

Output Parameter:
- `draw` - the drawing context.

Level: beginner

-seealso: `PetscDraw`, `PETSC_DRAW_IMAGE`, `PETSC_DRAW_X`, `PetscDrawSetSave()`, `PetscDrawSetFromOptions()`, `PetscDrawCreate()`, `PetscDrawDestroy()`

# External Links
$(_doc_external("Sys/PetscDrawOpenImage"))
"""
function PetscDrawOpenImage(petsclib::PetscLibType, comm::MPI_Comm, filename::String, w::Cint, h::Cint, draw::PetscDraw) end

@for_petsc function PetscDrawOpenImage(petsclib::$UnionPetscLib, comm::MPI_Comm, filename::String, w::Cint, h::Cint, draw::PetscDraw )

    @chk ccall(
               (:PetscDrawOpenImage, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, Cint, Cint, Ptr{PetscDraw}),
               comm, filename, w, h, draw,
              )


	return nothing
end 

"""
	PetscDrawOpenNull(petsclib::PetscLibType,comm::MPI_Comm, win::PetscDraw) 
Opens a null drawing context. All draw commands to
it are ignored.

Input Parameter:
- `comm` - MPI communicator

Output Parameter:
- `win` - the drawing context

Level: advanced

-seealso: `PetscDraw`, `PetscDrawIsNull()`, `PETSC_DRAW_NULL`, `PetscDrawOpenX()`

# External Links
$(_doc_external("Sys/PetscDrawOpenNull"))
"""
function PetscDrawOpenNull(petsclib::PetscLibType, comm::MPI_Comm, win::PetscDraw) end

@for_petsc function PetscDrawOpenNull(petsclib::$UnionPetscLib, comm::MPI_Comm, win::PetscDraw )

    @chk ccall(
               (:PetscDrawOpenNull, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PetscDraw}),
               comm, win,
              )


	return nothing
end 

"""
	yes::PetscBool = PetscDrawIsNull(petsclib::PetscLibType,draw::PetscDraw) 
Returns `PETSC_TRUE` if draw is a null draw object.

Not Collective

Input Parameter:
- `draw` - the draw context

Output Parameter:
- `yes` - `PETSC_TRUE` if it is a null draw object; otherwise `PETSC_FALSE`

Level: advanced

-seealso: `PetscDraw`, `PETSC_DRAW_NULL`, `PetscDrawOpenX()`

# External Links
$(_doc_external("Sys/PetscDrawIsNull"))
"""
function PetscDrawIsNull(petsclib::PetscLibType, draw::PetscDraw) end

@for_petsc function PetscDrawIsNull(petsclib::$UnionPetscLib, draw::PetscDraw )
	yes_ = Ref{PetscBool}()

    @chk ccall(
               (:PetscDrawIsNull, $petsc_library),
               PetscErrorCode,
               (PetscDraw, Ptr{PetscBool}),
               draw, yes_,
              )

	yes = yes_[]

	return yes
end 

"""
	PetscDrawOpenX(petsclib::PetscLibType,comm::MPI_Comm, display::String, title::String, x::Cint, y::Cint, w::Cint, h::Cint, draw::PetscDraw) 
Opens an X

Collective

Input Parameters:
- `comm`    - the communicator that will share X-window
- `display` - the X display on which to open, or `NULL` for the local machine
- `title`   - the title to put in the title bar, or `NULL` for no title
- `x`       - the x screen coordinates of the upper left corner of window (or `PETSC_DECIDE`)
- `y`       - the y screen coordinates of the upper left corner of window (or `PETSC_DECIDE`)
- `w`       - the screen width in pixels of (or `PETSC_DRAW_HALF_SIZE`, `PETSC_DRAW_FULL_SIZE`, or `PETSC_DRAW_THIRD_SIZE` or `PETSC_DRAW_QUARTER_SIZE`)
- `h`       - the screen height in pixels of (or `PETSC_DRAW_HALF_SIZE`, `PETSC_DRAW_FULL_SIZE`, or `PETSC_DRAW_THIRD_SIZE` or `PETSC_DRAW_QUARTER_SIZE`)

Output Parameter:
- `draw` - the drawing context.

Options Database Keys:
- `-nox`                    - Disables all x-windows output
- `-display <name>`         - Sets name of machine for the X display
- `-draw_pause <pause>`     - Sets time (in seconds) that the program pauses after `PetscDrawPause()` has been called
(0 is default, -1 implies until user input).
- `-draw_cmap <name>`       - Sets the colormap to use.
- `-draw_cmap_reverse`      - Reverses the colormap.
- `-draw_cmap_brighten`     - Brighten (0 < beta < 1) or darken (-1 < beta < 0) the colormap.
- `-draw_x_shared_colormap` - Causes PETSc to use a shared colormap. By default PETSc creates a separate color
for its windows, you must put the mouse into the graphics
window to see  the correct colors. This options forces
PETSc to use the default colormap which will usually result
in bad contour plots.
- `-draw_fast`              - Does not create colormap for contour plots.
- `-draw_double_buffer`     - Uses double buffering for smooth animation.
- `-geometry`               - Indicates location and size of window.

Level: beginner

-seealso: `PetscDrawFlush()`, `PetscDrawDestroy()`, `PetscDrawCreate()`, `PetscDrawOpnOpenGL()`

# External Links
$(_doc_external("Sys/PetscDrawOpenX"))
"""
function PetscDrawOpenX(petsclib::PetscLibType, comm::MPI_Comm, display::String, title::String, x::Cint, y::Cint, w::Cint, h::Cint, draw::PetscDraw) end

@for_petsc function PetscDrawOpenX(petsclib::$UnionPetscLib, comm::MPI_Comm, display::String, title::String, x::Cint, y::Cint, w::Cint, h::Cint, draw::PetscDraw )

    @chk ccall(
               (:PetscDrawOpenX, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}, Cint, Cint, Cint, Cint, Ptr{PetscDraw}),
               comm, display, title, x, y, w, h, draw,
              )


	return nothing
end 

"""
	PetscDrawZoom(petsclib::PetscLibType,draw::PetscDraw, func::external, ctx::Cvoid) 
Allows one to provide a function that gets called for zooming in on a drawing using the mouse buttons

Collective draw

Input Parameters:
- `draw` - the window where the graph will be made.
- `func` - users function that draws the graphic
- `ctx`  - pointer to any user required data

Level: advanced

-seealso: `PetscDraw`, `PetscDrawCreate()`

# External Links
$(_doc_external("Sys/PetscDrawZoom"))
"""
function PetscDrawZoom(petsclib::PetscLibType, draw::PetscDraw, func::external, ctx::Cvoid) end

@for_petsc function PetscDrawZoom(petsclib::$UnionPetscLib, draw::PetscDraw, func::external, ctx::Cvoid )

    @chk ccall(
               (:PetscDrawZoom, $petsc_library),
               PetscErrorCode,
               (PetscDraw, external, Ptr{Cvoid}),
               draw, func, ctx,
              )


	return nothing
end 

"""
	PetscDrawUtilitySetGamma(petsclib::PetscLibType,g::PetscReal) 

# External Links
$(_doc_external("Sys/PetscDrawUtilitySetGamma"))
"""
function PetscDrawUtilitySetGamma(petsclib::PetscLibType, g::PetscReal) end

@for_petsc function PetscDrawUtilitySetGamma(petsclib::$UnionPetscLib, g::$PetscReal )

    @chk ccall(
               (:PetscDrawUtilitySetGamma, $petsc_library),
               PetscErrorCode,
               ($PetscReal,),
               g,
              )


	return nothing
end 

#=
"""
	PetscDrawUtilitySetCmap(petsclib::PetscLibType,colormap::String, mapsize::Cint, char::Vector{unsigned}, char::Vector{unsigned}, char::Vector{unsigned}) 

# External Links
$(_doc_external("Sys/PetscDrawUtilitySetCmap"))
"""
function PetscDrawUtilitySetCmap(petsclib::PetscLibType, colormap::String, mapsize::Cint, char::Vector{unsigned}, char::Vector{unsigned}, char::Vector{unsigned}) end

@for_petsc function PetscDrawUtilitySetCmap(petsclib::$UnionPetscLib, colormap::String, mapsize::Cint, char::Vector{unsigned}, char::Vector{unsigned}, char::Vector{unsigned} )

    @chk ccall(
               (:PetscDrawUtilitySetCmap, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Cint, Ptr{unsigned}, Ptr{unsigned}, Ptr{unsigned}),
               colormap, mapsize, char, char, char,
              )


	return nothing
end 
=#

"""
	axis::PetscDrawAxis = PetscDrawAxisCreate(petsclib::PetscLibType,draw::PetscDraw) 
Generate the axis data structure.

Collective

Input Parameter:
- `draw` - `PetscDraw` object where axis to be made

Output Parameter:
- `axis` - the axis datastructure

-seealso: `PetscDrawLGCreate()`, `PetscDrawLG`, `PetscDrawSPCreate()`, `PetscDrawSP`, `PetscDrawHGCreate()`, `PetscDrawHG`, `PetscDrawBarCreate()`, `PetscDrawBar`, `PetscDrawLGGetAxis()`, `PetscDrawSPGetAxis()`,
`PetscDrawHGGetAxis()`, `PetscDrawBarGetAxis()`, `PetscDrawAxis`, `PetscDrawAxisDestroy()`, `PetscDrawAxisSetColors()`, `PetscDrawAxisSetLabels()`, `PetscDrawAxisSetLimits()`, `PetscDrawAxisGetLimits()`, `PetscDrawAxisSetHoldLimits()`,
`PetscDrawAxisDraw()`

# External Links
$(_doc_external("Sys/PetscDrawAxisCreate"))
"""
function PetscDrawAxisCreate(petsclib::PetscLibType, draw::PetscDraw) end

@for_petsc function PetscDrawAxisCreate(petsclib::$UnionPetscLib, draw::PetscDraw )
	axis_ = Ref{PetscDrawAxis}()

    @chk ccall(
               (:PetscDrawAxisCreate, $petsc_library),
               PetscErrorCode,
               (PetscDraw, Ptr{PetscDrawAxis}),
               draw, axis_,
              )

	axis = axis_[]

	return axis
end 

"""
	PetscDrawAxisDestroy(petsclib::PetscLibType,axis::PetscDrawAxis) 
Frees the space used by an axis structure.

Collective

Input Parameter:
- `axis` - the axis context

Level: advanced

-seealso: `PetscDraw`, `PetscDrawAxisCreate()`, `PetscDrawAxis`

# External Links
$(_doc_external("Sys/PetscDrawAxisDestroy"))
"""
function PetscDrawAxisDestroy(petsclib::PetscLibType, axis::PetscDrawAxis) end

@for_petsc function PetscDrawAxisDestroy(petsclib::$UnionPetscLib, axis::PetscDrawAxis )

    @chk ccall(
               (:PetscDrawAxisDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscDrawAxis},),
               axis,
              )


	return nothing
end 

"""
	PetscDrawAxisSetColors(petsclib::PetscLibType,axis::PetscDrawAxis, ac::Cint, tc::Cint, cc::Cint) 
Sets the colors to be used for the axis,
tickmarks, and text.

Logically Collective

Input Parameters:
- `axis` - the axis
- `ac`   - the color of the axis lines
- `tc`   - the color of the tick marks
- `cc`   - the color of the text strings

Level: advanced

-seealso: `PetscDraw`, `PetscDrawAxisCreate()`, `PetscDrawAxis`, `PetscDrawAxisSetLabels()`, `PetscDrawAxisDraw()`, `PetscDrawAxisSetLimits()`

# External Links
$(_doc_external("Sys/PetscDrawAxisSetColors"))
"""
function PetscDrawAxisSetColors(petsclib::PetscLibType, axis::PetscDrawAxis, ac::Cint, tc::Cint, cc::Cint) end

@for_petsc function PetscDrawAxisSetColors(petsclib::$UnionPetscLib, axis::PetscDrawAxis, ac::Cint, tc::Cint, cc::Cint )

    @chk ccall(
               (:PetscDrawAxisSetColors, $petsc_library),
               PetscErrorCode,
               (PetscDrawAxis, Cint, Cint, Cint),
               axis, ac, tc, cc,
              )


	return nothing
end 

"""
	PetscDrawAxisSetLabels(petsclib::PetscLibType,axis::PetscDrawAxis, top::String, xlabel::String, ylabel::String) 
Sets the x and y axis labels.

Logically Collective

Input Parameters:
- `axis`   - the axis
- `top`    - the label at the top of the image
- `xlabel` - the x axis label
- `ylabel` - the y axis label

Level: advanced

-seealso: `PetscDraw`, `PetscDrawAxisCreate()`, `PetscDrawAxis`, `PetscDrawAxisSetColors()`, `PetscDrawAxisDraw()`, `PetscDrawAxisSetLimits()`

# External Links
$(_doc_external("Sys/PetscDrawAxisSetLabels"))
"""
function PetscDrawAxisSetLabels(petsclib::PetscLibType, axis::PetscDrawAxis, top::String, xlabel::String, ylabel::String) end

@for_petsc function PetscDrawAxisSetLabels(petsclib::$UnionPetscLib, axis::PetscDrawAxis, top::String, xlabel::String, ylabel::String )

    @chk ccall(
               (:PetscDrawAxisSetLabels, $petsc_library),
               PetscErrorCode,
               (PetscDrawAxis, Ptr{Cchar}, Ptr{Cchar}, Ptr{Cchar}),
               axis, top, xlabel, ylabel,
              )


	return nothing
end 

"""
	PetscDrawAxisSetLimits(petsclib::PetscLibType,axis::PetscDrawAxis, xmin::PetscReal, xmax::PetscReal, ymin::PetscReal, ymax::PetscReal) 
Sets the limits (in user coords) of the axis

Logically Collective

Input Parameters:
- `axis` - the axis
- `xmin` - the lower x limit
- `xmax` - the upper x limit
- `ymin` - the lower y limit
- `ymax` - the upper y limit

Options Database Key:
- `-drawaxis_hold` - hold the initial set of axis limits for future plotting

Level: advanced

-seealso: `PetscDrawAxisSetHoldLimits()`, `PetscDrawAxisGetLimits()`, `PetscDrawAxisSetLabels()`, `PetscDrawAxisSetColors()`

# External Links
$(_doc_external("Sys/PetscDrawAxisSetLimits"))
"""
function PetscDrawAxisSetLimits(petsclib::PetscLibType, axis::PetscDrawAxis, xmin::PetscReal, xmax::PetscReal, ymin::PetscReal, ymax::PetscReal) end

@for_petsc function PetscDrawAxisSetLimits(petsclib::$UnionPetscLib, axis::PetscDrawAxis, xmin::$PetscReal, xmax::$PetscReal, ymin::$PetscReal, ymax::$PetscReal )

    @chk ccall(
               (:PetscDrawAxisSetLimits, $petsc_library),
               PetscErrorCode,
               (PetscDrawAxis, $PetscReal, $PetscReal, $PetscReal, $PetscReal),
               axis, xmin, xmax, ymin, ymax,
              )


	return nothing
end 

"""
	PetscDrawAxisGetLimits(petsclib::PetscLibType,axis::PetscDrawAxis, xmin::PetscReal, xmax::PetscReal, ymin::PetscReal, ymax::PetscReal) 
Gets the limits (in user coords) of the axis

Not Collective

Input Parameters:
- `axis` - the axis
- `xmin` - the lower x limit
- `xmax` - the upper x limit
- `ymin` - the lower y limit
- `ymax` - the upper y limit

Level: advanced

-seealso: `PetscDrawAxisCreate()`, `PetscDrawAxis`, `PetscDrawAxisSetHoldLimits()`, `PetscDrawAxisSetLimits()`, `PetscDrawAxisSetLabels()`, `PetscDrawAxisSetColors()`

# External Links
$(_doc_external("Sys/PetscDrawAxisGetLimits"))
"""
function PetscDrawAxisGetLimits(petsclib::PetscLibType, axis::PetscDrawAxis, xmin::PetscReal, xmax::PetscReal, ymin::PetscReal, ymax::PetscReal) end

@for_petsc function PetscDrawAxisGetLimits(petsclib::$UnionPetscLib, axis::PetscDrawAxis, xmin::$PetscReal, xmax::$PetscReal, ymin::$PetscReal, ymax::$PetscReal )

    @chk ccall(
               (:PetscDrawAxisGetLimits, $petsc_library),
               PetscErrorCode,
               (PetscDrawAxis, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
               axis, xmin, xmax, ymin, ymax,
              )


	return nothing
end 

"""
	PetscDrawAxisSetHoldLimits(petsclib::PetscLibType,axis::PetscDrawAxis, hold::PetscBool) 
Causes an axis to keep the same limits until this is called
again

Logically Collective

Input Parameters:
- `axis` - the axis
- `hold` - `PETSC_TRUE` - hold current limits, `PETSC_FALSE` allow limits to be changed

Level: advanced

-seealso: `PetscDrawAxisCreate()`, `PetscDrawAxis`, `PetscDrawAxisGetLimits()`, `PetscDrawAxisSetLimits()`, `PetscDrawAxisSetLabels()`, `PetscDrawAxisSetColors()`

# External Links
$(_doc_external("Sys/PetscDrawAxisSetHoldLimits"))
"""
function PetscDrawAxisSetHoldLimits(petsclib::PetscLibType, axis::PetscDrawAxis, hold::PetscBool) end

@for_petsc function PetscDrawAxisSetHoldLimits(petsclib::$UnionPetscLib, axis::PetscDrawAxis, hold::PetscBool )

    @chk ccall(
               (:PetscDrawAxisSetHoldLimits, $petsc_library),
               PetscErrorCode,
               (PetscDrawAxis, PetscBool),
               axis, hold,
              )


	return nothing
end 

"""
	PetscDrawAxisDraw(petsclib::PetscLibType,axis::PetscDrawAxis) 
draws an axis.

Collective

Input Parameter:
- `axis` - `PetscDrawAxis` structure

Level: advanced

-seealso: `PetscDrawAxisCreate()`, `PetscDrawAxis`, `PetscDrawAxisGetLimits()`, `PetscDrawAxisSetLimits()`, `PetscDrawAxisSetLabels()`, `PetscDrawAxisSetColors()`

# External Links
$(_doc_external("Sys/PetscDrawAxisDraw"))
"""
function PetscDrawAxisDraw(petsclib::PetscLibType, axis::PetscDrawAxis) end

@for_petsc function PetscDrawAxisDraw(petsclib::$UnionPetscLib, axis::PetscDrawAxis )

    @chk ccall(
               (:PetscDrawAxisDraw, $petsc_library),
               PetscErrorCode,
               (PetscDrawAxis,),
               axis,
              )


	return nothing
end 

"""
	PetscDrawLGGetAxis(petsclib::PetscLibType,lg::PetscDrawLG, axis::PetscDrawAxis) 
Gets the axis context associated with a line graph.
This is useful if one wants to change some axis property, such as
labels, color, etc. The axis context should not be destroyed by the
application code.

Not Collective, if lg is parallel then axis is parallel

Input Parameter:
- `lg` - the line graph context

Output Parameter:
- `axis` - the axis context

Level: advanced

-seealso: `PetscDrawLGCreate()`, `PetscDrawAxis`, `PetscDrawLG`

# External Links
$(_doc_external("Sys/PetscDrawLGGetAxis"))
"""
function PetscDrawLGGetAxis(petsclib::PetscLibType, lg::PetscDrawLG, axis::PetscDrawAxis) end

@for_petsc function PetscDrawLGGetAxis(petsclib::$UnionPetscLib, lg::PetscDrawLG, axis::PetscDrawAxis )

    @chk ccall(
               (:PetscDrawLGGetAxis, $petsc_library),
               PetscErrorCode,
               (PetscDrawLG, Ptr{PetscDrawAxis}),
               lg, axis,
              )


	return nothing
end 

"""
	PetscDrawLGGetDraw(petsclib::PetscLibType,lg::PetscDrawLG, draw::PetscDraw) 
Gets the draw context associated with a line graph.

Not Collective, if lg is parallel then draw is parallel

Input Parameter:
- `lg` - the line graph context

Output Parameter:
- `draw` - the draw context

Level: intermediate

-seealso: `PetscDrawLGCreate()`, `PetscDraw`, `PetscDrawLG`

# External Links
$(_doc_external("Sys/PetscDrawLGGetDraw"))
"""
function PetscDrawLGGetDraw(petsclib::PetscLibType, lg::PetscDrawLG, draw::PetscDraw) end

@for_petsc function PetscDrawLGGetDraw(petsclib::$UnionPetscLib, lg::PetscDrawLG, draw::PetscDraw )

    @chk ccall(
               (:PetscDrawLGGetDraw, $petsc_library),
               PetscErrorCode,
               (PetscDrawLG, Ptr{PetscDraw}),
               lg, draw,
              )


	return nothing
end 

"""
	PetscDrawLGSPDraw(petsclib::PetscLibType,lg::PetscDrawLG, spin::PetscDrawSP) 
Redraws a line graph and a scatter plot on the same `PetscDraw` they must share

Collective

Input Parameters:
- `lg`   - the line graph context
- `spin` - the scatter plot

Level: intermediate

-seealso: `PetscDrawLGDraw()`, `PetscDrawSPDraw()`

# External Links
$(_doc_external("Sys/PetscDrawLGSPDraw"))
"""
function PetscDrawLGSPDraw(petsclib::PetscLibType, lg::PetscDrawLG, spin::PetscDrawSP) end

@for_petsc function PetscDrawLGSPDraw(petsclib::$UnionPetscLib, lg::PetscDrawLG, spin::PetscDrawSP )

    @chk ccall(
               (:PetscDrawLGSPDraw, $petsc_library),
               PetscErrorCode,
               (PetscDrawLG, PetscDrawSP),
               lg, spin,
              )


	return nothing
end 

"""
	outlg::PetscDrawLG = PetscDrawLGCreate(petsclib::PetscLibType,draw::PetscDraw, dim::PetscInt) 
Creates a line graph data structure.

Collective

Input Parameters:
- `draw` - the window where the graph will be made.
- `dim`  - the number of curves which will be drawn

Output Parameter:
- `outlg` - the line graph context

Level: intermediate

-seealso: `PetscDrawLGDestroy()`, `PetscDrawLGAddPoint()`, `PetscDrawLGAddCommonPoint()`, `PetscDrawLGAddPoints()`, `PetscDrawLGDraw()`, `PetscDrawLGSave()`,
`PetscDrawLGView()`, `PetscDrawLGReset()`, `PetscDrawLGSetDimension()`, `PetscDrawLGGetDimension()`, `PetscDrawLGSetLegend()`, `PetscDrawLGGetAxis()`,
`PetscDrawLGGetDraw()`, `PetscDrawLGSetUseMarkers()`, `PetscDrawLGSetLimits()`, `PetscDrawLGSetColors()`, `PetscDrawLGSetOptionsPrefix()`, `PetscDrawLGSetFromOptions()`

# External Links
$(_doc_external("Sys/PetscDrawLGCreate"))
"""
function PetscDrawLGCreate(petsclib::PetscLibType, draw::PetscDraw, dim::PetscInt) end

@for_petsc function PetscDrawLGCreate(petsclib::$UnionPetscLib, draw::PetscDraw, dim::$PetscInt )
	outlg_ = Ref{PetscDrawLG}()

    @chk ccall(
               (:PetscDrawLGCreate, $petsc_library),
               PetscErrorCode,
               (PetscDraw, $PetscInt, Ptr{PetscDrawLG}),
               draw, dim, outlg_,
              )

	outlg = outlg_[]

	return outlg
end 

"""
	PetscDrawLGSetColors(petsclib::PetscLibType,lg::PetscDrawLG, colors::Vector{Cint}) 
Sets the color of each line graph drawn

Logically Collective

Input Parameters:
- `lg`     - the line graph context.
- `colors` - the colors, an array of length the value set with `PetscDrawLGSetDimension()`

Level: intermediate

-seealso: `PetscDrawLG`, `PetscDrawLGCreate()`, `PetscDrawLGSetDimension()`, `PetscDrawLGGetDimension()`

# External Links
$(_doc_external("Sys/PetscDrawLGSetColors"))
"""
function PetscDrawLGSetColors(petsclib::PetscLibType, lg::PetscDrawLG, colors::Vector{Cint}) end

@for_petsc function PetscDrawLGSetColors(petsclib::$UnionPetscLib, lg::PetscDrawLG, colors::Vector{Cint} )

    @chk ccall(
               (:PetscDrawLGSetColors, $petsc_library),
               PetscErrorCode,
               (PetscDrawLG, Ptr{Cint}),
               lg, colors,
              )


	return nothing
end 

"""
	PetscDrawLGSetLegend(petsclib::PetscLibType,lg::PetscDrawLG, names::String) 
sets the names of each curve plotted

Logically Collective

Input Parameters:
- `lg`    - the line graph context.
- `names` - the names for each curve

Level: intermediate

-seealso: `PetscDrawLGGetAxis()`, `PetscDrawAxis`, `PetscDrawAxisSetColors()`, `PetscDrawAxisSetLabels()`, `PetscDrawAxisSetHoldLimits()`

# External Links
$(_doc_external("Sys/PetscDrawLGSetLegend"))
"""
function PetscDrawLGSetLegend(petsclib::PetscLibType, lg::PetscDrawLG, names::String) end

@for_petsc function PetscDrawLGSetLegend(petsclib::$UnionPetscLib, lg::PetscDrawLG, names::String )
	names_ = Ref(pointer(names))

    @chk ccall(
               (:PetscDrawLGSetLegend, $petsc_library),
               PetscErrorCode,
               (PetscDrawLG, Ptr{Ptr{Cchar}}),
               lg, names_,
              )


	return nothing
end 

"""
	dim::PetscInt = PetscDrawLGGetDimension(petsclib::PetscLibType,lg::PetscDrawLG) 
Get the number of curves that are to be drawn.

Not Collective

Input Parameter:
- `lg` - the line graph context.

Output Parameter:
- `dim` - the number of curves.

Level: intermediate

-seealso: `PetscDrawLGC`, `PetscDrawLGCreate()`, `PetscDrawLGSetDimension()`

# External Links
$(_doc_external("Sys/PetscDrawLGGetDimension"))
"""
function PetscDrawLGGetDimension(petsclib::PetscLibType, lg::PetscDrawLG) end

@for_petsc function PetscDrawLGGetDimension(petsclib::$UnionPetscLib, lg::PetscDrawLG )
	dim_ = Ref{$PetscInt}()

    @chk ccall(
               (:PetscDrawLGGetDimension, $petsc_library),
               PetscErrorCode,
               (PetscDrawLG, Ptr{$PetscInt}),
               lg, dim_,
              )

	dim = dim_[]

	return dim
end 

"""
	PetscDrawLGSetDimension(petsclib::PetscLibType,lg::PetscDrawLG, dim::PetscInt) 
Change the number of curves that are to be drawn.

Logically Collective

Input Parameters:
- `lg`  - the line graph context.
- `dim` - the number of curves.

Level: intermediate

-seealso: `PetscDrawLGCreate()`, `PetscDrawLGGetDimension()`

# External Links
$(_doc_external("Sys/PetscDrawLGSetDimension"))
"""
function PetscDrawLGSetDimension(petsclib::PetscLibType, lg::PetscDrawLG, dim::PetscInt) end

@for_petsc function PetscDrawLGSetDimension(petsclib::$UnionPetscLib, lg::PetscDrawLG, dim::$PetscInt )

    @chk ccall(
               (:PetscDrawLGSetDimension, $petsc_library),
               PetscErrorCode,
               (PetscDrawLG, $PetscInt),
               lg, dim,
              )


	return nothing
end 

"""
	dim::PetscInt,n::PetscInt,x::Vector{PetscReal},y::Vector{PetscReal} = PetscDrawLGGetData(petsclib::PetscLibType,lg::PetscDrawLG) 
Get the data being plotted.

Not Collective

Input Parameter:
- `lg` - the line graph context

Output Parameters:
- `dim` - the number of curves
- `n`   - the number of points on each line
- `x`   - The x-value of each point, x[p * dim + c]
- `y`   - The y-value of each point, y[p * dim + c]

Level: intermediate

-seealso: `PetscDrawLGC`, `PetscDrawLGCreate()`, `PetscDrawLGGetDimension()`

# External Links
$(_doc_external("Sys/PetscDrawLGGetData"))
"""
function PetscDrawLGGetData(petsclib::PetscLibType, lg::PetscDrawLG) end

@for_petsc function PetscDrawLGGetData(petsclib::$UnionPetscLib, lg::PetscDrawLG )
	dim_ = Ref{$PetscInt}()
	n_ = Ref{$PetscInt}()
	x_ = Ref{Ptr{$PetscReal}}()
	y_ = Ref{Ptr{$PetscReal}}()

    @chk ccall(
               (:PetscDrawLGGetData, $petsc_library),
               PetscErrorCode,
               (PetscDrawLG, Ptr{$PetscInt}, Ptr{$PetscInt}, Ptr{Ptr{$PetscReal}}, Ptr{Ptr{$PetscReal}}),
               lg, dim_, n_, x_, y_,
              )

	dim = dim_[]
	n = n_[]
	x = unsafe_wrap(Array, x_[], VecGetLocalSize(petsclib, x); own = false)
	y = unsafe_wrap(Array, y_[], VecGetLocalSize(petsclib, x); own = false)

	return dim,n,x,y
end 

"""
	PetscDrawLGSetLimits(petsclib::PetscLibType,lg::PetscDrawLG, x_min::PetscReal, x_max::PetscReal, y_min::PetscReal, y_max::PetscReal) 
Sets the axis limits for a line graph. If more
points are added after this call, the limits will be adjusted to
include those additional points.

Logically Collective

Input Parameters:
- `lg`    - the line graph context
- `x_min` - the horizontal lower limit
- `x_max` - the horizontal upper limit
- `y_min` - the vertical lower limit
- `y_max` - the vertical upper limit

Level: intermediate

-seealso: `PetscDrawLGCreate()`, `PetscDrawLG`, `PetscDrawAxis`

# External Links
$(_doc_external("Sys/PetscDrawLGSetLimits"))
"""
function PetscDrawLGSetLimits(petsclib::PetscLibType, lg::PetscDrawLG, x_min::PetscReal, x_max::PetscReal, y_min::PetscReal, y_max::PetscReal) end

@for_petsc function PetscDrawLGSetLimits(petsclib::$UnionPetscLib, lg::PetscDrawLG, x_min::$PetscReal, x_max::$PetscReal, y_min::$PetscReal, y_max::$PetscReal )

    @chk ccall(
               (:PetscDrawLGSetLimits, $petsc_library),
               PetscErrorCode,
               (PetscDrawLG, $PetscReal, $PetscReal, $PetscReal, $PetscReal),
               lg, x_min, x_max, y_min, y_max,
              )


	return nothing
end 

"""
	PetscDrawLGReset(petsclib::PetscLibType,lg::PetscDrawLG) 
Clears line graph to allow for reuse with new data.

Logically Collective

Input Parameter:
- `lg` - the line graph context.

Level: intermediate

-seealso: `PetscDrawLG`, `PetscDrawLGCreate()`

# External Links
$(_doc_external("Sys/PetscDrawLGReset"))
"""
function PetscDrawLGReset(petsclib::PetscLibType, lg::PetscDrawLG) end

@for_petsc function PetscDrawLGReset(petsclib::$UnionPetscLib, lg::PetscDrawLG )

    @chk ccall(
               (:PetscDrawLGReset, $petsc_library),
               PetscErrorCode,
               (PetscDrawLG,),
               lg,
              )


	return nothing
end 

"""
	PetscDrawLGDestroy(petsclib::PetscLibType,lg::PetscDrawLG) 
Frees all space taken up by line graph data structure.

Collective

Input Parameter:
- `lg` - the line graph context

Level: intermediate

-seealso: `PetscDrawLG`, `PetscDrawLGCreate()`

# External Links
$(_doc_external("Sys/PetscDrawLGDestroy"))
"""
function PetscDrawLGDestroy(petsclib::PetscLibType, lg::PetscDrawLG) end

@for_petsc function PetscDrawLGDestroy(petsclib::$UnionPetscLib, lg::PetscDrawLG )

    @chk ccall(
               (:PetscDrawLGDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscDrawLG},),
               lg,
              )


	return nothing
end 

"""
	PetscDrawLGSetUseMarkers(petsclib::PetscLibType,lg::PetscDrawLG, flg::PetscBool) 
Causes the line graph object to draw a marker for each data

Logically Collective

Input Parameters:
- `lg`  - the linegraph context
- `flg` - should mark each data point

Options Database Key:
- `-lg_use_markers  <true,false>` - true means it draws a marker for each point

Level: intermediate

-seealso: `PetscDrawLG`, `PetscDrawLGCreate()`

# External Links
$(_doc_external("Sys/PetscDrawLGSetUseMarkers"))
"""
function PetscDrawLGSetUseMarkers(petsclib::PetscLibType, lg::PetscDrawLG, flg::PetscBool) end

@for_petsc function PetscDrawLGSetUseMarkers(petsclib::$UnionPetscLib, lg::PetscDrawLG, flg::PetscBool )

    @chk ccall(
               (:PetscDrawLGSetUseMarkers, $petsc_library),
               PetscErrorCode,
               (PetscDrawLG, PetscBool),
               lg, flg,
              )


	return nothing
end 

"""
	PetscDrawLGDraw(petsclib::PetscLibType,lg::PetscDrawLG) 
Redraws a line graph.

Collective

Input Parameter:
- `lg` - the line graph context

Level: intermediate

-seealso: `PetscDrawLG`, `PetscDrawSPDraw()`, `PetscDrawLGSPDraw()`, `PetscDrawLGReset()`

# External Links
$(_doc_external("Sys/PetscDrawLGDraw"))
"""
function PetscDrawLGDraw(petsclib::PetscLibType, lg::PetscDrawLG) end

@for_petsc function PetscDrawLGDraw(petsclib::$UnionPetscLib, lg::PetscDrawLG )

    @chk ccall(
               (:PetscDrawLGDraw, $petsc_library),
               PetscErrorCode,
               (PetscDrawLG,),
               lg,
              )


	return nothing
end 

"""
	PetscDrawLGSave(petsclib::PetscLibType,lg::PetscDrawLG) 
Saves a drawn image

Collective

Input Parameter:
- `lg` - The line graph context

Level: intermediate

-seealso: `PetscDrawLG`, `PetscDrawSave()`, `PetscDrawLGCreate()`, `PetscDrawLGGetDraw()`, `PetscDrawSetSave()`

# External Links
$(_doc_external("Sys/PetscDrawLGSave"))
"""
function PetscDrawLGSave(petsclib::PetscLibType, lg::PetscDrawLG) end

@for_petsc function PetscDrawLGSave(petsclib::$UnionPetscLib, lg::PetscDrawLG )

    @chk ccall(
               (:PetscDrawLGSave, $petsc_library),
               PetscErrorCode,
               (PetscDrawLG,),
               lg,
              )


	return nothing
end 

"""
	PetscDrawLGView(petsclib::PetscLibType,lg::PetscDrawLG, viewer::PetscViewer) 
Prints a line graph.

Collective

Input Parameters:
- `lg`     - the line graph context
- `viewer` - the viewer to view it with

Level: beginner

-seealso: `PetscDrawLG`, `PetscDrawLGCreate()`

# External Links
$(_doc_external("Sys/PetscDrawLGView"))
"""
function PetscDrawLGView(petsclib::PetscLibType, lg::PetscDrawLG, viewer::PetscViewer) end

@for_petsc function PetscDrawLGView(petsclib::$UnionPetscLib, lg::PetscDrawLG, viewer::PetscViewer )

    @chk ccall(
               (:PetscDrawLGView, $petsc_library),
               PetscErrorCode,
               (PetscDrawLG, PetscViewer),
               lg, viewer,
              )


	return nothing
end 

"""
	PetscDrawLGSetOptionsPrefix(petsclib::PetscLibType,lg::PetscDrawLG, prefix::String) 
Sets the prefix used for searching for all
`PetscDrawLG` options in the database.

Logically Collective

Input Parameters:
- `lg`     - the line graph context
- `prefix` - the prefix to prepend to all option names

Level: advanced

-seealso: `PetscDrawLG`, `PetscDrawLGSetFromOptions()`, `PetscDrawLGCreate()`

# External Links
$(_doc_external("Sys/PetscDrawLGSetOptionsPrefix"))
"""
function PetscDrawLGSetOptionsPrefix(petsclib::PetscLibType, lg::PetscDrawLG, prefix::String) end

@for_petsc function PetscDrawLGSetOptionsPrefix(petsclib::$UnionPetscLib, lg::PetscDrawLG, prefix::String )

    @chk ccall(
               (:PetscDrawLGSetOptionsPrefix, $petsc_library),
               PetscErrorCode,
               (PetscDrawLG, Ptr{Cchar}),
               lg, prefix,
              )


	return nothing
end 

"""
	PetscDrawLGSetFromOptions(petsclib::PetscLibType,lg::PetscDrawLG) 
Sets options related to the line graph object

Collective

Input Parameters:
- `lg` - the line graph context

Options Database Key:
- `-lg_use_markers  <true,false>` - true means it draws a marker for each point

Level: intermediate

-seealso: `PetscDrawLG`, `PetscDrawLGDestroy()`, `PetscDrawLGCreate()`

# External Links
$(_doc_external("Sys/PetscDrawLGSetFromOptions"))
"""
function PetscDrawLGSetFromOptions(petsclib::PetscLibType, lg::PetscDrawLG) end

@for_petsc function PetscDrawLGSetFromOptions(petsclib::$UnionPetscLib, lg::PetscDrawLG )

    @chk ccall(
               (:PetscDrawLGSetFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscDrawLG,),
               lg,
              )


	return nothing
end 

"""
	PetscDrawLGAddCommonPoint(petsclib::PetscLibType,lg::PetscDrawLG, x::PetscReal, y::PetscReal) 
Adds another point to each of the line graphs. All the points share
the same new X coordinate.  The new point must have an X coordinate larger than the old points.

Logically Collective

Input Parameters:
- `lg` - the line graph context
- `x`  - the common x coordinate point
- `y`  - the new y coordinate point for each curve.

Level: intermediate

-seealso: `PetscDrawLG`, `PetscDrawLGCreate()`, `PetscDrawLGAddPoints()`, `PetscDrawLGAddPoint()`, `PetscDrawLGReset()`, `PetscDrawLGDraw()`

# External Links
$(_doc_external("Sys/PetscDrawLGAddCommonPoint"))
"""
function PetscDrawLGAddCommonPoint(petsclib::PetscLibType, lg::PetscDrawLG, x::PetscReal, y::PetscReal) end

@for_petsc function PetscDrawLGAddCommonPoint(petsclib::$UnionPetscLib, lg::PetscDrawLG, x::$PetscReal, y::$PetscReal )

    @chk ccall(
               (:PetscDrawLGAddCommonPoint, $petsc_library),
               PetscErrorCode,
               (PetscDrawLG, $PetscReal, Ptr{$PetscReal}),
               lg, x, y,
              )


	return nothing
end 

"""
	PetscDrawLGAddPoint(petsclib::PetscLibType,lg::PetscDrawLG, x::PetscReal, y::PetscReal) 
Adds another point to each of the line graphs.
The new point must have an X coordinate larger than the old points.

Logically Collective

Input Parameters:
- `lg` - the line graph context
- `x`  - array containing the x coordinate for the point on each curve
- `y`  - array containing the y coordinate for the point on each curve

Level: intermediate

-seealso: `PetscDrawLG`, `PetscDrawLGCreate()`, `PetscDrawLGAddPoints()`, `PetscDrawLGAddCommonPoint()`, `PetscDrawLGReset()`, `PetscDrawLGDraw()`

# External Links
$(_doc_external("Sys/PetscDrawLGAddPoint"))
"""
function PetscDrawLGAddPoint(petsclib::PetscLibType, lg::PetscDrawLG, x::PetscReal, y::PetscReal) end

@for_petsc function PetscDrawLGAddPoint(petsclib::$UnionPetscLib, lg::PetscDrawLG, x::$PetscReal, y::$PetscReal )

    @chk ccall(
               (:PetscDrawLGAddPoint, $petsc_library),
               PetscErrorCode,
               (PetscDrawLG, Ptr{$PetscReal}, Ptr{$PetscReal}),
               lg, x, y,
              )


	return nothing
end 

"""
	PetscDrawLGAddPoints(petsclib::PetscLibType,lg::PetscDrawLG, n::PetscInt, xx::Vector{PetscReal}, yy::Vector{PetscReal}) 
Adds several points to each of the line graphs.
The new points must have an X coordinate larger than the old points.

Logically Collective

Input Parameters:
- `lg` - the line graph context
- `xx` - array of pointers that point to arrays containing the new x coordinates for each curve.
- `yy` - array of pointers that point to arrays containing the new y points for each curve.
- `n`  - number of points being added

Level: intermediate

-seealso: `PetscDrawLG`, `PetscDrawLGCreate()`, `PetscDrawLGAddPoint()`, `PetscDrawLGAddCommonPoint()`, `PetscDrawLGReset()`, `PetscDrawLGDraw()`

# External Links
$(_doc_external("Sys/PetscDrawLGAddPoints"))
"""
function PetscDrawLGAddPoints(petsclib::PetscLibType, lg::PetscDrawLG, n::PetscInt, xx::Vector{PetscReal}, yy::Vector{PetscReal}) end

@for_petsc function PetscDrawLGAddPoints(petsclib::$UnionPetscLib, lg::PetscDrawLG, n::$PetscInt, xx::Vector{$PetscReal}, yy::Vector{$PetscReal} )
	xx_ = Ref(pointer(xx))
	yy_ = Ref(pointer(yy))

    @chk ccall(
               (:PetscDrawLGAddPoints, $petsc_library),
               PetscErrorCode,
               (PetscDrawLG, $PetscInt, Ptr{Ptr{$PetscReal}}, Ptr{Ptr{$PetscReal}}),
               lg, n, xx_, yy_,
              )


	return nothing
end 

"""
	PetscDrawSplitViewPort(petsclib::PetscLibType,draw::PetscDraw) 
Splits a window shared by several processes into smaller
view ports. One for each process.

Collective

Input Parameter:
- `draw` - the drawing context

Level: advanced

-seealso: `PetscDrawDivideViewPort()`, `PetscDrawSetViewPort()`

# External Links
$(_doc_external("Sys/PetscDrawSplitViewPort"))
"""
function PetscDrawSplitViewPort(petsclib::PetscLibType, draw::PetscDraw) end

@for_petsc function PetscDrawSplitViewPort(petsclib::$UnionPetscLib, draw::PetscDraw )

    @chk ccall(
               (:PetscDrawSplitViewPort, $petsc_library),
               PetscErrorCode,
               (PetscDraw,),
               draw,
              )


	return nothing
end 

"""
	drawsp::PetscDrawSP = PetscDrawSPCreate(petsclib::PetscLibType,draw::PetscDraw, dim::Cint) 
Creates a scatter plot data structure.

Collective

Input Parameters:
- `draw` - the window where the graph will be made.
- `dim`  - the number of sets of points which will be drawn

Output Parameter:
- `drawsp` - the scatter plot context

Level: intermediate

-seealso: `PetscDrawLGCreate()`, `PetscDrawLG`, `PetscDrawBarCreate()`, `PetscDrawBar`, `PetscDrawHGCreate()`, `PetscDrawHG`, `PetscDrawSPDestroy()`, `PetscDraw`, `PetscDrawSP`, `PetscDrawSPSetDimension()`, `PetscDrawSPReset()`,
`PetscDrawSPAddPoint()`, `PetscDrawSPAddPoints()`, `PetscDrawSPDraw()`, `PetscDrawSPSave()`, `PetscDrawSPSetLimits()`, `PetscDrawSPGetAxis()`, `PetscDrawAxis`, `PetscDrawSPGetDraw()`

# External Links
$(_doc_external("Sys/PetscDrawSPCreate"))
"""
function PetscDrawSPCreate(petsclib::PetscLibType, draw::PetscDraw, dim::Cint) end

@for_petsc function PetscDrawSPCreate(petsclib::$UnionPetscLib, draw::PetscDraw, dim::Cint )
	drawsp_ = Ref{PetscDrawSP}()

    @chk ccall(
               (:PetscDrawSPCreate, $petsc_library),
               PetscErrorCode,
               (PetscDraw, Cint, Ptr{PetscDrawSP}),
               draw, dim, drawsp_,
              )

	drawsp = drawsp_[]

	return drawsp
end 

"""
	PetscDrawSPSetDimension(petsclib::PetscLibType,sp::PetscDrawSP, dim::Cint) 
Change the number of points that are added at each  `PetscDrawSPAddPoint()`

Not Collective

Input Parameters:
- `sp`  - the scatter plot context.
- `dim` - the number of point curves on this process

Level: intermediate

-seealso: `PetscDrawSP`, `PetscDrawSPCreate()`, `PetscDrawSPAddPoint()`, `PetscDrawSPAddPoints()`

# External Links
$(_doc_external("Sys/PetscDrawSPSetDimension"))
"""
function PetscDrawSPSetDimension(petsclib::PetscLibType, sp::PetscDrawSP, dim::Cint) end

@for_petsc function PetscDrawSPSetDimension(petsclib::$UnionPetscLib, sp::PetscDrawSP, dim::Cint )

    @chk ccall(
               (:PetscDrawSPSetDimension, $petsc_library),
               PetscErrorCode,
               (PetscDrawSP, Cint),
               sp, dim,
              )


	return nothing
end 

"""
	PetscDrawSPGetDimension(petsclib::PetscLibType,sp::PetscDrawSP, dim::Cint) 
Get the number of sets of points that are to be drawn at each `PetscDrawSPAddPoint()`

Not Collective

Input Parameter:
- `sp` - the scatter plot context.

Output Parameter:
- `dim` - the number of point curves on this process

Level: intermediate

-seealso: `PetscDrawSP`, `PetscDrawSPCreate()`, `PetscDrawSPAddPoint()`, `PetscDrawSPAddPoints()`

# External Links
$(_doc_external("Sys/PetscDrawSPGetDimension"))
"""
function PetscDrawSPGetDimension(petsclib::PetscLibType, sp::PetscDrawSP, dim::Cint) end

@for_petsc function PetscDrawSPGetDimension(petsclib::$UnionPetscLib, sp::PetscDrawSP, dim::Cint )

    @chk ccall(
               (:PetscDrawSPGetDimension, $petsc_library),
               PetscErrorCode,
               (PetscDrawSP, Ptr{Cint}),
               sp, dim,
              )


	return nothing
end 

"""
	PetscDrawSPReset(petsclib::PetscLibType,sp::PetscDrawSP) 
Clears scatter plot to allow for reuse with new data.

Not Collective

Input Parameter:
- `sp` - the scatter plot context.

Level: intermediate

-seealso: `PetscDrawSP`, `PetscDrawSPCreate()`, `PetscDrawSPAddPoint()`, `PetscDrawSPAddPoints()`, `PetscDrawSPDraw()`

# External Links
$(_doc_external("Sys/PetscDrawSPReset"))
"""
function PetscDrawSPReset(petsclib::PetscLibType, sp::PetscDrawSP) end

@for_petsc function PetscDrawSPReset(petsclib::$UnionPetscLib, sp::PetscDrawSP )

    @chk ccall(
               (:PetscDrawSPReset, $petsc_library),
               PetscErrorCode,
               (PetscDrawSP,),
               sp,
              )


	return nothing
end 

"""
	PetscDrawSPDestroy(petsclib::PetscLibType,sp::PetscDrawSP) 
Frees all space taken up by scatter plot data structure.

Collective

Input Parameter:
- `sp` - the scatter plot context

Level: intermediate

-seealso: `PetscDrawSPCreate()`, `PetscDrawSP`, `PetscDrawSPReset()`

# External Links
$(_doc_external("Sys/PetscDrawSPDestroy"))
"""
function PetscDrawSPDestroy(petsclib::PetscLibType, sp::PetscDrawSP) end

@for_petsc function PetscDrawSPDestroy(petsclib::$UnionPetscLib, sp::PetscDrawSP )

    @chk ccall(
               (:PetscDrawSPDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscDrawSP},),
               sp,
              )


	return nothing
end 

"""
	PetscDrawSPAddPoint(petsclib::PetscLibType,sp::PetscDrawSP, x::PetscReal, y::PetscReal) 
Adds another point to each of the scatter plot point curves.

Not Collective

Input Parameters:
- `sp` - the scatter plot data structure
- `x`  - the x coordinate values (of length dim) for the points of the curve
- `y`  - the y coordinate values (of length dim) for the points of the curve

Level: intermediate

-seealso: `PetscDrawSPAddPoints()`, `PetscDrawSP`, `PetscDrawSPCreate()`, `PetscDrawSPReset()`, `PetscDrawSPDraw()`, `PetscDrawSPAddPointColorized()`

# External Links
$(_doc_external("Sys/PetscDrawSPAddPoint"))
"""
function PetscDrawSPAddPoint(petsclib::PetscLibType, sp::PetscDrawSP, x::PetscReal, y::PetscReal) end

@for_petsc function PetscDrawSPAddPoint(petsclib::$UnionPetscLib, sp::PetscDrawSP, x::$PetscReal, y::$PetscReal )

    @chk ccall(
               (:PetscDrawSPAddPoint, $petsc_library),
               PetscErrorCode,
               (PetscDrawSP, Ptr{$PetscReal}, Ptr{$PetscReal}),
               sp, x, y,
              )


	return nothing
end 

"""
	PetscDrawSPAddPoints(petsclib::PetscLibType,sp::PetscDrawSP, n::Cint, xx::Vector{PetscReal}, yy::Vector{PetscReal}) 
Adds several points to each of the scatter plot point curves.

Not Collective

Input Parameters:
- `sp` - the scatter plot context
- `xx` - array of pointers that point to arrays containing the new x coordinates for each curve.
- `yy` - array of pointers that point to arrays containing the new y points for each curve.
- `n`  - number of points being added, each represents a subarray of length dim where dim is the value from `PetscDrawSPGetDimension()`

Level: intermediate

-seealso: `PetscDrawSPAddPoint()`, `PetscDrawSP`, `PetscDrawSPCreate()`, `PetscDrawSPReset()`, `PetscDrawSPDraw()`, `PetscDrawSPAddPointColorized()`

# External Links
$(_doc_external("Sys/PetscDrawSPAddPoints"))
"""
function PetscDrawSPAddPoints(petsclib::PetscLibType, sp::PetscDrawSP, n::Cint, xx::Vector{PetscReal}, yy::Vector{PetscReal}) end

@for_petsc function PetscDrawSPAddPoints(petsclib::$UnionPetscLib, sp::PetscDrawSP, n::Cint, xx::Vector{$PetscReal}, yy::Vector{$PetscReal} )
	xx_ = Ref(pointer(xx))
	yy_ = Ref(pointer(yy))

    @chk ccall(
               (:PetscDrawSPAddPoints, $petsc_library),
               PetscErrorCode,
               (PetscDrawSP, Cint, Ptr{Ptr{$PetscReal}}, Ptr{Ptr{$PetscReal}}),
               sp, n, xx_, yy_,
              )


	return nothing
end 

"""
	PetscDrawSPAddPointColorized(petsclib::PetscLibType,sp::PetscDrawSP, x::PetscReal, y::PetscReal, z::PetscReal) 
Adds another point to each of the scatter plots as well as a numeric value to be used to colorize the scatter point.

Not Collective

Input Parameters:
- `sp` - the scatter plot data structure
- `x`  - array of length dim containing the new x coordinate values for each of the point curves.
- `y`  - array of length dim containing the new y coordinate values for each of the point curves.
- `z`  - array of length dim containing the numeric values that will be mapped to [0,255] and used for scatter point colors.

Level: intermediate

-seealso: `PetscDrawSPAddPoints()`, `PetscDrawSP`, `PetscDrawSPCreate()`, `PetscDrawSPReset()`, `PetscDrawSPDraw()`, `PetscDrawSPAddPoint()`

# External Links
$(_doc_external("Sys/PetscDrawSPAddPointColorized"))
"""
function PetscDrawSPAddPointColorized(petsclib::PetscLibType, sp::PetscDrawSP, x::PetscReal, y::PetscReal, z::PetscReal) end

@for_petsc function PetscDrawSPAddPointColorized(petsclib::$UnionPetscLib, sp::PetscDrawSP, x::$PetscReal, y::$PetscReal, z::$PetscReal )

    @chk ccall(
               (:PetscDrawSPAddPointColorized, $petsc_library),
               PetscErrorCode,
               (PetscDrawSP, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
               sp, x, y, z,
              )


	return nothing
end 

"""
	PetscDrawSPDraw(petsclib::PetscLibType,sp::PetscDrawSP, clear::PetscBool) 
Redraws a scatter plot.

Collective

Input Parameters:
- `sp`    - the scatter plot context
- `clear` - clear the window before drawing the new plot

Level: intermediate

-seealso: `PetscDrawLGDraw()`, `PetscDrawLGSPDraw()`, `PetscDrawSP`, `PetscDrawSPCreate()`, `PetscDrawSPReset()`, `PetscDrawSPAddPoint()`, `PetscDrawSPAddPoints()`

# External Links
$(_doc_external("Sys/PetscDrawSPDraw"))
"""
function PetscDrawSPDraw(petsclib::PetscLibType, sp::PetscDrawSP, clear::PetscBool) end

@for_petsc function PetscDrawSPDraw(petsclib::$UnionPetscLib, sp::PetscDrawSP, clear::PetscBool )

    @chk ccall(
               (:PetscDrawSPDraw, $petsc_library),
               PetscErrorCode,
               (PetscDrawSP, PetscBool),
               sp, clear,
              )


	return nothing
end 

"""
	PetscDrawSPSave(petsclib::PetscLibType,sp::PetscDrawSP) 
Saves a drawn image

Collective

Input Parameter:
- `sp` - the scatter plot context

Level: intermediate

-seealso: `PetscDrawSPCreate()`, `PetscDrawSPGetDraw()`, `PetscDrawSetSave()`, `PetscDrawSave()`

# External Links
$(_doc_external("Sys/PetscDrawSPSave"))
"""
function PetscDrawSPSave(petsclib::PetscLibType, sp::PetscDrawSP) end

@for_petsc function PetscDrawSPSave(petsclib::$UnionPetscLib, sp::PetscDrawSP )

    @chk ccall(
               (:PetscDrawSPSave, $petsc_library),
               PetscErrorCode,
               (PetscDrawSP,),
               sp,
              )


	return nothing
end 

"""
	PetscDrawSPSetLimits(petsclib::PetscLibType,sp::PetscDrawSP, x_min::PetscReal, x_max::PetscReal, y_min::PetscReal, y_max::PetscReal) 
Sets the axis limits for a scatter plot. If more points are added after this call, the limits will be adjusted to include those additional points.

Not Collective

Input Parameters:
- `sp`    - the line graph context
- `x_min` - the horizontal lower limit
- `x_max` - the horizontal upper limit
- `y_min` - the vertical lower limit
- `y_max` - the vertical upper limit

Level: intermediate

-seealso: `PetscDrawSP`, `PetscDrawAxis`, `PetscDrawSPCreate()`, `PetscDrawSPDraw()`, `PetscDrawSPAddPoint()`, `PetscDrawSPAddPoints()`, `PetscDrawSPGetAxis()`

# External Links
$(_doc_external("Sys/PetscDrawSPSetLimits"))
"""
function PetscDrawSPSetLimits(petsclib::PetscLibType, sp::PetscDrawSP, x_min::PetscReal, x_max::PetscReal, y_min::PetscReal, y_max::PetscReal) end

@for_petsc function PetscDrawSPSetLimits(petsclib::$UnionPetscLib, sp::PetscDrawSP, x_min::$PetscReal, x_max::$PetscReal, y_min::$PetscReal, y_max::$PetscReal )

    @chk ccall(
               (:PetscDrawSPSetLimits, $petsc_library),
               PetscErrorCode,
               (PetscDrawSP, $PetscReal, $PetscReal, $PetscReal, $PetscReal),
               sp, x_min, x_max, y_min, y_max,
              )


	return nothing
end 

"""
	PetscDrawSPGetAxis(petsclib::PetscLibType,sp::PetscDrawSP, axis::PetscDrawAxis) 
Gets the axis context associated with a scatter plot

Not Collective

Input Parameter:
- `sp` - the scatter plot context

Output Parameter:
- `axis` - the axis context

Level: intermediate

-seealso: `PetscDrawSP`, `PetscDrawSPCreate()`, `PetscDrawSPDraw()`, `PetscDrawSPAddPoint()`, `PetscDrawSPAddPoints()`, `PetscDrawAxis`, `PetscDrawAxisCreate()`

# External Links
$(_doc_external("Sys/PetscDrawSPGetAxis"))
"""
function PetscDrawSPGetAxis(petsclib::PetscLibType, sp::PetscDrawSP, axis::PetscDrawAxis) end

@for_petsc function PetscDrawSPGetAxis(petsclib::$UnionPetscLib, sp::PetscDrawSP, axis::PetscDrawAxis )

    @chk ccall(
               (:PetscDrawSPGetAxis, $petsc_library),
               PetscErrorCode,
               (PetscDrawSP, Ptr{PetscDrawAxis}),
               sp, axis,
              )


	return nothing
end 

"""
	PetscDrawSPGetDraw(petsclib::PetscLibType,sp::PetscDrawSP, draw::PetscDraw) 
Gets the draw context associated with a scatter plot

Not Collective

Input Parameter:
- `sp` - the scatter plot context

Output Parameter:
- `draw` - the draw context

Level: intermediate

-seealso: `PetscDrawSP`, `PetscDrawSPCreate()`, `PetscDrawSPDraw()`, `PetscDraw`

# External Links
$(_doc_external("Sys/PetscDrawSPGetDraw"))
"""
function PetscDrawSPGetDraw(petsclib::PetscLibType, sp::PetscDrawSP, draw::PetscDraw) end

@for_petsc function PetscDrawSPGetDraw(petsclib::$UnionPetscLib, sp::PetscDrawSP, draw::PetscDraw )

    @chk ccall(
               (:PetscDrawSPGetDraw, $petsc_library),
               PetscErrorCode,
               (PetscDrawSP, Ptr{PetscDraw}),
               sp, draw,
              )


	return nothing
end 

"""
	hist::PetscDrawHG = PetscDrawHGCreate(petsclib::PetscLibType,draw::PetscDraw, bins::Cint) 
Creates a histogram data structure.

Collective

Input Parameters:
- `draw` - The window where the graph will be made
- `bins` - The number of bins to use

Output Parameter:
- `hist` - The histogram context

Level: intermediate

-seealso: `PetscDrawHGDestroy()`, `PetscDrawHG`, `PetscDrawBarCreate()`, `PetscDrawBar`, `PetscDrawLGCreate()`, `PetscDrawLG`, `PetscDrawSPCreate()`, `PetscDrawSP`,
`PetscDrawHGSetNumberBins()`, `PetscDrawHGReset()`, `PetscDrawHGAddValue()`, `PetscDrawHGDraw()`, `PetscDrawHGSave()`, `PetscDrawHGView()`, `PetscDrawHGSetColor()`,
`PetscDrawHGSetLimits()`, `PetscDrawHGCalcStats()`, `PetscDrawHGIntegerBins()`, `PetscDrawHGGetAxis()`, `PetscDrawAxis`, `PetscDrawHGGetDraw()`

# External Links
$(_doc_external("Sys/PetscDrawHGCreate"))
"""
function PetscDrawHGCreate(petsclib::PetscLibType, draw::PetscDraw, bins::Cint) end

@for_petsc function PetscDrawHGCreate(petsclib::$UnionPetscLib, draw::PetscDraw, bins::Cint )
	hist_ = Ref{PetscDrawHG}()

    @chk ccall(
               (:PetscDrawHGCreate, $petsc_library),
               PetscErrorCode,
               (PetscDraw, Cint, Ptr{PetscDrawHG}),
               draw, bins, hist_,
              )

	hist = hist_[]

	return hist
end 

"""
	PetscDrawHGSetNumberBins(petsclib::PetscLibType,hist::PetscDrawHG, bins::Cint) 
Change the number of bins that are to be drawn in the histogram

Logically Collective

Input Parameters:
- `hist` - The histogram context.
- `bins` - The number of bins.

Level: intermediate

-seealso: `PetscDrawHGCreate()`, `PetscDrawHG`, `PetscDrawHGDraw()`, `PetscDrawHGIntegerBins()`

# External Links
$(_doc_external("Sys/PetscDrawHGSetNumberBins"))
"""
function PetscDrawHGSetNumberBins(petsclib::PetscLibType, hist::PetscDrawHG, bins::Cint) end

@for_petsc function PetscDrawHGSetNumberBins(petsclib::$UnionPetscLib, hist::PetscDrawHG, bins::Cint )

    @chk ccall(
               (:PetscDrawHGSetNumberBins, $petsc_library),
               PetscErrorCode,
               (PetscDrawHG, Cint),
               hist, bins,
              )


	return nothing
end 

"""
	PetscDrawHGReset(petsclib::PetscLibType,hist::PetscDrawHG) 
Clears histogram to allow for reuse with new data.

Logically Collective

Input Parameter:
- `hist` - The histogram context.

Level: intermediate

-seealso: `PetscDrawHGCreate()`, `PetscDrawHG`, `PetscDrawHGDraw()`, `PetscDrawHGAddValue()`

# External Links
$(_doc_external("Sys/PetscDrawHGReset"))
"""
function PetscDrawHGReset(petsclib::PetscLibType, hist::PetscDrawHG) end

@for_petsc function PetscDrawHGReset(petsclib::$UnionPetscLib, hist::PetscDrawHG )

    @chk ccall(
               (:PetscDrawHGReset, $petsc_library),
               PetscErrorCode,
               (PetscDrawHG,),
               hist,
              )


	return nothing
end 

"""
	PetscDrawHGDestroy(petsclib::PetscLibType,hist::PetscDrawHG) 
Frees all space taken up by histogram data structure.

Collective

Input Parameter:
- `hist` - The histogram context

Level: intermediate

-seealso: `PetscDrawHGCreate()`, `PetscDrawHG`

# External Links
$(_doc_external("Sys/PetscDrawHGDestroy"))
"""
function PetscDrawHGDestroy(petsclib::PetscLibType, hist::PetscDrawHG) end

@for_petsc function PetscDrawHGDestroy(petsclib::$UnionPetscLib, hist::PetscDrawHG )

    @chk ccall(
               (:PetscDrawHGDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscDrawHG},),
               hist,
              )


	return nothing
end 

"""
	PetscDrawHGAddValue(petsclib::PetscLibType,hist::PetscDrawHG, value::PetscReal) 
Adds another value to the histogram.

Logically Collective

Input Parameters:
- `hist`  - The histogram
- `value` - The value

Level: intermediate

-seealso: `PetscDrawHGCreate()`, `PetscDrawHG`, `PetscDrawHGDraw()`, `PetscDrawHGReset()`, `PetscDrawHGAddWeightedValue()`

# External Links
$(_doc_external("Sys/PetscDrawHGAddValue"))
"""
function PetscDrawHGAddValue(petsclib::PetscLibType, hist::PetscDrawHG, value::PetscReal) end

@for_petsc function PetscDrawHGAddValue(petsclib::$UnionPetscLib, hist::PetscDrawHG, value::$PetscReal )

    @chk ccall(
               (:PetscDrawHGAddValue, $petsc_library),
               PetscErrorCode,
               (PetscDrawHG, $PetscReal),
               hist, value,
              )


	return nothing
end 

"""
	PetscDrawHGAddWeightedValue(petsclib::PetscLibType,hist::PetscDrawHG, value::PetscReal, weight::PetscReal) 
Adds another value to the histogram with a weight.

Logically Collective

Input Parameters:
- `hist`   - The histogram
- `value`  - The value
- `weight` - The value weight

Level: intermediate

-seealso: `PetscDrawHGCreate()`, `PetscDrawHG`, `PetscDrawHGDraw()`, `PetscDrawHGReset()`, `PetscDrawHGAddValue()`

# External Links
$(_doc_external("Sys/PetscDrawHGAddWeightedValue"))
"""
function PetscDrawHGAddWeightedValue(petsclib::PetscLibType, hist::PetscDrawHG, value::PetscReal, weight::PetscReal) end

@for_petsc function PetscDrawHGAddWeightedValue(petsclib::$UnionPetscLib, hist::PetscDrawHG, value::$PetscReal, weight::$PetscReal )

    @chk ccall(
               (:PetscDrawHGAddWeightedValue, $petsc_library),
               PetscErrorCode,
               (PetscDrawHG, $PetscReal, $PetscReal),
               hist, value, weight,
              )


	return nothing
end 

"""
	PetscDrawHGDraw(petsclib::PetscLibType,hist::PetscDrawHG) 
Redraws a histogram.

Collective

Input Parameter:
- `hist` - The histogram context

Level: intermediate

-seealso: `PetscDrawHGCreate()`, `PetscDrawHG`, `PetscDrawHGAddValue()`, `PetscDrawHGReset()`

# External Links
$(_doc_external("Sys/PetscDrawHGDraw"))
"""
function PetscDrawHGDraw(petsclib::PetscLibType, hist::PetscDrawHG) end

@for_petsc function PetscDrawHGDraw(petsclib::$UnionPetscLib, hist::PetscDrawHG )

    @chk ccall(
               (:PetscDrawHGDraw, $petsc_library),
               PetscErrorCode,
               (PetscDrawHG,),
               hist,
              )


	return nothing
end 

"""
	PetscDrawHGSave(petsclib::PetscLibType,hg::PetscDrawHG) 
Saves a drawn image

Collective

Input Parameter:
- `hg` - The histogram context

Level: intermediate

-seealso: `PetscDrawSave()`, `PetscDrawHGCreate()`, `PetscDrawHGGetDraw()`, `PetscDrawSetSave()`, `PetscDrawHGDraw()`

# External Links
$(_doc_external("Sys/PetscDrawHGSave"))
"""
function PetscDrawHGSave(petsclib::PetscLibType, hg::PetscDrawHG) end

@for_petsc function PetscDrawHGSave(petsclib::$UnionPetscLib, hg::PetscDrawHG )

    @chk ccall(
               (:PetscDrawHGSave, $petsc_library),
               PetscErrorCode,
               (PetscDrawHG,),
               hg,
              )


	return nothing
end 

"""
	PetscDrawHGView(petsclib::PetscLibType,hist::PetscDrawHG, viewer::PetscViewer) 
Prints the histogram information to a viewer

Not Collective

Input Parameters:
- `hist`   - The histogram context
- `viewer` - The viewer to view it with

Level: beginner

-seealso: `PetscDrawHG`, `PetscViewer`, `PetscDrawHGCreate()`, `PetscDrawHGGetDraw()`, `PetscDrawSetSave()`, `PetscDrawSave()`, `PetscDrawHGDraw()`

# External Links
$(_doc_external("Sys/PetscDrawHGView"))
"""
function PetscDrawHGView(petsclib::PetscLibType, hist::PetscDrawHG, viewer::PetscViewer) end

@for_petsc function PetscDrawHGView(petsclib::$UnionPetscLib, hist::PetscDrawHG, viewer::PetscViewer )

    @chk ccall(
               (:PetscDrawHGView, $petsc_library),
               PetscErrorCode,
               (PetscDrawHG, PetscViewer),
               hist, viewer,
              )


	return nothing
end 

"""
	PetscDrawHGSetColor(petsclib::PetscLibType,hist::PetscDrawHG, color::Cint) 
Sets the color the bars will be drawn with.

Logically Collective

Input Parameters:
- `hist`  - The histogram context
- `color` - one of the colors defined in petscdraw.h or `PETSC_DRAW_ROTATE` to make each bar a different color

Level: intermediate

-seealso: `PetscDrawHG`, `PetscDrawHGCreate()`, `PetscDrawHGGetDraw()`, `PetscDrawSetSave()`, `PetscDrawSave()`, `PetscDrawHGDraw()`, `PetscDrawHGGetAxis()`

# External Links
$(_doc_external("Sys/PetscDrawHGSetColor"))
"""
function PetscDrawHGSetColor(petsclib::PetscLibType, hist::PetscDrawHG, color::Cint) end

@for_petsc function PetscDrawHGSetColor(petsclib::$UnionPetscLib, hist::PetscDrawHG, color::Cint )

    @chk ccall(
               (:PetscDrawHGSetColor, $petsc_library),
               PetscErrorCode,
               (PetscDrawHG, Cint),
               hist, color,
              )


	return nothing
end 

"""
	PetscDrawHGSetLimits(petsclib::PetscLibType,hist::PetscDrawHG, x_min::PetscReal, x_max::PetscReal, y_min::Cint, y_max::Cint) 
Sets the axis limits for a histogram. If more
points are added after this call, the limits will be adjusted to
include those additional points.

Logically Collective

Input Parameters:
- `hist`  - The histogram context
- `x_min` - the horizontal lower limit
- `x_max` - the horizontal upper limit
- `y_min` - the vertical lower limit
- `y_max` - the vertical upper limit

Level: intermediate

-seealso: `PetscDrawHG`, `PetscDrawHGCreate()`, `PetscDrawHGGetDraw()`, `PetscDrawSetSave()`, `PetscDrawSave()`, `PetscDrawHGDraw()`, `PetscDrawHGGetAxis()`

# External Links
$(_doc_external("Sys/PetscDrawHGSetLimits"))
"""
function PetscDrawHGSetLimits(petsclib::PetscLibType, hist::PetscDrawHG, x_min::PetscReal, x_max::PetscReal, y_min::Cint, y_max::Cint) end

@for_petsc function PetscDrawHGSetLimits(petsclib::$UnionPetscLib, hist::PetscDrawHG, x_min::$PetscReal, x_max::$PetscReal, y_min::Cint, y_max::Cint )

    @chk ccall(
               (:PetscDrawHGSetLimits, $petsc_library),
               PetscErrorCode,
               (PetscDrawHG, $PetscReal, $PetscReal, Cint, Cint),
               hist, x_min, x_max, y_min, y_max,
              )


	return nothing
end 

"""
	PetscDrawHGCalcStats(petsclib::PetscLibType,hist::PetscDrawHG, calc::PetscBool) 
Turns on calculation of descriptive statistics associated with the histogram

Not Collective

Input Parameters:
- `hist` - The histogram context
- `calc` - Flag for calculation

Level: intermediate

-seealso: `PetscDrawHG`, `PetscDrawHGCreate()`, `PetscDrawHGAddValue()`, `PetscDrawHGView()`, `PetscDrawHGDraw()`

# External Links
$(_doc_external("Sys/PetscDrawHGCalcStats"))
"""
function PetscDrawHGCalcStats(petsclib::PetscLibType, hist::PetscDrawHG, calc::PetscBool) end

@for_petsc function PetscDrawHGCalcStats(petsclib::$UnionPetscLib, hist::PetscDrawHG, calc::PetscBool )

    @chk ccall(
               (:PetscDrawHGCalcStats, $petsc_library),
               PetscErrorCode,
               (PetscDrawHG, PetscBool),
               hist, calc,
              )


	return nothing
end 

"""
	PetscDrawHGIntegerBins(petsclib::PetscLibType,hist::PetscDrawHG, ints::PetscBool) 
Turns on integer width bins

Not Collective

Input Parameters:
- `hist` - The histogram context
- `ints` - Flag for integer width bins

Level: intermediate

-seealso: `PetscDrawHG`, `PetscDrawHGCreate()`, `PetscDrawHGAddValue()`, `PetscDrawHGView()`, `PetscDrawHGDraw()`, `PetscDrawHGSetColor()`

# External Links
$(_doc_external("Sys/PetscDrawHGIntegerBins"))
"""
function PetscDrawHGIntegerBins(petsclib::PetscLibType, hist::PetscDrawHG, ints::PetscBool) end

@for_petsc function PetscDrawHGIntegerBins(petsclib::$UnionPetscLib, hist::PetscDrawHG, ints::PetscBool )

    @chk ccall(
               (:PetscDrawHGIntegerBins, $petsc_library),
               PetscErrorCode,
               (PetscDrawHG, PetscBool),
               hist, ints,
              )


	return nothing
end 

"""
	PetscDrawHGGetAxis(petsclib::PetscLibType,hist::PetscDrawHG, axis::PetscDrawAxis) 
Gets the axis context associated with a histogram.
This is useful if one wants to change some axis property, such as
labels, color, etc. The axis context should not be destroyed by the
application code.

Not Collective, axis is parallel if hist is parallel

Input Parameter:
- `hist` - The histogram context

Output Parameter:
- `axis` - The axis context

Level: intermediate

-seealso: `PetscDrawHG`, `PetscDrawAxis`, `PetscDrawHGCreate()`, `PetscDrawHGAddValue()`, `PetscDrawHGView()`, `PetscDrawHGDraw()`, `PetscDrawHGSetColor()`, `PetscDrawHGSetLimits()`

# External Links
$(_doc_external("Sys/PetscDrawHGGetAxis"))
"""
function PetscDrawHGGetAxis(petsclib::PetscLibType, hist::PetscDrawHG, axis::PetscDrawAxis) end

@for_petsc function PetscDrawHGGetAxis(petsclib::$UnionPetscLib, hist::PetscDrawHG, axis::PetscDrawAxis )

    @chk ccall(
               (:PetscDrawHGGetAxis, $petsc_library),
               PetscErrorCode,
               (PetscDrawHG, Ptr{PetscDrawAxis}),
               hist, axis,
              )


	return nothing
end 

"""
	PetscDrawHGGetDraw(petsclib::PetscLibType,hist::PetscDrawHG, draw::PetscDraw) 
Gets the draw context associated with a histogram.

Not Collective, draw is parallel if hist is parallel

Input Parameter:
- `hist` - The histogram context

Output Parameter:
- `draw` - The draw context

Level: intermediate

-seealso: `PetscDraw`, `PetscDrawHG`, `PetscDrawHGCreate()`, `PetscDrawHGAddValue()`, `PetscDrawHGView()`, `PetscDrawHGDraw()`, `PetscDrawHGSetColor()`, `PetscDrawAxis`, `PetscDrawHGSetLimits()`

# External Links
$(_doc_external("Sys/PetscDrawHGGetDraw"))
"""
function PetscDrawHGGetDraw(petsclib::PetscLibType, hist::PetscDrawHG, draw::PetscDraw) end

@for_petsc function PetscDrawHGGetDraw(petsclib::$UnionPetscLib, hist::PetscDrawHG, draw::PetscDraw )

    @chk ccall(
               (:PetscDrawHGGetDraw, $petsc_library),
               PetscErrorCode,
               (PetscDrawHG, Ptr{PetscDraw}),
               hist, draw,
              )


	return nothing
end 

"""
	bar::PetscDrawBar = PetscDrawBarCreate(petsclib::PetscLibType,draw::PetscDraw) 
Creates a bar graph data structure.

Collective

Input Parameter:
- `draw` - The window where the graph will be made

Output Parameter:
- `bar` - The bar graph context

-seealso: `PetscDrawBar`, `PetscDrawLGCreate()`, `PetscDrawLG`, `PetscDrawSPCreate()`, `PetscDrawSP`, `PetscDrawHGCreate()`, `PetscDrawHG`, `PetscDrawBarDestroy()`, `PetscDrawBarSetData()`,
`PetscDrawBarDraw()`, `PetscDrawBarSave()`, `PetscDrawBarSetColor()`, `PetscDrawBarSort()`, `PetscDrawBarSetLimits()`, `PetscDrawBarGetAxis()`, `PetscDrawAxis`,
`PetscDrawBarGetDraw()`, `PetscDrawBarSetFromOptions()`

# External Links
$(_doc_external("Sys/PetscDrawBarCreate"))
"""
function PetscDrawBarCreate(petsclib::PetscLibType, draw::PetscDraw) end

@for_petsc function PetscDrawBarCreate(petsclib::$UnionPetscLib, draw::PetscDraw )
	bar_ = Ref{PetscDrawBar}()

    @chk ccall(
               (:PetscDrawBarCreate, $petsc_library),
               PetscErrorCode,
               (PetscDraw, Ptr{PetscDrawBar}),
               draw, bar_,
              )

	bar = bar_[]

	return bar
end 

"""
	PetscDrawBarSetData(petsclib::PetscLibType,bar::PetscDrawBar, bins::PetscInt, data::Vector{PetscReal}, labels::String) 
Set the data for a bar graph

Logically Collective

Input Parameters:
- `bar`    - The bar graph context.
- `bins`   - number of items
- `data`   - values of each item
- `labels` - optional label for each bar, `NULL` terminated array of strings

Level: intermediate

-seealso: `PetscDrawBar`, `PetscDrawBarCreate()`, `PetscDrawBarDraw()`

# External Links
$(_doc_external("Sys/PetscDrawBarSetData"))
"""
function PetscDrawBarSetData(petsclib::PetscLibType, bar::PetscDrawBar, bins::PetscInt, data::Vector{PetscReal}, labels::String) end

@for_petsc function PetscDrawBarSetData(petsclib::$UnionPetscLib, bar::PetscDrawBar, bins::$PetscInt, data::Vector{$PetscReal}, labels::String )
	labels_ = Ref(pointer(labels))

    @chk ccall(
               (:PetscDrawBarSetData, $petsc_library),
               PetscErrorCode,
               (PetscDrawBar, $PetscInt, Ptr{$PetscReal}, Ptr{Ptr{Cchar}}),
               bar, bins, data, labels_,
              )


	return nothing
end 

"""
	PetscDrawBarDestroy(petsclib::PetscLibType,bar::PetscDrawBar) 
Frees all space taken up by bar graph data structure.

Collective

Input Parameter:
- `bar` - The bar graph context

Level: intermediate

-seealso: `PetscDrawBar`, `PetscDrawBarCreate()`

# External Links
$(_doc_external("Sys/PetscDrawBarDestroy"))
"""
function PetscDrawBarDestroy(petsclib::PetscLibType, bar::PetscDrawBar) end

@for_petsc function PetscDrawBarDestroy(petsclib::$UnionPetscLib, bar::PetscDrawBar )

    @chk ccall(
               (:PetscDrawBarDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscDrawBar},),
               bar,
              )


	return nothing
end 

"""
	PetscDrawBarDraw(petsclib::PetscLibType,bar::PetscDrawBar) 
Redraws a bar graph.

Collective

Input Parameter:
- `bar` - The bar graph context

Level: intermediate

-seealso: `PetscDrawBar`, `PetscDrawBarCreate()`, `PetscDrawBarSetData()`

# External Links
$(_doc_external("Sys/PetscDrawBarDraw"))
"""
function PetscDrawBarDraw(petsclib::PetscLibType, bar::PetscDrawBar) end

@for_petsc function PetscDrawBarDraw(petsclib::$UnionPetscLib, bar::PetscDrawBar )

    @chk ccall(
               (:PetscDrawBarDraw, $petsc_library),
               PetscErrorCode,
               (PetscDrawBar,),
               bar,
              )


	return nothing
end 

"""
	PetscDrawBarSave(petsclib::PetscLibType,bar::PetscDrawBar) 
Saves a drawn bar graph

Collective

Input Parameter:
- `bar` - The bar graph context

Level: intermediate

-seealso: `PetscDrawSave()`, `PetscDrawBar`, `PetscDrawBarCreate()`, `PetscDrawBarGetDraw()`, `PetscDrawSetSave()`, `PetscDrawBarSetData()`

# External Links
$(_doc_external("Sys/PetscDrawBarSave"))
"""
function PetscDrawBarSave(petsclib::PetscLibType, bar::PetscDrawBar) end

@for_petsc function PetscDrawBarSave(petsclib::$UnionPetscLib, bar::PetscDrawBar )

    @chk ccall(
               (:PetscDrawBarSave, $petsc_library),
               PetscErrorCode,
               (PetscDrawBar,),
               bar,
              )


	return nothing
end 

"""
	PetscDrawBarSetColor(petsclib::PetscLibType,bar::PetscDrawBar, color::Cint) 
Sets the color the bars will be drawn with.

Logically Collective

Input Parameters:
- `bar`   - The bar graph context
- `color` - one of the colors defined in petscdraw.h or `PETSC_DRAW_ROTATE` to make each bar a
different color

Level: intermediate

-seealso: `PetscDrawBarCreate()`, `PetscDrawBar`, `PetscDrawBarSetData()`, `PetscDrawBarDraw()`, `PetscDrawBarGetAxis()`

# External Links
$(_doc_external("Sys/PetscDrawBarSetColor"))
"""
function PetscDrawBarSetColor(petsclib::PetscLibType, bar::PetscDrawBar, color::Cint) end

@for_petsc function PetscDrawBarSetColor(petsclib::$UnionPetscLib, bar::PetscDrawBar, color::Cint )

    @chk ccall(
               (:PetscDrawBarSetColor, $petsc_library),
               PetscErrorCode,
               (PetscDrawBar, Cint),
               bar, color,
              )


	return nothing
end 

"""
	PetscDrawBarSort(petsclib::PetscLibType,bar::PetscDrawBar, sort::PetscBool, tolerance::PetscReal) 
Sorts the values before drawing the bar chart, the bars will be in ascending order from left to right

Logically Collective

Input Parameters:
- `bar`       - The bar graph context
- `sort`      - `PETSC_TRUE` to sort the values
- `tolerance` - discard values less than tolerance

Level: intermediate

-seealso: `PetscDrawBar`, `PetscDrawBarCreate()`, `PetscDrawBarSetData()`, `PetscDrawBarSetColor()`, `PetscDrawBarDraw()`, `PetscDrawBarGetAxis()`

# External Links
$(_doc_external("Sys/PetscDrawBarSort"))
"""
function PetscDrawBarSort(petsclib::PetscLibType, bar::PetscDrawBar, sort::PetscBool, tolerance::PetscReal) end

@for_petsc function PetscDrawBarSort(petsclib::$UnionPetscLib, bar::PetscDrawBar, sort::PetscBool, tolerance::$PetscReal )

    @chk ccall(
               (:PetscDrawBarSort, $petsc_library),
               PetscErrorCode,
               (PetscDrawBar, PetscBool, $PetscReal),
               bar, sort, tolerance,
              )


	return nothing
end 

"""
	PetscDrawBarSetLimits(petsclib::PetscLibType,bar::PetscDrawBar, y_min::PetscReal, y_max::PetscReal) 
Sets the axis limits for a bar graph. If more
points are added after this call, the limits will be adjusted to
include those additional points.

Logically Collective

Input Parameters:
- `bar`   - The bar graph context
- `y_min` - The lower limit
- `y_max` - The upper limit

Level: intermediate

-seealso: `PetscDrawBar`, `PetscDrawBarCreate()`, `PetscDrawBarGetAxis()`, `PetscDrawBarSetData()`, `PetscDrawBarDraw()`

# External Links
$(_doc_external("Sys/PetscDrawBarSetLimits"))
"""
function PetscDrawBarSetLimits(petsclib::PetscLibType, bar::PetscDrawBar, y_min::PetscReal, y_max::PetscReal) end

@for_petsc function PetscDrawBarSetLimits(petsclib::$UnionPetscLib, bar::PetscDrawBar, y_min::$PetscReal, y_max::$PetscReal )

    @chk ccall(
               (:PetscDrawBarSetLimits, $petsc_library),
               PetscErrorCode,
               (PetscDrawBar, $PetscReal, $PetscReal),
               bar, y_min, y_max,
              )


	return nothing
end 

"""
	PetscDrawBarGetAxis(petsclib::PetscLibType,bar::PetscDrawBar, axis::PetscDrawAxis) 
Gets the axis context associated with a bar graph.
This is useful if one wants to change some axis property, such as
labels, color, etc. The axis context should not be destroyed by the
application code.

Not Collective, axis is parallel if bar is parallel

Input Parameter:
- `bar` - The bar graph context

Output Parameter:
- `axis` - The axis context

Level: intermediate

-seealso: `PetscDrawBar`, `PetscDrawBarCreate()`, `PetscDrawAxis`, `PetscDrawAxisCreate()`

# External Links
$(_doc_external("Sys/PetscDrawBarGetAxis"))
"""
function PetscDrawBarGetAxis(petsclib::PetscLibType, bar::PetscDrawBar, axis::PetscDrawAxis) end

@for_petsc function PetscDrawBarGetAxis(petsclib::$UnionPetscLib, bar::PetscDrawBar, axis::PetscDrawAxis )

    @chk ccall(
               (:PetscDrawBarGetAxis, $petsc_library),
               PetscErrorCode,
               (PetscDrawBar, Ptr{PetscDrawAxis}),
               bar, axis,
              )


	return nothing
end 

"""
	PetscDrawBarGetDraw(petsclib::PetscLibType,bar::PetscDrawBar, draw::PetscDraw) 
Gets the draw context associated with a bar graph.

Not Collective, draw is parallel if bar is parallel

Input Parameter:
- `bar` - The bar graph context

Output Parameter:
- `draw` - The draw context

Level: intermediate

-seealso: `PetscDrawBar`, `PetscDraw`, `PetscDrawBarCreate()`, `PetscDrawBarDraw()`

# External Links
$(_doc_external("Sys/PetscDrawBarGetDraw"))
"""
function PetscDrawBarGetDraw(petsclib::PetscLibType, bar::PetscDrawBar, draw::PetscDraw) end

@for_petsc function PetscDrawBarGetDraw(petsclib::$UnionPetscLib, bar::PetscDrawBar, draw::PetscDraw )

    @chk ccall(
               (:PetscDrawBarGetDraw, $petsc_library),
               PetscErrorCode,
               (PetscDrawBar, Ptr{PetscDraw}),
               bar, draw,
              )


	return nothing
end 

"""
	PetscDrawBarSetFromOptions(petsclib::PetscLibType,bar::PetscDrawBar) 
Sets options related to the display of the `PetscDrawBar`

Collective

Input Parameter:
- `bar` - the bar graph context

Options Database Key:
- `-bar_sort` - sort the entries before drawing the bar graph

Level: intermediate

-seealso: `PetscDrawBar`, `PetscDrawBarDestroy()`, `PetscDrawBarCreate()`, `PetscDrawBarSort()`

# External Links
$(_doc_external("Sys/PetscDrawBarSetFromOptions"))
"""
function PetscDrawBarSetFromOptions(petsclib::PetscLibType, bar::PetscDrawBar) end

@for_petsc function PetscDrawBarSetFromOptions(petsclib::$UnionPetscLib, bar::PetscDrawBar )

    @chk ccall(
               (:PetscDrawBarSetFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscDrawBar,),
               bar,
              )


	return nothing
end 

