"""
	PetscDrawAppendTitle(petsclib::PetscLibType, draw::PetscDraw, title::String) 
Appends to the title of a `PetscDraw` context.

Collective

Input Parameters:
- `draw`  - the graphics context
- `title` - the title

Level: advanced

See also: `PetscDraw`, `PetscDrawSetTitle()`, `PetscDrawGetTitle()`

# External Links
$(_doc_external("Draw/PetscDrawAppendTitle"))
"""
function PetscDrawAppendTitle(petsclib::PetscLibType, draw::PetscDraw, title::String)
    error("PetscDrawAppendTitle: no generated method for these argument types")
end

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
	PetscDrawArrow(petsclib::PetscLibType, draw::PetscDraw, xl::PetscReal, yl::PetscReal, xr::PetscReal, yr::PetscReal, cl::Cint) 
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

See also: `PetscDraw`, `PetscDrawLine()`, `PetscDrawLineSetWidth()`, `PetscDrawLineGetWidth()`, `PetscDrawRectangle()`, `PetscDrawTriangle()`, `PetscDrawEllipse()`,
`PetscDrawMarker()`, `PetscDrawPoint()`

# External Links
$(_doc_external("Draw/PetscDrawArrow"))
"""
function PetscDrawArrow(petsclib::PetscLibType, draw::PetscDraw, xl::Real, yl::Real, xr::Real, yr::Real, cl::Cint)
    error("PetscDrawArrow: no generated method for these argument types")
end

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
	axis::PetscDrawAxis = PetscDrawAxisCreate(petsclib::PetscLibType, draw::PetscDraw) 
Generate the axis data structure.

Collective

Input Parameter:
- `draw` - `PetscDraw` object where axis to be made

Output Parameter:
- `axis` - the axis datastructure

See also: `PetscDrawLGCreate()`, `PetscDrawLG`, `PetscDrawSPCreate()`, `PetscDrawSP`, `PetscDrawHGCreate()`, `PetscDrawHG`, `PetscDrawBarCreate()`, `PetscDrawBar`, `PetscDrawLGGetAxis()`, `PetscDrawSPGetAxis()`,
`PetscDrawHGGetAxis()`, `PetscDrawBarGetAxis()`, `PetscDrawAxis`, `PetscDrawAxisDestroy()`, `PetscDrawAxisSetColors()`, `PetscDrawAxisSetLabels()`, `PetscDrawAxisSetLimits()`, `PetscDrawAxisGetLimits()`, `PetscDrawAxisSetHoldLimits()`,
`PetscDrawAxisDraw()`

# External Links
$(_doc_external("Draw/PetscDrawAxisCreate"))
"""
function PetscDrawAxisCreate(petsclib::PetscLibType, draw::PetscDraw)
    error("PetscDrawAxisCreate: no generated method for these argument types")
end

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
	PetscDrawAxisDestroy(petsclib::PetscLibType, axis::Union{PetscDrawAxis, Ref{PetscDrawAxis}}) 
Frees the space used by an axis structure.

Collective

Input Parameter:
- `axis` - the axis context

Level: advanced

See also: `PetscDraw`, `PetscDrawAxisCreate()`, `PetscDrawAxis`

# External Links
$(_doc_external("Draw/PetscDrawAxisDestroy"))
"""
function PetscDrawAxisDestroy(petsclib::PetscLibType, axis::Union{PetscDrawAxis, Ref{PetscDrawAxis}})
    error("PetscDrawAxisDestroy: no generated method for these argument types")
end

@for_petsc function PetscDrawAxisDestroy(petsclib::$UnionPetscLib, axis::Union{PetscDrawAxis, Ref{PetscDrawAxis}} )
	axis_ = axis isa Base.RefValue ? axis : Ref{PetscDrawAxis}(axis)

    @chk ccall(
               (:PetscDrawAxisDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscDrawAxis},),
               axis_,
              )


	return nothing
end 

"""
	PetscDrawAxisDraw(petsclib::PetscLibType, axis::PetscDrawAxis) 
draws an axis.

Collective

Input Parameter:
- `axis` - `PetscDrawAxis` structure

Level: advanced

See also: `PetscDrawAxisCreate()`, `PetscDrawAxis`, `PetscDrawAxisGetLimits()`, `PetscDrawAxisSetLimits()`, `PetscDrawAxisSetLabels()`, `PetscDrawAxisSetColors()`

# External Links
$(_doc_external("Draw/PetscDrawAxisDraw"))
"""
function PetscDrawAxisDraw(petsclib::PetscLibType, axis::PetscDrawAxis)
    error("PetscDrawAxisDraw: no generated method for these argument types")
end

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
	xmin::PetscReal,xmax::PetscReal,ymin::PetscReal,ymax::PetscReal = PetscDrawAxisGetLimits(petsclib::PetscLibType, axis::PetscDrawAxis) 
Gets the limits (in user coords) of the axis

Not Collective

Input Parameters:
- `axis` - the axis
- `xmin` - the lower x limit
- `xmax` - the upper x limit
- `ymin` - the lower y limit
- `ymax` - the upper y limit

Level: advanced

See also: `PetscDrawAxisCreate()`, `PetscDrawAxis`, `PetscDrawAxisSetHoldLimits()`, `PetscDrawAxisSetLimits()`, `PetscDrawAxisSetLabels()`, `PetscDrawAxisSetColors()`

# External Links
$(_doc_external("Draw/PetscDrawAxisGetLimits"))
"""
function PetscDrawAxisGetLimits(petsclib::PetscLibType, axis::PetscDrawAxis)
    error("PetscDrawAxisGetLimits: no generated method for these argument types")
end

@for_petsc function PetscDrawAxisGetLimits(petsclib::$UnionPetscLib, axis::PetscDrawAxis )
	xmin_ = Ref{$PetscReal}()
	xmax_ = Ref{$PetscReal}()
	ymin_ = Ref{$PetscReal}()
	ymax_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscDrawAxisGetLimits, $petsc_library),
               PetscErrorCode,
               (PetscDrawAxis, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
               axis, xmin_, xmax_, ymin_, ymax_,
              )

	xmin = xmin_[]
	xmax = xmax_[]
	ymin = ymin_[]
	ymax = ymax_[]

	return xmin,xmax,ymin,ymax
end 

"""
	PetscDrawAxisSetColors(petsclib::PetscLibType, axis::PetscDrawAxis, ac::Cint, tc::Cint, cc::Cint) 
Sets the colors to be used for the axis,
tickmarks, and text.

Logically Collective

Input Parameters:
- `axis` - the axis
- `ac`   - the color of the axis lines
- `tc`   - the color of the tick marks
- `cc`   - the color of the text strings

Level: advanced

See also: `PetscDraw`, `PetscDrawAxisCreate()`, `PetscDrawAxis`, `PetscDrawAxisSetLabels()`, `PetscDrawAxisDraw()`, `PetscDrawAxisSetLimits()`

# External Links
$(_doc_external("Draw/PetscDrawAxisSetColors"))
"""
function PetscDrawAxisSetColors(petsclib::PetscLibType, axis::PetscDrawAxis, ac::Cint, tc::Cint, cc::Cint)
    error("PetscDrawAxisSetColors: no generated method for these argument types")
end

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
	PetscDrawAxisSetHoldLimits(petsclib::PetscLibType, axis::PetscDrawAxis, hold::PetscBool) 
Causes an axis to keep the same limits until this is called
again

Logically Collective

Input Parameters:
- `axis` - the axis
- `hold` - `PETSC_TRUE` - hold current limits, `PETSC_FALSE` allow limits to be changed

Level: advanced

See also: `PetscDrawAxisCreate()`, `PetscDrawAxis`, `PetscDrawAxisGetLimits()`, `PetscDrawAxisSetLimits()`, `PetscDrawAxisSetLabels()`, `PetscDrawAxisSetColors()`

# External Links
$(_doc_external("Draw/PetscDrawAxisSetHoldLimits"))
"""
function PetscDrawAxisSetHoldLimits(petsclib::PetscLibType, axis::PetscDrawAxis, hold::PetscBool)
    error("PetscDrawAxisSetHoldLimits: no generated method for these argument types")
end

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
	PetscDrawAxisSetLabels(petsclib::PetscLibType, axis::PetscDrawAxis, top::String, xlabel::String, ylabel::String) 
Sets the x and y axis labels.

Logically Collective

Input Parameters:
- `axis`   - the axis
- `top`    - the label at the top of the image
- `xlabel` - the x axis label
- `ylabel` - the y axis label

Level: advanced

See also: `PetscDraw`, `PetscDrawAxisCreate()`, `PetscDrawAxis`, `PetscDrawAxisSetColors()`, `PetscDrawAxisDraw()`, `PetscDrawAxisSetLimits()`

# External Links
$(_doc_external("Draw/PetscDrawAxisSetLabels"))
"""
function PetscDrawAxisSetLabels(petsclib::PetscLibType, axis::PetscDrawAxis, top::String, xlabel::String, ylabel::String)
    error("PetscDrawAxisSetLabels: no generated method for these argument types")
end

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
	PetscDrawAxisSetLimits(petsclib::PetscLibType, axis::PetscDrawAxis, xmin::PetscReal, xmax::PetscReal, ymin::PetscReal, ymax::PetscReal) 
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

See also: `PetscDrawAxisSetHoldLimits()`, `PetscDrawAxisGetLimits()`, `PetscDrawAxisSetLabels()`, `PetscDrawAxisSetColors()`

# External Links
$(_doc_external("Draw/PetscDrawAxisSetLimits"))
"""
function PetscDrawAxisSetLimits(petsclib::PetscLibType, axis::PetscDrawAxis, xmin::Real, xmax::Real, ymin::Real, ymax::Real)
    error("PetscDrawAxisSetLimits: no generated method for these argument types")
end

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
	PetscDrawBOP(petsclib::PetscLibType, draw::PetscDraw) 
Begins a new page or frame on the selected graphical device.

Logically Collective

Input Parameter:
- `draw` - the drawing context

Level: advanced

See also: `PetscDrawEOP()`, `PetscDrawClear()`

# External Links
$(_doc_external("Draw/PetscDrawBOP"))
"""
function PetscDrawBOP(petsclib::PetscLibType, draw::PetscDraw)
    error("PetscDrawBOP: no generated method for these argument types")
end

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
	bar::PetscDrawBar = PetscDrawBarCreate(petsclib::PetscLibType, draw::PetscDraw) 
Creates a bar graph data structure.

Collective

Input Parameter:
- `draw` - The window where the graph will be made

Output Parameter:
- `bar` - The bar graph context

See also: `PetscDrawBar`, `PetscDrawLGCreate()`, `PetscDrawLG`, `PetscDrawSPCreate()`, `PetscDrawSP`, `PetscDrawHGCreate()`, `PetscDrawHG`, `PetscDrawBarDestroy()`, `PetscDrawBarSetData()`,
`PetscDrawBarDraw()`, `PetscDrawBarSave()`, `PetscDrawBarSetColor()`, `PetscDrawBarSort()`, `PetscDrawBarSetLimits()`, `PetscDrawBarGetAxis()`, `PetscDrawAxis`,
`PetscDrawBarGetDraw()`, `PetscDrawBarSetFromOptions()`

# External Links
$(_doc_external("Draw/PetscDrawBarCreate"))
"""
function PetscDrawBarCreate(petsclib::PetscLibType, draw::PetscDraw)
    error("PetscDrawBarCreate: no generated method for these argument types")
end

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
	PetscDrawBarDestroy(petsclib::PetscLibType, bar::Union{PetscDrawBar, Ref{PetscDrawBar}}) 
Frees all space taken up by bar graph data structure.

Collective

Input Parameter:
- `bar` - The bar graph context

Level: intermediate

See also: `PetscDrawBar`, `PetscDrawBarCreate()`

# External Links
$(_doc_external("Draw/PetscDrawBarDestroy"))
"""
function PetscDrawBarDestroy(petsclib::PetscLibType, bar::Union{PetscDrawBar, Ref{PetscDrawBar}})
    error("PetscDrawBarDestroy: no generated method for these argument types")
end

@for_petsc function PetscDrawBarDestroy(petsclib::$UnionPetscLib, bar::Union{PetscDrawBar, Ref{PetscDrawBar}} )
	bar_ = bar isa Base.RefValue ? bar : Ref{PetscDrawBar}(bar)

    @chk ccall(
               (:PetscDrawBarDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscDrawBar},),
               bar_,
              )


	return nothing
end 

"""
	PetscDrawBarDraw(petsclib::PetscLibType, bar::PetscDrawBar) 
Redraws a bar graph.

Collective

Input Parameter:
- `bar` - The bar graph context

Level: intermediate

See also: `PetscDrawBar`, `PetscDrawBarCreate()`, `PetscDrawBarSetData()`

# External Links
$(_doc_external("Draw/PetscDrawBarDraw"))
"""
function PetscDrawBarDraw(petsclib::PetscLibType, bar::PetscDrawBar)
    error("PetscDrawBarDraw: no generated method for these argument types")
end

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
	axis::PetscDrawAxis = PetscDrawBarGetAxis(petsclib::PetscLibType, bar::PetscDrawBar) 
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

See also: `PetscDrawBar`, `PetscDrawBarCreate()`, `PetscDrawAxis`, `PetscDrawAxisCreate()`

# External Links
$(_doc_external("Draw/PetscDrawBarGetAxis"))
"""
function PetscDrawBarGetAxis(petsclib::PetscLibType, bar::PetscDrawBar)
    error("PetscDrawBarGetAxis: no generated method for these argument types")
end

@for_petsc function PetscDrawBarGetAxis(petsclib::$UnionPetscLib, bar::PetscDrawBar )
	axis_ = Ref{PetscDrawAxis}()

    @chk ccall(
               (:PetscDrawBarGetAxis, $petsc_library),
               PetscErrorCode,
               (PetscDrawBar, Ptr{PetscDrawAxis}),
               bar, axis_,
              )

	axis = axis_[]

	return axis
end 

"""
	draw::PetscDraw = PetscDrawBarGetDraw(petsclib::PetscLibType, bar::PetscDrawBar) 
Gets the draw context associated with a bar graph.

Not Collective, draw is parallel if bar is parallel

Input Parameter:
- `bar` - The bar graph context

Output Parameter:
- `draw` - The draw context

Level: intermediate

See also: `PetscDrawBar`, `PetscDraw`, `PetscDrawBarCreate()`, `PetscDrawBarDraw()`

# External Links
$(_doc_external("Draw/PetscDrawBarGetDraw"))
"""
function PetscDrawBarGetDraw(petsclib::PetscLibType, bar::PetscDrawBar)
    error("PetscDrawBarGetDraw: no generated method for these argument types")
end

@for_petsc function PetscDrawBarGetDraw(petsclib::$UnionPetscLib, bar::PetscDrawBar )
	draw_ = Ref{PetscDraw}()

    @chk ccall(
               (:PetscDrawBarGetDraw, $petsc_library),
               PetscErrorCode,
               (PetscDrawBar, Ptr{PetscDraw}),
               bar, draw_,
              )

	draw = draw_[]

	return draw
end 

"""
	PetscDrawBarSave(petsclib::PetscLibType, bar::PetscDrawBar) 
Saves a drawn bar graph

Collective

Input Parameter:
- `bar` - The bar graph context

Level: intermediate

See also: `PetscDrawSave()`, `PetscDrawBar`, `PetscDrawBarCreate()`, `PetscDrawBarGetDraw()`, `PetscDrawSetSave()`, `PetscDrawBarSetData()`

# External Links
$(_doc_external("Draw/PetscDrawBarSave"))
"""
function PetscDrawBarSave(petsclib::PetscLibType, bar::PetscDrawBar)
    error("PetscDrawBarSave: no generated method for these argument types")
end

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
	PetscDrawBarSetColor(petsclib::PetscLibType, bar::PetscDrawBar, color::Cint) 
Sets the color the bars will be drawn with.

Logically Collective

Input Parameters:
- `bar`   - The bar graph context
- `color` - one of the colors defined in petscdraw.h or `PETSC_DRAW_ROTATE` to make each bar a
different color

Level: intermediate

See also: `PetscDrawBarCreate()`, `PetscDrawBar`, `PetscDrawBarSetData()`, `PetscDrawBarDraw()`, `PetscDrawBarGetAxis()`

# External Links
$(_doc_external("Draw/PetscDrawBarSetColor"))
"""
function PetscDrawBarSetColor(petsclib::PetscLibType, bar::PetscDrawBar, color::Cint)
    error("PetscDrawBarSetColor: no generated method for these argument types")
end

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
	PetscDrawBarSetData(petsclib::PetscLibType, bar::PetscDrawBar, bins::PetscInt, data::Vector{PetscReal}, labels::String) 
Set the data for a bar graph

Logically Collective

Input Parameters:
- `bar`    - The bar graph context.
- `bins`   - number of items
- `data`   - values of each item
- `labels` - optional label for each bar, `NULL` terminated array of strings

Level: intermediate

See also: `PetscDrawBar`, `PetscDrawBarCreate()`, `PetscDrawBarDraw()`

# External Links
$(_doc_external("Draw/PetscDrawBarSetData"))
"""
function PetscDrawBarSetData(petsclib::PetscLibType, bar::PetscDrawBar, bins::Integer, data::AbstractVector{<:Number}, labels::String)
    error("PetscDrawBarSetData: no generated method for these argument types")
end

@for_petsc function PetscDrawBarSetData(petsclib::$UnionPetscLib, bar::PetscDrawBar, bins::$PetscInt, data::Vector{$PetscReal}, labels::String )
	labels_ = Ref{Ptr{Cchar}}(labels isa Ptr ? labels : pointer(labels))

    @chk ccall(
               (:PetscDrawBarSetData, $petsc_library),
               PetscErrorCode,
               (PetscDrawBar, $PetscInt, Ptr{$PetscReal}, Ptr{Ptr{Cchar}}),
               bar, bins, data, labels_,
              )


	return nothing
end 

"""
	PetscDrawBarSetFromOptions(petsclib::PetscLibType, bar::PetscDrawBar) 
Sets options related to the display of the `PetscDrawBar`

Collective

Input Parameter:
- `bar` - the bar graph context

Options Database Key:
- `-bar_sort` - sort the entries before drawing the bar graph

Level: intermediate

See also: `PetscDrawBar`, `PetscDrawBarDestroy()`, `PetscDrawBarCreate()`, `PetscDrawBarSort()`

# External Links
$(_doc_external("Draw/PetscDrawBarSetFromOptions"))
"""
function PetscDrawBarSetFromOptions(petsclib::PetscLibType, bar::PetscDrawBar)
    error("PetscDrawBarSetFromOptions: no generated method for these argument types")
end

@for_petsc function PetscDrawBarSetFromOptions(petsclib::$UnionPetscLib, bar::PetscDrawBar )

    @chk ccall(
               (:PetscDrawBarSetFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscDrawBar,),
               bar,
              )


	return nothing
end 

"""
	PetscDrawBarSetLimits(petsclib::PetscLibType, bar::PetscDrawBar, y_min::PetscReal, y_max::PetscReal) 
Sets the axis limits for a bar graph. If more
points are added after this call, the limits will be adjusted to
include those additional points.

Logically Collective

Input Parameters:
- `bar`   - The bar graph context
- `y_min` - The lower limit
- `y_max` - The upper limit

Level: intermediate

See also: `PetscDrawBar`, `PetscDrawBarCreate()`, `PetscDrawBarGetAxis()`, `PetscDrawBarSetData()`, `PetscDrawBarDraw()`

# External Links
$(_doc_external("Draw/PetscDrawBarSetLimits"))
"""
function PetscDrawBarSetLimits(petsclib::PetscLibType, bar::PetscDrawBar, y_min::Real, y_max::Real)
    error("PetscDrawBarSetLimits: no generated method for these argument types")
end

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
	PetscDrawBarSort(petsclib::PetscLibType, bar::PetscDrawBar, sort::PetscBool, tolerance::PetscReal) 
Sorts the values before drawing the bar chart, the bars will be in ascending order from left to right

Logically Collective

Input Parameters:
- `bar`       - The bar graph context
- `sort`      - `PETSC_TRUE` to sort the values
- `tolerance` - discard values less than tolerance

Level: intermediate

See also: `PetscDrawBar`, `PetscDrawBarCreate()`, `PetscDrawBarSetData()`, `PetscDrawBarSetColor()`, `PetscDrawBarDraw()`, `PetscDrawBarGetAxis()`

# External Links
$(_doc_external("Draw/PetscDrawBarSort"))
"""
function PetscDrawBarSort(petsclib::PetscLibType, bar::PetscDrawBar, sort::PetscBool, tolerance::Real)
    error("PetscDrawBarSort: no generated method for these argument types")
end

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
	PetscDrawCheckResizedWindow(petsclib::PetscLibType, draw::PetscDraw) 
Checks if the user has resized the window.

Collective

Input Parameter:
- `draw` - the window

Level: advanced

See also: `PetscDraw`, `PetscDrawResizeWindow()`

# External Links
$(_doc_external("Draw/PetscDrawCheckResizedWindow"))
"""
function PetscDrawCheckResizedWindow(petsclib::PetscLibType, draw::PetscDraw)
    error("PetscDrawCheckResizedWindow: no generated method for these argument types")
end

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
	PetscDrawClear(petsclib::PetscLibType, draw::PetscDraw) 
Clears graphical output. All processors must call this routine.
Does not return until the draw in context is clear.

Collective

Input Parameter:
- `draw` - the drawing context

Level: intermediate

See also: `PetscDrawBOP()`, `PetscDrawEOP()`

# External Links
$(_doc_external("Draw/PetscDrawClear"))
"""
function PetscDrawClear(petsclib::PetscLibType, draw::PetscDraw)
    error("PetscDrawClear: no generated method for these argument types")
end

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
	i::Cint,j::Cint = PetscDrawCoordinateToPixel(petsclib::PetscLibType, draw::PetscDraw, x::PetscReal, y::PetscReal) 
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

See also: `PetscDraw`

# External Links
$(_doc_external("Draw/PetscDrawCoordinateToPixel"))
"""
function PetscDrawCoordinateToPixel(petsclib::PetscLibType, draw::PetscDraw, x::Real, y::Real)
    error("PetscDrawCoordinateToPixel: no generated method for these argument types")
end

@for_petsc function PetscDrawCoordinateToPixel(petsclib::$UnionPetscLib, draw::PetscDraw, x::$PetscReal, y::$PetscReal )
	i_ = Ref{Cint}()
	j_ = Ref{Cint}()

    @chk ccall(
               (:PetscDrawCoordinateToPixel, $petsc_library),
               PetscErrorCode,
               (PetscDraw, $PetscReal, $PetscReal, Ptr{Cint}, Ptr{Cint}),
               draw, x, y, i_, j_,
              )

	i = i_[]
	j = j_[]

	return i,j
end 

"""
	indraw::PetscDraw = PetscDrawCreate(petsclib::PetscLibType, comm::MPI_Comm, display::String, title::String, x::Cint, y::Cint, w::Cint, h::Cint) 
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

See also: `PetscDrawSetType()`, `PetscDrawSetFromOptions()`, `PetscDrawDestroy()`, `PetscDrawLGCreate()`, `PetscDrawSPCreate()`,
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
$(_doc_external("Draw/PetscDrawCreate"))
"""
function PetscDrawCreate(petsclib::PetscLibType, comm::MPI_Comm, display::String, title::String, x::Cint, y::Cint, w::Cint, h::Cint)
    error("PetscDrawCreate: no generated method for these argument types")
end

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
	PetscDrawDestroy(petsclib::PetscLibType, draw::Union{PetscDraw, Ref{PetscDraw}}) 
Deletes a draw context.

Collective

Input Parameter:
- `draw` - the drawing context

Level: beginner

See also: `PetscDraw`, `PetscDrawCreate()`

# External Links
$(_doc_external("Draw/PetscDrawDestroy"))
"""
function PetscDrawDestroy(petsclib::PetscLibType, draw::Union{PetscDraw, Ref{PetscDraw}})
    error("PetscDrawDestroy: no generated method for these argument types")
end

@for_petsc function PetscDrawDestroy(petsclib::$UnionPetscLib, draw::Union{PetscDraw, Ref{PetscDraw}} )
	draw_ = draw isa Base.RefValue ? draw : Ref{PetscDraw}(draw)

    @chk ccall(
               (:PetscDrawDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscDraw},),
               draw_,
              )


	return nothing
end 

"""
	PetscDrawEOP(petsclib::PetscLibType, draw::PetscDraw) 
Ends a page or frame on the selected graphical device.

Logically Collective

Input Parameter:
- `draw` - the drawing context

Level: advanced

See also: `PetscDrawBOP()`, `PetscDrawClear()`

# External Links
$(_doc_external("Draw/PetscDrawEOP"))
"""
function PetscDrawEOP(petsclib::PetscLibType, draw::PetscDraw)
    error("PetscDrawEOP: no generated method for these argument types")
end

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
	PetscDrawEllipse(petsclib::PetscLibType, draw::PetscDraw, x::PetscReal, y::PetscReal, a::PetscReal, b::PetscReal, c::Cint) 
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

See also: `PetscDraw`, `PetscDrawRectangle()`, `PetscDrawTriangle()`, `PetscDrawMarker()`, `PetscDrawPoint()`, `PetscDrawString()`, `PetscDrawArrow()`

# External Links
$(_doc_external("Draw/PetscDrawEllipse"))
"""
function PetscDrawEllipse(petsclib::PetscLibType, draw::PetscDraw, x::Real, y::Real, a::Real, b::Real, c::Cint)
    error("PetscDrawEllipse: no generated method for these argument types")
end

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
	PetscDrawFinalizePackage(petsclib::PetscLibType) 
This function destroys everything in the PETSc interface to the `PetscDraw` package. It is
called from `PetscFinalize()`.

Level: developer

See also: `PetscDraw`, `PetscFinalize()`

# External Links
$(_doc_external("Draw/PetscDrawFinalizePackage"))
"""
function PetscDrawFinalizePackage(petsclib::PetscLibType)
    error("PetscDrawFinalizePackage: no generated method for these argument types")
end

@for_petsc function PetscDrawFinalizePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscDrawFinalizePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	PetscDrawFlush(petsclib::PetscLibType, draw::PetscDraw) 
Flushes graphical output.

Collective

Input Parameter:
- `draw` - the drawing context

Level: beginner

See also: `PetscDraw`, `PetscDrawClear()`

# External Links
$(_doc_external("Draw/PetscDrawFlush"))
"""
function PetscDrawFlush(petsclib::PetscLibType, draw::PetscDraw)
    error("PetscDrawFlush: no generated method for these argument types")
end

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
	xl::PetscReal,yl::PetscReal,xr::PetscReal,yr::PetscReal = PetscDrawGetBoundingBox(petsclib::PetscLibType, draw::PetscDraw) 
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

See also: `PetscDraw`, `PetscDrawPushCurrentPoint()`, `PetscDrawPopCurrentPoint()`, `PetscDrawSetCurrentPoint()`

# External Links
$(_doc_external("Draw/PetscDrawGetBoundingBox"))
"""
function PetscDrawGetBoundingBox(petsclib::PetscLibType, draw::PetscDraw)
    error("PetscDrawGetBoundingBox: no generated method for these argument types")
end

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
	xl::PetscReal,yl::PetscReal,xr::PetscReal,yr::PetscReal = PetscDrawGetCoordinates(petsclib::PetscLibType, draw::PetscDraw) 
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

See also: `PetscDraw`, `PetscDrawSetCoordinates()`

# External Links
$(_doc_external("Draw/PetscDrawGetCoordinates"))
"""
function PetscDrawGetCoordinates(petsclib::PetscLibType, draw::PetscDraw)
    error("PetscDrawGetCoordinates: no generated method for these argument types")
end

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
	x::PetscReal,y::PetscReal = PetscDrawGetCurrentPoint(petsclib::PetscLibType, draw::PetscDraw) 
Gets the current draw point, some codes use this point to determine where to draw next

Not Collective

Input Parameter:
- `draw` - the drawing context

Output Parameters:
- `x` - horizontal coordinate of the current point
- `y` - vertical coordinate of the current point

Level: intermediate

See also: `PetscDraw`, `PetscDrawPushCurrentPoint()`, `PetscDrawPopCurrentPoint()`, `PetscDrawSetCurrentPoint()`

# External Links
$(_doc_external("Draw/PetscDrawGetCurrentPoint"))
"""
function PetscDrawGetCurrentPoint(petsclib::PetscLibType, draw::PetscDraw)
    error("PetscDrawGetCurrentPoint: no generated method for these argument types")
end

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
	mtype::PetscDrawMarkerType = PetscDrawGetMarkerType(petsclib::PetscLibType, draw::PetscDraw) 
gets the type of marker to display with `PetscDrawMarker()`

Not Collective

Input Parameters:
- `draw`  - the drawing context
- `mtype` - either `PETSC_DRAW_MARKER_CROSS` (default) or `PETSC_DRAW_MARKER_POINT`

Level: beginner

See also: `PetscDraw`, `PetscDrawPoint()`, `PetscDrawMarker()`, `PetscDrawSetMarkerType()`, `PetscDrawMarkerType`

# External Links
$(_doc_external("Draw/PetscDrawGetMarkerType"))
"""
function PetscDrawGetMarkerType(petsclib::PetscLibType, draw::PetscDraw)
    error("PetscDrawGetMarkerType: no generated method for these argument types")
end

@for_petsc function PetscDrawGetMarkerType(petsclib::$UnionPetscLib, draw::PetscDraw )
	mtype_ = Ref{PetscDrawMarkerType}()

    @chk ccall(
               (:PetscDrawGetMarkerType, $petsc_library),
               PetscErrorCode,
               (PetscDraw, Ptr{PetscDrawMarkerType}),
               draw, mtype_,
              )

	mtype = mtype_[]

	return mtype
end 

"""
	button::PetscDrawButton,x_user::PetscReal,y_user::PetscReal,x_phys::PetscReal,y_phys::PetscReal = PetscDrawGetMouseButton(petsclib::PetscLibType, draw::PetscDraw) 
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

See also: `PetscDraw`, `PetscDrawButton`

# External Links
$(_doc_external("Draw/PetscDrawGetMouseButton"))
"""
function PetscDrawGetMouseButton(petsclib::PetscLibType, draw::PetscDraw)
    error("PetscDrawGetMouseButton: no generated method for these argument types")
end

@for_petsc function PetscDrawGetMouseButton(petsclib::$UnionPetscLib, draw::PetscDraw )
	button_ = Ref{PetscDrawButton}()
	x_user_ = Ref{$PetscReal}()
	y_user_ = Ref{$PetscReal}()
	x_phys_ = Ref{$PetscReal}()
	y_phys_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscDrawGetMouseButton, $petsc_library),
               PetscErrorCode,
               (PetscDraw, Ptr{PetscDrawButton}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
               draw, button_, x_user_, y_user_, x_phys_, y_phys_,
              )

	button = button_[]
	x_user = x_user_[]
	y_user = y_user_[]
	x_phys = x_phys_[]
	y_phys = y_phys_[]

	return button,x_user,y_user,x_phys,y_phys
end 

"""
	lpause::PetscReal = PetscDrawGetPause(petsclib::PetscLibType, draw::PetscDraw) 
Gets the amount of time that program pauses after
a `PetscDrawPause()` is called.

Not Collective

Input Parameters:
- `draw`   - the drawing object
- `lpause` - number of seconds to pause, -1 implies until user input

Level: intermediate

See also: `PetscDraw`, `PetscDrawSetPause()`, `PetscDrawPause()`

# External Links
$(_doc_external("Draw/PetscDrawGetPause"))
"""
function PetscDrawGetPause(petsclib::PetscLibType, draw::PetscDraw)
    error("PetscDrawGetPause: no generated method for these argument types")
end

@for_petsc function PetscDrawGetPause(petsclib::$UnionPetscLib, draw::PetscDraw )
	lpause_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscDrawGetPause, $petsc_library),
               PetscErrorCode,
               (PetscDraw, Ptr{$PetscReal}),
               draw, lpause_,
              )

	lpause = lpause_[]

	return lpause
end 

"""
	popup::PetscDraw = PetscDrawGetPopup(petsclib::PetscLibType, draw::PetscDraw) 
Creates a popup window associated with a `PetscDraw` window.

Collective

Input Parameter:
- `draw` - the original window

Output Parameter:
- `popup` - the new popup window

Level: advanced

See also: `PetscDraw`, `PetscDrawScalePopup()`, `PetscDrawCreate()`

# External Links
$(_doc_external("Draw/PetscDrawGetPopup"))
"""
function PetscDrawGetPopup(petsclib::PetscLibType, draw::PetscDraw)
    error("PetscDrawGetPopup: no generated method for these argument types")
end

@for_petsc function PetscDrawGetPopup(petsclib::$UnionPetscLib, draw::PetscDraw )
	popup_ = Ref{PetscDraw}()

    @chk ccall(
               (:PetscDrawGetPopup, $petsc_library),
               PetscErrorCode,
               (PetscDraw, Ptr{PetscDraw}),
               draw, popup_,
              )

	popup = popup_[]

	return popup
end 

"""
	sdraw::PetscDraw = PetscDrawGetSingleton(petsclib::PetscLibType, draw::PetscDraw) 
Gain access to a `PetscDraw` object as if it were owned
by the one process.

Collective

Input Parameter:
- `draw` - the original window

Output Parameter:
- `sdraw` - the singleton window

Level: advanced

See also: `PetscDraw`, `PetscDrawRestoreSingleton()`, `PetscViewerGetSingleton()`, `PetscViewerRestoreSingleton()`

# External Links
$(_doc_external("Draw/PetscDrawGetSingleton"))
"""
function PetscDrawGetSingleton(petsclib::PetscLibType, draw::PetscDraw)
    error("PetscDrawGetSingleton: no generated method for these argument types")
end

@for_petsc function PetscDrawGetSingleton(petsclib::$UnionPetscLib, draw::PetscDraw )
	sdraw_ = Ref{PetscDraw}()

    @chk ccall(
               (:PetscDrawGetSingleton, $petsc_library),
               PetscErrorCode,
               (PetscDraw, Ptr{PetscDraw}),
               draw, sdraw_,
              )

	sdraw = sdraw_[]

	return sdraw
end 

"""
	title::Ptr{Cchar} = PetscDrawGetTitle(petsclib::PetscLibType, draw::PetscDraw) 
Gets pointer to title of a `PetscDraw` context.

Not Collective

Input Parameter:
- `draw` - the graphics context

Output Parameter:
- `title` - the title

Level: intermediate

See also: `PetscDraw`, `PetscDrawSetTitle()`

# External Links
$(_doc_external("Draw/PetscDrawGetTitle"))
"""
function PetscDrawGetTitle(petsclib::PetscLibType, draw::PetscDraw)
    error("PetscDrawGetTitle: no generated method for these argument types")
end

@for_petsc function PetscDrawGetTitle(petsclib::$UnionPetscLib, draw::PetscDraw )
	title_ = Ref{Ptr{Cchar}}()

    @chk ccall(
               (:PetscDrawGetTitle, $petsc_library),
               PetscErrorCode,
               (PetscDraw, Ptr{Ptr{Cchar}}),
               draw, title_,
              )

	title = title_[]

	return title
end 

"""
	type::PetscDrawType = PetscDrawGetType(petsclib::PetscLibType, draw::PetscDraw) 
Gets the `PetscDraw` type as a string from the `PetscDraw` object.

Not Collective

Input Parameter:
- `draw` - Krylov context

Output Parameter:
- `type` - name of PetscDraw method

Level: advanced

See also: `PetscDraw`, `PetscDrawType`, `PetscDrawSetType()`, `PetscDrawCreate()`, `PetscObjectTypeCompare()`, `PetscObjectTypeCompareAny()`

# External Links
$(_doc_external("Draw/PetscDrawGetType"))
"""
function PetscDrawGetType(petsclib::PetscLibType, draw::PetscDraw)
    error("PetscDrawGetType: no generated method for these argument types")
end

@for_petsc function PetscDrawGetType(petsclib::$UnionPetscLib, draw::PetscDraw )
	type_ = Ref{PetscDrawType}()

    @chk ccall(
               (:PetscDrawGetType, $petsc_library),
               PetscErrorCode,
               (PetscDraw, Ptr{PetscDrawType}),
               draw, type_,
              )

	type = type_[] == C_NULL ? "" : unsafe_string(type_[])

	return type
end 

"""
	xl::PetscReal,yl::PetscReal,xr::PetscReal,yr::PetscReal = PetscDrawGetViewPort(petsclib::PetscLibType, draw::PetscDraw) 
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

See also: `PetscDraw`, `PetscDrawSplitViewPort()`, `PetscDrawSetViewPort()`

# External Links
$(_doc_external("Draw/PetscDrawGetViewPort"))
"""
function PetscDrawGetViewPort(petsclib::PetscLibType, draw::PetscDraw)
    error("PetscDrawGetViewPort: no generated method for these argument types")
end

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
	w::Cint,h::Cint = PetscDrawGetWindowSize(petsclib::PetscLibType, draw::PetscDraw) 
Gets the size of the window.

Not Collective

Input Parameter:
- `draw` - the window

Output Parameters:
- `w` - the window width
- `h` - the window height

Level: intermediate

See also: `PetscDraw`, `PetscDrawResizeWindow()`, `PetscDrawCheckResizedWindow()`

# External Links
$(_doc_external("Draw/PetscDrawGetWindowSize"))
"""
function PetscDrawGetWindowSize(petsclib::PetscLibType, draw::PetscDraw)
    error("PetscDrawGetWindowSize: no generated method for these argument types")
end

@for_petsc function PetscDrawGetWindowSize(petsclib::$UnionPetscLib, draw::PetscDraw )
	w_ = Ref{Cint}()
	h_ = Ref{Cint}()

    @chk ccall(
               (:PetscDrawGetWindowSize, $petsc_library),
               PetscErrorCode,
               (PetscDraw, Ptr{Cint}, Ptr{Cint}),
               draw, w_, h_,
              )

	w = w_[]
	h = h_[]

	return w,h
end 

"""
	PetscDrawHGAddValue(petsclib::PetscLibType, hist::PetscDrawHG, value::PetscReal) 
Adds another value to the histogram.

Logically Collective

Input Parameters:
- `hist`  - The histogram
- `value` - The value

Level: intermediate

See also: `PetscDrawHGCreate()`, `PetscDrawHG`, `PetscDrawHGDraw()`, `PetscDrawHGReset()`, `PetscDrawHGAddWeightedValue()`

# External Links
$(_doc_external("Draw/PetscDrawHGAddValue"))
"""
function PetscDrawHGAddValue(petsclib::PetscLibType, hist::PetscDrawHG, value::Real)
    error("PetscDrawHGAddValue: no generated method for these argument types")
end

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
	PetscDrawHGAddWeightedValue(petsclib::PetscLibType, hist::PetscDrawHG, value::PetscReal, weight::PetscReal) 
Adds another value to the histogram with a weight.

Logically Collective

Input Parameters:
- `hist`   - The histogram
- `value`  - The value
- `weight` - The value weight

Level: intermediate

See also: `PetscDrawHGCreate()`, `PetscDrawHG`, `PetscDrawHGDraw()`, `PetscDrawHGReset()`, `PetscDrawHGAddValue()`

# External Links
$(_doc_external("Draw/PetscDrawHGAddWeightedValue"))
"""
function PetscDrawHGAddWeightedValue(petsclib::PetscLibType, hist::PetscDrawHG, value::Real, weight::Real)
    error("PetscDrawHGAddWeightedValue: no generated method for these argument types")
end

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
	PetscDrawHGCalcStats(petsclib::PetscLibType, hist::PetscDrawHG, calc::PetscBool) 
Turns on calculation of descriptive statistics associated with the histogram

Not Collective

Input Parameters:
- `hist` - The histogram context
- `calc` - Flag for calculation

Level: intermediate

See also: `PetscDrawHG`, `PetscDrawHGCreate()`, `PetscDrawHGAddValue()`, `PetscDrawHGView()`, `PetscDrawHGDraw()`

# External Links
$(_doc_external("Draw/PetscDrawHGCalcStats"))
"""
function PetscDrawHGCalcStats(petsclib::PetscLibType, hist::PetscDrawHG, calc::PetscBool)
    error("PetscDrawHGCalcStats: no generated method for these argument types")
end

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
	hist::PetscDrawHG = PetscDrawHGCreate(petsclib::PetscLibType, draw::PetscDraw, bins::Cint) 
Creates a histogram data structure.

Collective

Input Parameters:
- `draw` - The window where the graph will be made
- `bins` - The number of bins to use

Output Parameter:
- `hist` - The histogram context

Level: intermediate

See also: `PetscDrawHGDestroy()`, `PetscDrawHG`, `PetscDrawBarCreate()`, `PetscDrawBar`, `PetscDrawLGCreate()`, `PetscDrawLG`, `PetscDrawSPCreate()`, `PetscDrawSP`,
`PetscDrawHGSetNumberBins()`, `PetscDrawHGReset()`, `PetscDrawHGAddValue()`, `PetscDrawHGDraw()`, `PetscDrawHGSave()`, `PetscDrawHGView()`, `PetscDrawHGSetColor()`,
`PetscDrawHGSetLimits()`, `PetscDrawHGCalcStats()`, `PetscDrawHGIntegerBins()`, `PetscDrawHGGetAxis()`, `PetscDrawAxis`, `PetscDrawHGGetDraw()`

# External Links
$(_doc_external("Draw/PetscDrawHGCreate"))
"""
function PetscDrawHGCreate(petsclib::PetscLibType, draw::PetscDraw, bins::Cint)
    error("PetscDrawHGCreate: no generated method for these argument types")
end

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
	PetscDrawHGDestroy(petsclib::PetscLibType, hist::Union{PetscDrawHG, Ref{PetscDrawHG}}) 
Frees all space taken up by histogram data structure.

Collective

Input Parameter:
- `hist` - The histogram context

Level: intermediate

See also: `PetscDrawHGCreate()`, `PetscDrawHG`

# External Links
$(_doc_external("Draw/PetscDrawHGDestroy"))
"""
function PetscDrawHGDestroy(petsclib::PetscLibType, hist::Union{PetscDrawHG, Ref{PetscDrawHG}})
    error("PetscDrawHGDestroy: no generated method for these argument types")
end

@for_petsc function PetscDrawHGDestroy(petsclib::$UnionPetscLib, hist::Union{PetscDrawHG, Ref{PetscDrawHG}} )
	hist_ = hist isa Base.RefValue ? hist : Ref{PetscDrawHG}(hist)

    @chk ccall(
               (:PetscDrawHGDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscDrawHG},),
               hist_,
              )


	return nothing
end 

"""
	PetscDrawHGDraw(petsclib::PetscLibType, hist::PetscDrawHG) 
Redraws a histogram.

Collective

Input Parameter:
- `hist` - The histogram context

Level: intermediate

See also: `PetscDrawHGCreate()`, `PetscDrawHG`, `PetscDrawHGAddValue()`, `PetscDrawHGReset()`

# External Links
$(_doc_external("Draw/PetscDrawHGDraw"))
"""
function PetscDrawHGDraw(petsclib::PetscLibType, hist::PetscDrawHG)
    error("PetscDrawHGDraw: no generated method for these argument types")
end

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
	axis::PetscDrawAxis = PetscDrawHGGetAxis(petsclib::PetscLibType, hist::PetscDrawHG) 
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

See also: `PetscDrawHG`, `PetscDrawAxis`, `PetscDrawHGCreate()`, `PetscDrawHGAddValue()`, `PetscDrawHGView()`, `PetscDrawHGDraw()`, `PetscDrawHGSetColor()`, `PetscDrawHGSetLimits()`

# External Links
$(_doc_external("Draw/PetscDrawHGGetAxis"))
"""
function PetscDrawHGGetAxis(petsclib::PetscLibType, hist::PetscDrawHG)
    error("PetscDrawHGGetAxis: no generated method for these argument types")
end

@for_petsc function PetscDrawHGGetAxis(petsclib::$UnionPetscLib, hist::PetscDrawHG )
	axis_ = Ref{PetscDrawAxis}()

    @chk ccall(
               (:PetscDrawHGGetAxis, $petsc_library),
               PetscErrorCode,
               (PetscDrawHG, Ptr{PetscDrawAxis}),
               hist, axis_,
              )

	axis = axis_[]

	return axis
end 

"""
	draw::PetscDraw = PetscDrawHGGetDraw(petsclib::PetscLibType, hist::PetscDrawHG) 
Gets the draw context associated with a histogram.

Not Collective, draw is parallel if hist is parallel

Input Parameter:
- `hist` - The histogram context

Output Parameter:
- `draw` - The draw context

Level: intermediate

See also: `PetscDraw`, `PetscDrawHG`, `PetscDrawHGCreate()`, `PetscDrawHGAddValue()`, `PetscDrawHGView()`, `PetscDrawHGDraw()`, `PetscDrawHGSetColor()`, `PetscDrawAxis`, `PetscDrawHGSetLimits()`

# External Links
$(_doc_external("Draw/PetscDrawHGGetDraw"))
"""
function PetscDrawHGGetDraw(petsclib::PetscLibType, hist::PetscDrawHG)
    error("PetscDrawHGGetDraw: no generated method for these argument types")
end

@for_petsc function PetscDrawHGGetDraw(petsclib::$UnionPetscLib, hist::PetscDrawHG )
	draw_ = Ref{PetscDraw}()

    @chk ccall(
               (:PetscDrawHGGetDraw, $petsc_library),
               PetscErrorCode,
               (PetscDrawHG, Ptr{PetscDraw}),
               hist, draw_,
              )

	draw = draw_[]

	return draw
end 

"""
	PetscDrawHGIntegerBins(petsclib::PetscLibType, hist::PetscDrawHG, ints::PetscBool) 
Turns on integer width bins

Not Collective

Input Parameters:
- `hist` - The histogram context
- `ints` - Flag for integer width bins

Level: intermediate

See also: `PetscDrawHG`, `PetscDrawHGCreate()`, `PetscDrawHGAddValue()`, `PetscDrawHGView()`, `PetscDrawHGDraw()`, `PetscDrawHGSetColor()`

# External Links
$(_doc_external("Draw/PetscDrawHGIntegerBins"))
"""
function PetscDrawHGIntegerBins(petsclib::PetscLibType, hist::PetscDrawHG, ints::PetscBool)
    error("PetscDrawHGIntegerBins: no generated method for these argument types")
end

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
	PetscDrawHGReset(petsclib::PetscLibType, hist::PetscDrawHG) 
Clears histogram to allow for reuse with new data.

Logically Collective

Input Parameter:
- `hist` - The histogram context.

Level: intermediate

See also: `PetscDrawHGCreate()`, `PetscDrawHG`, `PetscDrawHGDraw()`, `PetscDrawHGAddValue()`

# External Links
$(_doc_external("Draw/PetscDrawHGReset"))
"""
function PetscDrawHGReset(petsclib::PetscLibType, hist::PetscDrawHG)
    error("PetscDrawHGReset: no generated method for these argument types")
end

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
	PetscDrawHGSave(petsclib::PetscLibType, hg::PetscDrawHG) 
Saves a drawn image

Collective

Input Parameter:
- `hg` - The histogram context

Level: intermediate

See also: `PetscDrawSave()`, `PetscDrawHGCreate()`, `PetscDrawHGGetDraw()`, `PetscDrawSetSave()`, `PetscDrawHGDraw()`

# External Links
$(_doc_external("Draw/PetscDrawHGSave"))
"""
function PetscDrawHGSave(petsclib::PetscLibType, hg::PetscDrawHG)
    error("PetscDrawHGSave: no generated method for these argument types")
end

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
	PetscDrawHGSetColor(petsclib::PetscLibType, hist::PetscDrawHG, color::Cint) 
Sets the color the bars will be drawn with.

Logically Collective

Input Parameters:
- `hist`  - The histogram context
- `color` - one of the colors defined in petscdraw.h or `PETSC_DRAW_ROTATE` to make each bar a different color

Level: intermediate

See also: `PetscDrawHG`, `PetscDrawHGCreate()`, `PetscDrawHGGetDraw()`, `PetscDrawSetSave()`, `PetscDrawSave()`, `PetscDrawHGDraw()`, `PetscDrawHGGetAxis()`

# External Links
$(_doc_external("Draw/PetscDrawHGSetColor"))
"""
function PetscDrawHGSetColor(petsclib::PetscLibType, hist::PetscDrawHG, color::Cint)
    error("PetscDrawHGSetColor: no generated method for these argument types")
end

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
	PetscDrawHGSetLimits(petsclib::PetscLibType, hist::PetscDrawHG, x_min::PetscReal, x_max::PetscReal, y_min::Cint, y_max::Cint) 
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

See also: `PetscDrawHG`, `PetscDrawHGCreate()`, `PetscDrawHGGetDraw()`, `PetscDrawSetSave()`, `PetscDrawSave()`, `PetscDrawHGDraw()`, `PetscDrawHGGetAxis()`

# External Links
$(_doc_external("Draw/PetscDrawHGSetLimits"))
"""
function PetscDrawHGSetLimits(petsclib::PetscLibType, hist::PetscDrawHG, x_min::Real, x_max::Real, y_min::Cint, y_max::Cint)
    error("PetscDrawHGSetLimits: no generated method for these argument types")
end

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
	PetscDrawHGSetNumberBins(petsclib::PetscLibType, hist::PetscDrawHG, bins::Cint) 
Change the number of bins that are to be drawn in the histogram

Logically Collective

Input Parameters:
- `hist` - The histogram context.
- `bins` - The number of bins.

Level: intermediate

See also: `PetscDrawHGCreate()`, `PetscDrawHG`, `PetscDrawHGDraw()`, `PetscDrawHGIntegerBins()`

# External Links
$(_doc_external("Draw/PetscDrawHGSetNumberBins"))
"""
function PetscDrawHGSetNumberBins(petsclib::PetscLibType, hist::PetscDrawHG, bins::Cint)
    error("PetscDrawHGSetNumberBins: no generated method for these argument types")
end

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
	PetscDrawHGView(petsclib::PetscLibType, hist::PetscDrawHG, viewer::PetscViewer) 
Prints the histogram information to a viewer

Not Collective

Input Parameters:
- `hist`   - The histogram context
- `viewer` - The viewer to view it with

Level: beginner

See also: `PetscDrawHG`, `PetscViewer`, `PetscDrawHGCreate()`, `PetscDrawHGGetDraw()`, `PetscDrawSetSave()`, `PetscDrawSave()`, `PetscDrawHGDraw()`

# External Links
$(_doc_external("Draw/PetscDrawHGView"))
"""
function PetscDrawHGView(petsclib::PetscLibType, hist::PetscDrawHG, viewer::PetscViewer)
    error("PetscDrawHGView: no generated method for these argument types")
end

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
	PetscDrawIndicatorFunction(petsclib::PetscLibType, draw::PetscDraw, xmin::PetscReal, xmax::PetscReal, ymin::PetscReal, ymax::PetscReal, c::Cint, indicator::external, ctx::Ptr{Cvoid}) 
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

See also: `PetscDraw`

# External Links
$(_doc_external("Draw/PetscDrawIndicatorFunction"))
"""
function PetscDrawIndicatorFunction(petsclib::PetscLibType, draw::PetscDraw, xmin::Real, xmax::Real, ymin::Real, ymax::Real, c::Cint, indicator::external, ctx::Ptr{Cvoid})
    error("PetscDrawIndicatorFunction: no generated method for these argument types")
end

@for_petsc function PetscDrawIndicatorFunction(petsclib::$UnionPetscLib, draw::PetscDraw, xmin::$PetscReal, xmax::$PetscReal, ymin::$PetscReal, ymax::$PetscReal, c::Cint, indicator::external, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:PetscDrawIndicatorFunction, $petsc_library),
               PetscErrorCode,
               (PetscDraw, $PetscReal, $PetscReal, $PetscReal, $PetscReal, Cint, external, Ptr{Cvoid}),
               draw, xmin, xmax, ymin, ymax, c, indicator, ctx,
              )


	return nothing
end 

"""
	PetscDrawInitializePackage(petsclib::PetscLibType) 
This function initializes everything in the `PetscDraw` package. It is called
from PetscDLLibraryRegister_petsc() when using dynamic libraries, and on the call to `PetscInitialize()`
when using shared or static libraries.

Level: developer

See also: `PetscDraw`, `PetscInitialize()`

# External Links
$(_doc_external("Draw/PetscDrawInitializePackage"))
"""
function PetscDrawInitializePackage(petsclib::PetscLibType)
    error("PetscDrawInitializePackage: no generated method for these argument types")
end

@for_petsc function PetscDrawInitializePackage(petsclib::$UnionPetscLib)

    @chk ccall(
               (:PetscDrawInitializePackage, $petsc_library),
               PetscErrorCode,
               (),
              )


	return nothing
end 

"""
	yes::PetscBool = PetscDrawIsNull(petsclib::PetscLibType, draw::PetscDraw) 
Returns `PETSC_TRUE` if draw is a null draw object.

Not Collective

Input Parameter:
- `draw` - the draw context

Output Parameter:
- `yes` - `PETSC_TRUE` if it is a null draw object; otherwise `PETSC_FALSE`

Level: advanced

See also: `PetscDraw`, `PETSC_DRAW_NULL`, `PetscDrawOpenX()`

# External Links
$(_doc_external("Draw/PetscDrawIsNull"))
"""
function PetscDrawIsNull(petsclib::PetscLibType, draw::PetscDraw)
    error("PetscDrawIsNull: no generated method for these argument types")
end

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
	PetscDrawLGAddCommonPoint(petsclib::PetscLibType, lg::PetscDrawLG, x::PetscReal, y::Vector{PetscReal}) 
Adds another point to each of the line graphs. All the points share
the same new X coordinate.  The new point must have an X coordinate larger than the old points.

Logically Collective

Input Parameters:
- `lg` - the line graph context
- `x`  - the common x coordinate point
- `y`  - the new y coordinate point for each curve.

Level: intermediate

See also: `PetscDrawLG`, `PetscDrawLGCreate()`, `PetscDrawLGAddPoints()`, `PetscDrawLGAddPoint()`, `PetscDrawLGReset()`, `PetscDrawLGDraw()`

# External Links
$(_doc_external("Draw/PetscDrawLGAddCommonPoint"))
"""
function PetscDrawLGAddCommonPoint(petsclib::PetscLibType, lg::PetscDrawLG, x::Real, y::AbstractVector{<:Number})
    error("PetscDrawLGAddCommonPoint: no generated method for these argument types")
end

@for_petsc function PetscDrawLGAddCommonPoint(petsclib::$UnionPetscLib, lg::PetscDrawLG, x::$PetscReal, y::Vector{$PetscReal} )

    @chk ccall(
               (:PetscDrawLGAddCommonPoint, $petsc_library),
               PetscErrorCode,
               (PetscDrawLG, $PetscReal, Ptr{$PetscReal}),
               lg, x, y,
              )


	return nothing
end 

"""
	PetscDrawLGAddPoint(petsclib::PetscLibType, lg::PetscDrawLG, x::Vector{PetscReal}, y::Vector{PetscReal}) 
Adds another point to each of the line graphs.
The new point must have an X coordinate larger than the old points.

Logically Collective

Input Parameters:
- `lg` - the line graph context
- `x`  - array containing the x coordinate for the point on each curve
- `y`  - array containing the y coordinate for the point on each curve

Level: intermediate

See also: `PetscDrawLG`, `PetscDrawLGCreate()`, `PetscDrawLGAddPoints()`, `PetscDrawLGAddCommonPoint()`, `PetscDrawLGReset()`, `PetscDrawLGDraw()`

# External Links
$(_doc_external("Draw/PetscDrawLGAddPoint"))
"""
function PetscDrawLGAddPoint(petsclib::PetscLibType, lg::PetscDrawLG, x::AbstractVector{<:Number}, y::AbstractVector{<:Number})
    error("PetscDrawLGAddPoint: no generated method for these argument types")
end

@for_petsc function PetscDrawLGAddPoint(petsclib::$UnionPetscLib, lg::PetscDrawLG, x::Vector{$PetscReal}, y::Vector{$PetscReal} )

    @chk ccall(
               (:PetscDrawLGAddPoint, $petsc_library),
               PetscErrorCode,
               (PetscDrawLG, Ptr{$PetscReal}, Ptr{$PetscReal}),
               lg, x, y,
              )


	return nothing
end 

"""
	PetscDrawLGAddPoints(petsclib::PetscLibType, lg::PetscDrawLG, n::PetscInt, xx::Union{Ptr, AbstractArray{PetscReal}}, yy::Union{Ptr, AbstractArray{PetscReal}}) 
Adds several points to each of the line graphs.
The new points must have an X coordinate larger than the old points.

Logically Collective

Input Parameters:
- `lg` - the line graph context
- `xx` - array of pointers that point to arrays containing the new x coordinates for each curve.
- `yy` - array of pointers that point to arrays containing the new y points for each curve.
- `n`  - number of points being added

Level: intermediate

See also: `PetscDrawLG`, `PetscDrawLGCreate()`, `PetscDrawLGAddPoint()`, `PetscDrawLGAddCommonPoint()`, `PetscDrawLGReset()`, `PetscDrawLGDraw()`

# External Links
$(_doc_external("Draw/PetscDrawLGAddPoints"))
"""
function PetscDrawLGAddPoints(petsclib::PetscLibType, lg::PetscDrawLG, n::Integer, xx::Union{Ptr, AbstractArray{<:Number}}, yy::Union{Ptr, AbstractArray{<:Number}})
    error("PetscDrawLGAddPoints: no generated method for these argument types")
end

@for_petsc function PetscDrawLGAddPoints(petsclib::$UnionPetscLib, lg::PetscDrawLG, n::$PetscInt, xx::Union{Ptr, AbstractArray{$PetscReal}}, yy::Union{Ptr, AbstractArray{$PetscReal}} )
	xx_ = Ref{Ptr{$PetscReal}}(xx isa Ptr ? xx : pointer(xx))
	yy_ = Ref{Ptr{$PetscReal}}(yy isa Ptr ? yy : pointer(yy))

    @chk ccall(
               (:PetscDrawLGAddPoints, $petsc_library),
               PetscErrorCode,
               (PetscDrawLG, $PetscInt, Ptr{Ptr{$PetscReal}}, Ptr{Ptr{$PetscReal}}),
               lg, n, xx_, yy_,
              )


	return nothing
end 

"""
	outlg::PetscDrawLG = PetscDrawLGCreate(petsclib::PetscLibType, draw::PetscDraw, dim::PetscInt) 
Creates a line graph data structure.

Collective

Input Parameters:
- `draw` - the window where the graph will be made.
- `dim`  - the number of curves which will be drawn

Output Parameter:
- `outlg` - the line graph context

Level: intermediate

See also: `PetscDrawLGDestroy()`, `PetscDrawLGAddPoint()`, `PetscDrawLGAddCommonPoint()`, `PetscDrawLGAddPoints()`, `PetscDrawLGDraw()`, `PetscDrawLGSave()`,
`PetscDrawLGView()`, `PetscDrawLGReset()`, `PetscDrawLGSetDimension()`, `PetscDrawLGGetDimension()`, `PetscDrawLGSetLegend()`, `PetscDrawLGGetAxis()`,
`PetscDrawLGGetDraw()`, `PetscDrawLGSetUseMarkers()`, `PetscDrawLGSetLimits()`, `PetscDrawLGSetColors()`, `PetscDrawLGSetOptionsPrefix()`, `PetscDrawLGSetFromOptions()`

# External Links
$(_doc_external("Draw/PetscDrawLGCreate"))
"""
function PetscDrawLGCreate(petsclib::PetscLibType, draw::PetscDraw, dim::Integer)
    error("PetscDrawLGCreate: no generated method for these argument types")
end

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
	PetscDrawLGDestroy(petsclib::PetscLibType, lg::Union{PetscDrawLG, Ref{PetscDrawLG}}) 
Frees all space taken up by line graph data structure.

Collective

Input Parameter:
- `lg` - the line graph context

Level: intermediate

See also: `PetscDrawLG`, `PetscDrawLGCreate()`

# External Links
$(_doc_external("Draw/PetscDrawLGDestroy"))
"""
function PetscDrawLGDestroy(petsclib::PetscLibType, lg::Union{PetscDrawLG, Ref{PetscDrawLG}})
    error("PetscDrawLGDestroy: no generated method for these argument types")
end

@for_petsc function PetscDrawLGDestroy(petsclib::$UnionPetscLib, lg::Union{PetscDrawLG, Ref{PetscDrawLG}} )
	lg_ = lg isa Base.RefValue ? lg : Ref{PetscDrawLG}(lg)

    @chk ccall(
               (:PetscDrawLGDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscDrawLG},),
               lg_,
              )


	return nothing
end 

"""
	PetscDrawLGDraw(petsclib::PetscLibType, lg::PetscDrawLG) 
Redraws a line graph.

Collective

Input Parameter:
- `lg` - the line graph context

Level: intermediate

See also: `PetscDrawLG`, `PetscDrawSPDraw()`, `PetscDrawLGSPDraw()`, `PetscDrawLGReset()`

# External Links
$(_doc_external("Draw/PetscDrawLGDraw"))
"""
function PetscDrawLGDraw(petsclib::PetscLibType, lg::PetscDrawLG)
    error("PetscDrawLGDraw: no generated method for these argument types")
end

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
	axis::PetscDrawAxis = PetscDrawLGGetAxis(petsclib::PetscLibType, lg::PetscDrawLG) 
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

See also: `PetscDrawLGCreate()`, `PetscDrawAxis`, `PetscDrawLG`

# External Links
$(_doc_external("Draw/PetscDrawLGGetAxis"))
"""
function PetscDrawLGGetAxis(petsclib::PetscLibType, lg::PetscDrawLG)
    error("PetscDrawLGGetAxis: no generated method for these argument types")
end

@for_petsc function PetscDrawLGGetAxis(petsclib::$UnionPetscLib, lg::PetscDrawLG )
	axis_ = Ref{PetscDrawAxis}()

    @chk ccall(
               (:PetscDrawLGGetAxis, $petsc_library),
               PetscErrorCode,
               (PetscDrawLG, Ptr{PetscDrawAxis}),
               lg, axis_,
              )

	axis = axis_[]

	return axis
end 

"""
	dim::PetscInt,n::PetscInt,x::Vector{PetscReal},y::Vector{PetscReal} = PetscDrawLGGetData(petsclib::PetscLibType, lg::PetscDrawLG) 
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

See also: `PetscDrawLGC`, `PetscDrawLGCreate()`, `PetscDrawLGGetDimension()`

# External Links
$(_doc_external("Draw/PetscDrawLGGetData"))
"""
function PetscDrawLGGetData(petsclib::PetscLibType, lg::PetscDrawLG)
    error("PetscDrawLGGetData: no generated method for these argument types")
end

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
	x = x_[] == C_NULL ? $PetscReal[] : unsafe_wrap(Array, x_[], n * dim; own = false)
	y = y_[] == C_NULL ? $PetscReal[] : unsafe_wrap(Array, y_[], n * dim; own = false)

	return dim,n,x,y
end 

"""
	dim::PetscInt = PetscDrawLGGetDimension(petsclib::PetscLibType, lg::PetscDrawLG) 
Get the number of curves that are to be drawn.

Not Collective

Input Parameter:
- `lg` - the line graph context.

Output Parameter:
- `dim` - the number of curves.

Level: intermediate

See also: `PetscDrawLGC`, `PetscDrawLGCreate()`, `PetscDrawLGSetDimension()`

# External Links
$(_doc_external("Draw/PetscDrawLGGetDimension"))
"""
function PetscDrawLGGetDimension(petsclib::PetscLibType, lg::PetscDrawLG)
    error("PetscDrawLGGetDimension: no generated method for these argument types")
end

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
	draw::PetscDraw = PetscDrawLGGetDraw(petsclib::PetscLibType, lg::PetscDrawLG) 
Gets the draw context associated with a line graph.

Not Collective, if lg is parallel then draw is parallel

Input Parameter:
- `lg` - the line graph context

Output Parameter:
- `draw` - the draw context

Level: intermediate

See also: `PetscDrawLGCreate()`, `PetscDraw`, `PetscDrawLG`

# External Links
$(_doc_external("Draw/PetscDrawLGGetDraw"))
"""
function PetscDrawLGGetDraw(petsclib::PetscLibType, lg::PetscDrawLG)
    error("PetscDrawLGGetDraw: no generated method for these argument types")
end

@for_petsc function PetscDrawLGGetDraw(petsclib::$UnionPetscLib, lg::PetscDrawLG )
	draw_ = Ref{PetscDraw}()

    @chk ccall(
               (:PetscDrawLGGetDraw, $petsc_library),
               PetscErrorCode,
               (PetscDrawLG, Ptr{PetscDraw}),
               lg, draw_,
              )

	draw = draw_[]

	return draw
end 

"""
	PetscDrawLGReset(petsclib::PetscLibType, lg::PetscDrawLG) 
Clears line graph to allow for reuse with new data.

Logically Collective

Input Parameter:
- `lg` - the line graph context.

Level: intermediate

See also: `PetscDrawLG`, `PetscDrawLGCreate()`

# External Links
$(_doc_external("Draw/PetscDrawLGReset"))
"""
function PetscDrawLGReset(petsclib::PetscLibType, lg::PetscDrawLG)
    error("PetscDrawLGReset: no generated method for these argument types")
end

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
	PetscDrawLGSPDraw(petsclib::PetscLibType, lg::PetscDrawLG, spin::PetscDrawSP) 
Redraws a line graph and a scatter plot on the same `PetscDraw` they must share

Collective

Input Parameters:
- `lg`   - the line graph context
- `spin` - the scatter plot

Level: intermediate

See also: `PetscDrawLGDraw()`, `PetscDrawSPDraw()`

# External Links
$(_doc_external("Draw/PetscDrawLGSPDraw"))
"""
function PetscDrawLGSPDraw(petsclib::PetscLibType, lg::PetscDrawLG, spin::PetscDrawSP)
    error("PetscDrawLGSPDraw: no generated method for these argument types")
end

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
	PetscDrawLGSave(petsclib::PetscLibType, lg::PetscDrawLG) 
Saves a drawn image

Collective

Input Parameter:
- `lg` - The line graph context

Level: intermediate

See also: `PetscDrawLG`, `PetscDrawSave()`, `PetscDrawLGCreate()`, `PetscDrawLGGetDraw()`, `PetscDrawSetSave()`

# External Links
$(_doc_external("Draw/PetscDrawLGSave"))
"""
function PetscDrawLGSave(petsclib::PetscLibType, lg::PetscDrawLG)
    error("PetscDrawLGSave: no generated method for these argument types")
end

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
	PetscDrawLGSetColors(petsclib::PetscLibType, lg::PetscDrawLG, colors::Vector{Cint}) 
Sets the color of each line graph drawn

Logically Collective

Input Parameters:
- `lg`     - the line graph context.
- `colors` - the colors, an array of length the value set with `PetscDrawLGSetDimension()`

Level: intermediate

See also: `PetscDrawLG`, `PetscDrawLGCreate()`, `PetscDrawLGSetDimension()`, `PetscDrawLGGetDimension()`

# External Links
$(_doc_external("Draw/PetscDrawLGSetColors"))
"""
function PetscDrawLGSetColors(petsclib::PetscLibType, lg::PetscDrawLG, colors::Vector{Cint})
    error("PetscDrawLGSetColors: no generated method for these argument types")
end

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
	PetscDrawLGSetDimension(petsclib::PetscLibType, lg::PetscDrawLG, dim::PetscInt) 
Change the number of curves that are to be drawn.

Logically Collective

Input Parameters:
- `lg`  - the line graph context.
- `dim` - the number of curves.

Level: intermediate

See also: `PetscDrawLGCreate()`, `PetscDrawLGGetDimension()`

# External Links
$(_doc_external("Draw/PetscDrawLGSetDimension"))
"""
function PetscDrawLGSetDimension(petsclib::PetscLibType, lg::PetscDrawLG, dim::Integer)
    error("PetscDrawLGSetDimension: no generated method for these argument types")
end

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
	PetscDrawLGSetFromOptions(petsclib::PetscLibType, lg::PetscDrawLG) 
Sets options related to the line graph object

Collective

Input Parameters:
- `lg` - the line graph context

Options Database Key:
- `-lg_use_markers (true|false)` - true means it draws a marker for each point

Level: intermediate

See also: `PetscDrawLG`, `PetscDrawLGDestroy()`, `PetscDrawLGCreate()`

# External Links
$(_doc_external("Draw/PetscDrawLGSetFromOptions"))
"""
function PetscDrawLGSetFromOptions(petsclib::PetscLibType, lg::PetscDrawLG)
    error("PetscDrawLGSetFromOptions: no generated method for these argument types")
end

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
	PetscDrawLGSetLegend(petsclib::PetscLibType, lg::PetscDrawLG, names::String) 
sets the names of each curve plotted

Logically Collective

Input Parameters:
- `lg`    - the line graph context.
- `names` - the names for each curve

Level: intermediate

See also: `PetscDrawLGGetAxis()`, `PetscDrawAxis`, `PetscDrawAxisSetColors()`, `PetscDrawAxisSetLabels()`, `PetscDrawAxisSetHoldLimits()`

# External Links
$(_doc_external("Draw/PetscDrawLGSetLegend"))
"""
function PetscDrawLGSetLegend(petsclib::PetscLibType, lg::PetscDrawLG, names::String)
    error("PetscDrawLGSetLegend: no generated method for these argument types")
end

@for_petsc function PetscDrawLGSetLegend(petsclib::$UnionPetscLib, lg::PetscDrawLG, names::String )
	names_ = Ref{Ptr{Cchar}}(names isa Ptr ? names : pointer(names))

    @chk ccall(
               (:PetscDrawLGSetLegend, $petsc_library),
               PetscErrorCode,
               (PetscDrawLG, Ptr{Ptr{Cchar}}),
               lg, names_,
              )


	return nothing
end 

"""
	PetscDrawLGSetLimits(petsclib::PetscLibType, lg::PetscDrawLG, x_min::PetscReal, x_max::PetscReal, y_min::PetscReal, y_max::PetscReal) 
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

See also: `PetscDrawLGCreate()`, `PetscDrawLG`, `PetscDrawAxis`

# External Links
$(_doc_external("Draw/PetscDrawLGSetLimits"))
"""
function PetscDrawLGSetLimits(petsclib::PetscLibType, lg::PetscDrawLG, x_min::Real, x_max::Real, y_min::Real, y_max::Real)
    error("PetscDrawLGSetLimits: no generated method for these argument types")
end

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
	PetscDrawLGSetOptionsPrefix(petsclib::PetscLibType, lg::PetscDrawLG, prefix::String) 
Sets the prefix used for searching for all
`PetscDrawLG` options in the database.

Logically Collective

Input Parameters:
- `lg`     - the line graph context
- `prefix` - the prefix to prepend to all option names

Level: advanced

See also: `PetscDrawLG`, `PetscDrawLGSetFromOptions()`, `PetscDrawLGCreate()`

# External Links
$(_doc_external("Draw/PetscDrawLGSetOptionsPrefix"))
"""
function PetscDrawLGSetOptionsPrefix(petsclib::PetscLibType, lg::PetscDrawLG, prefix::String)
    error("PetscDrawLGSetOptionsPrefix: no generated method for these argument types")
end

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
	PetscDrawLGSetUseMarkers(petsclib::PetscLibType, lg::PetscDrawLG, flg::PetscBool) 
Causes the line graph object to draw a marker for each data-point.

Logically Collective

Input Parameters:
- `lg`  - the linegraph context
- `flg` - should mark each data point

Options Database Key:
- `-lg_use_markers (true|false)` - true means it draws a marker for each point

Level: intermediate

See also: `PetscDrawLG`, `PetscDrawLGCreate()`

# External Links
$(_doc_external("Draw/PetscDrawLGSetUseMarkers"))
"""
function PetscDrawLGSetUseMarkers(petsclib::PetscLibType, lg::PetscDrawLG, flg::PetscBool)
    error("PetscDrawLGSetUseMarkers: no generated method for these argument types")
end

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
	PetscDrawLGView(petsclib::PetscLibType, lg::PetscDrawLG, viewer::PetscViewer) 
Prints a line graph.

Collective

Input Parameters:
- `lg`     - the line graph context
- `viewer` - the viewer to view it with

Level: beginner

See also: `PetscDrawLG`, `PetscDrawLGCreate()`

# External Links
$(_doc_external("Draw/PetscDrawLGView"))
"""
function PetscDrawLGView(petsclib::PetscLibType, lg::PetscDrawLG, viewer::PetscViewer)
    error("PetscDrawLGView: no generated method for these argument types")
end

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
	PetscDrawLine(petsclib::PetscLibType, draw::PetscDraw, xl::PetscReal, yl::PetscReal, xr::PetscReal, yr::PetscReal, cl::Cint) 
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

See also: `PetscDraw`, `PetscDrawArrow()`, `PetscDrawLineSetWidth()`, `PetscDrawLineGetWidth()`, `PetscDrawRectangle()`, `PetscDrawTriangle()`, `PetscDrawEllipse()`,
`PetscDrawMarker()`, `PetscDrawPoint()`

# External Links
$(_doc_external("Draw/PetscDrawLine"))
"""
function PetscDrawLine(petsclib::PetscLibType, draw::PetscDraw, xl::Real, yl::Real, xr::Real, yr::Real, cl::Cint)
    error("PetscDrawLine: no generated method for these argument types")
end

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
	width::PetscReal = PetscDrawLineGetWidth(petsclib::PetscLibType, draw::PetscDraw) 
Gets the line width for future draws.  The width is
relative to the user coordinates of the window; 0.0 denotes the natural
width; 1.0 denotes the interior viewport.

Not Collective

Input Parameter:
- `draw` - the drawing context

Output Parameter:
- `width` - the width in user coordinates

Level: advanced

See also: `PetscDraw`, `PetscDrawLineSetWidth()`, `PetscDrawLine()`, `PetscDrawArrow()`

# External Links
$(_doc_external("Draw/PetscDrawLineGetWidth"))
"""
function PetscDrawLineGetWidth(petsclib::PetscLibType, draw::PetscDraw)
    error("PetscDrawLineGetWidth: no generated method for these argument types")
end

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
	PetscDrawLineSetWidth(petsclib::PetscLibType, draw::PetscDraw, width::PetscReal) 
Sets the line width for future draws.  The width is
relative to the user coordinates of the window; 0.0 denotes the natural
width; 1.0 denotes the entire viewport.

Not Collective

Input Parameters:
- `draw`  - the drawing context
- `width` - the width in user coordinates

Level: advanced

See also: `PetscDraw`, `PetscDrawLineGetWidth()`, `PetscDrawLine()`, `PetscDrawArrow()`

# External Links
$(_doc_external("Draw/PetscDrawLineSetWidth"))
"""
function PetscDrawLineSetWidth(petsclib::PetscLibType, draw::PetscDraw, width::Real)
    error("PetscDrawLineSetWidth: no generated method for these argument types")
end

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
	PetscDrawMarker(petsclib::PetscLibType, draw::PetscDraw, xl::PetscReal, yl::PetscReal, cl::Cint) 
draws a marker onto a drawable.

Not Collective

Input Parameters:
- `draw` - the drawing context
- `xl`   - horizontal coordinate of the marker
- `yl`   - vertical coordinate of the marker
- `cl`   - the color of the marker

Level: beginner

See also: `PetscDraw`, `PetscDrawPoint()`, `PetscDrawString()`, `PetscDrawSetMarkerType()`, `PetscDrawGetMarkerType()`

# External Links
$(_doc_external("Draw/PetscDrawMarker"))
"""
function PetscDrawMarker(petsclib::PetscLibType, draw::PetscDraw, xl::Real, yl::Real, cl::Cint)
    error("PetscDrawMarker: no generated method for these argument types")
end

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
	draw::PetscDraw = PetscDrawOpenImage(petsclib::PetscLibType, comm::MPI_Comm, filename::String, w::Cint, h::Cint) 
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

See also: `PetscDraw`, `PETSC_DRAW_IMAGE`, `PETSC_DRAW_X`, `PetscDrawSetSave()`, `PetscDrawSetFromOptions()`, `PetscDrawCreate()`, `PetscDrawDestroy()`

# External Links
$(_doc_external("Draw/PetscDrawOpenImage"))
"""
function PetscDrawOpenImage(petsclib::PetscLibType, comm::MPI_Comm, filename::String, w::Cint, h::Cint)
    error("PetscDrawOpenImage: no generated method for these argument types")
end

@for_petsc function PetscDrawOpenImage(petsclib::$UnionPetscLib, comm::MPI_Comm, filename::String, w::Cint, h::Cint )
	draw_ = Ref{PetscDraw}()

    @chk ccall(
               (:PetscDrawOpenImage, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, Cint, Cint, Ptr{PetscDraw}),
               comm, filename, w, h, draw_,
              )

	draw = draw_[]

	return draw
end 

"""
	win::PetscDraw = PetscDrawOpenNull(petsclib::PetscLibType, comm::MPI_Comm) 
Opens a null drawing context. All draw commands to
it are ignored.

Input Parameter:
- `comm` - MPI communicator

Output Parameter:
- `win` - the drawing context

Level: advanced

See also: `PetscDraw`, `PetscDrawIsNull()`, `PETSC_DRAW_NULL`, `PetscDrawOpenX()`

# External Links
$(_doc_external("Draw/PetscDrawOpenNull"))
"""
function PetscDrawOpenNull(petsclib::PetscLibType, comm::MPI_Comm)
    error("PetscDrawOpenNull: no generated method for these argument types")
end

@for_petsc function PetscDrawOpenNull(petsclib::$UnionPetscLib, comm::MPI_Comm )
	win_ = Ref{PetscDraw}()

    @chk ccall(
               (:PetscDrawOpenNull, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{PetscDraw}),
               comm, win_,
              )

	win = win_[]

	return win
end 

"""
	draw::PetscDraw = PetscDrawOpenX(petsclib::PetscLibType, comm::MPI_Comm, display::String, title::String, x::Cint, y::Cint, w::Cint, h::Cint) 
Opens an X-window for use with the `PetscDraw` routines.

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
- `-display name`           - Sets name of machine for the X display
- `-draw_pause pause`       - Sets time (in seconds) that the program pauses after `PetscDrawPause()` has been called
(0 is default, -1 implies until user input).
- `-draw_cmap name`         - Sets the colormap to use.
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

See also: `PetscDrawFlush()`, `PetscDrawDestroy()`, `PetscDrawCreate()`

# External Links
$(_doc_external("Draw/PetscDrawOpenX"))
"""
function PetscDrawOpenX(petsclib::PetscLibType, comm::MPI_Comm, display::String, title::String, x::Cint, y::Cint, w::Cint, h::Cint)
    error("PetscDrawOpenX: no generated method for these argument types")
end

@for_petsc function PetscDrawOpenX(petsclib::$UnionPetscLib, comm::MPI_Comm, display::String, title::String, x::Cint, y::Cint, w::Cint, h::Cint )
	draw_ = Ref{PetscDraw}()

    @chk ccall(
               (:PetscDrawOpenX, $petsc_library),
               PetscErrorCode,
               (MPI_Comm, Ptr{Cchar}, Ptr{Cchar}, Cint, Cint, Cint, Cint, Ptr{PetscDraw}),
               comm, display, title, x, y, w, h, draw_,
              )

	draw = draw_[]

	return draw
end 

"""
	PetscDrawPause(petsclib::PetscLibType, draw::PetscDraw) 
Waits n seconds or until user input, depending on input
to `PetscDrawSetPause()`.

Collective

Input Parameter:
- `draw` - the drawing context

Level: beginner

See also: `PetscDraw`, `PetscDrawSetPause()`, `PetscDrawGetPause()`

# External Links
$(_doc_external("Draw/PetscDrawPause"))
"""
function PetscDrawPause(petsclib::PetscLibType, draw::PetscDraw)
    error("PetscDrawPause: no generated method for these argument types")
end

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
	x::PetscReal,y::PetscReal = PetscDrawPixelToCoordinate(petsclib::PetscLibType, draw::PetscDraw, i::Cint, j::Cint) 
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

See also: `PetscDraw`

# External Links
$(_doc_external("Draw/PetscDrawPixelToCoordinate"))
"""
function PetscDrawPixelToCoordinate(petsclib::PetscLibType, draw::PetscDraw, i::Cint, j::Cint)
    error("PetscDrawPixelToCoordinate: no generated method for these argument types")
end

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
	PetscDrawPoint(petsclib::PetscLibType, draw::PetscDraw, xl::PetscReal, yl::PetscReal, cl::Cint) 
draws a point onto a drawable.

Not Collective

Input Parameters:
- `draw` - the drawing context
- `xl`   - horizatonal coordinate of the point
- `yl`   - vertical coordinate of the point
- `cl`   - the color of the point

Level: beginner

See also: `PetscDraw`, `PetscDrawPointPixel()`, `PetscDrawPointSetSize()`, `PetscDrawLine()`, `PetscDrawRectangle()`, `PetscDrawTriangle()`, `PetscDrawEllipse()`,
`PetscDrawMarker()`, `PetscDrawString()`, `PetscDrawArrow()`

# External Links
$(_doc_external("Draw/PetscDrawPoint"))
"""
function PetscDrawPoint(petsclib::PetscLibType, draw::PetscDraw, xl::Real, yl::Real, cl::Cint)
    error("PetscDrawPoint: no generated method for these argument types")
end

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
	PetscDrawPointPixel(petsclib::PetscLibType, draw::PetscDraw, x::Cint, y::Cint, c::Cint) 
draws a point onto a drawable, in pixel coordinates

Not Collective

Input Parameters:
- `draw` - the drawing context
- `x`    - horizontal pixel coordinates of the point
- `y`    - vertical pixel coordinates of the point
- `c`    - the color of the point

Level: beginner

See also: `PetscDraw`, `PetscDrawPoint()`, `PetscDrawPointSetSize()`

# External Links
$(_doc_external("Draw/PetscDrawPointPixel"))
"""
function PetscDrawPointPixel(petsclib::PetscLibType, draw::PetscDraw, x::Cint, y::Cint, c::Cint)
    error("PetscDrawPointPixel: no generated method for these argument types")
end

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
	PetscDrawPointSetSize(petsclib::PetscLibType, draw::PetscDraw, width::PetscReal) 
Sets the point size for future draws.  The size is
relative to the user coordinates of the window; 0.0 denotes the natural
width, 1.0 denotes the entire viewport.

Not Collective

Input Parameters:
- `draw`  - the drawing context
- `width` - the width in user coordinates

Level: advanced

See also: `PetscDraw`, `PetscDrawPoint()`, `PetscDrawMarker()`

# External Links
$(_doc_external("Draw/PetscDrawPointSetSize"))
"""
function PetscDrawPointSetSize(petsclib::PetscLibType, draw::PetscDraw, width::Real)
    error("PetscDrawPointSetSize: no generated method for these argument types")
end

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
	PetscDrawPopCurrentPoint(petsclib::PetscLibType, draw::PetscDraw) 
Pops a current draw point (discarding it)

Not Collective

Input Parameter:
- `draw` - the drawing context

Level: intermediate

See also: `PetscDraw`, `PetscDrawPushCurrentPoint()`, `PetscDrawSetCurrentPoint()`, `PetscDrawGetCurrentPoint()`

# External Links
$(_doc_external("Draw/PetscDrawPopCurrentPoint"))
"""
function PetscDrawPopCurrentPoint(petsclib::PetscLibType, draw::PetscDraw)
    error("PetscDrawPopCurrentPoint: no generated method for these argument types")
end

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
	PetscDrawPushCurrentPoint(petsclib::PetscLibType, draw::PetscDraw, x::PetscReal, y::PetscReal) 
Pushes a new current draw point, retaining the old one, some codes use this point to determine where to draw next

Not Collective

Input Parameters:
- `draw` - the drawing context
- `x`    - horizontal coordinate of the current point
- `y`    - vertical coordinate of the current point

Level: intermediate

See also: `PetscDraw`, `PetscDrawPopCurrentPoint()`, `PetscDrawGetCurrentPoint()`

# External Links
$(_doc_external("Draw/PetscDrawPushCurrentPoint"))
"""
function PetscDrawPushCurrentPoint(petsclib::PetscLibType, draw::PetscDraw, x::Real, y::Real)
    error("PetscDrawPushCurrentPoint: no generated method for these argument types")
end

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
	PetscDrawRectangle(petsclib::PetscLibType, draw::PetscDraw, xl::PetscReal, yl::PetscReal, xr::PetscReal, yr::PetscReal, c1::Cint, c2::Cint, c3::Cint, c4::Cint) 
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

See also: `PetscDraw`, `PetscDrawLine()`, `PetscDrawTriangle()`, `PetscDrawEllipse()`,
`PetscDrawMarker()`, `PetscDrawPoint()`, `PetscDrawString()`, `PetscDrawArrow()`

# External Links
$(_doc_external("Draw/PetscDrawRectangle"))
"""
function PetscDrawRectangle(petsclib::PetscLibType, draw::PetscDraw, xl::Real, yl::Real, xr::Real, yr::Real, c1::Cint, c2::Cint, c3::Cint, c4::Cint)
    error("PetscDrawRectangle: no generated method for these argument types")
end

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
	PetscDrawRegister(petsclib::PetscLibType, sname::String, fnc::external) 
Adds a method to the graphics package.

Not Collective, No Fortran Support

Input Parameters:
- `sname`    - name of a new user-defined graphics class
- `function` - routine to create method context

Level: developer

See also: `PetscDraw`, `PetscDrawRegisterAll()`, `PetscDrawRegisterDestroy()`, `PetscDrawType`, `PetscDrawSetType()`

# External Links
$(_doc_external("Draw/PetscDrawRegister"))
"""
function PetscDrawRegister(petsclib::PetscLibType, sname::String, fnc::external)
    error("PetscDrawRegister: no generated method for these argument types")
end

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
	PetscDrawResizeWindow(petsclib::PetscLibType, draw::PetscDraw, w::Cint, h::Cint) 
Allows one to resize a window from a program.

Collective

Input Parameters:
- `draw` - the window
- `w`    - the new width of the window
- `h`    - the new height of the window

Level: intermediate

See also: `PetscDraw`, `PetscDrawCheckResizedWindow()`

# External Links
$(_doc_external("Draw/PetscDrawResizeWindow"))
"""
function PetscDrawResizeWindow(petsclib::PetscLibType, draw::PetscDraw, w::Cint, h::Cint)
    error("PetscDrawResizeWindow: no generated method for these argument types")
end

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
	PetscDrawRestoreSingleton(petsclib::PetscLibType, draw::PetscDraw, sdraw::PetscDraw) 
Remove access to a `PetscDraw` object obtained with `PetscDrawGetSingleton()`
by the one process.

Collective

Input Parameters:
- `draw`  - the original window
- `sdraw` - the singleton window

Level: advanced

See also: `PetscDraw`, `PetscDrawGetSingleton()`, `PetscViewerGetSingleton()`, `PetscViewerRestoreSingleton()`

# External Links
$(_doc_external("Draw/PetscDrawRestoreSingleton"))
"""
function PetscDrawRestoreSingleton(petsclib::PetscLibType, draw::PetscDraw, sdraw::PetscDraw)
    error("PetscDrawRestoreSingleton: no generated method for these argument types")
end

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
	x::PetscReal,y::PetscReal = PetscDrawSPAddPoint(petsclib::PetscLibType, sp::PetscDrawSP) 
Adds another point to each of the scatter plot point curves.

Not Collective

Input Parameters:
- `sp` - the scatter plot data structure
- `x`  - the x coordinate values (of length dim) for the points of the curve
- `y`  - the y coordinate values (of length dim) for the points of the curve

Level: intermediate

See also: `PetscDrawSPAddPoints()`, `PetscDrawSP`, `PetscDrawSPCreate()`, `PetscDrawSPReset()`, `PetscDrawSPDraw()`, `PetscDrawSPAddPointColorized()`

# External Links
$(_doc_external("Draw/PetscDrawSPAddPoint"))
"""
function PetscDrawSPAddPoint(petsclib::PetscLibType, sp::PetscDrawSP)
    error("PetscDrawSPAddPoint: no generated method for these argument types")
end

@for_petsc function PetscDrawSPAddPoint(petsclib::$UnionPetscLib, sp::PetscDrawSP )
	x_ = Ref{$PetscReal}()
	y_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscDrawSPAddPoint, $petsc_library),
               PetscErrorCode,
               (PetscDrawSP, Ptr{$PetscReal}, Ptr{$PetscReal}),
               sp, x_, y_,
              )

	x = x_[]
	y = y_[]

	return x,y
end 

"""
	x::PetscReal,y::PetscReal,z::PetscReal = PetscDrawSPAddPointColorized(petsclib::PetscLibType, sp::PetscDrawSP) 
Adds another point to each of the scatter plots as well as a numeric value to be used to colorize the scatter point.

Not Collective

Input Parameters:
- `sp` - the scatter plot data structure
- `x`  - array of length dim containing the new x coordinate values for each of the point curves.
- `y`  - array of length dim containing the new y coordinate values for each of the point curves.
- `z`  - array of length dim containing the numeric values that will be mapped to [0,255] and used for scatter point colors.

Level: intermediate

See also: `PetscDrawSPAddPoints()`, `PetscDrawSP`, `PetscDrawSPCreate()`, `PetscDrawSPReset()`, `PetscDrawSPDraw()`, `PetscDrawSPAddPoint()`

# External Links
$(_doc_external("Draw/PetscDrawSPAddPointColorized"))
"""
function PetscDrawSPAddPointColorized(petsclib::PetscLibType, sp::PetscDrawSP)
    error("PetscDrawSPAddPointColorized: no generated method for these argument types")
end

@for_petsc function PetscDrawSPAddPointColorized(petsclib::$UnionPetscLib, sp::PetscDrawSP )
	x_ = Ref{$PetscReal}()
	y_ = Ref{$PetscReal}()
	z_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscDrawSPAddPointColorized, $petsc_library),
               PetscErrorCode,
               (PetscDrawSP, Ptr{$PetscReal}, Ptr{$PetscReal}, Ptr{$PetscReal}),
               sp, x_, y_, z_,
              )

	x = x_[]
	y = y_[]
	z = z_[]

	return x,y,z
end 

"""
	PetscDrawSPAddPoints(petsclib::PetscLibType, sp::PetscDrawSP, n::Cint, xx::Union{Ptr, AbstractArray{PetscReal}}, yy::Union{Ptr, AbstractArray{PetscReal}}) 
Adds several points to each of the scatter plot point curves.

Not Collective

Input Parameters:
- `sp` - the scatter plot context
- `xx` - array of pointers that point to arrays containing the new x coordinates for each curve.
- `yy` - array of pointers that point to arrays containing the new y points for each curve.
- `n`  - number of points being added, each represents a subarray of length dim where dim is the value from `PetscDrawSPGetDimension()`

Level: intermediate

See also: `PetscDrawSPAddPoint()`, `PetscDrawSP`, `PetscDrawSPCreate()`, `PetscDrawSPReset()`, `PetscDrawSPDraw()`, `PetscDrawSPAddPointColorized()`

# External Links
$(_doc_external("Draw/PetscDrawSPAddPoints"))
"""
function PetscDrawSPAddPoints(petsclib::PetscLibType, sp::PetscDrawSP, n::Cint, xx::Union{Ptr, AbstractArray{<:Number}}, yy::Union{Ptr, AbstractArray{<:Number}})
    error("PetscDrawSPAddPoints: no generated method for these argument types")
end

@for_petsc function PetscDrawSPAddPoints(petsclib::$UnionPetscLib, sp::PetscDrawSP, n::Cint, xx::Union{Ptr, AbstractArray{$PetscReal}}, yy::Union{Ptr, AbstractArray{$PetscReal}} )
	xx_ = Ref{Ptr{$PetscReal}}(xx isa Ptr ? xx : pointer(xx))
	yy_ = Ref{Ptr{$PetscReal}}(yy isa Ptr ? yy : pointer(yy))

    @chk ccall(
               (:PetscDrawSPAddPoints, $petsc_library),
               PetscErrorCode,
               (PetscDrawSP, Cint, Ptr{Ptr{$PetscReal}}, Ptr{Ptr{$PetscReal}}),
               sp, n, xx_, yy_,
              )


	return nothing
end 

"""
	drawsp::PetscDrawSP = PetscDrawSPCreate(petsclib::PetscLibType, draw::PetscDraw, dim::Cint) 
Creates a scatter plot data structure.

Collective

Input Parameters:
- `draw` - the window where the graph will be made.
- `dim`  - the number of sets of points which will be drawn

Output Parameter:
- `drawsp` - the scatter plot context

Level: intermediate

See also: `PetscDrawLGCreate()`, `PetscDrawLG`, `PetscDrawBarCreate()`, `PetscDrawBar`, `PetscDrawHGCreate()`, `PetscDrawHG`, `PetscDrawSPDestroy()`, `PetscDraw`, `PetscDrawSP`, `PetscDrawSPSetDimension()`, `PetscDrawSPReset()`,
`PetscDrawSPAddPoint()`, `PetscDrawSPAddPoints()`, `PetscDrawSPDraw()`, `PetscDrawSPSave()`, `PetscDrawSPSetLimits()`, `PetscDrawSPGetAxis()`, `PetscDrawAxis`, `PetscDrawSPGetDraw()`

# External Links
$(_doc_external("Draw/PetscDrawSPCreate"))
"""
function PetscDrawSPCreate(petsclib::PetscLibType, draw::PetscDraw, dim::Cint)
    error("PetscDrawSPCreate: no generated method for these argument types")
end

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
	PetscDrawSPDestroy(petsclib::PetscLibType, sp::Union{PetscDrawSP, Ref{PetscDrawSP}}) 
Frees all space taken up by scatter plot data structure.

Collective

Input Parameter:
- `sp` - the scatter plot context

Level: intermediate

See also: `PetscDrawSPCreate()`, `PetscDrawSP`, `PetscDrawSPReset()`

# External Links
$(_doc_external("Draw/PetscDrawSPDestroy"))
"""
function PetscDrawSPDestroy(petsclib::PetscLibType, sp::Union{PetscDrawSP, Ref{PetscDrawSP}})
    error("PetscDrawSPDestroy: no generated method for these argument types")
end

@for_petsc function PetscDrawSPDestroy(petsclib::$UnionPetscLib, sp::Union{PetscDrawSP, Ref{PetscDrawSP}} )
	sp_ = sp isa Base.RefValue ? sp : Ref{PetscDrawSP}(sp)

    @chk ccall(
               (:PetscDrawSPDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscDrawSP},),
               sp_,
              )


	return nothing
end 

"""
	PetscDrawSPDraw(petsclib::PetscLibType, sp::PetscDrawSP, clear::PetscBool) 
Redraws a scatter plot.

Collective

Input Parameters:
- `sp`    - the scatter plot context
- `clear` - clear the window before drawing the new plot

Level: intermediate

See also: `PetscDrawLGDraw()`, `PetscDrawLGSPDraw()`, `PetscDrawSP`, `PetscDrawSPCreate()`, `PetscDrawSPReset()`, `PetscDrawSPAddPoint()`, `PetscDrawSPAddPoints()`

# External Links
$(_doc_external("Draw/PetscDrawSPDraw"))
"""
function PetscDrawSPDraw(petsclib::PetscLibType, sp::PetscDrawSP, clear::PetscBool)
    error("PetscDrawSPDraw: no generated method for these argument types")
end

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
	axis::PetscDrawAxis = PetscDrawSPGetAxis(petsclib::PetscLibType, sp::PetscDrawSP) 
Gets the axis context associated with a scatter plot

Not Collective

Input Parameter:
- `sp` - the scatter plot context

Output Parameter:
- `axis` - the axis context

Level: intermediate

See also: `PetscDrawSP`, `PetscDrawSPCreate()`, `PetscDrawSPDraw()`, `PetscDrawSPAddPoint()`, `PetscDrawSPAddPoints()`, `PetscDrawAxis`, `PetscDrawAxisCreate()`

# External Links
$(_doc_external("Draw/PetscDrawSPGetAxis"))
"""
function PetscDrawSPGetAxis(petsclib::PetscLibType, sp::PetscDrawSP)
    error("PetscDrawSPGetAxis: no generated method for these argument types")
end

@for_petsc function PetscDrawSPGetAxis(petsclib::$UnionPetscLib, sp::PetscDrawSP )
	axis_ = Ref{PetscDrawAxis}()

    @chk ccall(
               (:PetscDrawSPGetAxis, $petsc_library),
               PetscErrorCode,
               (PetscDrawSP, Ptr{PetscDrawAxis}),
               sp, axis_,
              )

	axis = axis_[]

	return axis
end 

"""
	dim::Cint = PetscDrawSPGetDimension(petsclib::PetscLibType, sp::PetscDrawSP) 
Get the number of sets of points that are to be drawn at each `PetscDrawSPAddPoint()`

Not Collective

Input Parameter:
- `sp` - the scatter plot context.

Output Parameter:
- `dim` - the number of point curves on this process

Level: intermediate

See also: `PetscDrawSP`, `PetscDrawSPCreate()`, `PetscDrawSPAddPoint()`, `PetscDrawSPAddPoints()`

# External Links
$(_doc_external("Draw/PetscDrawSPGetDimension"))
"""
function PetscDrawSPGetDimension(petsclib::PetscLibType, sp::PetscDrawSP)
    error("PetscDrawSPGetDimension: no generated method for these argument types")
end

@for_petsc function PetscDrawSPGetDimension(petsclib::$UnionPetscLib, sp::PetscDrawSP )
	dim_ = Ref{Cint}()

    @chk ccall(
               (:PetscDrawSPGetDimension, $petsc_library),
               PetscErrorCode,
               (PetscDrawSP, Ptr{Cint}),
               sp, dim_,
              )

	dim = dim_[]

	return dim
end 

"""
	draw::PetscDraw = PetscDrawSPGetDraw(petsclib::PetscLibType, sp::PetscDrawSP) 
Gets the draw context associated with a scatter plot

Not Collective

Input Parameter:
- `sp` - the scatter plot context

Output Parameter:
- `draw` - the draw context

Level: intermediate

See also: `PetscDrawSP`, `PetscDrawSPCreate()`, `PetscDrawSPDraw()`, `PetscDraw`

# External Links
$(_doc_external("Draw/PetscDrawSPGetDraw"))
"""
function PetscDrawSPGetDraw(petsclib::PetscLibType, sp::PetscDrawSP)
    error("PetscDrawSPGetDraw: no generated method for these argument types")
end

@for_petsc function PetscDrawSPGetDraw(petsclib::$UnionPetscLib, sp::PetscDrawSP )
	draw_ = Ref{PetscDraw}()

    @chk ccall(
               (:PetscDrawSPGetDraw, $petsc_library),
               PetscErrorCode,
               (PetscDrawSP, Ptr{PetscDraw}),
               sp, draw_,
              )

	draw = draw_[]

	return draw
end 

"""
	PetscDrawSPReset(petsclib::PetscLibType, sp::PetscDrawSP) 
Clears scatter plot to allow for reuse with new data.

Not Collective

Input Parameter:
- `sp` - the scatter plot context.

Level: intermediate

See also: `PetscDrawSP`, `PetscDrawSPCreate()`, `PetscDrawSPAddPoint()`, `PetscDrawSPAddPoints()`, `PetscDrawSPDraw()`

# External Links
$(_doc_external("Draw/PetscDrawSPReset"))
"""
function PetscDrawSPReset(petsclib::PetscLibType, sp::PetscDrawSP)
    error("PetscDrawSPReset: no generated method for these argument types")
end

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
	PetscDrawSPSave(petsclib::PetscLibType, sp::PetscDrawSP) 
Saves a drawn image

Collective

Input Parameter:
- `sp` - the scatter plot context

Level: intermediate

See also: `PetscDrawSPCreate()`, `PetscDrawSPGetDraw()`, `PetscDrawSetSave()`, `PetscDrawSave()`

# External Links
$(_doc_external("Draw/PetscDrawSPSave"))
"""
function PetscDrawSPSave(petsclib::PetscLibType, sp::PetscDrawSP)
    error("PetscDrawSPSave: no generated method for these argument types")
end

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
	PetscDrawSPSetDimension(petsclib::PetscLibType, sp::PetscDrawSP, dim::Cint) 
Change the number of points that are added at each  `PetscDrawSPAddPoint()`

Not Collective

Input Parameters:
- `sp`  - the scatter plot context.
- `dim` - the number of point curves on this process

Level: intermediate

See also: `PetscDrawSP`, `PetscDrawSPCreate()`, `PetscDrawSPAddPoint()`, `PetscDrawSPAddPoints()`

# External Links
$(_doc_external("Draw/PetscDrawSPSetDimension"))
"""
function PetscDrawSPSetDimension(petsclib::PetscLibType, sp::PetscDrawSP, dim::Cint)
    error("PetscDrawSPSetDimension: no generated method for these argument types")
end

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
	PetscDrawSPSetLimits(petsclib::PetscLibType, sp::PetscDrawSP, x_min::PetscReal, x_max::PetscReal, y_min::PetscReal, y_max::PetscReal) 
Sets the axis limits for a scatter plot. If more points are added after this call, the limits will be adjusted to include those additional points.

Not Collective

Input Parameters:
- `sp`    - the line graph context
- `x_min` - the horizontal lower limit
- `x_max` - the horizontal upper limit
- `y_min` - the vertical lower limit
- `y_max` - the vertical upper limit

Level: intermediate

See also: `PetscDrawSP`, `PetscDrawAxis`, `PetscDrawSPCreate()`, `PetscDrawSPDraw()`, `PetscDrawSPAddPoint()`, `PetscDrawSPAddPoints()`, `PetscDrawSPGetAxis()`

# External Links
$(_doc_external("Draw/PetscDrawSPSetLimits"))
"""
function PetscDrawSPSetLimits(petsclib::PetscLibType, sp::PetscDrawSP, x_min::Real, x_max::Real, y_min::Real, y_max::Real)
    error("PetscDrawSPSetLimits: no generated method for these argument types")
end

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
	PetscDrawSave(petsclib::PetscLibType, draw::PetscDraw) 
Saves a drawn image

Collective

Input Parameter:
- `draw` - the drawing context

Level: advanced

See also: `PetscDraw`, `PetscDrawSetSave()`

# External Links
$(_doc_external("Draw/PetscDrawSave"))
"""
function PetscDrawSave(petsclib::PetscLibType, draw::PetscDraw)
    error("PetscDrawSave: no generated method for these argument types")
end

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
	PetscDrawSaveMovie(petsclib::PetscLibType, draw::PetscDraw) 
Saves a movie from previously saved images

Collective

Input Parameter:
- `draw` - the drawing context

Level: advanced

See also: `PetscDraw`, `PetscDrawSetSave()`, `PetscDrawSetSaveMovie()`

# External Links
$(_doc_external("Draw/PetscDrawSaveMovie"))
"""
function PetscDrawSaveMovie(petsclib::PetscLibType, draw::PetscDraw)
    error("PetscDrawSaveMovie: no generated method for these argument types")
end

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
	PetscDrawScalePopup(petsclib::PetscLibType, popup::PetscDraw, min::PetscReal, max::PetscReal) 
draws a contour scale window.

Collective

Input Parameters:
- `popup` - the window (often a window obtained via `PetscDrawGetPopup()`
- `min`   - minimum value being plotted
- `max`   - maximum value being plotted

Level: intermediate

See also: `PetscDraw`, `PetscDrawGetPopup()`, `PetscDrawTensorContour()`

# External Links
$(_doc_external("Draw/PetscDrawScalePopup"))
"""
function PetscDrawScalePopup(petsclib::PetscLibType, popup::PetscDraw, min::Real, max::Real)
    error("PetscDrawScalePopup: no generated method for these argument types")
end

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
	PetscDrawSetCoordinates(petsclib::PetscLibType, draw::PetscDraw, xl::PetscReal, yl::PetscReal, xr::PetscReal, yr::PetscReal) 
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

See also: `PetscDraw`, `PetscDrawGetCoordinates()`

# External Links
$(_doc_external("Draw/PetscDrawSetCoordinates"))
"""
function PetscDrawSetCoordinates(petsclib::PetscLibType, draw::PetscDraw, xl::Real, yl::Real, xr::Real, yr::Real)
    error("PetscDrawSetCoordinates: no generated method for these argument types")
end

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
	PetscDrawSetCurrentPoint(petsclib::PetscLibType, draw::PetscDraw, x::PetscReal, y::PetscReal) 
Sets the current draw point, some codes use this point to determine where to draw next

Not Collective

Input Parameters:
- `draw` - the drawing context
- `x`    - horizontal coordinate of the current point
- `y`    - vertical coordinate of the current point

Level: intermediate

See also: `PetscDraw`, `PetscDrawPushCurrentPoint()`, `PetscDrawPopCurrentPoint()`, `PetscDrawGetCurrentPoint()`

# External Links
$(_doc_external("Draw/PetscDrawSetCurrentPoint"))
"""
function PetscDrawSetCurrentPoint(petsclib::PetscLibType, draw::PetscDraw, x::Real, y::Real)
    error("PetscDrawSetCurrentPoint: no generated method for these argument types")
end

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
	PetscDrawSetDisplay(petsclib::PetscLibType, draw::PetscDraw, display::String) 
Sets the display where a `PetscDraw` object will be displayed

Input Parameters:
- `draw`    - the drawing context
- `display` - the X windows display

Level: advanced

See also: `PetscDraw`, `PetscDrawOpenX()`, `PetscDrawCreate()`

# External Links
$(_doc_external("Draw/PetscDrawSetDisplay"))
"""
function PetscDrawSetDisplay(petsclib::PetscLibType, draw::PetscDraw, display::String)
    error("PetscDrawSetDisplay: no generated method for these argument types")
end

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
	PetscDrawSetDoubleBuffer(petsclib::PetscLibType, draw::PetscDraw) 
Sets a window to be double buffered.

Logically Collective

Input Parameter:
- `draw` - the drawing context

Level: intermediate

See also: `PetscDraw`, `PetscDrawOpenX()`, `PetscDrawCreate()`

# External Links
$(_doc_external("Draw/PetscDrawSetDoubleBuffer"))
"""
function PetscDrawSetDoubleBuffer(petsclib::PetscLibType, draw::PetscDraw)
    error("PetscDrawSetDoubleBuffer: no generated method for these argument types")
end

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
	PetscDrawSetFromOptions(petsclib::PetscLibType, draw::PetscDraw) 
Sets the graphics type from the options database.
Defaults to a PETSc X Windows graphics.

Collective

Input Parameter:
- `draw` - the graphics context

Options Database Keys:
- `-nox`                              - do not use X graphics (ignore graphics calls, but run program correctly)
- `-nox_warning`                      - when X Windows support is not installed this prevents the warning message from being printed
- `-draw_pause seconds`               - -1 indicates wait for mouse input, -2 indicates pause when window is to be destroyed
- `-draw_marker_type (x|point)`       - set the marker type
- `-draw_save [filename]`             - (X Windows only) saves each image before it is cleared to a file
- `-draw_save_final_image [filename]` - (X Windows only) saves the final image displayed in a window
- `-draw_save_movie`                  - converts image files to a movie  at the end of the run. See `PetscDrawSetSave()`
- `-draw_save_single_file`            - saves each new image in the same file, normally each new image is saved in a new file with 'filename/filename_%d.ext'
- `-draw_save_on_clear`               - saves an image on each clear, mainly for debugging
- `-draw_save_on_flush`               - saves an image on each flush, mainly for debugging

Level: intermediate

See also: `PetscDraw`, `PetscDrawCreate()`, `PetscDrawSetType()`, `PetscDrawSetSave()`, `PetscDrawSetSaveFinalImage()`, `PetscDrawPause()`, `PetscDrawSetPause()`

# External Links
$(_doc_external("Draw/PetscDrawSetFromOptions"))
"""
function PetscDrawSetFromOptions(petsclib::PetscLibType, draw::PetscDraw)
    error("PetscDrawSetFromOptions: no generated method for these argument types")
end

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
	PetscDrawSetMarkerType(petsclib::PetscLibType, draw::PetscDraw, mtype::PetscDrawMarkerType) 
sets the type of marker to display with `PetscDrawMarker()`

Not Collective

Input Parameters:
- `draw`  - the drawing context
- `mtype` - either `PETSC_DRAW_MARKER_CROSS` (default) or `PETSC_DRAW_MARKER_POINT`

Options Database Key:
- `-draw_marker_type` - x or point

Level: beginner

See also: `PetscDraw`, `PetscDrawPoint()`, `PetscDrawMarker()`, `PetscDrawGetMarkerType()`, `PetscDrawMarkerType`

# External Links
$(_doc_external("Draw/PetscDrawSetMarkerType"))
"""
function PetscDrawSetMarkerType(petsclib::PetscLibType, draw::PetscDraw, mtype::PetscDrawMarkerType)
    error("PetscDrawSetMarkerType: no generated method for these argument types")
end

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
	PetscDrawSetOptionsPrefix(petsclib::PetscLibType, draw::PetscDraw, prefix::String) 
Sets the prefix used for searching for all
`PetscDraw` options in the database.

Logically Collective

Input Parameters:
- `draw`   - the draw context
- `prefix` - the prefix to prepend to all option names

Level: advanced

See also: `PetscDraw`, `PetscDrawSetFromOptions()`, `PetscDrawCreate()`

# External Links
$(_doc_external("Draw/PetscDrawSetOptionsPrefix"))
"""
function PetscDrawSetOptionsPrefix(petsclib::PetscLibType, draw::PetscDraw, prefix::String)
    error("PetscDrawSetOptionsPrefix: no generated method for these argument types")
end

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
	PetscDrawSetPause(petsclib::PetscLibType, draw::PetscDraw, lpause::PetscReal) 
Sets the amount of time that program pauses after
a `PetscDrawPause()` is called.

Logically Collective

Input Parameters:
- `draw`   - the drawing object
- `lpause` - number of seconds to pause, -1 implies until user input, -2 pauses only on the `PetscDrawDestroy()`

Options Database Key:
- `-draw_pause value` - set the time to pause

Level: intermediate

See also: `PetscDraw`, `PetscDrawGetPause()`, `PetscDrawPause()`

# External Links
$(_doc_external("Draw/PetscDrawSetPause"))
"""
function PetscDrawSetPause(petsclib::PetscLibType, draw::PetscDraw, lpause::Real)
    error("PetscDrawSetPause: no generated method for these argument types")
end

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
	PetscDrawSetSave(petsclib::PetscLibType, draw::PetscDraw, filename::String) 
Saves images produced in a `PetscDraw` into a file

Collective

Input Parameters:
- `draw`     - the graphics context
- `filename` - name of the file, if .ext then uses name of draw object plus .ext using .ext to determine the image type

Options Database Keys:
- `-draw_save filename filename`      - `filename` could be `name.ext` or `.ext` (where .ext determines the type of graphics file to save, for example .png)
- `-draw_save_final_image [filename]` - saves the final image displayed in a window
- `-draw_save_single_file`            - saves each new image in the same file, normally each new image is saved in a new file with filename/filename_%d.ext

Level: intermediate

See also: `PetscDraw`, `PetscDrawOpenX()`, `PetscDrawOpenImage()`, `PetscDrawSetFromOptions()`, `PetscDrawCreate()`, `PetscDrawDestroy()`, `PetscDrawSetSaveFinalImage()`

# External Links
$(_doc_external("Draw/PetscDrawSetSave"))
"""
function PetscDrawSetSave(petsclib::PetscLibType, draw::PetscDraw, filename::String)
    error("PetscDrawSetSave: no generated method for these argument types")
end

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
	PetscDrawSetSaveFinalImage(petsclib::PetscLibType, draw::PetscDraw, filename::String) 
Saves the final image produced in a `PetscDraw` into a file

Collective

Input Parameters:
- `draw`     - the graphics context
- `filename` - name of the file, if NULL or empty uses name set with `PetscDrawSetSave()` or the name of the draw object

Options Database Key:
- `-draw_save_final_image filename` - filename could be name.ext or .ext (where .ext determines the type of graphics file to save, for example .png)

Level: intermediate

See also: `PetscDraw`, `PetscDrawSetSave()`, `PetscDrawSetFromOptions()`, `PetscDrawCreate()`, `PetscDrawDestroy()`

# External Links
$(_doc_external("Draw/PetscDrawSetSaveFinalImage"))
"""
function PetscDrawSetSaveFinalImage(petsclib::PetscLibType, draw::PetscDraw, filename::String)
    error("PetscDrawSetSaveFinalImage: no generated method for these argument types")
end

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
	PetscDrawSetSaveMovie(petsclib::PetscLibType, draw::PetscDraw, movieext::String) 
Saves a movie produced from a `PetscDraw` into a file

Collective

Input Parameters:
- `draw`     - the graphics context
- `movieext` - optional extension defining the movie format

Options Database Key:
- `-draw_save_movie .ext` - saves a movie with extension .ext

Level: intermediate

See also: `PetscDraw`, `PetscDrawSetSave()`, `PetscDrawSetFromOptions()`, `PetscDrawCreate()`, `PetscDrawDestroy()`

# External Links
$(_doc_external("Draw/PetscDrawSetSaveMovie"))
"""
function PetscDrawSetSaveMovie(petsclib::PetscLibType, draw::PetscDraw, movieext::String)
    error("PetscDrawSetSaveMovie: no generated method for these argument types")
end

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
	PetscDrawSetTitle(petsclib::PetscLibType, draw::PetscDraw, title::String) 
Sets the title of a `PetscDraw` context.

Collective

Input Parameters:
- `draw`  - the graphics context
- `title` - the title

Level: intermediate

See also: `PetscDraw`, `PetscDrawGetTitle()`, `PetscDrawAppendTitle()`

# External Links
$(_doc_external("Draw/PetscDrawSetTitle"))
"""
function PetscDrawSetTitle(petsclib::PetscLibType, draw::PetscDraw, title::String)
    error("PetscDrawSetTitle: no generated method for these argument types")
end

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
	PetscDrawSetType(petsclib::PetscLibType, draw::PetscDraw, type::PetscDrawType) 
Builds graphics object for a particular implementation

Collective

Input Parameters:
- `draw` - the graphics context
- `type` - for example, `PETSC_DRAW_X`

Options Database Key:
- `-draw_type (x|null|win32|tikz|image)` - Sets the type; see `PetscDrawType`

Level: intermediate

See also: `PetscDraw`, `PETSC_DRAW_X`, `PETSC_DRAW_TIKZ`, `PETSC_DRAW_IMAGE`, `PetscDrawSetFromOptions()`, `PetscDrawCreate()`, `PetscDrawDestroy()`, `PetscDrawType`

# External Links
$(_doc_external("Draw/PetscDrawSetType"))
"""
function PetscDrawSetType(petsclib::PetscLibType, draw::PetscDraw, type::PetscDrawType)
    error("PetscDrawSetType: no generated method for these argument types")
end

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
	PetscDrawSetViewPort(petsclib::PetscLibType, draw::PetscDraw, xl::PetscReal, yl::PetscReal, xr::PetscReal, yr::PetscReal) 
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

See also: `PetscDrawGetViewPort()`, `PetscDraw`, `PetscDrawSplitViewPort()`, `PetscDrawViewPortsCreate()`

# External Links
$(_doc_external("Draw/PetscDrawSetViewPort"))
"""
function PetscDrawSetViewPort(petsclib::PetscLibType, draw::PetscDraw, xl::Real, yl::Real, xr::Real, yr::Real)
    error("PetscDrawSetViewPort: no generated method for these argument types")
end

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
	PetscDrawSetVisible(petsclib::PetscLibType, draw::PetscDraw, visible::PetscBool) 
Sets if the drawing surface (the 'window') is visible on its display.

Input Parameters:
- `draw`    - the drawing window
- `visible` - if the surface should be visible

Level: intermediate

See also: `PetscDraw`

# External Links
$(_doc_external("Draw/PetscDrawSetVisible"))
"""
function PetscDrawSetVisible(petsclib::PetscLibType, draw::PetscDraw, visible::PetscBool)
    error("PetscDrawSetVisible: no generated method for these argument types")
end

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
	PetscDrawSplitViewPort(petsclib::PetscLibType, draw::PetscDraw) 
Splits a window shared by several processes into smaller
view ports. One for each process.

Collective

Input Parameter:
- `draw` - the drawing context

Level: advanced

See also: `PetscDrawDivideViewPort()`, `PetscDrawSetViewPort()`

# External Links
$(_doc_external("Draw/PetscDrawSplitViewPort"))
"""
function PetscDrawSplitViewPort(petsclib::PetscLibType, draw::PetscDraw)
    error("PetscDrawSplitViewPort: no generated method for these argument types")
end

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
	PetscDrawString(petsclib::PetscLibType, draw::PetscDraw, xl::PetscReal, yl::PetscReal, cl::Cint, text::String) 
draws text onto a drawable.

Not Collective

Input Parameters:
- `draw` - the drawing context
- `xl`   - coordinate of lower left corner of text
- `yl`   - coordinate of lower left corner of text
- `cl`   - the color of the text
- `text` - the text to draw

Level: beginner

See also: `PetscDraw`, `PetscDrawStringVertical()`, `PetscDrawStringCentered()`, `PetscDrawStringBoxed()`, `PetscDrawStringSetSize()`,
`PetscDrawStringGetSize()`, `PetscDrawLine()`, `PetscDrawRectangle()`, `PetscDrawTriangle()`, `PetscDrawEllipse()`,
`PetscDrawMarker()`, `PetscDrawPoint()`

# External Links
$(_doc_external("Draw/PetscDrawString"))
"""
function PetscDrawString(petsclib::PetscLibType, draw::PetscDraw, xl::Real, yl::Real, cl::Cint, text::String)
    error("PetscDrawString: no generated method for these argument types")
end

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
	w::PetscReal,h::PetscReal = PetscDrawStringBoxed(petsclib::PetscLibType, draw::PetscDraw, sxl::PetscReal, syl::PetscReal, sc::Cint, bc::Cint, text::String) 
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

See also: `PetscDraw`, `PetscDrawStringVertical()`, `PetscDrawString()`, `PetscDrawStringCentered()`, `PetscDrawStringSetSize()`,
`PetscDrawStringGetSize()`

# External Links
$(_doc_external("Draw/PetscDrawStringBoxed"))
"""
function PetscDrawStringBoxed(petsclib::PetscLibType, draw::PetscDraw, sxl::Real, syl::Real, sc::Cint, bc::Cint, text::String)
    error("PetscDrawStringBoxed: no generated method for these argument types")
end

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
	PetscDrawStringCentered(petsclib::PetscLibType, draw::PetscDraw, xc::PetscReal, yl::PetscReal, cl::Cint, text::String) 
draws text onto a drawable centered at a point

Not Collective

Input Parameters:
- `draw` - the drawing context
- `xc`   - the coordinates of right-left center of text
- `yl`   - the coordinates of lower edge of text
- `cl`   - the color of the text
- `text` - the text to draw

Level: beginner

See also: `PetscDraw`, `PetscDrawStringVertical()`, `PetscDrawString()`, `PetscDrawStringBoxed()`, `PetscDrawStringSetSize()`,
`PetscDrawStringGetSize()`

# External Links
$(_doc_external("Draw/PetscDrawStringCentered"))
"""
function PetscDrawStringCentered(petsclib::PetscLibType, draw::PetscDraw, xc::Real, yl::Real, cl::Cint, text::String)
    error("PetscDrawStringCentered: no generated method for these argument types")
end

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
	width::PetscReal,height::PetscReal = PetscDrawStringGetSize(petsclib::PetscLibType, draw::PetscDraw) 
Gets the size for character text.  The width is
relative to the user coordinates of the window.

Not Collective

Input Parameters:
- `draw`   - the drawing context
- `width`  - the width in user coordinates
- `height` - the character height

Level: advanced

See also: `PetscDraw`, `PetscDrawStringVertical()`, `PetscDrawString()`, `PetscDrawStringCentered()`, `PetscDrawStringBoxed()`,
`PetscDrawStringSetSize()`

# External Links
$(_doc_external("Draw/PetscDrawStringGetSize"))
"""
function PetscDrawStringGetSize(petsclib::PetscLibType, draw::PetscDraw)
    error("PetscDrawStringGetSize: no generated method for these argument types")
end

@for_petsc function PetscDrawStringGetSize(petsclib::$UnionPetscLib, draw::PetscDraw )
	width_ = Ref{$PetscReal}()
	height_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscDrawStringGetSize, $petsc_library),
               PetscErrorCode,
               (PetscDraw, Ptr{$PetscReal}, Ptr{$PetscReal}),
               draw, width_, height_,
              )

	width = width_[]
	height = height_[]

	return width,height
end 

"""
	PetscDrawStringSetSize(petsclib::PetscLibType, draw::PetscDraw, width::PetscReal, height::PetscReal) 
Sets the size for character text.

Not Collective

Input Parameters:
- `draw`   - the drawing context
- `width`  - the width in user coordinates
- `height` - the character height in user coordinates

Level: advanced

See also: `PetscDraw`, `PetscDrawStringVertical()`, `PetscDrawString()`, `PetscDrawStringCentered()`, `PetscDrawStringBoxed()`,
`PetscDrawStringGetSize()`

# External Links
$(_doc_external("Draw/PetscDrawStringSetSize"))
"""
function PetscDrawStringSetSize(petsclib::PetscLibType, draw::PetscDraw, width::Real, height::Real)
    error("PetscDrawStringSetSize: no generated method for these argument types")
end

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
	PetscDrawStringVertical(petsclib::PetscLibType, draw::PetscDraw, xl::PetscReal, yl::PetscReal, cl::Cint, text::String) 
draws text onto a drawable.

Not Collective

Input Parameters:
- `draw` - the drawing context
- `xl`   - coordinate of upper left corner of text
- `yl`   - coordinate of upper left corner of text
- `cl`   - the color of the text
- `text` - the text to draw

Level: beginner

See also: `PetscDraw`, `PetscDrawString()`, `PetscDrawStringCentered()`, `PetscDrawStringBoxed()`, `PetscDrawStringSetSize()`,
`PetscDrawStringGetSize()`

# External Links
$(_doc_external("Draw/PetscDrawStringVertical"))
"""
function PetscDrawStringVertical(petsclib::PetscLibType, draw::PetscDraw, xl::Real, yl::Real, cl::Cint, text::String)
    error("PetscDrawStringVertical: no generated method for these argument types")
end

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
	PetscDrawTensorContour(petsclib::PetscLibType, draw::PetscDraw, m::Cint, n::Cint, xi::Vector{PetscReal}, yi::Vector{PetscReal}, v::Vector{PetscReal}) 
draws a contour plot for a two-dimensional array

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

See also: `PetscDraw`, `PetscDrawTensorContourPatch()`, `PetscDrawScalePopup()`

# External Links
$(_doc_external("Draw/PetscDrawTensorContour"))
"""
function PetscDrawTensorContour(petsclib::PetscLibType, draw::PetscDraw, m::Cint, n::Cint, xi::AbstractVector{<:Number}, yi::AbstractVector{<:Number}, v::AbstractVector{<:Number})
    error("PetscDrawTensorContour: no generated method for these argument types")
end

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
	x::PetscReal,y::PetscReal,v::PetscReal = PetscDrawTensorContourPatch(petsclib::PetscLibType, draw::PetscDraw, m::Cint, n::Cint, min::PetscReal, max::PetscReal) 
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

See also: `PetscDraw`, `PetscDrawTensorContour()`

# External Links
$(_doc_external("Draw/PetscDrawTensorContourPatch"))
"""
function PetscDrawTensorContourPatch(petsclib::PetscLibType, draw::PetscDraw, m::Cint, n::Cint, min::Real, max::Real)
    error("PetscDrawTensorContourPatch: no generated method for these argument types")
end

@for_petsc function PetscDrawTensorContourPatch(petsclib::$UnionPetscLib, draw::PetscDraw, m::Cint, n::Cint, min::$PetscReal, max::$PetscReal )
	x_ = Ref{$PetscReal}()
	y_ = Ref{$PetscReal}()
	v_ = Ref{$PetscReal}()

    @chk ccall(
               (:PetscDrawTensorContourPatch, $petsc_library),
               PetscErrorCode,
               (PetscDraw, Cint, Cint, Ptr{$PetscReal}, Ptr{$PetscReal}, $PetscReal, $PetscReal, Ptr{$PetscReal}),
               draw, m, n, x_, y_, min, max, v_,
              )

	x = x_[]
	y = y_[]
	v = v_[]

	return x,y,v
end 

"""
	PetscDrawTriangle(petsclib::PetscLibType, draw::PetscDraw, x1::PetscReal, y_1::PetscReal, x2::PetscReal, y2::PetscReal, x3::PetscReal, y3::PetscReal, c1::Cint, c2::Cint, c3::Cint) 
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

See also: `PetscDraw`, `PetscDrawLine()`, `PetscDrawRectangle()`, `PetscDrawEllipse()`, `PetscDrawMarker()`, `PetscDrawPoint()`, `PetscDrawArrow()`

# External Links
$(_doc_external("Draw/PetscDrawTriangle"))
"""
function PetscDrawTriangle(petsclib::PetscLibType, draw::PetscDraw, x1::Real, y_1::Real, x2::Real, y2::Real, x3::Real, y3::Real, c1::Cint, c2::Cint, c3::Cint)
    error("PetscDrawTriangle: no generated method for these argument types")
end

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
	PetscDrawUtilitySetCmap(petsclib::PetscLibType, colormap::String, mapsize::Cint, char::Vector{Cuchar}, M_char::Vector{Cuchar}, M_char_2::Vector{Cuchar}) 
Populate the RGB entries of a colormap from a named palette, honoring options-database
overrides for the colormap name, reversal, and brightness.

Not Collective

Input Parameters:
- `colormap` - the name of the colormap (e.g. `"hue"`, `"gray"`, `"jet"`, `"viridis"`), or `NULL`/empty for the default
- `mapsize`  - the number of colormap entries to fill

Output Parameters:
- `R` - the red channel of length `mapsize`
- `G` - the green channel of length `mapsize`
- `B` - the blue channel of length `mapsize`

Options Database Keys:
- `-draw_cmap name`           - select the colormap by name
- `-draw_cmap_reverse`        - reverse the colormap
- `-draw_cmap_brighten value` - brighten (positive) or darken (negative) the colormap; value must be in `(-1, 1)`

Level: developer

See also: `PetscDraw`, `PetscDrawUtilitySetGamma()`

# External Links
$(_doc_external("Draw/PetscDrawUtilitySetCmap"))
"""
function PetscDrawUtilitySetCmap(petsclib::PetscLibType, colormap::String, mapsize::Cint, char::Vector{Cuchar}, M_char::Vector{Cuchar}, M_char_2::Vector{Cuchar})
    error("PetscDrawUtilitySetCmap: no generated method for these argument types")
end

@for_petsc function PetscDrawUtilitySetCmap(petsclib::$UnionPetscLib, colormap::String, mapsize::Cint, char::Vector{Cuchar}, M_char::Vector{Cuchar}, M_char_2::Vector{Cuchar} )

    @chk ccall(
               (:PetscDrawUtilitySetCmap, $petsc_library),
               PetscErrorCode,
               (Ptr{Cchar}, Cint, Ptr{Cuchar}, Ptr{Cuchar}, Ptr{Cuchar}),
               colormap, mapsize, char, M_char, M_char_2,
              )


	return nothing
end 

"""
	PetscDrawUtilitySetGamma(petsclib::PetscLibType, g::PetscReal) 
Set the monitor gamma-correction value used by the drawing colormap utilities.

Not Collective

Input Parameter:
- `g` - the gamma value; a typical value is 2.0

Level: developer

See also: `PetscDraw`, `PetscDrawUtilitySetCmap()`

# External Links
$(_doc_external("Draw/PetscDrawUtilitySetGamma"))
"""
function PetscDrawUtilitySetGamma(petsclib::PetscLibType, g::Real)
    error("PetscDrawUtilitySetGamma: no generated method for these argument types")
end

@for_petsc function PetscDrawUtilitySetGamma(petsclib::$UnionPetscLib, g::$PetscReal )

    @chk ccall(
               (:PetscDrawUtilitySetGamma, $petsc_library),
               PetscErrorCode,
               ($PetscReal,),
               g,
              )


	return nothing
end 

"""
	PetscDrawView(petsclib::PetscLibType, indraw::PetscDraw, viewer::PetscViewer) 
Prints the `PetscDraw` data structure.

Collective

Input Parameters:
- `indraw` - the `PetscDraw` context
- `viewer` - visualization context

See PetscDrawSetFromOptions() for options database keys

See also: `PetscDraw`, `PetscViewerASCIIOpen()`, `PetscViewer`

# External Links
$(_doc_external("Draw/PetscDrawView"))
"""
function PetscDrawView(petsclib::PetscLibType, indraw::PetscDraw, viewer::PetscViewer)
    error("PetscDrawView: no generated method for these argument types")
end

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
	PetscDrawViewFromOptions(petsclib::PetscLibType, A::PetscDraw, obj, name::String) 
View a `PetscDraw` from the option database

Collective

Input Parameters:
- `A`    - the `PetscDraw` context
- `obj`  - Optional object
- `name` - command line option

Options Database Key:
- `-name [viewertype][:...]` - option name and values. See `PetscObjectViewFromOptions()` for the possible arguments

Level: intermediate

See also: `PetscDraw`, `PetscDrawView`, `PetscObjectViewFromOptions()`, `PetscDrawCreate()`

# External Links
$(_doc_external("Draw/PetscDrawViewFromOptions"))
"""
function PetscDrawViewFromOptions(petsclib::PetscLibType, A::PetscDraw, obj, name::String)
    error("PetscDrawViewFromOptions: no generated method for these argument types")
end

@for_petsc function PetscDrawViewFromOptions(petsclib::$UnionPetscLib, A::PetscDraw, obj, name::String )

    @chk ccall(
               (:PetscDrawViewFromOptions, $petsc_library),
               PetscErrorCode,
               (PetscDraw, PetscObject, Ptr{Cchar}),
               A, obj, name,
              )


	return nothing
end 

"""
	newports::Ptr{PetscDrawViewPorts} = PetscDrawViewPortsCreate(petsclib::PetscLibType, draw::PetscDraw, nports::PetscInt) 
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

See also: `PetscDrawSplitViewPort()`, `PetscDrawSetViewPort()`, `PetscDrawViewPortsSet()`, `PetscDrawViewPortsDestroy()`

# External Links
$(_doc_external("Draw/PetscDrawViewPortsCreate"))
"""
function PetscDrawViewPortsCreate(petsclib::PetscLibType, draw::PetscDraw, nports::Integer)
    error("PetscDrawViewPortsCreate: no generated method for these argument types")
end

@for_petsc function PetscDrawViewPortsCreate(petsclib::$UnionPetscLib, draw::PetscDraw, nports::$PetscInt )
	newports_ = Ref{Ptr{PetscDrawViewPorts}}()

    @chk ccall(
               (:PetscDrawViewPortsCreate, $petsc_library),
               PetscErrorCode,
               (PetscDraw, $PetscInt, Ptr{Ptr{PetscDrawViewPorts}}),
               draw, nports, newports_,
              )

	newports = newports_[]

	return newports
end 

"""
	newports::Ptr{PetscDrawViewPorts} = PetscDrawViewPortsCreateRect(petsclib::PetscLibType, draw::PetscDraw, nx::PetscInt, ny::PetscInt) 
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

See also: `PetscDrawSplitViewPort()`, `PetscDrawSetViewPort()`, `PetscDrawViewPortsSet()`, `PetscDrawViewPortsDestroy()`, `PetscDrawViewPorts`

# External Links
$(_doc_external("Draw/PetscDrawViewPortsCreateRect"))
"""
function PetscDrawViewPortsCreateRect(petsclib::PetscLibType, draw::PetscDraw, nx::Integer, ny::Integer)
    error("PetscDrawViewPortsCreateRect: no generated method for these argument types")
end

@for_petsc function PetscDrawViewPortsCreateRect(petsclib::$UnionPetscLib, draw::PetscDraw, nx::$PetscInt, ny::$PetscInt )
	newports_ = Ref{Ptr{PetscDrawViewPorts}}()

    @chk ccall(
               (:PetscDrawViewPortsCreateRect, $petsc_library),
               PetscErrorCode,
               (PetscDraw, $PetscInt, $PetscInt, Ptr{Ptr{PetscDrawViewPorts}}),
               draw, nx, ny, newports_,
              )

	newports = newports_[]

	return newports
end 

"""
	PetscDrawViewPortsDestroy(petsclib::PetscLibType, ports::Vector{PetscDrawViewPorts}) 
frees a `PetscDrawViewPorts` object

Collective on the `PetscDraw` inside `ports`

Input Parameter:
- `ports` - the `PetscDrawViewPorts` object

Level: advanced

See also: `PetscDrawViewPorts`, `PetscDrawSplitViewPort()`, `PetscDrawSetViewPort()`, `PetscDrawViewPortsSet()`, `PetscDrawViewPortsCreate()`

# External Links
$(_doc_external("Draw/PetscDrawViewPortsDestroy"))
"""
function PetscDrawViewPortsDestroy(petsclib::PetscLibType, ports::Vector{PetscDrawViewPorts})
    error("PetscDrawViewPortsDestroy: no generated method for these argument types")
end

@for_petsc function PetscDrawViewPortsDestroy(petsclib::$UnionPetscLib, ports::Vector{PetscDrawViewPorts} )

    @chk ccall(
               (:PetscDrawViewPortsDestroy, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscDrawViewPorts},),
               ports,
              )


	return nothing
end 

"""
	PetscDrawViewPortsSet(petsclib::PetscLibType, ports::Vector{PetscDrawViewPorts}, port::PetscInt) 
sets a draw object to use a particular subport

Logically Collective on the `PetscDraw` inside `ports`

Input Parameters:
- `ports` - the `PetscDrawViewPorts` object
- `port`  - the port number, from 0 to nports-1

Level: advanced

See also: `PetscDrawViewPorts`, `PetscDrawSplitViewPort()`, `PetscDrawSetViewPort()`, `PetscDrawViewPortsDestroy()`, `PetscDrawViewPortsCreate()`

# External Links
$(_doc_external("Draw/PetscDrawViewPortsSet"))
"""
function PetscDrawViewPortsSet(petsclib::PetscLibType, ports::Vector{PetscDrawViewPorts}, port::Integer)
    error("PetscDrawViewPortsSet: no generated method for these argument types")
end

@for_petsc function PetscDrawViewPortsSet(petsclib::$UnionPetscLib, ports::Vector{PetscDrawViewPorts}, port::$PetscInt )

    @chk ccall(
               (:PetscDrawViewPortsSet, $petsc_library),
               PetscErrorCode,
               (Ptr{PetscDrawViewPorts}, $PetscInt),
               ports, port,
              )


	return nothing
end 

"""
	PetscDrawZoom(petsclib::PetscLibType, draw::PetscDraw, func::external, ctx::Ptr{Cvoid}) 
Allows one to provide a function that gets called for zooming in on a drawing using the mouse buttons

Collective draw

Input Parameters:
- `draw` - the window where the graph will be made.
- `func` - users function that draws the graphic
- `ctx`  - pointer to any application required data

Calling sequence of func:
- `draw` - the `PetscDraw` object to zoom on
- `ctx`  - the context for the zooming operation

Level: advanced

See also: `PetscDraw`, `PetscDrawCreate()`

# External Links
$(_doc_external("Draw/PetscDrawZoom"))
"""
function PetscDrawZoom(petsclib::PetscLibType, draw::PetscDraw, func::external, ctx::Ptr{Cvoid})
    error("PetscDrawZoom: no generated method for these argument types")
end

@for_petsc function PetscDrawZoom(petsclib::$UnionPetscLib, draw::PetscDraw, func::external, ctx::Ptr{Cvoid} )

    @chk ccall(
               (:PetscDrawZoom, $petsc_library),
               PetscErrorCode,
               (PetscDraw, external, Ptr{Cvoid}),
               draw, func, ctx,
              )


	return nothing
end 

