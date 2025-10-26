mutable struct PetscFEGeom
    mode::PetscFEGeomMode
    isAffine::PetscBool
    dim::PetscInt
    dimEmbed::PetscInt
    numCells::PetscInt
    numPoints::PetscInt
    xi::Ptr{PetscReal}
    v::Ptr{PetscReal}
    J::Ptr{PetscReal}
    invJ::Ptr{PetscReal}
    detJ::Ptr{PetscReal}
    n::Ptr{PetscReal}
    face::Ptr{NTuple{4, PetscInt}}
    suppJ::Ptr{NTuple{2, PetscReal}}
    suppInvJ::Ptr{NTuple{2, PetscReal}}
    suppDetJ::Ptr{NTuple{2, PetscReal}}
    PetscFEGeom() = new()
end 

mutable struct VecTaggerBox
    min::PetscScalar
    max::PetscScalar
    VecTaggerBox() = new()
end 

mutable struct PetscDrawViewPorts
    nports::PetscInt
    xl::Ptr{PetscReal}
    xr::Ptr{PetscReal}
    yl::Ptr{PetscReal}
    yr::Ptr{PetscReal}
    draw::PetscDraw
    port_xl::PetscReal 
    port_yl::PetscReal 
    port_xr::PetscReal 
    port_yr::PetscReal 
    
    PetscDrawViewPorts() = new()
end 

mutable struct PetscFormKey
    label::DMLabel
    value::PetscInt
    field::PetscInt
    part::PetscInt
    PetscFormKey() = new()
end 

mutable struct DMStagStencil
    loc::DMStagStencilLocation
    i::PetscInt 
    j::PetscInt 
    k::PetscInt 
    c::PetscInt 
    
    DMStagStencil() = new()
end 

mutable struct MatStencil
    k::PetscInt 
    j::PetscInt 
    i::PetscInt 
    c::PetscInt 
    
    MatStencil() = new()
end 

mutable struct MatInfo
    block_size::PetscLogDouble
    nz_allocated::PetscLogDouble 
    nz_used::PetscLogDouble 
    nz_unneeded::PetscLogDouble 
    
    memory::PetscLogDouble
    assemblies::PetscLogDouble
    mallocs::PetscLogDouble
    fill_ratio_given::PetscLogDouble 
    fill_ratio_needed::PetscLogDouble 
    
    factor_mallocs::PetscLogDouble
    MatInfo() = new()
end 

mutable struct MatFactorInfo
    diagonal_fill::PetscReal
    usedt::PetscReal
    dt::PetscReal
    dtcol::PetscReal
    dtcount::PetscReal
    fill::PetscReal
    levels::PetscReal
    pivotinblocks::PetscReal
    zeropivot::PetscReal
    shifttype::PetscReal
    shiftamount::PetscReal
    factoronhost::PetscBool
    solveonhost::PetscBool
    MatFactorInfo() = new()
end 

mutable struct LandauStaticData
    invJ::Ptr{Cvoid}
    D::Ptr{Cvoid}
    B::Ptr{Cvoid}
    alpha::Ptr{Cvoid}
    beta::Ptr{Cvoid}
    invMass::Ptr{Cvoid}
    w::Ptr{Cvoid}
    x::Ptr{Cvoid}
    y::Ptr{Cvoid}
    z::Ptr{Cvoid}
    Eq_m::Ptr{Cvoid}
    f::Ptr{Cvoid}
    dfdx::Ptr{Cvoid}
    dfdy::Ptr{Cvoid}
    dfdz::Ptr{Cvoid}
    dim_::Cint 
    ns_::Cint 
    nip_::Cint 
    nq_::Cint 
    nb_::Cint 
    
    NCells::Ptr{Cvoid}
    species_offset::Ptr{Cvoid}
    mat_offset::Ptr{Cvoid}
    elem_offset::Ptr{Cvoid}
    ip_offset::Ptr{Cvoid}
    ipf_offset::Ptr{Cvoid}
    ipfdf_data::Ptr{Cvoid}
    maps::Ptr{Cvoid}
    coo_elem_offsets::Ptr{Cvoid}
    coo_elem_point_offsets::Ptr{Cvoid}
    coo_elem_fullNb::Ptr{Cvoid}
    coo_vals::Ptr{Cvoid}
    lambdas::Ptr{Cvoid}
    coo_n_cellsTot::LandauIdx
    coo_size::LandauIdx
    coo_max_fullnb::LandauIdx
    LandauStaticData() = new()
end 

mutable struct pointInterpolationP
    scale::PetscReal
    gid::LandauIdx
    pointInterpolationP() = new()
end 

mutable struct PetscStack
    _function::Ptr{NTuple{PETSCSTACKSIZE, Cchar}}
    file::Ptr{NTuple{PETSCSTACKSIZE, Cchar}}
    line::NTuple{PETSCSTACKSIZE, Cint}
    petscroutine::NTuple{PETSCSTACKSIZE, Cint}
    currentsize::Cint
    hotdepth::Cint
    check::PetscBool
    PetscStack() = new()
end 

mutable struct JacActionCtx
    dm::PetscDM
    u::PetscVec
    J::PetscMat
    user::Ptr{Cvoid}
    JacActionCtx() = new()
end 

mutable struct TSMonitorDMDARayCtx
    ray::PetscVec
    scatter::VecScatter
    viewer::PetscViewer
    lgctx::TSMonitorLGCtx
    TSMonitorDMDARayCtx() = new()
end 

mutable struct PCMPIServerAddresses
    n::PetscInt
    addr::Ptr{NTuple{3, Cvoid}}
    PCMPIServerAddresses() = new()
end 

mutable struct PetscFVFaceGeom
    normal::NTuple{3, PetscReal}
    centroid::NTuple{3, PetscReal}
    grad::NTuple{2, PetscScalar}
    PetscFVFaceGeom() = new()
end 

mutable struct PetscFVCellGeom
    centroid::NTuple{3, PetscReal}
    volume::PetscReal
    PetscFVCellGeom() = new()
end 

mutable struct PetscViewerAndFormat
    viewer::PetscViewer
    format::PetscViewerFormat
    view_interval::PetscInt
    data::Ptr{Cvoid}
    data_destroy::Ptr{PetscCtxDestroyFn}
    PetscViewerAndFormat() = new()
end 

mutable struct PetscSFNode
    rank::PetscInt
    index::PetscInt
    PetscSFNode() = new()
end 

mutable struct DMDALocalInfo
    da::PetscDM
    dim::PetscInt 
    dof::PetscInt 
    sw::PetscInt 
    
    mx::PetscInt 
    my::PetscInt 
    mz::PetscInt 
    
    xs::PetscInt 
    ys::PetscInt 
    zs::PetscInt 
    
    xm::PetscInt 
    ym::PetscInt 
    zm::PetscInt 
    
    gxs::PetscInt 
    gys::PetscInt 
    gzs::PetscInt 
    
    gxm::PetscInt 
    gym::PetscInt 
    gzm::PetscInt 
    
    bx::DMBoundaryType 
    by::DMBoundaryType 
    bz::DMBoundaryType 
    
    st::DMDAStencilType
    DMDALocalInfo() = new()
end 

mutable struct PetscEventPerfInfo
    id::Cint
    active::PetscBool
    visible::PetscBool
    depth::Cint
    count::Cint
    flops::PetscLogDouble
    flops2::PetscLogDouble
    flopsTmp::PetscLogDouble
    time::PetscLogDouble
    time2::PetscLogDouble
    timeTmp::PetscLogDouble
    syncTime::PetscLogDouble
    dof::NTuple{8, PetscLogDouble}
    errors::NTuple{8, PetscLogDouble}
    numMessages::PetscLogDouble
    messageLength::PetscLogDouble
    numReductions::PetscLogDouble
    memIncrease::PetscLogDouble
    mallocIncrease::PetscLogDouble
    mallocSpace::PetscLogDouble
    mallocIncreaseEvent::PetscLogDouble
    PetscEventPerfInfo() = new()
end 

mutable struct PetscLogEventInfo
    name::Ptr{Cchar}
    classid::PetscClassId
    collective::PetscBool
    PetscLogEventInfo() = new()
end 

mutable struct PetscLogClassInfo
    name::Ptr{Cchar}
    classid::PetscClassId
    PetscLogClassInfo() = new()
end 

mutable struct PetscLogStageInfo
    name::Ptr{Cchar}
    PetscLogStageInfo() = new()
end 

mutable struct DMDACoor
    x::PetscScalar 
    y::PetscScalar 
    z::PetscScalar 
    
    DMDACoor() = new()
end 

