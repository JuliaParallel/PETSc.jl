struct PetscFEGeom
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

    PetscFEGeom(mode,isAffine,dim,dimEmbed,numCells,numPoints,xi,v,J,invJ,detJ,n,face,suppJ,suppInvJ,suppDetJ) = new(mode::PetscFEGeomMode,isAffine::PetscBool,dim::PetscInt,dimEmbed::PetscInt,numCells::PetscInt,numPoints::PetscInt,xi::Ptr{PetscReal},v::Ptr{PetscReal},J::Ptr{PetscReal},invJ::Ptr{PetscReal},detJ::Ptr{PetscReal},n::Ptr{PetscReal},face::Ptr{NTuple{4, PetscInt}},suppJ::Ptr{NTuple{2, PetscReal}},suppInvJ::Ptr{NTuple{2, PetscReal}},suppDetJ::Ptr{NTuple{2, PetscReal}})

end 

struct VecTaggerBox
    min::PetscScalar
    max::PetscScalar

    VecTaggerBox() = new()

    VecTaggerBox(min,max) = new(min::PetscScalar,max::PetscScalar)

end 

struct PetscDrawViewPorts
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

    PetscDrawViewPorts(nports,xl,xr,yl,yr,draw,port_xl) = new(nports::PetscInt,xl::Ptr{PetscReal},xr::Ptr{PetscReal},yl::Ptr{PetscReal},yr::Ptr{PetscReal},draw::PetscDraw,port_xl::PetscReal ,    port_yl::PetscReal ,    port_xr::PetscReal ,    port_yr::PetscReal)

end 

struct PetscFormKey
    label::DMLabel
    value::PetscInt
    field::PetscInt
    part::PetscInt

    PetscFormKey() = new()

    PetscFormKey(label,value,field,part) = new(label::DMLabel,value::PetscInt,field::PetscInt,part::PetscInt)

end 

struct DMStagStencil
    loc::DMStagStencilLocation
    i::PetscInt 
    j::PetscInt 
    k::PetscInt 
    c::PetscInt 
    

    DMStagStencil() = new()

    DMStagStencil(loc,i,j,k,c) = new(loc::DMStagStencilLocation,i::PetscInt ,    j::PetscInt ,    k::PetscInt ,    c::PetscInt)

end 

struct MatStencil
    # Must match PETSc C struct order: struct MatStencil { PetscInt i,j,k,c; };
    i::PetscInt
    j::PetscInt
    k::PetscInt
    c::PetscInt

    MatStencil() = new()
        
    # Convenience constructor accepting integers in i,j,k,c order
    MatStencil(i,j,k,c) = new(i::PetscInt, j::PetscInt, k::PetscInt, c::PetscInt)
end 

struct MatInfo
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

    MatInfo(block_size,nz_allocated,memory,assemblies,mallocs,fill_ratio_given,factor_mallocs) = new(block_size::PetscLogDouble,nz_allocated::PetscLogDouble,nz_used::PetscLogDouble ,nz_unneeded::PetscLogDouble,memory::PetscLogDouble,assemblies::PetscLogDouble,mallocs::PetscLogDouble,fill_ratio_given::PetscLogDouble ,fill_ratio_needed::PetscLogDouble,factor_mallocs::PetscLogDouble)

end 

struct MatFactorInfo
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

    MatFactorInfo(diagonal_fill,usedt,dt,dtcol,dtcount,fill,levels,pivotinblocks,zeropivot,shifttype,shiftamount,factoronhost,solveonhost) = new(diagonal_fill::PetscReal,usedt::PetscReal,dt::PetscReal,dtcol::PetscReal,dtcount::PetscReal,fill::PetscReal,levels::PetscReal,pivotinblocks::PetscReal,zeropivot::PetscReal,shifttype::PetscReal,shiftamount::PetscReal,factoronhost::PetscBool,solveonhost::PetscBool)

end 

struct LandauStaticData
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

    LandauStaticData(invJ,D,B,alpha,beta,invMass,w,x,y,z,Eq_m,f,dfdx,dfdy,dfdz,dim_,NCells,species_offset,mat_offset,elem_offset,ip_offset,ipf_offset,ipfdf_data,maps,coo_elem_offsets,coo_elem_point_offsets,coo_elem_fullNb,coo_vals,lambdas,coo_n_cellsTot,coo_size,coo_max_fullnb) = new(invJ::Ptr{Cvoid},D::Ptr{Cvoid},B::Ptr{Cvoid},alpha::Ptr{Cvoid},beta::Ptr{Cvoid},invMass::Ptr{Cvoid},w::Ptr{Cvoid},x::Ptr{Cvoid},y::Ptr{Cvoid},z::Ptr{Cvoid},Eq_m::Ptr{Cvoid},f::Ptr{Cvoid},dfdx::Ptr{Cvoid},dfdy::Ptr{Cvoid},dfdz::Ptr{Cvoid},dim_::Cint,ns_::Cint ,    nip_::Cint ,    nq_::Cint ,    nb_::Cint ,NCells::Ptr{Cvoid},species_offset::Ptr{Cvoid},mat_offset::Ptr{Cvoid},elem_offset::Ptr{Cvoid},ip_offset::Ptr{Cvoid},ipf_offset::Ptr{Cvoid},ipfdf_data::Ptr{Cvoid},maps::Ptr{Cvoid},coo_elem_offsets::Ptr{Cvoid},coo_elem_point_offsets::Ptr{Cvoid},coo_elem_fullNb::Ptr{Cvoid},coo_vals::Ptr{Cvoid},lambdas::Ptr{Cvoid},coo_n_cellsTot::LandauIdx,coo_size::LandauIdx,coo_max_fullnb::LandauIdx)

end 

struct PetscStack
    _function::Ptr{NTuple{PETSCSTACKSIZE, Cchar}}
    file::Ptr{NTuple{PETSCSTACKSIZE, Cchar}}
    line::NTuple{PETSCSTACKSIZE, Cint}
    petscroutine::NTuple{PETSCSTACKSIZE, Cint}
    currentsize::Cint
    hotdepth::Cint
    check::PetscBool

    PetscStack() = new()

    PetscStack(_function,file,line,petscroutine,currentsize,hotdepth,check) = new(_function::Ptr{NTuple{PETSCSTACKSIZE, Cchar}},file::Ptr{NTuple{PETSCSTACKSIZE, Cchar}},line::NTuple{PETSCSTACKSIZE, Cint},petscroutine::NTuple{PETSCSTACKSIZE, Cint},currentsize::Cint,hotdepth::Cint,check::PetscBool)

end 

struct JacActionCtx
    dm::PetscDM
    u::PetscVec
    J::PetscMat
    user::Ptr{Cvoid}

    JacActionCtx() = new()

    JacActionCtx(dm,u,J,user) = new(dm::PetscDM,u::PetscVec,J::PetscMat,user::Ptr{Cvoid})

end 

struct TSMonitorDMDARayCtx
    ray::PetscVec
    scatter::VecScatter
    viewer::PetscViewer
    lgctx::TSMonitorLGCtx

    TSMonitorDMDARayCtx() = new()

    TSMonitorDMDARayCtx(ray,scatter,viewer,lgctx) = new(ray::PetscVec,scatter::VecScatter,viewer::PetscViewer,lgctx::TSMonitorLGCtx)

end 

struct PCMPIServerAddresses
    n::PetscInt
    addr::Ptr{NTuple{3, Cvoid}}

    PCMPIServerAddresses() = new()

    PCMPIServerAddresses(n,addr) = new(n::PetscInt,addr::Ptr{NTuple{3, Cvoid}})

end 

struct PetscFVFaceGeom
    normal::NTuple{3, PetscReal}
    centroid::NTuple{3, PetscReal}
    grad::NTuple{2, PetscScalar}

    PetscFVFaceGeom() = new()

    PetscFVFaceGeom(normal,centroid,grad) = new(normal::NTuple{3, PetscReal},centroid::NTuple{3, PetscReal},grad::NTuple{2, PetscScalar})

end 

struct PetscFVCellGeom
    centroid::NTuple{3, PetscReal}
    volume::PetscReal

    PetscFVCellGeom() = new()

    PetscFVCellGeom(centroid,volume) = new(centroid::NTuple{3, PetscReal},volume::PetscReal)

end 

struct PetscViewerAndFormat
    viewer::PetscViewer
    format::PetscViewerFormat
    view_interval::PetscInt
    data::Ptr{Cvoid}
    data_destroy::Ptr{PetscCtxDestroyFn}

    PetscViewerAndFormat() = new()

    PetscViewerAndFormat(viewer,format,view_interval,data,data_destroy) = new(viewer::PetscViewer,format::PetscViewerFormat,view_interval::PetscInt,data::Ptr{Cvoid},data_destroy::Ptr{PetscCtxDestroyFn})

end 

struct PetscSFNode
    rank::PetscInt
    index::PetscInt

    PetscSFNode() = new()

    PetscSFNode(rank,index) = new(rank::PetscInt,index::PetscInt)

end 

struct DMDALocalInfo
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

    DMDALocalInfo(da,dim,mx,xs,xm,gxs,gxm,bx,st) = new(da::PetscDM,dim::PetscInt ,    dof::PetscInt ,    sw::PetscInt ,   mx::PetscInt ,    my::PetscInt ,    mz::PetscInt ,xs::PetscInt ,    ys::PetscInt ,    zs::PetscInt ,    xm::PetscInt ,    ym::PetscInt ,    zm::PetscInt ,    gxs::PetscInt ,    gys::PetscInt ,    gzs::PetscInt ,    gxm::PetscInt ,    gym::PetscInt ,    gzm::PetscInt ,    bx::DMBoundaryType ,    by::DMBoundaryType ,    bz::DMBoundaryType ,    st::DMDAStencilType)

end 

struct PetscEventPerfInfo
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

    PetscEventPerfInfo(id,active,visible,depth,count,flops,flops2,flopsTmp,time,time2,timeTmp,syncTime,dof,errors,numMessages,messageLength,numReductions,memIncrease,mallocIncrease,mallocSpace,mallocIncreaseEvent) = new(id::Cint,active::PetscBool,visible::PetscBool,depth::Cint,count::Cint,flops::PetscLogDouble,flops2::PetscLogDouble,flopsTmp::PetscLogDouble,time::PetscLogDouble,time2::PetscLogDouble,timeTmp::PetscLogDouble,syncTime::PetscLogDouble,dof::NTuple{8, PetscLogDouble},errors::NTuple{8, PetscLogDouble},numMessages::PetscLogDouble,messageLength::PetscLogDouble,numReductions::PetscLogDouble,memIncrease::PetscLogDouble,mallocIncrease::PetscLogDouble,mallocSpace::PetscLogDouble,mallocIncreaseEvent::PetscLogDouble)

end 

struct PetscLogEventInfo
    name::Ptr{Cchar}
    classid::PetscClassId
    collective::PetscBool

    PetscLogEventInfo() = new()

    PetscLogEventInfo(name,classid,collective) = new(name::Ptr{Cchar},classid::PetscClassId,collective::PetscBool)

end 

struct PetscLogClassInfo
    name::Ptr{Cchar}
    classid::PetscClassId

    PetscLogClassInfo() = new()

    PetscLogClassInfo(name,classid) = new(name::Ptr{Cchar},classid::PetscClassId)

end 

struct PetscLogStageInfo
    name::Ptr{Cchar}

    PetscLogStageInfo() = new()

    PetscLogStageInfo(name) = new(name::Ptr{Cchar})

end 

struct DMDACoor
    x::PetscScalar 
    y::PetscScalar 
    z::PetscScalar 
    

    DMDACoor() = new()

    DMDACoor(x) = new(x::PetscScalar ,    y::PetscScalar ,    z::PetscScalar)

end 

