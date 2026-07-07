mutable struct PETScTSIntegrator{
    algType <: PETScTSAlgorithm,
    uType,
    tType,
    pType,
    solType,
    optType,
    cbCacheType,
    libType,
    tsType,
    vecType,
    cbCtxType,
} <: SciMLBase.AbstractODEIntegrator{algType, true, uType, tType}
    alg::algType
    u::uType
    uprev::uType
    t::tType
    tprev::tType
    dt::tType
    p::pType
    opts::optType
    u_modified::Bool
    derivative_discontinuity::Bool
    tdir::tType
    sizeu::Tuple
    sol::solType
    callback_cache::cbCacheType
    petsclib::libType
    ts::tsType
    u_petsc::vecType
    cb_ctx::cbCtxType
    initialized::Bool
    done::Bool
    retcode::SciMLBase.ReturnCode.T
end
