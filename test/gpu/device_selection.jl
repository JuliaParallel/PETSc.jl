# test/gpu/device_selection.jl
# PETScCUDAExt and PETScAMDGPUExt: PETSc uses the device the Julia GPU package
# has selected, and an explicit -device_select_* option still overrides it.
#
# Not part of runtests.jl: it needs a GPU. With a PETSc built with CUDA or HIP
# it checks the selection; with two or more devices, a device other than PETSc's
# own default. With any other PETSc it checks that the option is left unset,
# since that build never reads it. Run it like test/gpu/device_arrays.jl:
#
#   julia --project=ENV test/gpu/device_selection.jl cuda    # or hip

using Test, Libdl
using PETSc, MPI

const gpu = get(ARGS, 1, "")
if gpu == "cuda"
    using CUDA
    const option = "-device_select_cuda"
    const device_function = :VecCreateSeqCUDA
    devices() = collect(CUDA.devices())
    select!(d) = CUDA.device!(d)
    selected_id() = CUDA.deviceid(CUDA.device())
elseif gpu == "hip"
    using AMDGPU
    const option = "-device_select_hip"
    const device_function = :VecCreateSeqHIP
    devices() = AMDGPU.devices()
    select!(d) = AMDGPU.device!(d)
    selected_id() = AMDGPU.device_id(AMDGPU.device()) - 1
else
    error("pass the GPU backend as the argument: cuda or hip, got $gpu")
end

MPI.Initialized() || MPI.Init()
petsclib = PETSc.petsclibs[1]

# The id of the device behind PETSc's current device context. Called through
# dlsym: the generated wrappers of these two functions pass the handles by value.
function petsc_device_id(petsclib::PETSc.LibPETSc.PetscLibType{ST, IT}) where {ST, IT}
    h = Libdl.dlopen(petsclib.petsc_library)
    dctx, dev = Ref{Ptr{Cvoid}}(), Ref{Ptr{Cvoid}}()
    id = Ref{IT}()
    errs = (
        ccall(dlsym(h, :PetscDeviceContextGetCurrentContext), Cint, (Ptr{Ptr{Cvoid}},), dctx),
        ccall(dlsym(h, :PetscDeviceContextGetDevice), Cint, (Ptr{Cvoid}, Ptr{Ptr{Cvoid}}), dctx[], dev),
        ccall(dlsym(h, :PetscDeviceGetDeviceId), Cint, (Ptr{Cvoid}, Ptr{IT}), dev[], id),
    )
    all(iszero, errs) || error("PETSc device query failed with error codes $errs")
    return Int(id[])
end

@testset "device selection ($gpu)" begin
    if !PETSc.has_symbol(petsclib, device_function)
        @testset "a PETSc without $gpu leaves $option unset" begin
            PETSc.initialize(petsclib)
            global_options = PETSc.LibPETSc.PetscOptions{typeof(petsclib)}(C_NULL, petsclib.age; own = false)
            @test !PETSc.LibPETSc.PetscOptionsHasName(petsclib, global_options, "", option)
            PETSc.finalize(petsclib)
        end
    else
        # PETSc's own default on one rank is device 0: pick the last one instead
        select!(last(devices()))
        @testset "PETSc follows the Julia GPU package" begin
            PETSc.initialize(petsclib)
            @test petsc_device_id(petsclib) == selected_id()
            PETSc.finalize(petsclib)
        end

        @testset "an explicit option overrides it" begin
            PETSc.initialize(petsclib; options = [option, "0"])
            @test petsc_device_id(petsclib) == 0
            PETSc.finalize(petsclib)
        end
    end
end
