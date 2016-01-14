using PETSc
if VERSION >= v"0.5.0-dev+7720"
  using Base.Test
else
  using BaseTestNext
end

# determine scalar type of current run
global ST = Float64  # scalar type

function RC(x::Number)
# used to test real, complex
  if ST == Float64
    return Float64(real(x))
  elseif ST == Float32
    return Float32(real(x))
  else  # scalar_type == 3
    return complex(x)
  end
end

function RC(x::AbstractArray)
# used to test real, complex
  if ST == Float64
    tmp = similar(x, ST)
    for i=1:length(x)
      tmp[i] = Float64(real(x[i]))
    end
    return tmp
  elseif ST == Float32
    tmp = similar(x, ST)
    for i=1:length(x)
      tmp[i] = ST(real(x[i]))
    end
    return tmp
  else  # scalar_type == 3
    return x
  end
end

# convert to PetscReal
function RT(x::Number)
  if ST == Float64 || ST == Complex128
    return Float64(x)
  else
    return Float32(x)
  end

end

