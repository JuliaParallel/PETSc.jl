# provide some most commonly used options, leave rest as low level
# common options: Orthogonilization type, KSP type, PC type
# use PC context created as part of KSP

export KSP, setoptions!

type KSP{T, MType}
  pksp::C.KSP{T}
  ppc::C.PC{T}
  own_pc::Bool  # is the PC owned by the ksp context
  A::Mat{T, MType}
end

function KSP{T, MType}(A::Mat{T, MType}, pc_mat::Mat{T, MType}=A)
  ksp_arr = Array(C.KSP{T}, 1)
  pc_arr = Array(C.PC{T}, 1)

  chk(C.KSPCreate(A.comm, ksp_arr))
  ksp = ksp_arr[1]
  C.KSPGetPC(ksp, pc_arr)
  pc = pc_arr[1]

  chk(C.KSPSetOperators(ksp, A.p, pc_mat.p))
 # call KSPSetOptions from here 

  return KSP{T, MType}(ksp, pc, true, A)
end  # add finalizer


function KSPDestroy(ksp::KSP)

  tmp = Array(PetscBool, 1)
  C.PetscFinalized(eltype(mat), tmp)
   
  if tmp[1] == 0  # if petsc has not been finalized yet
    if !ksp.own_pc
      C.PCDestroy([ksp.ppc])
    end

    C.KSPDestroy([ksp.pksp])
  end

   # if Petsc has been finalized, let the OS deallocate the memory
end

function setoptions!{T, MType}(ksp::KSP{T, MType}; opts...)
# sets some options in the options database, call KSPSetFromOptions, then reset the
# opts becomes an Array of tuples
# options database to the original state
# this prevents unexpected dynamic phenomena like setting an option for one KSP 
# contex and having it still be set for another
# the keys in the dictionary should have the prepended -
  # to get options we have to provide a string buffer of sufficient length
  # to be populated with the returned string
  # what is a 'sufficient length'? no way to know, so make it 256 characters
  # for now and enlarge if needed later

  # set the options in the (global) options databse
  opts_orig, opts_unset = setoptions!(T, opts)

  # copy options into ksp object
  chk(C.KSPSetFromOptions(ksp.pksp))
 
  # reset the options in the databse
  unset_options!(T, opts_orig, opts_unset)
 
  return nothing
end


function setoptions!(T::DataType, opts)
# opts is any iterable container of tuples containing a key and a value, both
# strings

  println("setting options")
  println("typeof(opts) = ", typeof(opts))

  len = Csize_t(256)
  string_buff = (Array(UInt8, len))
  string_buff2 = string(string_buff)
  null_str = string(UInt8[0])  # null prefix
  opts_orig = Dict{UTF8String, UTF8String}()  # options with existing values
  opts_unset = Set{UTF8String}()  # options without existing values
  isset = Ref{PetscBool}()

  println("opts = ", opts)
  for i in opts
    # record the original option value
    println("i = ", i)
    println("typeof(i) = ", typeof(i))
    chk(C.PetscOptionsGetString(T, null_str, addPrefix(i[1]), string_buff2, len, isset))

    if isset[] != 0  # if an option with the specified name was found
      str = bytestring(pointer(string_buff2))
      opts_orig[string(i[1])] = str
      println("value of option $i is ", str)
    else  # option has not previously been set
      push!(opts_unset, string(i[1]))
    end

    # now set the the option
    chk(C.PetscOptionsSetValue(T, addPrefix(i[1]), string(i[2])))
  end  # end loop over opts

  return opts_orig, opts_unset

end

function addPrefix(val)
  if !startswith(string(val), '-')
  return string("-", val)
  else
    return string(val)
  end
end


function unset_options!(T::DataType, opts_orig, opts_unset::Set)
# opts_orig is an iterable container of tuples containing keys and values

  # reset options that had original values
  for i in opts_orig
    chk(C.PetscOptionsSetValues(T, addPrefix(i[1]), string(i[2])))
  end
  # now reset the options that were not previously set
  for i in opts_unset
    chk(C.PetscOptionsClearValue(T, addPrefix(i)))
  end

end


# can use options databse instead
function settolerances{T}(ksp::KSP{T}; rtol=1e-8, abstol=1e-12, dtol=1e5, maxits=size(ksp.A, 1))

  C.KSPSetTolerances(ksp, rtol, abstol, dtol, maxits)
end

# A ldiv B 

import Base: A_ldiv_B!
function A_ldiv_B!{T}(ksp::KSP{T}, b::Vec{T}, x::Vec{T})
# perform the solve
# users should specify all the options they want to use
# before calling this function
# if solving multiple rhs with the same matrix A,
# the preconditioner is resued automatically
# if A changes, the preconditioner is recomputed

#  C.KSPSetUp(ksp.pksp)   # this is called by KSPSolve if needed
                   # decreases logging accurace of setup operations


  # assemble the matrix
  AssemblyBegin(ksp.A, C.MAT_FINAL_ASSEMBLY)
  AssemblyEnd(ksp.A, C.MAT_FINAL_ASSEMBLY)

  # assemble the vector
  AssemblyBegin(x, C.MAT_FINAL_ASSEMBLY)
  AssemblyEnd(x, C.MAT_FINAL_ASSEMBLY)
 
#  chk(C.KSPSetFromOptions(ksp.pksp))
  chk(C.KSPSolve(ksp.pksp, b.p, x.p))

  reason_arr = Array(Cint, 1)
  chk(C.KSPGetConvergedReason(ksp.pksp, reason_arr))
  reason = reason_arr[1]



  if reason < 0
    println(STDERR, "Warning: KSP Solve did not converge")
  end

end

import Base: \
function (\){T}(ksp::KSP{T}, b::Vec{T})
# this only works for square systems
  x = similar(b)
  A_ldiv_B!(ksp, b, x)
  return x
end



  
