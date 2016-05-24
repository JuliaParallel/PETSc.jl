using PETSc


function write_names(f)
nms = names(PETSc, true, true)

for i in nms
  val = @eval(PETSc.$i) 
  println("symbol = ", i)
  println("typeof(val) = ", typeof(val))
  if !(typeof(val) <: Dict || typeof(val) <: PETSc.Vec || typeof(val) <: ObjectIdDict)
    println("val = ", val)
  end

#=
  if typeof(val) <: ObjectIdDict
    println("keys(val) = ", keys(val))
  end
=#
#=  
  if typeof(val) == Function
    println("methods = \n", methods(val))
    print("\n")
  end
=#
  println(f, string(i))
  print("\n")
end

return nothing
end

f = open("src/index.md", "w")
println(f, "```@docs")
write_names(f)
println(f, "```")
close(f)
