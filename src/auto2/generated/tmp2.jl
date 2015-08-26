function replace_symbol(ex, sym_old::Symbol, sym_new)
# do a recursive descent replace one symbol with another

#  @assert typeof(ex) == Expr
  println("receiving expression ", ex)
  println("typeof(ex) = ", typeof(ex))
  println("typeof(sym_old) = ", typeof(sym_old))
  println("typeof(sym_new) = ", typeof(sym_new))

  if typeof(ex) == Symbol  # if this is a symbol
      println("  found symbol ", ex)
      if ex == sym_old
        println("  performing replacement ", ex, " with ", sym_new)
        return sym_new
      else
        println("  returning original symbol")
        return ex
      end
        
  elseif typeof(ex)  == Expr  # keep recursing
  
    for i=1:length(ex.args)
#        println("  processing sub expression ", ex.args[i])
         println("  recursing expression ", ex.args[i]) 
         ex.args[i] =  replace_symbol(ex.args[i], sym_old, sym_new)
    end  # end loop over args
  else # we don't know/care what this expression is
    println("  not modify unknown expression ", ex)
    return ex
  end  # end if ... elseifa

  println("Warning, got to end of replace_symbol")
  println("ex = ", ex)
  println("typeof(ex) = ", typeof(ex))

  return ex

end

ex = :(PetscViewer{S <: PetscScalar})
sym_old = :PetscViewer
sym_new = :(PetscViewer{S <: PetscScalar})

println("old expression = ", ex)
replace_symbol(ex, sym_old, sym_new)

println("new expression = ", ex)
