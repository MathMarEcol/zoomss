# This function checks the model output against a previous run to see where any differences are.

fZooMSS_CheckIdent <- function(out,out_old){
  nm <- names(out_old$model$param)
  for(a in 1:length(nm)){
    b <- eval(parse(text = paste0("identical(out_old$model$param$",nm[a],",out$model$param$",nm[a],")")))
    if (b == FALSE &
        eval(parse(text = paste0("!is.null(out_old$model$param$",nm[a],")"))) &
        eval(parse(text = paste0("!is.null(out$model$param$",nm[a],")")))
    ){
      print(paste0("Param ",nm[a]," is different."))
    }#else{print(paste0("Param ",nm[a]," is equal."))}
  }

  nm2 <- names(out_old$model)
  for(a in 1:length(nm2)){
    b <- eval(parse(text = paste0("identical(out_old$model$",nm2[a],",out$model$",nm2[a],")")))
    if (b == FALSE &
        eval(parse(text = paste0("!is.null(out_old$model$",nm2[a],")"))) &
        eval(parse(text = paste0("!is.null(out$model$",nm2[a],")")))
    ){
      print(paste0("Model ",nm2[a]," is different."))
    }#else{print(paste0("Model Setup ",nm2[a]," is equal."))}
  }

}