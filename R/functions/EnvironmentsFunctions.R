DatLongToEnvironment <- function(dat.long){
  # split dat.long into environment with index corresponding to gene
  my.env <- new.env()
  dat.long %>%
    group_by(gene) %>%
    do(AddToEnviron(., env = my.env, cname = "gene"))
  return(my.env)
}

AddToEnviron <- function(dat.sub, env, cname="gene"){
  # add data split into environment, name by column name 
#   print(head(dat.sub))
  jname <- as.character(dat.sub[[cname]][[1]])
#   print(jname)
  env[[jname]] <- dat.sub
}