## sim_ril.R

# these functions aren't very good
ri2 <-
function(L,n.gen=20,m=10,obligate.chiasma=TRUE, selfing=FALSE)
{
  f1 <- create_parent(L,c(1,2))
  par1 <- cross(f1,f1,m,obligate.chiasma)

  if(selfing) par2 <- par1
  else par2 <- cross(f1,f1,m,obligate.chiasma)

  for(i in 1:n.gen) {
    c1 <- cross(par1,par2,m,obligate.chiasma)
    if(selfing)
        par1 <- par2 <- c1
    else {
        c2 <- cross(par1,par2,m,obligate.chiasma)
        par1 <- c1
        par2 <- c2
    }
  }
  par1
}

ri8 <-
function(L,n.gen=20,m=10,obligate.chiasma=TRUE, selfing=FALSE)
{
  f1a <- create_parent(L,c(1,2))
  f1b <- create_parent(L,c(3,4))
  f1c <- create_parent(L,c(5,6))
  f1d <- create_parent(L,c(7,8))
  par1 <- cross(f1a,f1b,m,obligate.chiasma)

  if(selfing) par2 <- par1
  else par2 <- cross(f1c,f1d,m,obligate.chiasma)

  if(length(n.gen)==1) {
    for(i in 1:(n.gen+1)) {
      c1 <- cross(par1,par2,m,obligate.chiasma)
      if(selfing)
          par1 <- par2 <- c1
      else {
          c2 <- cross(par1,par2,m,obligate.chiasma)
          par1 <- c1
          par2 <- c2
      }
    }
    return(par1)
  }
  else {
    result <- vector("list",length(n.gen))
    names(result) <- n.gen
    n.gen <- c(-1,n.gen)
    for(j in 2:length(n.gen)) {
      for(i in (n.gen[j-1]+2):(n.gen[j]+1)) {
        c1 <- cross(par1,par2,m,obligate.chiasma)
        if(selfing)
            par1 <- par2 <- c1
        else {
            c2 <- cross(par1,par2,m,obligate.chiasma)
            par1 <- c1
            par2 <- c2
        }
      }
      result[[j-1]] <- par1
    }
    return(result)
  }

}

