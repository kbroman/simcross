## sim_cross.R

# meiosis_sub <- 
#
# Simulate the locations of crossovers on a meiotic product
# via the chi-square model (m=0 corresponds to no interference)
# and potentially with an obligate chiasma.
#
#' @importFrom stats uniroot rpois qpois dpois rbinom runif
#
meiosis_sub <-
function(L, m=10, obligate.chiasma=TRUE)
{
  if(obligate.chiasma) { # adjust mean no. chiasmata
    if(L <= 50) stop("L must be > 50 cM")
    if(m==0) f <- function(Lstar,f.L,f.m=0) f.L-Lstar/(1-exp(-Lstar/50)) 
    else {
      f <- function(Lstar,f.L,f.m=0)
        {
          lambdastar <- Lstar/50*(f.m+1)
          temp <- lambdastar
          for(i in 1:length(temp))
            temp[i] <- sum(exp(-lambdastar[i] + (0:f.m)*log(lambdastar[i])-
                               lgamma((0:f.m)+1)) * (f.m+1-(0:f.m))/(f.m+1))
          f.L - Lstar/(1-temp)
        }
    }

    Lstar <- uniroot(f,c(1e-5,L+1),f.m=m,f.L=L)$root
  }
  else Lstar <- L

  if(m==0) { # no interference
    if(!obligate.chiasma)  # no obligate chiasma
      n.xo <- rpois(1,Lstar/100)
    else {
      up <- qpois(1e-14,Lstar/50,lower.tail=FALSE)
      p <- dpois(1:up,Lstar/50)/ppois(0,Lstar/50)
      n.chi <- sample(1:up,1,prob=p)
      n.xo <- rbinom(1,n.chi,0.5)
    }
    if(n.xo==0) xo <- NULL
    else xo <- sort(runif(n.xo,0,L))
  }
  else { # chi-square model
    n.chi <- 0
    while(n.chi == 0) {
      n.pts <- rpois(1,Lstar/50*(m+1))
      first <- sample(1:(m+1),1)
      if(first <= n.pts || !obligate.chiasma) n.chi <- 1
    }
    if(first > n.pts)
      xo <- NULL
    else {
      pt.loc <- sort(runif(n.pts,0,L))
      chi.loc <- pt.loc[seq(first,length(pt.loc),by=m+1)]
      n.xo <- rbinom(1,length(chi.loc),0.5)
      if(n.xo==0) xo <- NULL
      else if(length(chi.loc)==1) xo <- chi.loc
      else xo <- sort(sample(chi.loc,n.xo,replace=FALSE))
    }
  }
    
  if(length(xo) == 0) xo <- NULL
  xo
}

# create_parent
#
#' Create a parent object
#'
#' Create a parent object
#'
#' @param L chromosome length in cM
#' @param allele vector of alleles, of length 1 or 2
#' @return a list with two components, for the two chromosomes.
#' Each is a matrix with two rows: locations of crossovers, and
#' the alleles in the segments
#'
#' @keywords datagen
#' @export
#'
#' @examples
#' create_parent(100, 1)
#' create_parent(100, 1:2)
create_parent <-
function(L, allele=1)
{
  if(length(allele) == 1) allele <- rep(allele,2)
  if(length(allele) != 2)
    stop("allele should be of length 1 or 2")
  
  list(mat=rbind(c(0,L),allele[1]),
       pat=rbind(c(0,L),allele[2]))
}


# meiosis
#
#' Simulate meiosis
#'
#' Output a random meiotic product from an input individual.
#'
#' @param parent An individual object, as output by
#' \code{\link{create_parent}} or \code{\link{cross}}
#' @param m interference parameter for chi-square model
#' @param obligate.chiasma If TRUE, simulate meiosis with an
#' obligate chiasma on the four-strand bundle
#' @return A matrix with two rows: locations of crossovers, and the
#' allele in each segment
#'
#' @keywords datagen
#' @export
#'
#' @examples
#' ind <- create_parent(100, 1:2)
#' prod <- meiosis(ind)
meiosis <-
function(parent, m=10, obligate.chiasma=TRUE)
{
  L <- parent$mat[1,ncol(parent$mat)]
  if(abs(parent$pat[1,ncol(parent$pat)] - L) > 1e-13)
    stop("There is a problem with the parent's data structure.")

  product <- meiosis_sub(L, m, obligate.chiasma)
  a <- sample(1:2,1)

  if(length(product)==0) return(parent[[a]])

  else {
    for(i in 1:length(product)) {
      if(i == 1) 
        result <- parent[[a]][,parent[[a]][1,]<product[1],drop=FALSE]
      else {
        temp <- parent[[a]][1,]>=product[i-1] & parent[[a]][1,]<product[i]
        result <- cbind(result,parent[[a]][,temp])
      }
      u <- parent[[a]][2,parent[[a]][1,]>=product[i]]
      result <- cbind(result,c(product[i],u[1]))
      a <- 3-a
    }
    temp <- parent[[a]][1,]>=product[length(product)]
    result <- cbind(result,parent[[a]][,temp])
  }

  # clean out excess stuff in the result
  if(ncol(result)>2) {
    keep <- rep(TRUE,ncol(result))
    for(i in 2:(ncol(result)-1)) 
      if(result[2,i] == result[2,i+1])
        keep[i] <- FALSE
  }
#  print(result)
  result[,keep,drop=FALSE]
}


# cross
#
#' Cross two individuals
#'
#' Simulate the cross of two individuals to create a
#' single progeny
#'
#' @param mom An individual object, as produced by \code{\link{create_parent}}
#' or \code{\link{cross}}
#' @param dad An individual object, as produced by \code{\link{create_parent}}
#' or \code{\link{cross}}
#' @param m interference parameter for chi-square model
#' @param obligate.chiasma If TRUE, simulate meiosis with an
#' obligate chiasma on the four-strand bundle
#' @param xchr If TRUE, simulate X chromosome
#' @param male If TRUE, simulate a male (matters only if \code{xchr=TRUE})
#' @return A matrix with two rows: locations of crossovers, and the
#' allele in each segment
#' 
#' @keywords datagen
#' @export
#'
#' @examples
#' mom <- create_parent(100, 1:2)
#' dad <- create_parent(100, 1:2)
#' child <- cross(mom, dad)
cross <-
function(mom, dad, m=10, obligate.chiasma=TRUE, xchr=FALSE, male=FALSE)
{
  if(!xchr) {
    return(list(mat=meiosis(mom,m,obligate.chiasma),
                pat=meiosis(dad,m,obligate.chiasma)))
  }
  else {
    if(male)
      return(list(mat=meiosis(mom,m,obligate.chiasma),
                  pat=dad$pat))
    else
      return(list(mat=meiosis(mom,m,obligate.chiasma),
                  pat=dad$mat))
  }
}

# these functions aren't very good
ri2 <-
function(L,n.gen=20,m=10,obligate.chiasma=TRUE)
{
  f1 <- create_parent(L,c(1,2))
  par1 <- cross(f1,f1,m,obligate.chiasma)
  par2 <- cross(f1,f1,m,obligate.chiasma)
  for(i in 1:n.gen) {
    c1 <- cross(par1,par2,m,obligate.chiasma)
    c2 <- cross(par1,par2,m,obligate.chiasma)
    par1 <- c1
    par2 <- c2
  }
  par1
}

ri8 <-
function(L,n.gen=20,m=10,obligate.chiasma=TRUE)
{
  f1a <- create_parent(L,c(1,2))
  f1b <- create_parent(L,c(3,4))
  f1c <- create_parent(L,c(5,6))
  f1d <- create_parent(L,c(7,8))
  par1 <- cross(f1a,f1b,m,obligate.chiasma)
  par2 <- cross(f1c,f1d,m,obligate.chiasma)
  if(length(n.gen)==1) {
    for(i in 1:(n.gen+1)) {
      c1 <- cross(par1,par2,m,obligate.chiasma)
      c2 <- cross(par1,par2,m,obligate.chiasma)
      par1 <- c1
      par2 <- c2
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
        c2 <- cross(par1,par2,m,obligate.chiasma)
        par1 <- c1
        par2 <- c2
      }
      result[[j-1]] <- par1
    }
    return(result)
  }
        
}
    
    
# where_het
#
#' Find heterozygous regions
#'
#' Find regions of heterozygosity in an individual
#'
#' @param ind An individual object, as output be \code{\link{create_parent}}
#' or \code{\link{cross}}
#' @return A matrix with two columns; each row indicates the start
#' and end of a region where the individual is heterozygous
#'
#' @export
#' @examples
#' mom <- create_parent(100, 1:2)
#' dad <- create_parent(100, 1:2)
#' child <- cross(mom, dad)
#' where_het(child)
where_het <-
function(ind)
{
  if(ncol(ind$mat)==ncol(ind$pat) && all(ind$mat == ind$pat)) {
    return(NULL)
  }
  u <- sort(unique(c(ind$mat[1,],ind$pat[1,])))
  het <- NULL
  for(i in 2:length(u)) {
    mat <- ind$mat[,ind$mat[1,] >= u[i],drop=FALSE]
    mat <- mat[2,1]

    pat <- ind$pat[,ind$pat[1,] >= u[i],drop=FALSE]
    pat <- pat[2,1]

    if(mat!=pat) { # heterozygous
      if(is.null(het)) het <- cbind(u[i-1],u[i])
      else het <- rbind(het,c(u[i-1],u[i]))
    }
  }

  # clean up
  if(nrow(het) > 1) {
    keep <- rep(TRUE,nrow(het))
    for(j in 2:nrow(het)) {
      if(het[j,1] == het[j-1,2]) {
        het[j,1] <- het[j-1,1]
        keep[j-1] <- FALSE
      }
    }
    het <- het[keep,,drop=FALSE]
  }
  het
}
