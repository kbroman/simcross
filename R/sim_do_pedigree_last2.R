sim_do_pedigree_last2 <- function(ngen=12, npairs=30, nkids_per=5, n.ccgen=15,
                                  design=c("nosib", "random"),
                                  nsample_ngen=150)
{
  design <- match.arg(design)
  npairs_la2 <- ceiling(nsample_ngen/nkids_per)
  nkids_la2 <- ceiling(npairs_la2*2/npairs)
  
  ## last but 2 generation
  ped <- sim_do_pedigree(ngen = ngen-1, npairs = npairs,
                         ccgen = rep(n.ccgen, npairs),
                         nkids_per=nkids_la2, design=design)
  id <- ped[,"id"]
  mom <- ped[,"mom"]
  dad <- ped[,"dad"]
  sex <- ped[,"sex"]
  gen <- ped[,"gen"]
  do <- ped[,"do"]
  
  ## last generation
  n.last <- npairs_la2*nkids_per
  kids <- 1:(n.last)+max(id)

  wh <- which(ped[,"gen"] == ngen-1 & ped[, "sex"] == 0 & ped[, "do"] == 1)
  moms <- ped[wh, "mom"]
  wh <- which(ped[,"gen"] == ngen-1 & ped[, "sex"] == 1 & ped[, "do"] == 1)
  dads <- ped[wh, "dad"]
  rownames(ped) <- ped[, "id"]

  while(design=="nosib") { # sample until no sibs
    dads <- sample(dads)
    if(all(ped[dads, "dad"] != ped[moms, "dad"]))
        break
  }

  wh <- sort(sample(n.last, size=nsample_ngen))
  id <- c(id, kids[wh])
  mom <- c(mom, rep(moms, each=nkids_per)[wh])
  dad <- c(dad, rep(dads, each=nkids_per)[wh])
  sex <- c(sex, rep_len(c(0,1), length.out=n.last)[wh])
  gen <- c(gen, rep(ngen, n.last)[wh])
  do <- c(do, rep(1, length(wh)))
  
  result <- cbind(id=id, mom=mom, dad=dad, sex=sex, gen=gen, do=do)
  storage.mode(result) <- "integer"
  result
}
