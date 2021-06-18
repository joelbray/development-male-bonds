#Model Code

#Title: Immature Male Chimpanzees’ (Pan troglodytes schweinfurthii) Social Relationships with Adult Males, but Not Peers, Persist into Adulthood
#Authors: Joel Bray, Carson Murray, Ian Gilby, Maggie Stanton

############################################################
#Immature Associations with Adult Males
############################################################

#Dyadic association during immaturity as a predictor of adult-adult association

m.strength.fm.dai <- map2stan(
  alist(
    
    dai ~ dzagamma2( p, mu , scale ),
    
    logit(p) ~ ap + ap_id[id1] + ap_id[id2] + ap_dyad[dyad] + ap_year[year] +
      (bp_rank + bp_rank_id[id1])*rank_id1 + (bp_rank + bp_rank_id[id2])*rank_id2 + 
      (bp_age + bp_age_id[id1])*age_id1 + (bp_age + bp_age_id[id2])*age_id2 + 
      (bp_rankDiff + bp_rankDiff_dyad[dyad])*rank_diff +
      bp_ageDiff*age_diff +
      bp_kin*maternal_kin + 
      bp_kinXageDiff*maternal_kin*age_diff +
      bp_kinXrankDiff*maternal_kin*rank_diff +
      bp_dev_dai4*dev_dai4 +
      bp_days*days,
    
    log(mu) ~ am + am_id[id1] + am_id[id2] + am_dyad[dyad] + am_year[year] +
      (bm_rank + bm_rank_id[id1])*rank_id1 + (bm_rank + bm_rank_id[id2])*rank_id2 + 
      (bm_age + bm_age_id[id1])*age_id1 + (bm_age + bm_age_id[id2])*age_id2 + 
      (bm_rankDiff + bm_rankDiff_dyad[dyad])*rank_diff +
      bm_ageDiff*age_diff +
      bm_kin*maternal_kin + 
      bm_kinXageDiff*maternal_kin*age_diff +
      bm_kinXrankDiff*maternal_kin*rank_diff +
      bm_dev_dai4*dev_dai4 +
      bm_days*days,
    
    c(ap_year, am_year)[year] ~ dmvnormNC( sigma_year, Rho_year ),
    c(ap_id, am_id, bp_rank_id, bm_rank_id, bp_age_id, bm_age_id)[id1] ~ dmvnormNC( sigma_id, Rho_id ),
    c(ap_dyad, am_dyad, bp_rankDiff_dyad, bm_rankDiff_dyad)[dyad] ~ dmvnormNC( sigma_dyad, Rho_dyad ),
    ap ~ dnorm(0,2),
    am ~ dnorm(0,2),
    c(bp_rank, bm_rank, bp_age, bm_age, bp_kin, bm_kin, bp_rankDiff, bm_rankDiff, bp_ageDiff, bm_ageDiff, bp_kinXageDiff, bm_kinXageDiff, bp_kinXrankDiff, bm_kinXrankDiff, bp_dev_dai4, bm_dev_dai4, bp_days, bm_days) ~ dnorm(0,2),
    c(sigma_year, sigma_id, sigma_dyad) ~ dexp(1), 
    c(Rho_year, Rho_id, Rho_dyad) ~ dlkjcorr(3),
    scale ~ dexp(1)
    
  ),
  
  data=data.fm.dai, cores=3 , chains=3 , warmup=2000, iter=4000, constraints=list(scale="lower=0") , control=list(adapt_delta=0.9), WAIC=TRUE
  
)

#Dyadic association during immaturity as a predictor of adult-adult grooming

m.strength.fm.grm <- map2stan(
  alist(
    
    grm ~ dzagamma2( p, mu , scale ),
    
    logit(p) ~ ap + ap_id[id1] + ap_id[id2] + ap_dyad[dyad] + ap_year[year] +
      (bp_rank + bp_rank_id[id1])*rank_id1 + (bp_rank + bp_rank_id[id2])*rank_id2 + 
      (bp_age + bp_age_id[id1])*age_id1 + (bp_age + bp_age_id[id2])*age_id2 + 
      (bp_rankDiff + bp_rankDiff_dyad[dyad])*rank_diff +
      bp_ageDiff*age_diff +
      bp_kin*maternal_kin + 
      bp_kinXageDiff*maternal_kin*age_diff +
      bp_kinXrankDiff*maternal_kin*rank_diff +
      bp_dev_dai4*dev_dai4 +
      bp_days*days,
    
    log(mu) ~ am + am_id[id1] + am_id[id2] + am_dyad[dyad] + am_year[year] +
      (bm_rank + bm_rank_id[id1])*rank_id1 + (bm_rank + bm_rank_id[id2])*rank_id2 + 
      (bm_age + bm_age_id[id1])*age_id1 + (bm_age + bm_age_id[id2])*age_id2 + 
      (bm_rankDiff + bm_rankDiff_dyad[dyad])*rank_diff +
      bm_ageDiff*age_diff +
      bm_kin*maternal_kin + 
      bm_kinXageDiff*maternal_kin*age_diff +
      bm_kinXrankDiff*maternal_kin*rank_diff +
      bm_dev_dai4*dev_dai4 +
      bm_days*days,
    
    c(ap_year, am_year)[year] ~ dmvnormNC( sigma_year, Rho_year ),
    c(ap_id, am_id, bp_rank_id, bm_rank_id, bp_age_id, bm_age_id)[id1] ~ dmvnormNC( sigma_id, Rho_id ),
    c(ap_dyad, am_dyad, bp_rankDiff_dyad, bm_rankDiff_dyad)[dyad] ~ dmvnormNC( sigma_dyad, Rho_dyad ),
    ap ~ dnorm(0,2),
    am ~ dnorm(0,2),
    c(bp_rank, bm_rank, bp_age, bm_age, bp_kin, bm_kin, bp_rankDiff, bm_rankDiff, bp_ageDiff, bm_ageDiff, bp_kinXageDiff, bm_kinXageDiff, bp_kinXrankDiff, bm_kinXrankDiff, bp_dev_dai4, bm_dev_dai4, bp_days, bm_days) ~ dnorm(0,2),
    c(sigma_year, sigma_id, sigma_dyad) ~ dexp(1), 
    c(Rho_year, Rho_id, Rho_dyad) ~ dlkjcorr(3),
    scale ~ dexp(1)
    
  ),
  
  data=data.fm.grm, cores=3 , chains=3 , warmup=2000, iter=4000, constraints=list(scale="lower=0") , control=list(adapt_delta=0.9), WAIC=TRUE
  
)

############################################################
#Associations Among Immature Peers
############################################################

#Dyadic association during immaturity as a predictor of adult-adult association

m.strength.ff.dai <- map2stan(
  alist(
    
    dai ~ dzagamma2( p, mu , scale ),
    
    logit(p) ~ ap + ap_id[id1] + ap_id[id2] + ap_dyad[dyad] + ap_year[year] +
      (bp_rank + bp_rank_id[id1])*rank_id1 + (bp_rank + bp_rank_id[id2])*rank_id2 + 
      (bp_age + bp_age_id[id1])*age_id1 + (bp_age + bp_age_id[id2])*age_id2 + 
      (bp_rankDiff + bp_rankDiff_dyad[dyad])*rank_diff +
      bp_ageDiff*age_diff +
      bp_dev_dai4*dev_dai4 +
      bp_days*days,
    
    log(mu) ~ am + am_id[id1] + am_id[id2] + am_dyad[dyad] + am_year[year] +
      (bm_rank + bm_rank_id[id1])*rank_id1 + (bm_rank + bm_rank_id[id2])*rank_id2 + 
      (bm_age + bm_age_id[id1])*age_id1 + (bm_age + bm_age_id[id2])*age_id2 + 
      (bm_rankDiff + bm_rankDiff_dyad[dyad])*rank_diff +
      bm_ageDiff*age_diff +
      bm_dev_dai4*dev_dai4 +
      bm_days*days,
    
    c(ap_year, am_year)[year] ~ dmvnormNC( sigma_year, Rho_year ),
    c(ap_id, am_id, bp_rank_id, bm_rank_id, bp_age_id, bm_age_id)[id1] ~ dmvnormNC( sigma_id, Rho_id ),
    c(ap_dyad, am_dyad, bp_rankDiff_dyad, bm_rankDiff_dyad)[dyad] ~ dmvnormNC( sigma_dyad, Rho_dyad ),
    ap ~ dnorm(0,2),
    am ~ dnorm(0,2),
    c(bp_rank, bm_rank, bp_age, bm_age, bp_rankDiff, bm_rankDiff, bp_ageDiff, bm_ageDiff, bp_dev_dai4, bm_dev_dai4, bp_days, bm_days) ~ dnorm(0,2),
    c(sigma_year, sigma_id, sigma_dyad) ~ dexp(1), 
    c(Rho_year, Rho_id, Rho_dyad) ~ dlkjcorr(3),
    scale ~ dexp(1)
    
  ),
  
  data=data.ff.dai, cores=3 , chains=3 , warmup=2000, iter=4000, constraints=list(scale="lower=0") , control=list(adapt_delta=0.9), WAIC=TRUE
  
)

#Dyadic association during immaturity as a predictor of adult-adult grooming

m.strength.ff.grm <- map2stan(
  alist(
    
    grm ~ dzagamma2( p, mu , scale ),
    
    logit(p) ~ ap + ap_id[id1] + ap_id[id2] + ap_dyad[dyad] + ap_year[year] +
      (bp_rank + bp_rank_id[id1])*rank_id1 + (bp_rank + bp_rank_id[id2])*rank_id2 + 
      (bp_age + bp_age_id[id1])*age_id1 + (bp_age + bp_age_id[id2])*age_id2 + 
      (bp_rankDiff + bp_rankDiff_dyad[dyad])*rank_diff +
      bp_ageDiff*age_diff +
      bp_dev_dai4*dev_dai4 +
      bp_days*days,
    
    log(mu) ~ am + am_id[id1] + am_id[id2] + am_dyad[dyad] + am_year[year] +
      (bm_rank + bm_rank_id[id1])*rank_id1 + (bm_rank + bm_rank_id[id2])*rank_id2 + 
      (bm_age + bm_age_id[id1])*age_id1 + (bm_age + bm_age_id[id2])*age_id2 + 
      (bm_rankDiff + bm_rankDiff_dyad[dyad])*rank_diff +
      bm_ageDiff*age_diff +
      bm_dev_dai4*dev_dai4 +
      bm_days*days,
    
    c(ap_year, am_year)[year] ~ dmvnormNC( sigma_year, Rho_year ),
    c(ap_id, am_id, bp_rank_id, bm_rank_id, bp_age_id, bm_age_id)[id1] ~ dmvnormNC( sigma_id, Rho_id ),
    c(ap_dyad, am_dyad, bp_rankDiff_dyad, bm_rankDiff_dyad)[dyad] ~ dmvnormNC( sigma_dyad, Rho_dyad ),
    ap ~ dnorm(0,2),
    am ~ dnorm(0,2),
    c(bp_rank, bm_rank, bp_age, bm_age, bp_rankDiff, bm_rankDiff, bp_ageDiff, bm_ageDiff, bp_dev_dai4, bm_dev_dai4, bp_days, bm_days) ~ dnorm(0,2),
    c(sigma_year, sigma_id, sigma_dyad) ~ dexp(1), 
    c(Rho_year, Rho_id, Rho_dyad) ~ dlkjcorr(3),
    scale ~ dexp(1)
    
  ),
  
  data=data.ff.grm, cores=3 , chains=3 , warmup=2000, iter=4000, constraints=list(scale="lower=0") , control=list(adapt_delta=0.9), WAIC=TRUE
  
)