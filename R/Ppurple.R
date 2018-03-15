#' @import GenomicRanges
#' @import data.table
#' @import Rsamtools
#' @import gUtils
#' @import bamUtils



#' @name ssegment
#' @rdname internal
#' @title Internal function utilizing DNAcopy to segment a coverage profile
#' @description 
#' 
#' Internal function utilizing DNAcopy to segment a coverage profile
#' @param tcov GRanges of binned genome-wide coverage
#' @return GRanges of genomewise segments of piecewise constant coverage
#' @keywords internal
#' @author Marcin Imielinski
ssegment = function(cov, verbose = verbose){
    new.sl = seqlengths(cov)
    ix = which(!is.na(cov$y))
    if (verbose)
      pmessage('sending ', length(ix), ' segments to DNAcopy')
    cna = CNA(log(cov$y[ix]), as.character(seqnames(cov))[ix], start(cov)[ix], data.type = 'logratio')
    gc()
    seg = segment(smooth.CNA(cna), alpha = 1e-5, verbose = 0)
    out = seg2gr(seg$out, new.sl) ## remove seqlengths that have not been segmented
    out = gr.fix(out, new.sl, drop = T)
    if (verbose)
      pmessage('\t ..finished segmentation')
    return(out)
}


#' @name log.sum.exp
#' @rdname internal
#' @title Internal function doing log sum exp
#' @description 
#' 
#' Internal function
#' @param x vector of log probabilities
#' @return log.sum.exp of input
#' @keywords internal
#' @author Marcin Imielinski
log.sum.exp<- function(x){
    offset <- max(x)
    log(sum(exp(x - offset))) + offset
}

#' @name llpois
#' @rdname internal
#' @title Internal function doing product of poisson log likelihoods
#' @description 
#' 
#' Internal function
#' @param slfx numeric vector of sum log factorials of data
#' @param sx numeric vector of sum of data
#' @param n numeric vector of data ns
#' @param lambda numeric vector of lambda parameters
#' @return log poisson joint likelihood
#' @keywords internal
#' @author Marcin Imielinski
llpois = function(slfx, sx, n, lambda) -n*lambda - slfx + sx*log(lambda)

#' @name lfactorial
#' @rdname internal
#' @title Internal function doing approximate log factorial
#' @description 
#' 
#' Internal function
#' @param x integer vector for whiich to compute approximate log factorial
#' @return log factorial
#' @keywords internal
#' @author Marcin Imielinski
lfactorial = function(x) {
  y = suppressWarnings(log(factorial(x)));
  if (any(ix <- is.infinite(y)))
    y[ix] = x[ix]*log(x[ix]) - x[ix] + 1   ## use stirling's approx
  return(y)
}

#' @name llnorml
#' @rdname internal
#' @title Internal function doing product of normal log likelihoods from sufficient stats
#' @description 
#' 
#' Internal function
#' @param x numeric vector observed mean
#' @param sos numeric vector observed sum of squares
#' @param n interger vector number of observations
#' @param mu mean param of normal distribution
#' @param sd sigma param of normal distribution
#' 
#' @author Marcin Imielinski
llnorm = function(x, sos, n, mu, sd)
  -n*log(sqrt(2*3.141593)*sd)-(n*(mu-x)^2 + sos)/(2*sd^2)


#' @name hapseg
#' @rdname internal
#' @title Internal function doing simple implementation of Carter, Getz 2011 hapseg
#' @description 
#'
#' Phases grouped (i.e. presegmented) hets and computes cluster centers (lambdas)
#' and expected sufficient statistics for joint poisson
#' per group (i.e. segment)
#' 
#' Internal function
#' @param hets data.table of genome wide het sites with $alt, $ref and field $j that groups the hets into segments 
#' @keywords internal
#' @author Marcin Imielinski
hapseg = function(hets, verbose = FALSE)
{
  ll = data.table(j = unique(hets$j), ll = -1e100, diff = 1e100)
  setkey(ll, j)

  ## replicating hets to include two phases (TRUE = alt high, ref low, FALSE = alt low, ref high)
  dat = hets[, .(j = rep(j,2), alt = rep(alt,2), ref = rep(ref,2), phase = rep(c(TRUE, FALSE), each = .N))][, i := rep(1:(.N/2),2), by = j]

  ## populate initial value for posterior probability of cluster membership, which we set to probability 1 based on ? alt>ref
  ## i.e. p is 1 iff phase is TRUE and alt>ref or phase is FALSE and ref>=alt, otherwise 0
#  dat[, p := ifelse(phase, as.numeric(alt>ref), as.numeric(alt<=ref))]  
  p0 = 1
  dat[, p := ifelse(phase, ifelse(alt>ref, p0, 1-p0), ifelse(as.numeric(alt<=ref), p0, 1-p0))]
  setkeyv(dat, c("j", "i"))
  tol = 1
  iteration = 1
  while (any(ll$diff > tol))
  {
    if (verbose)
      pmessage('Hapseg iteration ', iteration, ':')

    ##
    ## M STEP: estimating poisson cluster centers
    ##
    lambda = dat[, .(
      lhigh = sum(p[phase]*alt[phase] + p[!phase]*ref[!phase])/(.N/2), ## dividing by number of hets in each j slice (i.e. which .N/2 since we "replicated")
      llow = sum(p[!phase]*alt[!phase] + p[phase]*ref[phase])/(.N/2)
    ), keyby = .(j)]
    dat$lhigh = lambda[.(dat$j), lhigh] ## merge lhigh / llow estimates back to dat which tracks het phases across all hets
    dat$llow = lambda[.(dat$j), llow]

    ##
    ## E STEP
    ##
    dat[phase == TRUE, ll := log(0.5) + dpois(alt, lhigh, log = TRUE) + dpois(ref, llow, log = TRUE)] ## assume 0.5 prior of + vs - phase 
    dat[phase == FALSE, ll := log(0.5) + dpois(alt, llow, log = TRUE) + dpois(ref, lhigh, log = TRUE)]
    dat[, lsei := log.sum.exp(ll), by = .(j, i)] ## sum over phases (duplicating lsei across phases)
    dat[, lsej := sum(lsei[phase]), by = .(j)] ## sum over segment (don't overcount hence select for phase == TRUE)
    dat[, p := exp(ll-lsei), by = .(j, i)] ## p is the new posterior probability of phase in each het i 

    ll.new = dat[i == 1 & phase == TRUE, .(ll = lsej), keyby = j]
    ll.new$ll.old =  ll[.(ll.new$j), ll]
    ll.new[, diff := ll-ll.old]
    ll = ll.new
    iteration = iteration + 1
    if (verbose)
      pmessage('\tLL diff: ', max(ll$diff), ' tol:', tol)
  }

  ## final go around compute our expected sufficient statistics
  ## these summary stats will be used by llpois
  segs = dat[, .(    
    high.lambda = sum(p[phase]*alt[phase] + p[!phase]*ref[!phase])/(.N/2),
    low.lambda = sum(p[!phase]*alt[!phase] + p[phase]*ref[phase])/(.N/2),
    high.sx = sum(p[phase]*alt[phase] + p[!phase]*ref[!phase]),
    low.sx = sum(p[!phase]*alt[!phase] + p[phase]*ref[phase]),
    high.slfx = sum(p[phase]*lfactorial(alt[phase]) + p[!phase]*lfactorial(ref[!phase])),
    low.slfx = sum(p[!phase]*lfactorial(alt[!phase]) + p[phase]*lfactorial(ref[phase])),
    n = .N/2
  ), keyby = .(j)]

  return(list(segs = segs, phased = dat[, .(j, i, alt, ref, phase, p)]))
}
                                                                                                                                       
#' @name ppurple
#' @title Probabilistic purity ploidy estimation
#' @description
#' Computes posterior probability of purity ploidy for data with coverage cov (granges with value $y specifying coverage)
#' hets with $alt, $ref specifying tumor and normal alt and ref allelic counts, and (optional) segments via EM algorithm. 
#' 
#' If segments not specified, then it is inferred from the segmentation of cov via CBS (using DNA copy)
#'
#' EM begins with initial (user specified) grid of values and then refines at 10 fold (or value of refine arg, >1)
#' around the modes, meaning that it reruns the calculation with a finer grid completely surrounding the modes 
#' to output the final solution which is a data.table mapping $purity and $ploidy values
#' to a posterior probability $p
#' @author Marcin Imielinski
#' @param cov GRanges or data.table of genome wide coverage tiles with field $y specifying normalized coverage
#' @param hets GRanges or data.table of germline hets with fields $alt, $ref specifying alt and ref counts of hets in tumor
#' @param segs GRanges of pre-computed segments (optional, if NULL these will be computed via DNA copy of cov)
#' @param purities numeric vector of ploidies to sweep in grid (default from 1 to 5, 0.2 increment)
#' @param ploidies numeric vector of purities to sweep in grid (default from 0 to 1, 0.1 increment)
#' @param refine integer scalar of how many fold refinement of purities x ploidies grid to perform after initial run (default = 10, careful to not make too big)
#' @param ignore.sex logical flag whether to throw out sex chromosomes (X, Y, chrX, chrY)
#' @param min.bins minimum number of coverage bins in a seg for processing (default = 0)
#' @param min.hets minimum number of hets in a seg for downstream processing (default = 0)
#' @param K integer scalar number of copy states to model (default 20)
#' @param verbose logical flag
#' @param mc.cores integer scalar to parallelize (default = 1)
#' @param numchunks how many chunks to parallelize over (default = mc.cores)
#' @export
#' @examples
#'
#' cov = fread(system.file("extdata", "coverage.csv", package = "Ppurple"))
#' hets = fread(system.file("extdata", "hets.csv", package = "Ppurple"))
#' segs = fread(system.file("extdata", "segs.csv", package = "Ppurple"))
#'
#' pp = ppurple(cov = cov, hets = hets, verbose = TRUE)
#' pp = ppurple(cov = sample(dt2gr(cov), 10000), hets = sample(dt2gr(hets), 10000), segs = segs, verbose = TRUE)
#' 
ppurple = function(cov, hets = NULL, segs = NULL, purities = seq(0, 1.0, 0.1), ploidies = seq(1, 5, 0.2), refine = 10, K = 20, min.bins = 0, min.hets = 0, ignore.sex = FALSE, verbose = FALSE, min.p = 0.0001, mc.cores = 1, binsize = NULL, numchunks = mc.cores)
{
  if (!is.null(hets))
  {
      if (is(hets, 'GRanges'))
        {
          if (length(hets)>0)
            hets = gr2dt(hets)
          else
            hets = NULL
        }
      else if (nrow(hets)==0)      
        hets = NULL

      if (!is.null(hets))
        {
          if (is.null(hets$alt))
            stop('hets should have a column $alt specifying the alt count at each het site')
          
          if (is.null(hets$ref))
            stop('hets should have a column $ref specifying the ref count at each het site')
        }
    }

  if (is.null(cov$y))
      stop('cov should have a column $y specifying the normalized coverage value at each interval')

  if (!inherits(cov, 'GRanges'))
    cov = dt2gr(cov)
  
  if (!inherits(cov, 'GRanges'))
    cov = dt2gr(cov)
  
  if (verbose)
    pmessage('Fitting initial grid of ', length(purities), ' purity and ', length(ploidies), ' ploidy combinations.')

  if (any(iix <- is.infinite(cov$y)))
    cov$y[iix]= NA

  cov = cov[!is.na(cov$y), ]

  ## collapse coverage to "rough" binsize prior to segmentation
  if (!is.null(binsize))
    {
      pmessage('Collapsing coverage to binsize of ', binsize, '.')
      tmp.cov = gr2dt(cov)[!is.infinite(y),][which(y<quantile(y, 0.999, na.rm = TRUE)), ]
      collapsed.cov = tmp.cov[ ,.(y = mean(y, na.rm = TRUE)), by = .(seqnames, start = floor(start/binsize)* binsize+1)]
      collapsed.cov[, end := start + binsize -1]
      cov = dt2gr(collapsed.cov)
    }

  
  if (length(cov)==0)
    stop("No non empty cov's remaining")

  if (is.null(segs))
  {
    if (verbose)
      pmessage('Segments not provided so doing internal segmentation via DNAcopy')
    segs = ssegment(cov, verbose)
  }

  if (!inherits(segs, 'GRanges'))
    segs = dt2gr(segs)
  
  if (!is.null(hets))
    hets$j = gr.match(dt2gr(hets), segs)

  if (ignore.sex)
  {
    if (verbose)
      pmessage('Removing sex chromosomes')
    segs = segs[!(gsub('chr', '', seqnames(segs)) %in% c('X', 'Y')), ]
    if (!is.null(hets))
      hets = hets[!(gsub('chr', '', seqnames)) %in% c('X', 'Y'), ]
    cov = cov[!(gsub('chr', '', seqnames(cov)) %in% c('X', 'Y')), ]
  }

  ## aggregate coverage around segs
  segs.gr = segs;
  cov$j = gr.match(cov, segs.gr)

  segs = gr2dt(cov)[, .(y = mean(y, na.rm = TRUE), sos = sum((y-mean(y, na.rm = TRUE))^2, na.rm = TRUE), nbins = .N),
                    keyby = .(j = j)][nbins>min.bins, ]

  if (verbose & min.bins > 0)
    pmessage('Removing segs with fewer than ', min.bins, ' bins')

  if (!is.null(hets))
  {
    if (verbose)
      pmessage('Running hapseg on hets')
    segs.h = hapseg(hets[!is.na(j),], verbose)$seg
    segs.h$y.high = segs.h$high.lambda
    segs.h$y.low = segs.h$low.lambda
    segs.h$nbins.h = segs.h$n
    if (verbose & min.hets > 0)
      pmessage('Removing het segs with fewer than ', min.hets, ' hets')
    
    segs.h = segs.h[nbins.h>min.hets, ]
    rho.h = mean(c(hets$alt, hets$ref), na.rm = TRUE)
  }
  else
  {
    segs.h = data.table(j = as.numeric(NA), y.high = as.numeric(NA), y.low = as.numeric(NA), nbins.h = as.numeric(NA))
    rho.h = NA
  }

  rho = mean(cov$y, na.rm = TRUE)

  purities = round(sort(unique(purities)), 5) ## gets rid of weird numeric errors with seq
  ploidies = round(sort(unique(ploidies)), 5)

  pp = as.data.table(expand.grid(alpha = purities, tau = ploidies))

  pp.l = split(pp, rep(1:numchunks, nrow(pp)/numchunks*2)[1:nrow(pp)])

  pta = rbindlist(mclapply(pp.l, function(pp) ppemgrid(pp = pp, segs = segs, segs.h = segs.h, rho = rho, rho.h = rho.h, K = K, verbose = verbose)$pta, mc.cores = mc.cores))

  ## renormalize posterior probabilities 
  ll = log.sum.exp(pta$log.px_spat)
  pta[, log.pat_xsp := log.px_spat - ll]
  pta[, pta := exp(log.pat_xsp)]

  if (!is.na(refine))
    if (round(refine)>1)
    {
      refine = round(refine)      
      gridw.pu = min(diff(purities))
      gridw.pl = min(diff(ploidies))

      ##
      ## FIRST refine purities using all ploidies
      ##
      if (verbose)
        pmessage('Refining purities around initially provided ploidies')
      ## draw rectangles around all regions above min.p
      modes = pta[pta>min.p, ]
      modes[, pli := match(tau, ploidies)]
      modes[, pui := match(alpha, purities)]
      pp  = rbindlist(lapply(1:nrow(modes), function(x)
        as.data.table(expand.grid(
          tau = ploidies,
          alpha = modes[x, seq(purities[pmax(pui-1, 1)], purities[pmin(length(purities), pui+1)], gridw.pu/refine)]))))

      pp = pp[!duplicated(cbind(alpha, tau)), ]
      pp.l = split(pp, rep(1:numchunks, nrow(pp)/numchunks*2)[1:nrow(pp)])

      if (verbose)
        pmessage('Refining ', nrow(pp), ' solutions around ', nrow(modes), ' pixels from first round')

      pta = rbindlist(mclapply(pp.l, function(pp) ppemgrid(pp = pp, segs = segs, segs.h = segs.h, rho = rho, rho.h = rho.h, K = K, verbose = verbose)$pta, mc.cores = mc.cores))
      
      ## renormalize posteiror probabilities
      ll = log.sum.exp(pta$log.px_spat)
      pta[, log.pat_xsp := log.px_spat - ll]
      pta[, pta := exp(log.pat_xsp)]

      if (verbose)
        pmessage('Refining purities and ploidies together')

      ## draw rectangles around all regions above min.p
      purities = sort(unique(pta$alpha))
      ploidies = sort(unique(pta$tau))
      gridw.pu = min(diff(purities))
      gridw.pl = min(diff(ploidies))
      modes = pta[rev(order(log.pat_xsp)), ][pta>min.p | (1:.N %in% 1:5), ]
      modes[, pli := match(tau, ploidies)]
      modes[, pui := match(alpha, purities)]
      pp  = rbindlist(lapply(1:nrow(modes), function(x)
        as.data.table(expand.grid(
          tau = modes[x, seq(ploidies[pmax(pli-1, 1)], ploidies[pmin(length(ploidies), pli+1)], gridw.pl/refine)],
          alpha = modes[x, seq(purities[pmax(pui-1, 1)], purities[pmin(length(purities), pui+1)], gridw.pu)]))))

      pp = pp[!duplicated(cbind(alpha, tau)), ]
      pp.l = split(pp, rep(1:numchunks, nrow(pp)/numchunks*2)[1:nrow(pp)])

      if (verbose)
        pmessage('Refining ', nrow(pp), ' solutions around ', nrow(modes), ' pixels from second round')

      pta = rbindlist(mclapply(pp.l, function(pp) ppemgrid(pp = pp, segs = segs, segs.h = segs.h, rho = rho, rho.h = rho.h, K = K, verbose = verbose)$pta, mc.cores = mc.cores))
      
      ## renormalize posteiror probabilities
      ll = log.sum.exp(pta$log.px_spat)
      pta[, log.pat_xsp := log.px_spat - ll]
      pta[, pta := exp(log.pat_xsp)]
    }

  res = pta[rev(order(log.pat_xsp)), .(purity = alpha, ploidy = tau, prob = pta)]
  return(res)
}

#' @name ppemgrid
#' @rdname internal
#' @title ppemgrid internal function
#' @description
#' 
#' Given data.tlable of segs and segs.h populated with summary stats
#' sets up grid across purities and ploidies and K copy states and deploys ppem on it
#' returns posterior probability across purities x ploidy combos given data and optimal
#' params
#' 
#' @author Marcin Imielinski
#' @param purities vector of purities to sweep
#' @param ploidies vector of ploidies to sweep
#' @param pp  (optional) data.table of purity x ploidy x alpha x tau for explicit combos to sweep
#' @param segs data.table of per segment mean coverage data, ns,  sos's across segs
#' @param segs.h data.table of phased cluster centers across segs
#' @param rho  numeric scalar mean total coverage 
#' @param rho.h  numeric scalar mean het count
#' @param k.dist integer scalar optimization constant of how many copy states to explicitly compute probabilities for relative to the MLE value
#' @param K integer scalar total number of copy number states
#' @param verbose logical flag 
#' @keywords internal
ppemgrid = function(purities = NULL, ploidies = NULL, pp = NULL, segs, segs.h, rho = 1, rho.h = 1, k.dist = 3, K = 20, verbose = TRUE)
{

  if (!is.null(pp))
  {
    betagammas = pp
    purities = unique(pp$alpha)
    ploidies = unique(pp$tau)
  }
  else
  {    
    betagammas = as.data.table(expand.grid(alpha = purities, tau = ploidies))
  }

  segs = segs[!is.na(j), ]

  sd0 = sqrt(var(segs$y, na.rm = TRUE))

  if (verbose)
    pmessage('Running ppemgrid with ', length(purities), ' purities [ ', min(purities), ' .. ', max(purities), '] and ', length(ploidies), ' ploidies [ ', min(ploidies), ' .. ', max(ploidies), ' ], rho = ', signif(rho, 3), ', het rho = ', signif(rho.h,3), '.')

  betagammas[, beta := rho*alpha / (2*(1-alpha) + alpha*tau)]
  betagammas[, gamma :=  rho*(2*(1-alpha))/(2*(1-alpha)+alpha*tau)]
  ijk = as.data.table(expand.grid(j = segs$j, i = 1:nrow(betagammas), k = (-1):K))
  segs.grid = cbind(segs[.(ijk$j)], betagammas[ijk$i, ], data.table(k = ijk$k))
  
  set.seed(42)
  parts = unique(c(-Inf, quantile(sample(segs$y, 1e4, replace = TRUE, prob = segs$nbins), na.rm = TRUE, probs = seq(0, 1, length.out = 10)), Inf))

  ## first set total copy number, prune a priori unlikely states to save compute
  segs.grid[, mu := gamma + beta * k]
  segs.grid[, pname := cut(mu, parts)]
  segs.grid[, partition := as.integer(pname)]
  segs.grid[, sd := sd0]
  segs.grid[, kdist := abs(y-mu)/beta]
  setkeyv(segs.grid, c("alpha", "tau", "j", "kdist"))
  segs.grid[, kord := 1:.N, .(alpha, tau, j)]
  ##  segs.grid = segs.grid[k<0 | kdist<=k.dist, ]
  segs.grid = segs.grid[k<0 | kord<k.dist, ] ## pick k.dist y that are closest to mu 
  segs.grid[, sos.k := nbins*(mu-y)^2 + sos] ## we recompute sum of sq (sos.k) around the various cluster centers (mu) using the intra segment sos (sos) and segment mean (y) and num bins (nbins)
  segs.grid[, first := (k == -1)]
  
  ## data.table to expand k by k.high and k low --> note every segment will have a k.low = 0 
  khl = rbind(as.data.table(expand.grid(k = 0:K, k.high = 0:K))[k.high<=k, ][, k.low := k-k.high][k.low<=k.high,],
              data.table(k = -1, k.high = -1, k.low = -1)) ## enumerating all combinations that add up to k
  setkey(khl, k)
  segs.grid = merge(segs.grid, segs.h, by = 'j', all.x = TRUE) ## adding het data to the total coverage
  segs.grid = merge(segs.grid, khl, by = 'k', all.x = TRUE, allow.cartesian = TRUE) ## expanding the copy states by khl (ie to add rows corresponding to high low combos)

  ## formula for beta gamma slightly different for hets 
  segs.grid[, gamma.h := gamma*rho.h/rho]
  segs.grid[, beta.h := 2*beta*rho.h/rho] ## het beta is twice regular beta after sloidy correction
  segs.grid[, mu.high := gamma.h + beta.h * k.high]
  segs.grid[, mu.low := gamma.h + beta.h * k.low]
  segs.grid[, first := (k == -1)]
  segs.grid[mu ==0, mu := mu+1e-6] ## fix the mu == 0 to avoid Nan and -Inf issues
  segs.grid[mu.high ==0, mu.high := mu.high+1e-6] ## fix the mu == 0 to avoid Nan and -Inf issues
  segs.grid[mu.low ==0, mu.low := mu.low+1e-6] ## fix the mu == 0 to avoid Nan and -Inf issues
  setkey(segs.grid, first)

  res = ppem(segs.grid, segs, segs.h, verbose = verbose, sd0 = sd0)
  return(list(pta = res$pta, sd = res$sd, pk = res$pk, pkh = res$pkh))
}


#' @name ppem
#' @rdname internal
#' @title ppem 
#' @description
#'
#' Utility function taking in
#' segs = grid of total copy number segments x alpha x tau x k
#'         --> y, mu, sd, sos, sos.k, nbins
#' segs.h = data.table storing grid of "high" and "low" haplotype segments x alpha x tau x k
#'         --> y.high, y.low, high.slfx, high.sx, low.slfx, low.sx, nbins
#'
#' and returns locally optimal solution with tolerance <= tol or after max_iter iterations
#' initialized with sd0 and copy number prior with tau.sd = 1
#' the output is list with fields, each containing the following data.tables 
#' $pp = posterior probability of purity ploidy combo alpha tau 
#' $segs = grid of segments with final output
#' $segs.h = grid of haplotype segments with final output
#' $pi.k = mixing probability across total copy states for each alpha tau 
#' $pi.kh = mixing probability across haplotype specific copy states for each alpha tau
#' 
#' @author Marcin Imielinski
#' @param segs.grid data.table collating phased het count, mean coverage, sos for every seg x purity x ploidy x copy state combos
#' @param segs data.table of total coverage segs inputted to ppemgrid
#' @param segs.h data.table of het segs inputted to ppemgrid
#' @param sd0 numeric scalar initial sd guess for total coverage clusters
#' @param sd0.k numeric scalar initial sd guess for distribution of integer copy number
#' @param use.tot logical flag of whether to use total coverage in inference
#' @param use.uniform logical flag of whether to use het coverage in inferendce
#' @param tol numeric tolerance
#' @param fix.sd logical flag whether to fix.sd to initial value in inference
#' @param max_iter integer scalar max iteration
#' @param verbose logical flag 
#' @keywords internal
ppem = function(segs.grid, segs = NULL, segs.h = NULL,
                 sd0 = segs.grid[!duplicated(j), sqrt(var(y, na.rm = TRUE))],
                 sd0.k = 3,
                 use.tot = TRUE,
                 use.het = TRUE, 
                 use.uniform = FALSE, 
                 tol = 1,
                 fix.sd = NULL,
                 max_iter = 100, verbose = FALSE)
{
  pta.h = pta = NULL  
  
  ## YYYYY
  ll.init = -1e100
  iteration = 0
  convergence = FALSE

  if (verbose)
    pmessage('Starting Ppurple EM with total segment grid containing ', nrow(segs.grid), ' rows')

  ## pi.k is our ploidy condition prior probability of copy states
  ## (which ends up softly specifying / constraining ploidy)
  ## we initialize to norm distribution centered at that ploidy

  pi.k0 = as.data.table(expand.grid(tau = unique(segs.grid$tau), k = unique(segs.grid$k)))
  pi.k0[k<0, pk := 0.1]
  pi.k0 = pi.k0[k>=0, pk := dnorm(k, tau, sd0.k)+0.0001]
  pi.k0[, pk := pk/sum(pk), by = tau]
  pi.k = copy(pi.k0)
  setkeyv(pi.k, c('tau', 'k'))


  segs.grid$sd = sd0

  segs.grid$pk =  pi.k[.(segs.grid$tau, segs.grid$k), pk]
  
  datarange = 100*quantile(segs.grid$y, na.rm = TRUE, c(0.99)) ## param for uniform dist

  segs.grid[, log.pxj_ktasp.high := 0]
  segs.grid[, log.pxj_ktasp.low := 0]


  ## set up uniform components likelihoods if we will be including this in the model
  if (!use.uniform)
    segs.grid[k<0, log.pxj_ktasp := -1e100] ## essentially placeholder

  use.tot = as.numeric(use.tot)
  use.het = as.numeric(use.het)

  start = Sys.time()
  while (convergence == FALSE && iteration < max_iter){
    iteration = iteration + 1

    if (verbose)
      pmessage(paste('Ppurple EM iteration', iteration, ':'))
#    print(Sys.time() - start)

    #################
    ## E step
    #################

    ## for total cn
    ## compute total copy number log likelihood using sos, sd, and pk
    segs.grid[k>=0, log.pxj_ktasp.tot := llnorm(y[1], sos[1], nbins[1], mu[1], sd[1]) + log(pk[1]), by = .(j, k, tau, alpha)]
    ## computing poisson likelihood for those segments that actually have hets 
    segs.grid[k>=0 & !is.na(y.high), log.pxj_ktasp.high := llpois(high.slfx, high.sx, nbins.h, mu.high), by = .(j, k, tau, alpha)]
    segs.grid[k>=0 & !is.na(y.low), log.pxj_ktasp.low := llpois(low.slfx, low.sx, nbins.h, mu.low), by = .(j, k, tau, alpha)]
    if (use.uniform)
      segs.grid[k<0, log.pxj_ktasp := nbins*log(1/datarange) + log(pk), by = .(j, k, tau, alpha)]

    ## we are summing across all joint high low states to get the overall likelihood of that total copy state
    segs.grid[k>=0, log.pxj_ktasp := (log.pxj_ktasp.tot[1]*use.tot
      + use.het*log.sum.exp(log.pxj_ktasp.high + log.pxj_ktasp.low)),  ## this log.sum.exp is over allelic copy states that add to a given k
      by = .(j, k, tau, alpha)]
    ## segs.grid[k>=0, log.pxj_ktasp.het := (0*log.pxj_ktasp.tot[1]*use.tot
    ##   + use.het*log.sum.exp(log.pxj_ktasp.high + log.pxj_ktasp.low)), 
    ##   by = .(j, k, tau, alpha)]
    segs.grid[k.low<=0, log.pxj_spat.norm := log.sum.exp(log.pxj_ktasp), ## note: every segment x k is guaranteed to have a k.low = 0, so we marginalizing over copy states in each seg
              by = .(j,alpha,tau)] ## log likelihood of segment j given sigma, pi (after integrating out k, k.high, k.low)
    segs.grid[, log.pk_xjspat := log.pxj_ktasp-log.pxj_spat.norm]  ## this is log posterior prob of k, tau, alpha aka r_jkat
    segs.grid[, r_jkat := exp(log.pk_xjspat)]  ## this is posterior prob of k given tau, alpha, sigma, pi
 
    if (!is.null(pta))
      ll.old = pta$log.px_spat
    else
      ll.old = ll.init

    ## combine seg and het likelihoods into a single log likelihood for alpha tau combo
    pta = segs.grid[.(TRUE), .(log.px_spat = sum(log.pxj_spat.norm)), keyby = .(tau, alpha)]  ## .(TRUE) gives one row per segment 
    pta[, log.pat_xsp := log.px_spat - log.sum.exp(pta$log.px_spat)]
    pta[, log.px_spat_old := ll.old]
    pta[, ll.diff := log.px_spat - log.px_spat_old]
    pta[, pta := exp(log.pat_xsp)]

    #################
    ## M step
    #################

    ## this part (in absence of uniform distribution) will estimate sd for copy number prior (which is discretized normal distribution
    ## centered at tau)  that updated prior to segs.grid

    ## re-compute pk ie copy number prior ie by finding new sd param for "discrete normal"
    pi.k = segs.grid[k.low<=0 , .(n = sum(r_jkat*nbins)), keyby = .(alpha, tau, k)] ## one segment per k, a, t

    ## remember - this is sd for the copy number prior which is centered at tau
    sd.k = segs.grid[k.low<=0, .(sd = sqrt((sum((r_jkat*nbins+1)*(k - tau)^2))/(sum((r_jkat*nbins+1))))), keyby = .(alpha, tau)] ## note: different sd.k for every alpha, tau

    ## merge with prior probabilities (only relevant for uniform component, i.e. if we are using them)
    pi.k = merge(merge(pi.k, pi.k0, by = c('tau', 'k')), sd.k, by = c('alpha', 'tau')) ##
    if (use.uniform)
      pi.k[, pk0 := (n[k==-1]+pk[k==-1])/sum(n+pk), by = .(alpha, tau)] ## we'll only keep this estimate for k = -1
    else
      pi.k[, pk0 := 0]
    pi.k[k<0, pk := pk0]
    pi.k[k>=0, pk := dnorm(k, tau, sd), by = .(alpha, tau)] ## reapply prior for non-uniform components using sd 
    pi.k[k>=0, pk := (pk/sum(pk))*(1-pk0), by = .(alpha, tau)] ## uniform component gets pk0, and the rest is to make sure we sum to 1
    setkeyv(pi.k, c('alpha', 'tau', 'k'))

    segs.grid$pk = pi.k[.(segs.grid$alpha, segs.grid$tau, segs.grid$k), pk]

    ## compute sd
    if (is.null(fix.sd))
    {
      ## k>=0 & k.low<=0 ensures that we have non duplicate segments per alpha, tau
      sd = segs.grid[k>=0 & k.low<=0, .(sd = sqrt(sum(r_jkat*sos.k)/sum(r_jkat*nbins))), keyby = .(alpha,tau)] ## estimate different sd for every alpha tau
      segs.grid$sd = sd[.(segs.grid[, .(alpha,tau)]), sd]
    }

    diff = max(pta$log.px_spat - ll.old)

    if (verbose)
    {
      pmessage('\tLL diff:', diff, ' tol: ', tol)
#      print('purity ploidy')
#      pta[, p := round(pta,3)]
#  print(pta[order(log.pat_xsp, decreasing = T)][1:10, ])
    }

    ## print(pta[order(log.pat_xsp, decreasing = T)][1:10, ])

    ## print('sd')
    ## print(sd[alpha == 0.1 & tau %in% c(3.8, 5), ])

    ## print('sd.k')
    ## print(sd.k[alpha == 0.1 & tau %in% c(3.8, 4), ])

    ## print('pi.k')
    ## print(dcast.data.table(pi.k[alpha == 0.1 & tau %in% c(3.8, 6), .(alpha, tau, k, p = round(pk, 2))], tau ~ k, value.var = 'p'))


    if(all(diff < tol)){convergence = TRUE}
  }
  return(list(pta = pta, segs = segs.grid, pk = pi.k))
}


pmessage = function(..., pre = 'Ppurple')
  message(pre, ' ', paste0(as.character(Sys.time()), ': '), ...)
