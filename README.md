[![Build Status](https://travis-ci.org/mskilab/gUtils.svg?branch=master)](https://travis-ci.org/mskilab/Ppurple)
[![codecov.io](https://img.shields.io/codecov/c/github/mskilab/gUtils.svg)](https://codecov.io/github/mskilab/Ppurple?branch=master)


Ppurple
=======
Probabilistic purity ploidy estimation


Installation
------------

1. Install dependent packages and latest Bioconductor (if you haven't already)

```{r}
source("https://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")
```

2. Install devtools from CRAN (if you don't have it already)

```{r}
install.packages('devtools')
```

3. Install mskilab R dependencies (gUtils and bamUtils)


```{r}
devtools::install_github('mskilab/gUtils)
devtools::install_github('mskilab/bamUtils)
```

4. Install Ppurple


```{r}
devtools::install_github('mskilab/Ppurple')
```


Tutorial
------------

Ppurple uses EM to infer purity and ploidy from total coverage and heterozygote
counts, provided as data.frames, data.tables, or GRanges.  It returns a
data.table of ranked solutions, associated with a posterior probability. 

### Load Ppurple

```R
library(Ppurple)
```
### Load coverage, hets, and segs


```R
> cov = fread(system.file("extdata", "coverage.csv", package = "Ppurple"))
> head(cov)
```


<table>
<thead><tr><th scope=col>seqnames</th><th scope=col>start</th><th
scope=col>end</th><th scope=col>width</th><th scope=col>strand</th><th
scope=col>y</th></tr></thead>
<tbody>
<tr><td>1        </td><td> 79401   </td><td> 79600   </td><td>200
</td><td>*        </td><td>1.5395624</td></tr>
<tr><td>1        </td><td>531201   </td><td>531400   </td><td>200
</td><td>*        </td><td>1.9518525</td></tr>
<tr><td>1        </td><td>555401   </td><td>555600   </td><td>200
</td><td>*        </td><td>0.8810031</td></tr>
<tr><td>1        </td><td>571601   </td><td>571800   </td><td>200
</td><td>*        </td><td>0.8270474</td></tr>
<tr><td>1        </td><td>573601   </td><td>573800   </td><td>200
</td><td>*        </td><td>1.0702993</td></tr>
<tr><td>1        </td><td>617801   </td><td>618000   </td><td>200
</td><td>*        </td><td>1.3891481</td></tr>
</tbody>
</table>

```R
> hets = fread(system.file("extdata", "hets.csv", package = "Ppurple"))
> head(hets)
```

<table>
<thead><tr><th scope=col>seqnames</th><th scope=col>start</th><th
scope=col>end</th><th scope=col>width</th><th scope=col>strand</th><th
scope=col>ALT</th><th scope=col>REF</th><th scope=col>alt</th><th
scope=col>ref</th></tr></thead>
<tbody>
<tr><td>1      </td><td> 779322</td><td> 779322</td><td>1      </td><td>*
</td><td>G      </td><td>A      </td><td>36     </td><td>37     </td></tr>
<tr><td>1      </td><td> 998395</td><td> 998395</td><td>1      </td><td>*
</td><td>G      </td><td>A      </td><td>27     </td><td>33     </td></tr>
<tr><td>1      </td><td> 998501</td><td> 998501</td><td>1      </td><td>*
</td><td>C      </td><td>G      </td><td>21     </td><td>29     </td></tr>
<tr><td>1      </td><td>1158277</td><td>1158277</td><td>1      </td><td>*
</td><td>A      </td><td>G      </td><td>29     </td><td>23     </td></tr>
<tr><td>1      </td><td>1160665</td><td>1160665</td><td>1      </td><td>*
</td><td>A      </td><td>G      </td><td>27     </td><td>24     </td></tr>
<tr><td>1      </td><td>1206619</td><td>1206619</td><td>1      </td><td>*
</td><td>A      </td><td>C      </td><td> 1     </td><td>52     </td></tr>
</tbody>
</table>



### Run Ppurple without precomputed segs


```R
> pp = ppurple(cov = cov, hets = hets, verbose = TRUE)
```

    Segments not provided so doing internal segmentation via DNAcopy
    sending  92654  segments to DNAcopy
    ... ...
    Ppurple EM iteration 3 :
    		LL diff:58.2742695834022 tol: 1
    Ppurple EM iteration 4 :
    		LL diff:0.912960716173984 tol: 1


#### output is a data.table of solutions and their probabilities, showing the most likely solution in the first row


```R
> pp[1,]
```

<table>
<thead><tr><th scope=col>purity</th><th scope=col>ploidy</th><th scope=col>prob</th></tr></thead>
<tbody>
	<tr><td>0.49</td><td>3.86</td><td>1   </td></tr>
</tbody>
</table>

### Run Ppurple with pre-computed segs


```R
> segs = fread(system.file("extdata", "segs.csv", package = "Ppurple"))
> head(segs)
```

<table style = 'font-size:50%'>
<thead><tr><th scope=col>seqnames</th><th scope=col>start</th><th scope=col>end</th><th scope=col>width</th><th scope=col>strand</th><th scope=col>ID</th><th scope=col>num.mark</th><th scope=col>seg.mean</th></tr></thead>
<tbody>
	<tr><td>1       </td><td>   79401</td><td> 6376801</td><td> 6297401</td><td>*       </td><td>Sample.1</td><td>183     </td><td> 0.0683 </td></tr>
	<tr><td>1       </td><td> 6437801</td><td> 8702401</td><td> 2264601</td><td>*       </td><td>Sample.1</td><td> 76     </td><td>-0.0967 </td></tr>
	<tr><td>1       </td><td> 8758201</td><td> 9043001</td><td>  284801</td><td>*       </td><td>Sample.1</td><td>  9     </td><td> 0.1808 </td></tr>
	<tr><td>1       </td><td> 9080201</td><td>20046801</td><td>10966601</td><td>*       </td><td>Sample.1</td><td>347     </td><td>-0.1132 </td></tr>
	<tr><td>1       </td><td>20113201</td><td>23387801</td><td> 3274601</td><td>*       </td><td>Sample.1</td><td>102     </td><td> 0.0931 </td></tr>
	<tr><td>1       </td><td>23421201</td><td>48154401</td><td>24733201</td><td>*       </td><td>Sample.1</td><td>795     </td><td>-0.1152 </td></tr>
</tbody>
</table>

```R
> pp = ppurple(cov = cov, hets = hets, segs = segs, verbose = TRUE)
```

    Fitting initial grid of 11 purity and 21 ploidy combinations.
    Hapseg iteration 1:
    		LL diff: 1e+100 tol:1
    Hapseg iteration 2:
    		LL diff: 5.17692387802526e-06 tol:1
    Running ppemgrid with 11 purities ranging from 0 to 1 and 21 ploidies ranging from 1 to 5 with rho of 1.00978885796089 het rho of 22.283605.
    ... ...
    Ppurple EM iteration 3 :
    		LL diff:58.6840663431212 tol: 1
    Ppurple EM iteration 4 :
    		LL diff:0.864216329413466 tol: 1


```R
>  pp[1,]
```


<table>
<thead><tr><th scope=col>purity</th><th scope=col>ploidy</th><th scope=col>prob</th></tr></thead>
<tbody>
	<tr><td>0.49</td><td>3.86</td><td>1   </td></tr>
</tbody>
</table>


Attributions
------------
> Marcin Imielinski - Weill Cornell Medicine; Core Member, New York Genome Center

> Aditya Desphande - Tri-I CBM PhD candidate, Weill Cornell Medicine

[license]: https://github.com/mskilab/gUtils/blob/master/LICENSE
[docs]: http://gutils.readthedocs.org/en/latest/index.html

