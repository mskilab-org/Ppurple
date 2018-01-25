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
<thead><tr><th scope=col>seqnames</th><th scope=col>start</th><th scope=col>end</th><th scope=col>width</th><th scope=col>strand</th><th scope=col>ratio</th><th scope=col>tum.counts</th><th scope=col>norm.counts</th><th scope=col>y</th></tr></thead>
<tbody>
	<tr><td>9        </td><td> 57584801</td><td> 57585000</td><td>200      </td><td>*        </td><td>       NA</td><td>       NA</td><td>       NA</td><td>       NA</td></tr>
	<tr><td>12       </td><td>  3410201</td><td>  3410400</td><td>200      </td><td>*        </td><td>0.8453500</td><td>0.8596121</td><td>0.9817843</td><td>0.8453500</td></tr>
	<tr><td>17       </td><td> 27437401</td><td> 27437600</td><td>200      </td><td>*        </td><td>0.8585241</td><td>1.0653026</td><td>1.1980380</td><td>0.8585241</td></tr>
	<tr><td>5        </td><td>163569001</td><td>163569200</td><td>200      </td><td>*        </td><td>1.0817177</td><td>1.1726743</td><td>1.0466792</td><td>1.0817177</td></tr>
	<tr><td>11       </td><td> 63494201</td><td> 63494400</td><td>200      </td><td>*        </td><td>0.7240008</td><td>0.7616872</td><td>1.0157520</td><td>0.7240008</td></tr>
	<tr><td>12       </td><td> 54107601</td><td> 54107800</td><td>200      </td><td>*        </td><td>0.8295073</td><td>0.9813019</td><td>1.1421747</td><td>0.8295073</td></tr>
</tbody>
</table>


```R
> hets = fread(system.file("extdata", "hets.csv", package = "Ppurple"))
> head(hets)
```

<table>
<thead><tr><th scope=col>seqnames</th><th scope=col>start</th><th scope=col>end</th><th scope=col>Tumor_Seq_Allele1</th><th scope=col>Reference_Allele</th><th scope=col>ref.count.t</th><th scope=col>alt.count.t</th><th scope=col>ref.count.n</th><th scope=col>alt.count.n</th><th scope=col>alt.frac.t</th><th scope=col>ref.frac.t</th><th scope=col>alt.frac.n</th><th scope=col>ref.frac.n</th><th scope=col>alt</th><th scope=col>ref</th></tr></thead>
<tbody>
	<tr><td>12       </td><td> 14628069</td><td> 14628069</td><td>G        </td><td>A        </td><td> 6       </td><td>31       </td><td>11       </td><td>20       </td><td>0.8378378</td><td>0.1621622</td><td>0.6451613</td><td>0.3548387</td><td>31       </td><td> 6       </td></tr>
	<tr><td>10       </td><td> 96954298</td><td> 96954298</td><td>G        </td><td>A        </td><td>25       </td><td>22       </td><td>27       </td><td>26       </td><td>0.4680851</td><td>0.5319149</td><td>0.4905660</td><td>0.5094340</td><td>22       </td><td>25       </td></tr>
	<tr><td>2        </td><td> 51229355</td><td> 51229355</td><td>A        </td><td>C        </td><td>27       </td><td>25       </td><td>24       </td><td>42       </td><td>0.4807692</td><td>0.5192308</td><td>0.6363636</td><td>0.3636364</td><td>25       </td><td>27       </td></tr>
	<tr><td>X        </td><td> 43498623</td><td> 43498623</td><td>G        </td><td>A        </td><td>10       </td><td>10       </td><td>22       </td><td>25       </td><td>0.5000000</td><td>0.5000000</td><td>0.5319149</td><td>0.4680851</td><td>10       </td><td>10       </td></tr>
	<tr><td>11       </td><td> 97932192</td><td> 97932192</td><td>G        </td><td>A        </td><td>12       </td><td> 4       </td><td>15       </td><td>17       </td><td>0.2500000</td><td>0.7500000</td><td>0.5312500</td><td>0.4687500</td><td> 4       </td><td>12       </td></tr>
	<tr><td>5        </td><td>145511921</td><td>145511921</td><td>A        </td><td>G        </td><td> 4       </td><td>27       </td><td>13       </td><td>19       </td><td>0.8709677</td><td>0.1290323</td><td>0.5937500</td><td>0.4062500</td><td>27       </td><td> 4       </td></tr>
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
    Warning message in dt2gr(hets):
    “coercing to GRanges via non-standard columns”Warning message in seg2gr(dt, seqlengths, seqinfo):
    “some seqnames in seg object were not included in provided seqlengths: MT”Running hapseg on hets
    Hapseg iteration 1:
    		LL diff: 1e+100 tol:1
    Hapseg iteration 2:
    		LL diff: 5.17692387802526e-06 tol:1
    Running ppemgrid with 11 purities ranging from 0 to 1 and 21 ploidies ranging from 1 to 5 with rho of 1.00978885796089 het rho of 22.283605.
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

