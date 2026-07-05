# proxymix: KL-Optimal Gaussian Mixture Proxies for Arbitrary Target Densities

Fits multivariate Gaussian-mixture proxies that are Kullback-Leibler
optimal to user-supplied target densities on real Euclidean space. Three
fitting regimes are unified under one verb: (i) closed-form moment
matching for a single component, (ii) classical expectation-maximisation
when independent samples are available, and (iii) importance-sampled
KLD-EM when the target can be evaluated point-wise but not (cheaply)
sampled. Closed-form Gaussian-mixture operators (density, sampling,
marginalisation, conditioning, divergence) round out the toolkit. The
conditioning operator drives multiple imputation of data missing at
random, covering the multimodal and heteroscedastic cases a
single-Gaussian model cannot represent. Implements the regime hierarchy
of Hoek and Elliott (2024)
[doi:10.1080/07362994.2024.2372605](https://doi.org/10.1080/07362994.2024.2372605)
.

## See also

Useful links:

- <https://github.com/max578/proxymix>

- <https://max578.github.io/proxymix/>

- Report bugs at <https://github.com/max578/proxymix/issues>

## Author

**Maintainer**: Max Moldovan <max.moldovan@gmail.com>
([ORCID](https://orcid.org/0000-0001-9680-8474))

Other contributors:

- Johannes van der Hoek (Foundational theory in Hoek & Elliott (2024).)
  \[contributor\]

- Robert J. Elliott (Foundational theory in Hoek & Elliott (2024).)
  \[contributor\]
