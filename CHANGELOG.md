

Version 0.1-0
=============

* Added C functions for calculating PDF, CDF, and pmax when
    doing predictions. Currently `predict()` allows for an
    undocumented input argument `useC = TRUE` (defaults to `FALSE`)
    for testing. Roughly 3-5 times faster than plain R.
* Several speed/performance improvements.
