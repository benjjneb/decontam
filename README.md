## Documentation

**An [introductory vignette](https://benjjneb.github.io/decontam/vignettes/decontam_intro.html) demonstrating how to use the decontam package to identify contaminants.**

**The [manuscript introducing decontam](https://doi.org/10.1186/s40168-018-0605-2) with benchmarking demonstrating how `decontam`-inating your data removes reagent sequences, improves accuracy, reduces batch effects, and prevents false-positive associations.**

More documentation available through the R help interface (e.g. `?isContaminant`) and at [the decontam web site](https://benjjneb.github.io/decontam).

## Installation

The decontam R package is most easily installed from the Bioconductor repository. To install the most recent version of the decontam package, start R (version "3.6") and enter:

```S
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("decontam")
```

See [the Bioconductor documentation for decontam](https://bioconductor.org/packages/release/bioc/html/decontam.html) for additional information.

The decontam R package is also available as a source package through github. Installation requires the ability to compile R packages. This means that R and the R tool-chain must be installed, which requires the [Xcode command-line tools](http://railsapps.github.io/xcode-command-line-tools.html) on Mac and [Rtools](https://cran.r-project.org/bin/windows/Rtools/) on Windows.

The easiest source installation method uses the `devtools` package:

```S
library(devtools)
devtools::install_github("benjjneb/decontam")
```

## Other Resources

Bugs and difficulties in using decontam are welcome on [the issue tracker](https://github.com/benjjneb/decontam/issues).

Planned feature improvements are also publicly catalogued at on the "Issues" page for decontam: https://github.com/benjjneb/decontam/issues
