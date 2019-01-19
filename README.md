## Documentation

**An [introductory vignette](https://benjjneb.github.io/decontam/vignettes/decontam_intro.html) demonstrating how to use the decontam package to identify contaminants.**

**The [preprint introducing decontam](https://doi.org/10.1101/221499) with benchmarking demonstrating how `decontam`-inating your data removes reagent sequences, improves accuracy, reduces batch effects, and prevents false-positive assocations.**

More documentation available through the R help interface (e.g. `?isContaminant`) and at [the decontam web site](https://benjjneb.github.io/decontam).

## Installation

The decontam R package is currently available as a source package through github, so installation requires the ability to compile R packages. This means that R and the R tool-chain must be installed, which requires the [Xcode command-line tools](http://railsapps.github.io/xcode-command-line-tools.html) on Mac and [Rtools](https://cran.r-project.org/bin/windows/Rtools/) on Windows.

The easiest installation method uses the `devtools` package:

```S
library(devtools)
devtools::install_github("benjjneb/decontam")
```

Alternatively you can install from source by hand. First download the zipped package and unzip it (or expand the tarball). Start a fresh R session, and enter the following.

```S
install.packages("path/to/decontam", repos = NULL, type = "source",
                 dependencies = c("Depends", "Suggests","Imports"))
```

## Other Resources

Bugs and difficulties in using decontam are welcome on [the issue tracker](https://github.com/benjjneb/decontam/issues).

Planned feature improvements are also publicly catalogued at on the "Issues" page for decontam: https://github.com/benjjneb/decontam/issues
