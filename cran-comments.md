## Re-submission of Archived Package mcmcplots (version 0.5.0)

This is a re-submission of the previously archived package `mcmcplots` (archived on 2025-07-10 as version 0.4.4).
Version 0.5.0 addresses all CRAN reviewer feedback, modernizes the HTML reporting interface, and expands compatibility across the modern MCMC ecosystem (posterior draws, 3D arrays, Stan/brms objects).

### Reason for Archival and Fixes

The package was archived due to uncorrected issues following reminders. The following fixes have been implemented:

* **Rd cross-reference package anchors:** Added package anchor `\link[denstrip]{denstrip}` in `man/caterplot.Rd` and `\link[utils]{browseURL}` in `man/mcmcplot.Rd`.
* **LazyData specification:** Removed `LazyData: yes` from `DESCRIPTION` since the package contains no `data/` directory.
* **Obsolete LazyLoad directive:** Removed obsolete `LazyLoad: yes` from `DESCRIPTION`.
* **Codoc and Rd usage:** Added the missing `browse = TRUE` argument to the `\usage` section of `man/mcmcplot.Rd`.
* **Metadata & Description:** Expanded `DESCRIPTION` description field, updated author metadata, added `URL` and `BugReports` links, and updated release date.

### Responses to CRAN Team Feedback

* **Missing `\value` tags and documentation of return values**: Added explicit `\value` tag in `man/caterpoints.Rd` documenting that the function has no return value and is called for side effects. Additionally, audited all package `.Rd` files to replace generic "Creates a plot" text with precise return values (e.g. invisibly returning plotted parameter names in `denplot`, `traplot`, `rmeanplot`) or explicit "No return value, called for side effects" descriptions. Added `\value` and `\keyword{internal}` tags to unexported internal Rd files (`.graypr.Rd`, `.html.end.Rd`, `.html.img.Rd`, `.html.begin.Rd`).
* **Unexported function call with `:::` in `man/to.greek.Rd`**: Removed the executable call `mcmcplots:::.to.greek(...)` and replaced the example with a descriptive comment, marking the file with `\keyword{internal}` so `.to.greek` remains cleanly unexported.
* **Examples structure and `\dontrun` vs `\donttest`**:
  - In `man/mcmcplot.Rd`, unwrapped a fast non-browser example (`mcmcplot(..., browse = FALSE)`) that executes in < 1 second, and wrapped interactive browser-launching calls inside `\donttest{}`.
  - In `man/caterplot.Rd`, `man/denplot.Rd`, `man/traplot.Rd`, and `man/mcmcplots-package.Rd`, converted optional legacy `R2WinBUGS` demonstration code into documented comments, eliminating unstated dependency warnings while preserving usage guidance. In `man/as.mcmc.bugs.Rd`, replaced `\dontrun{}` with `\donttest{}` guarded by `requireNamespace("R2WinBUGS", quietly = TRUE)`.
* **Restoration of graphical parameters (`par`)**: Added `oldpar <- par(mfrow = c(2, 2))` and trailing `par(oldpar)` in all examples that modify graphical parameters (`man/caterpoints.Rd`, `man/caterplot.Rd`, `man/mcmcplots-package.Rd`).


### Test environments
* Local: Windows 11, R 4.6.1
* Win-builder: Windows Server 2022, R Under development (unstable) (2026-10-14 r90673 ucrt)

### R CMD check results

There were no ERRORs or WARNINGs on either platform.

1 NOTE:
* Checking CRAN incoming feasibility: Package was archived on CRAN (expected for re-submission of an archived package).
