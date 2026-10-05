## Re-submission of Archived Package mcmcplots

This is a re-submission of the previously archived package `mcmcplots` (archived on 2025-07-10).

### Reason for Archival and Fixes

The package was archived due to uncorrected issues following reminders. The following fixes have been implemented:

* **Rd cross-reference package anchors:** Added package anchor `\link[denstrip]{denstrip}` in `man/caterplot.Rd` and `\link[utils]{browseURL}` in `man/mcmcplot.Rd`.
* **LazyData specification:** Removed `LazyData: yes` from `DESCRIPTION` since the package contains no `data/` directory.
* **Obsolete LazyLoad directive:** Removed obsolete `LazyLoad: yes` from `DESCRIPTION`.
* **Codoc and Rd usage:** Added the missing `browse = TRUE` argument to the `\usage` section of `man/mcmcplot.Rd`.
* **Metadata & Description:** Expanded `DESCRIPTION` description field, updated author metadata, added `URL` and `BugReports` links, and updated release date.

### Test environments
* Local: Windows 11, R 4.6.1
* Win-builder: Windows Server 2022, R Under development (unstable) (2026-09-30 r90605 ucrt)

### R CMD check results

There were no ERRORs or WARNINGs on either platform.

1 NOTE:
* Checking CRAN incoming feasibility: Package was archived on CRAN (expected for re-submission of an archived package).
