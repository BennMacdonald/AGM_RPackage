## Test environments
* local linux mint 17.3 64-bit, R 3.4.2
* win-builder (devel and release)

## R CMD check results
There were no ERRORs or WARNINGs. There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Frank Dondelinger <fdondelinger.work@gmail.com>'
  
  New submission
  Package was archived on CRAN
  
  
  Possibly mis-spelled words in DESCRIPTION:
    AGM (8:6)
    Dondelinger (8:11)
    Macdonald (9:5)
    al (8:26)
    et (8:23)
  
  CRAN repository db overrides:
    X-CRAN-Comment: Archived on 2019-12-21 as check problems were not
      corrected in time.

The issues with the package that caused it to fail the check have been addressed (see below). The words in the DESCRIPTION are spelled correctly.

## Downstream dependencies
No known dependencies.

## Further comments
Resubmission after it was removed from CRAN for throwing an error:
* Fixed code that was using the `class()` function incorrectly and throwing an error.