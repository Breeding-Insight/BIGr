## R CMD check results

0 errors \| 0 warnings \| 1 note

-   This is a new release.

## Resubmission

In this version I have:

-   Updated DESCRIPTION
-   Added return value for merge_MADCs
-   Added optional seed for check_ped
-   Added verbose option for the packages specified or converted cat() to message()/warning()

We couldn't find the source of the issue "You write information messages to the console that cannot be easily suppressed" in script R/utils.R. We replaced the if(verbose) cat(..) by if(verbose) message(..).
