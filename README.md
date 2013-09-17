radalytics
==========
assay results visualization

Radalytics is an R package for analyzing sample results and inventory data.

Rap on Inventory and Results:
Specimen Inventory and Results data really should be considered together
when attempting to discover new things about disease. 
if you're analyzing a specimin and discover something that leads you to want to do some more tests, it's important that you can actually grab some more specimens, or know that you just don't have any available for further analysis.  so we're going to try to do it all together.   wish us luck!

see http://www.stat.ufl.edu/system/r-pkg-tut/ for a good tutorial on getting started with r packages

to test everything out:
R CMD check radalytics

to build this package:
R CMD build radalytics

to install the package once you've built it:
R CMD INSTALL  radalytics_0.1.tar.gz
