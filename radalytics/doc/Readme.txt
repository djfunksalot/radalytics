This package implements many of the approaches described in the
paper "Monte Carlo Methods for Bayesian Analysis of Survival Data
Using Mixtures of Dirichlet Process Priors."

-----------------------
Setting Up the Programs
-----------------------

You should have three files: this README file, and the files combo.f and
Splus_code.  The file combo.f contains Fortran programs which are
designed to be called from within Splus.  The file Splus_code contains
Splus functions.  To use these programs, carry out the following three
steps.

1. Compile the program combo.f with the -c option (compile, but do not
   link) via the command:

   f77 -c combo.f
   or
   f90 -c combo.f

   This should create an object file named combo.o

2. Start Splus and then read in the required Splus code using:

   source("Splus_code")

3. Load the Fortran program into Splus using either dyn.open(),
   dyn.load(), dyn.load2(), dyn.load.shared(), or static loading with
   LOAD.  The method used depends on what is available or preferred on
   your system.  These are illustrated below.  Our experience is with
   Unix-type systems.

   In the new version of S (version 5) you execute:

   dyn.open("combo.o")

   In older versions of Splus/S, you do one of the following.

   If dyn.load() or dyn.load2() are available, then after starting Splus
   you execute either:

   dyn.load("combo.o") or dyn.load2("combo.o")

   If dyn.load2() does not work for you, you may have to use the size
   option and do something like this:

   dyn.load2("combo.o", size=250000)

   If dyn.load.shared() is available, you issue the Unix command:

   Splus SHLIB -o combo.so combo.o

   Then you start Splus and execute the command:

   dyn.load.shared("./combo.so")

   If none of the above are available, you must use static loading as
   follows.  Issue the Unix command:

   Splus LOAD combo.o

   This creates a local version of Splus named local.Sqpe in the current
   directory.  This local version contains the Fortran subroutines in
   combo.  If you execute Splus in this directory, the local version
   will be used.


-------------------------------
General Description of Programs
-------------------------------

After completing the above three steps, you are ready to work with the
programs.

The available algorithms are gibbs1, sis1, sis2, and ritcen which are
executed using Splus functions of the same names.  These functions
compute the posterior mean and variance of the survival function, and
also give upper and lower probability limits for the survival function.
All of these functions do the same thing, but they use different
algorithms. The function gibbs1 uses a Gibbs sampler; sis1 and sis2 use
sequential importance sampling.  ritcen is a specialized program
intended only for right-censored data; it works by generating an i.i.d.
sample from the posterior distribution.  The programs gibbs1, sis1, and
sis2 allow for general interval censoring.

Here is a sample job using these functions.  

out <- gibbs1(bcdat.rad,c(2,88,10),10000)
summary(out)
plot(out)

All our functions have the same required arguments and use the same
format for the data.  Thus, in the sample job one can replace gibbs1 by
sis1 or sis2.  In the sample job, bcdat.rad is a data set which happens
to involve interval censoring, but for right-censored data, one can also
use ritcen.

In all our functions, the three required arguments are: (1) a two-column
matrix containing the data to be analyzed, (2) a triple specifying the
prior distribution, and (3) the simulation sample size.  In the example,
the estimates produced by gibbs1 are based on 10,000 observations taken
periodically from the output of a Gibbs sampler.

The functions summary and plot are used to display the output of the
gibbs1 procedure.  The data set bcdat.rad contains the breast cancer
data used in our technical report (the data for the radiation only
group).  The second group of the breast cancer data (those receiving
both radiation and chemotherapy) is contained in bcdat.radchem.


-----------------------------
The Format of the Data Matrix
-----------------------------

The data matrix must have two columns.  For easy comparison of ritcen,
gibbs1, sis1, and sis2 with each other and the survfit program in Splus,
these programs can accept a data matrix in either of two formats:

  Format 1: Each row gives the lower and upper endpoints of the time
    interval for that observation.  An observed death time (an
    uncensored) observation) is represented as a degenerate interval
    with both endpoints equal.  For right-censored observations, the
    time value infinity may be coded as NA, Inf, or -99.

  Format 2: Column 1 gives the times, and column 2 gives an indicator of
    status with 1 = dead (uncensored) and 0 = alive (right-censored).
    This format is suitable for passing to the Splus survfit function.

Here is an example of interval-censored data coded in Format 1.  This is
the bcdat.rad data used above.  (Note: this data is already available in
the file Splus_code.  Executing the Splus command source("Splus_code")
automatically creates the matrix bcdat.rad.)

 45   Inf
  6    10
  0     7
 46   Inf
 46   Inf
  7    16
 17   Inf
  7    14
 37    44
  0     8
  4    11
 15   Inf
 11    15
 22   Inf
 46   Inf
 46   Inf
 25    37
 46   Inf
 26    40
 46   Inf
 27    34
 36    44
 46   Inf
 36    48
 37   Inf
 40   Inf
 17    25
 46   Inf
 11    18
 38   Inf
  5    12
 37   Inf
  0     5
 18   Inf
 24   Inf
 36   Inf
  5    11
 19    35
 17    25
 24   Inf
 32   Inf
 33   Inf
 19    26
 37   Inf
 34   Inf
 36   Inf

Here is an example of right-censored data coded in Format 2.  This is
the treatment group from the well known gehan data.

 10    1
  7    1
 32    0
 23    1
 22    1
  6    1
 16    1
 34    0
 32    0
 25    0
 11    0
 20    0
 19    0
  6    1
 17    0
 35    0
  6    1
 13    1
  9    0
  6    0
 10    0

Suppose you have already carried out the dynamic loading of the fortran
programs and the "sourcing" of the file Splus_code.  If you have a data
file named "input" (for example) which is entered in one of the two
formats described above, then you can read in the data and run one of
our programs (say, gibbs1) with the following commands:

data <- read.table("input")
out <- gibbs1(data,c(2,88,10),10000)


------------------------------------------------
Optional Arguments Common to All the Algorithms:
------------------------------------------------

The options adtim, ci, do.mcse, npass.ci, nrep.ci, seed, and reuse are
available in all our algorithms.  The default values of these options
are listed below.  Following each is a brief description of the option.

   ------------
   adtim = NULL
   ------------

In our functions, the default is to compute estimates only for the times
which occur in the data; for interval-censored data, estimates are
supplied only for times which are end points of the censoring intervals.
The optional argument adtim allows us to specify additional times at
which estimates are to be computed.  For example, setting
adtim=c(2,30,42,60) produces output at the specified times in addition
to the times which occur in the data.

   ------
   ci = 0
   ------

The computation of upper and lower probability limits for the survival
function is fairly time consuming, so the default (ci = 0) is NOT to
compute them.  If you set ci = .95 (or, say, ci = .90), then the program
will compute limits which give .95 (or .90) pointwise posterior
probability intervals for the survival function.  In general, setting ci
= beta produces estimates of the alpha and 1-alpha quantiles (with alpha
= (1-beta)/2) of the posterior distribution of the survival function
S(t) for all those times t included in the analysis.

   -----------
   do.mcse = T
   -----------

When do.mcse = T, the programs compute rough estimates of the Monte
Carlo standard errors in the estimates of the posterior mean and
standard deviation of the survival function.  When do.mcse = F, these
are not computed.  The method of computing the standard errors varies
according to the algorithm.  For gibbs1, the estimates are produced by
breaking the gibbs sampler output into batches and computing the
variability between batches.  The number of batches is controlled by the
optional argument nbatch which has the default value 20.  This nbatch
option is only used by gibbs1.  The functions sis1 and sis2 compute
Monte Carlo standard errors using the delta-method.  The function ritcen
uses simple sample variances to estimate the standard errors.  This is
legitimate since the estimates produced by ritcen are based upon i.i.d.
observations from the posterior.

   -----------
   seed = NULL
   -----------

If you want your simulation results to be exactly reproducible at later
times, then you need to specify seeds for the random number generator.
The seed is a vector of three positive integers; for example, seed =
c(27, 10, 49).

[There is a bug in the current version of the program.  Setting the seed
does not always guarantee reproducible results.  We are trying to
correct this problem.]

   --------------
   npass.ci = 3
   nrep.ci = NULL
   reuse = F
   --------------

The probability intervals are obtained by a stochastic Newton-Raphson
procedure.  These three options give more detailed control over this
procedure.  Every iteration of the Newton-Raphson procedure requires
another sample (from the Gibbs or SIS scheme).  npass.ci is the number
of samples used in all.  The sample sizes are given by nrep.ci.  If the
default settings are not producing accurate enough probability
intervals, you can choose your own values for these options.  For
example, if you set npass.ci = 5 and nrep.ci = 5000, the program will
use 5 samples of size 5000 in the Newton-Raphson procedure.  If you set
reuse = T, the SAME sample will be used in each iteration of the
Newton-Raphson procedure.


-------------
Other options
-------------

Individual functions have their own particular options which are
explained in the comment statements inside the functions.


------------------
Another Sample Job
------------------

Here is another sample job which uses two of the more frequently used
options.

out <- sis2(bcdat.rad,c(2,88,10),20000,adtim=c(2,30,42,60),ci=.95)
summary(out)
plot(out)


-----------------
The plot function
-----------------

The plot function displays the survival function and the probability
intervals (if the ci option has been used).  If no intervals were
computed, the plot function displays a band around the survival function
which has a half-width which is simply a multiple of the posterior
standard deviation of the survival function.  The default multiplier is
1, but the mult option allows you to specify an arbitrary multiple.  If
you do not want to display any band at all, set band = F.  If you do not
want to create a new plot, but want to add to an existing plot, then set
new = F.


-------------------
Superimposing Plots
-------------------

The plot function allows one to superimpose plots (by setting new=F).
Here is an example of this:

out1 <- gibbs1(bcdat.rad,c(2,88,10),20000,adtim=c(2,30,42,60))
out2 <- gibbs1(bcdat.radchem,c(2,88,10),20000,adtim=c(2,30,42,60))
plot(out1,lty1=1,lty2=1)
plot(out2,new=F,lty1=4,lty2=4)

The arguments lty1 and lty2 control the line types which are used for
plotting the survival function and the surrounding band respectively.
The defaults are lty1 = 1, lty2 = 2.


------
parmod
------

For comparison purposes, we have included the function parmod.  The
function parmod fits the standard parametric Bayesian model in which
theta has a gamma prior and, conditional on theta, the lifetimes are
i.i.d. exponential(theta).  The function calls are the same as for
gibbs1, sis1, etc.  The ci option is not available.  Here is an example
of the use of this function.

out <- parmod(bcdat.rad,c(2,88,10),20000,adtim=c(2,30,42,60))
summary(out)
plot(out)

Note: The value of totm (the third entry in c(2,88,10)) is ignored.

Note: The function parmod is designed for interval censored data.  If
ALL of the censoring is right-censoring, the current version of the
program will crash.  There must be at least one observation in the data
set with a finite interval (the upper end-point is finite).


---------------------------------------------
Comparing our estimates with the Kaplan-Meier
---------------------------------------------

Here are some examples running the various programs on survival data
(leukemia and ovarian) available in S-plus.  The functions Surv and
survfit are in the S-plus survival analysis package.  These examples all
use very large sample sizes.  We can get away with this only because we
are not computing the probability intervals (ci=0).

# Create the leukemia data matrix.
data <- Surv(leukemia$time,leukemia$status)
fit1 <- survfit(data)  # Fits a Kaplan-Meier.
fit2 <- ritcen(data,c(1,30,5),50000,adtim=c(0,80,120))
fit3 <- gibbs1(data,c(1,30,5),50000,adtim=c(0,80,120))
fit4 <- sis1(data,c(1,30,5),50000,adtim=c(0,80,120))
fit5 <- sis2(data,c(1,30,5),50000,adtim=c(0,80,120))
summary(fit1)
summary(fit2)
summary(fit3)
summary(fit4)
summary(fit5)
plot(fit1)
plot(fit2)
plot(fit3)
plot(fit4)
plot(fit5)

# Create the ovarian data matrix.
data <- Surv(ovarian$futime,ovarian$fustat)
fit1 <- survfit(data)
fit2 <- ritcen(data,c(2,1600,10),50000,adtim=0)
fit3 <- gibbs1(data,c(2,1600,10),50000,adtim=0)
fit4 <- sis1(data,c(2,1600,10),50000,adtim=0)
fit5 <- sis2(data,c(2,1600,10),50000,adtim=0)
summary(fit1)
summary(fit2)
summary(fit3)
summary(fit4)
summary(fit5)
plot(fit1)
plot(fit2)
plot(fit3)
plot(fit4)
plot(fit5)

# Choosing a prior (in ritcen) to get an answer close 
# to the Kaplan-Meier estimate.
data <- Surv(leukemia$time,leukemia$status)
fit1 <- survfit(data)  # Fits a Kaplan-Meier.
fit2 <- ritcen(data,c(.1,.1,.0001),100000,adtim=c(0,80,120))


-------------------------------
Limits on the Size of Data Sets
-------------------------------

The programs are currently limited to data sets having no more than 210
observations.  It is easy to modify the programs to deal with larger
data sets.  Just do the following.  In the file combo.f, search for all
occurrences of the line

      parameter(nmax=210, ...

and change nmax to the desired size.  Also, the function llfn and the
subroutines rptset and llsolv (which are used by the ritcen subroutine)
currently limit the maximum number of censored observations to no more
than 100.  So, if you wish to use the ritcen algorithm on data sets
containing more than 100 censored observations, you must search for all
occurrences of the line

      parameter(mxcen=100)

in the file combo.f, and replace 100 by the desired value.  After making
these changes, you re-compile the programs (as described earlier) and
repeat the dynamic loading.


--------------------------
Credit where Credit is Due
--------------------------

In our programs we use many routines written by others.  In particular,
we use some of the ranlib routines (available from StatLib) compiled and
written by Barry W. Brown and James Lovato.  Also, we use programs
written by Alfred H. Morris (available from netlib).  Further
information is available in the comment statements of the Fortran
programs.

