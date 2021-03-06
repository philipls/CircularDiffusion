Matlab/C routines to fit the circular diffusion model to noisy color patch data. To accompany
Smith, P. L., Saber, S., Corbett, E. A., & Lilburn, S. D. (2020). Modeling continuous outcome color decisions with the circular diffusion model: Metric and categorical properties, Psychological Review, xx, xxx-xxx.

I. Resources (Matlab/C)
---------------------

1.  ColData.mat - raw and summary data for the individual participants
2.  paicircle11 - encoding failure model with 3 color categories (S1, S2) Table 1 fit
3.  paicircle13 - encoding failure model with 5 color categories (S3, S4) Table 1 fit
4.  colplot11x - plots distributions of decision outcomes and RTs (e.g., Figure 5)
5.  qplot11x - plots Q x Q plot (e.g., Figure 6)
6.  angplotmd - predictions of color category model (Figure 13). 

Routines for the Jones-Pewsey phase-angle model are not included in this library. They are available on request.


II. Dependencies
-----------------
paicircle11 --> vdcircle300cls.c
paicircle13 --> vdcircle300cls.c

The C-code was compiled under Linux using the GNU scientific library, gsl, and the gslcblas library. These libraries provide Bessel functions and their roots. Other libraries may provide these functions in other environments. The C-code is compiled under Matlab with the command:

mex vdcircle300cls.c -lgsl -lgslcblas  -lm

III. Data Structures
--------------------

The individual participant data is contained in structures S1, S2, S3, S4 in ColData.mat. There are also raw data files S1raw, S2raw, S3raw, S4raw, which is are needed for plotting with angplotmd. The S1...S4 data structures are 3-element cell arrays, one for each level of stimulus discriminability. Each cell contains a ntrials x 3 matrix, with columns: stimulus angle (-pi - pi), response error (-pi - pi) and RT (s).

The .mat files S1fit, S2fit, S3fit, S4fit contain the parameters and fit statistics for the encoding failure model in Table 1 of the article.

IV. Calling Conventions
-----------------------

Help calls: "help paicircle11", etc. gives the calling conventions:

>> help paicircle11
  ==========================================================================
  Circular diffusion with drift anisotropies for color circle task.
  Rectangular variability in criterion.
    [ll,qaic,qbic,Pred,Gstuff] =  paicircle11(Pvar, Pfix, Sel, Data, trace)
     P = [v1...v3b, eta1a....eta3b, a, Ter, b1...b3, alpha, a1...a3, sa, pi1]
           1...6         7...12    13, 14,  15..17     18   19...21  22  23
   'Data' is a 3-element cell array
   3-category version (S1 and S2)
   Overdispersion set internally.
  ===========================================================================

For paicircle11, the call is:

   [ll,qaic,qbic,Pred,Gstuff] =  paicircle11(Pvar, Pfix, Sel, Data, trace)

The function takes a vector of variable parameters, Pvar, which are estimated during fitting, a vector of fixed parameters, Pfix, which remain fixed, a selector vector, Sel, used to select and assemble parameter vectors, a data matrix, Data, and an optional trace switch, which gives information about parameter bounds. How the elements of Pfix are treated depends on the internal logic of paicircle11.

S1fit contains the fit for S1. This was carried out with
Sel =[1     0     1     0     1     0     1     0     1     0     1     0     1     1     1     1 1     1     1     1     1     0     1];

The fitting routines generate predictions by rotating the stimulus into a canonical orientation (theta = 0) so that the first component of drift rate points along the positive x-axis, and the second component points along the positive y-axis. In canonical oriention these components can be interpreted as the radial and tangential components of drift rate, respectively. Smith (Journal of Mathematical Psychology, 91, 145-158, 2019) provides expressions for the radial and tangential components of drift rate for arbitrary orientations but they were not used in this article. The first and second eta components are the radial and tangential components of drift rate variability. The Sel vector above fixes all of the tangential components of drift rate and drift variability to their starting values. In the fits in Table 1, these parameters were set to small (effectively negligible) values. The second-to-last element of Sel sets the criterion variability parameter to zero. 

V. Fitting and Plotting
-----------------------

The functions paicircle11 and paicircle13 are minimized using Simplex by

>> setopt
>> pest = fminsearch(@paicircle11, P(Sel==1), options, P(Sel==0), Sel, S1)
>> P(Sel==1) = pest
>> [ll,qaic,qbic,Pred,Gstuff]=paicircle11(P(Sel==1), P(Sel==0), Sel, S1)
 
Plots of distributions of decision outcomes and decision times are obtained by

>> colplot11x(S1, Pred)

Q x Q plots are generated by

>> qlot11x(@paicircle11, P(Sel==1), P(Sel==0), Sel, S1)

(This routine takes a function parameter to make it general.)

paicircle11 and paicircle13 fit a color category model, whose parameters are b1...b3, alpha, a1...a3.  The b's are the norms of the category vectors; the a's are their locations, and alpha is the exponential decay parameter. The color category model can be removed by setting b1...b3 = 0 and setting the entries of Sel for the b's, the a's, and alpha to zero.

The color category predictions are plotted by

>> angplotmd(S1raw, Pred)

The raw data are needed here, because the signed response errors is calculated using the true stimulus and response angles. 

VI. Miscellaneous Notes and Cautions
------------------------------------

There are a number of context-sensitive flags and variables that are set internally inside the functions. Overdispersion is set manually on l. 22 of paicircle11 and paicircle13 by selecting the correct element of the Overdispersion array for the subject. The maximum time index is currently set to 4.0. There is also a combination of hard and soft (quadratically penalized) constraints on the parameters to keep Simplex out of parts of the space where things may go bad numerically. These are very ad hoc, but work well enough in practice, so long as you're aware of what the settings are. Setting the trace flag to 1 in paicircle11 and paicircle13 will give an indication if you're on the edge of the penalty region and may need to adjust things. 

The number of time steps on the range [0, tmax] is currently set to 300 and the number of hitting point bins on the range [0, 2pi] is currently set to 50. Both paicircle11 and paicircle13 save a copy of the working parameter vector to disk each time they are called. These calls add some disk-access overheads but are convenient if you want to exit a fit prematurely or to diagnose where something went bad. Because of the context-sensitive nature of these flags, these routines should not be treated as black-boxes, but instead as starting points that you should adapt to your own requirements. 

The infinite-series expression for the first-passage time density for the Bessel process (function dhamana) truncates at 50 terms. This can sometimes result in a small spike of artifactual probably mass near zero. There is a parameter "badix," which forces the density to zero in the region of t = 0 to control for this. For most participants badix = 5 works well, but if the best decision criterion is large, then badix may need to be increased. This is the case for S4, for whom the best criterion was a ~ 4.0. The fits for S4 in the article were obtained with badix = 25, which is set on l. 260 of paicircle13. Any problems of this kind are very apparent in the colplot11x plot. 

As usual, this code is supplied as is, for noncommercial, academic research purposes only, without any explicit or implied warranties or liabilities. 










