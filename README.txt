The m-files in this directory are associated with the following three papers:

[1] A. Bespalov, D. Silvester and F. Xu,
    Error estimation and adaptivity for stochastic collocation finite elements.
    Part I: single-level approximation.
    SIAM Journal on Scientific Computing, Vol. 44 (2022), Issue 5, pp. A3393-A3412.
    https://doi.org/10.1137/21M1446745
    Preprint
    http://arxiv.org/abs/2109.07320

[2] A. Bespalov and D. Silvester,
    Error estimation and adaptivity for stochastic collocation finite elements.
    Part II: multilevel approximation.
    SIAM Journal on Scientific Computing (to appear).
    Preprint
    https://arxiv.org/abs/2202.08902

[3] A. Bespalov, Dirk Praetorius, Thomas Round and Andrey Savinov,
    Goal-oriented error estimation and adaptivity for stochastic collocation FEM.
    Preprint
    https://arxiv.org/abs/2406.05028


The driver for generating adaptive single-level SC-FEM approximations is
singlelevelSC.m

The driver for generating adaptive multilevel SC-FEM approximations is
multilevelSC.m

The driver for generating adaptive single-level SC-FEM approximations for goal-oriented problems is
goafem_singlelevelSC.m

The default implementation is for a Unix architecture.
On a Windows machine, run the script-file
install_adaptive_scfem.m 
before running the above main drivers for the first time.


The diary files included in this directory were generated using MATLAB R2021a
running under Windows 10 Enterprise x64 Version 10.0.19042.

The test runs reproducing the numerical results in [1, Section 5] are saved
in the following diary files:
SCtest_sl_tp1.txt
SCtest_sl_tp2a.txt
SCtest_sl_tp2b.txt
SCtest_sl_tp2c.txt

The test runs reproducing the numerical results in [2, Section 4] are saved
in the following diary files:
SCtest_ml_tp1.txt
SCtest_ml_tp2c.txt
SCtest_sl_tp3.txt
SCtest_ml_tp3.txt

The T-IFISS* driver for the SG comparison run is
stoch_adapt_testproblem
The associated diary file is
SGtest5.1.txt
* Stochastic T-IFISS can be downloaded from
https://github.com/albespalov/Stochastic_T-IFISS

