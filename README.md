# robust_tail_numerics

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!! Under construction !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
This repositary is a companion to the paper Tail Analysis without Parametric Models: A Worst-case Perspective by Lam H.,
and Mottet C.. Here are made available ALL the tables, the datasets, and the figures given in the section discussing
numerical findings (Section 7). These results were generated using python. For
transparency and reproducibility purposes, the R-codes used in the computations are also made available.

We point out that these codes rely on functions available in the python packages robust_tail and distribution_pty (see
installation detail below), that we also developed. The motivation for building separate packages is to have readable,
and coherently organized codes. More specifically,

* robust_tail is a package that focus on solving Equation (5) and (EC.19) in the specific case where only 2 point masses
distribution functions are considered. In addition, RobustTail provides functions to estimate the parameters
$\underline \eta$, $\overline \eta$, $\underline \beta$, $\overline \beta$, and $\overline \nu$ of program (11).
* distribution_pty, as its name suggests, is a package containing functions related to the properties of some probability
distributions functions. As an example, the function Dlnorm and Dpareto gives the derivatives of the log-Normal and the
Pareto distribution.
