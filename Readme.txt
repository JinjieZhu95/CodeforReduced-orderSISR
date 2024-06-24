The code files for 'Reduced-order dynamics of coupled self-induced stochastic resonance quasiperiodic oscillators driven by varying noise intensities'
by Jinjie Zhu
e-mail: jinjiezhu[at]nuaa.edu.cn 
***************************************************
***************************************************

The main code:
Main4CoupledSISR.m
----Seven parts for the comparison between the results of the reduced model and the original stochastic model.
----The data loaded therein are mainly obtained from the following other codes.
***************************************************

Other codes:
Weibull2MLE.m
----Two-parameter Weibull estimation of the FPTD of the left and right branches.
----The data of 'ypositionsigma0d0*.mat' are from Monte Carlo simulations of the original stochastic system.
----paramsuml and paramsumr store the results as those in Table 1 in the main text.

Estimateofabversussigma.m
----Polynomial interpolations of the scale and shape parameters.

ThetaJumpleft2right.m
----Calculation of K(\theta).
