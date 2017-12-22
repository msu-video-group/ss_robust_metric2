Robust Saliency Map Comparison, version 2
=========================================

ss_robust_metric2 MATLAB function performs comparison of saliency video sequence with ground-truth saliency video
sequence using the method described in "EVALUATION METHODOLOGY" section of

> V. Lyudvichenko, M. Erofeev, Y. Gitman and D. Vatolin, "A semiautomatic saliency model and its application to video compression," 2017 13th IEEE International Conference on Intelligent Computer Communication and Processing (ICCP), Cluj-Napoca, 2017, pp. 403-410.
> doi: 10.1109/ICCP.2017.8117038

The comparison carried out by the function is invariant to all common transformation, thus it allows to carry out fare comparisons of saliency models.

Note that the [previous version](https://github.com/merofeev/ss_robust_metric) of this algorithm does not guarantee optimal transformation and the optimization algorithm has significantly more computational complexity.

Features
---------

 - Metric is invariant to
   - Any brightness correction function
   - Mixing method with Center Prior
 - Finds the optimal transformation of a saliency model to fit GT saliency according to PSNR measure

Links
-------
The paper is available [here](http://compression.ru/video/savam/pdf/A_semiautomatic_saliency_model_and_its_application_to_video_compression_ICCP_2017_0.pdf).

More details are available at [project page](http://compression.ru/video/savam/#downloads).

Robust Saliency Map Comparison, version 1: https://github.com/merofeev/ss_robust_metric
