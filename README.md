# MatFiber version of Fiber3 algorithm for automated fiber orientation analysis
In 1998, the Cardiac Mechanics Research Group at UCSD (CMRG) published an algorithm for automated analysis of muscle fiber orientations in histologic slides: [Karlon WJ, Covell JW, McCulloch AD, Hunter JJ, Omens JH. Automated measurement of myofiber disarray in transgenic mice with ventricular expression of ras. Anat Rec 252(4):612-625, 1998](http://www.ncbi.nlm.nih.gov/pubmed/9845212). The CMRG named their implementation Fiber3, and a version of Fiber 3 is available as part of their [freely available modeling software Continuity](http://continuity.ucsd.edu/Continuity).

Our group has used the algorithm primarily for analyzing collagen fiber orientations in picrosirius red-stained sections of myocardial scar tissue. We developed a MATLAB implementation of the algorithm named MatFiber, and verified that MatFiber gives the same results as Fiber3 when analyzing the same image using the same settings.

This repository contains the MatFiber code as well as some sample images to practice using the analysis.

Please note that the results of your analysis will depend heavily on your choice of subregion size (the algorithm produces one average orientation vector for each subregion) and your choice of threshold (the threshold excludes subregions with weaker overall image gradient information)! It is very important to test a range of subregion sizes and thresholds to establish parameters that produce an accurate analysis of your images. More information about these parameters is provided in comments within the MatFiber.m file.

If you use our MatFiber code, please cite the paper where we first employed it: [Fomovsky GM, Holmes JW. Evolution of scar structure, mechanics, and ventricular function after myocardial infarction in the rat. Am J Physiol Heart Circ Physiol 298:H221-H228, 2010](http://www.ncbi.nlm.nih.gov/pubmed/19897714).
