# context_models
Models and data files for our 2016 Cell Systems publication: 

This small repository of files documents the experimental data and mathematical models used in the paper titled "The Effect of Compositional Context in Synthetic Gene Networks" by Enoch Yeung et al. 
The model files are in the subfolder /models and coded as MATLAB scripts, i.e. .m files.  

To get started, open an instance of MATLAB and run the script RFPCFPmain.m or mSpinachMG.m.   You should see a series of plots, in each instance:

*two plots showing the expression of the gene system against time, 
*three plots showing the temporal dynamics of the underlying supercoiling states. 

The core functions are ConvModelRFPCFP_Journal.m, DivModelRFPCFP_Journal.m and TandModelRFPCFP_Journal.m.  These scripts define ODE models that capture the feedback between supercoiling and gene expression.  The derivation of the model is provided in the Methods (Supplement) section of the paper.  Feel free to contact me at eyeung@caltech.edu if you would like a walk-through or have any questions on the derivation.  

All three of these core functions rely on kcatmaxfunc.m, kTLfunc.m, kfmaxfunc.m, and m_gyrtopfunc.m to function.  Only kcatmaxfunc, kfmaxfunc, and m_gyrtopfunc are required for the supercoiling model. The function kTLfunc.m models the slow-down in translation in cell-free batch reactions, as amino acids, ATP and metabolites deplete from synthetic biocircuit expression.   More details on this process can be found in Z. Tuza et al. 2014.  A full model of resource depletion in cell-free batch reactions is quite complex, for this work we used an exponential decay model to capture the first-order effects.  

The plasmids for these experiments are also available on AddGene.  The data for the experiments is posted in separate folders for each set of experiments; if you need help understanding or processing the data, I can help. Feel free to email me directly if you have any questions at enoch.yeung@pnnl.gov or eyeung@caltech.edu. 

Happy coding!

Enoch 







