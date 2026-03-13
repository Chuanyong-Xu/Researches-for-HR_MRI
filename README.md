# Resarches-for-HR_MRI, the code and data for main results has been upload here
  
  (1) for behavioral analysis and Bayesian modelling
  ------------------ folder: MRI_behavioural
  
  1) run 'm_1_psychometric_switch.m' for Bayesian modelling;
  2) run the .py scripts in each of the subfolder in Model_fitting_Tom folder, for preparing the data and alterative model fitting;
  3) run 'm_5_CollectData_for_ModelCompare.m' for model comparing (LL, AIC, BIC);
  
  (2) for RNN simulation
  ------------------ folder: RNN/DeepRL_RNN/script
  1) run 'trian_valid.py' for training;
  2) run 'test.py' for testing;
  
  ------------------ folder: RNN
  1) run 'S0_CollectData_ForRNN.m' for organizing data;
  2) run 'S1_prSwi_MDS.m' for plotting the RNN performance compared to human;
  3) run 'S2_prepare_RNNforRDM.m' for extracting the hidden layer activations for RSA analysis;
  4) run 'S03_1_Decoding_mvpa_ners_rois_individuals.m' for Gaussian fitting for errors;
  5) run 'S03_2_Decoding_mvpa_diffs_rois_individuals.m' for Gaussian fitting for difficulty;
  6) run 'S04_Decoding_mvpa_NersDiffs_rois_individuals.m' for 3-D visualization for combination of the errors and difficulty.

  (3) GLM analysis for fMRI, MVPA for fMRI, RSA for fMRI and RNN
  ------------------ folder: fMRI
  Due to the limitation of upload, in this folder, only the 2nd results of each GLM were upload.

  
  1) GLM analysis code for MRI in folder of ---spmGLM---, the order was marker by 'Sxxxxx.m';
  2) MVPA code for MRI in folder of ---MVPA---, the order was marker by 'Sxxxxx.py';
  3) RSA or Geometry analysis code for fMRI and RNN in folder of ---RSA_Geometry---.
