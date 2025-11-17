# MVPA; python3.10
# 2024/09/10 by Chuanyong Xu (refer to ning mei. @author: nmei; Created on Mon Apr 17 11:46:05 2023)
# rule_list, rule_response, difficulties, confidence, rule_swirch
import nilearn
import os,gc
import numpy as np
import pandas as pd
import nibabel as nib
from glob import glob
from nilearn.maskers import NiftiMasker
from nilearn.image import load_img
from nilearn import plotting, image
from matplotlib import pyplot as plt
#from utils import get_sphere_data,concept_mapping
from sklearn.model_selection import check_cv,cross_validate
from sklearn.dummy import DummyClassifier
from sklearn.svm import LinearSVC
from sklearn.calibration import CalibratedClassifierCV
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import VarianceThreshold
#from sklearn.pipeline import make_pipeline # ******for non upsampling only******
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedShuffleSplit
from scipy.stats import ttest_1samp
try:
    from sklearn.utils.testing import ignore_warnings
    from sklearn.exceptions import ConvergenceWarning
except:
    from sklearn.utils._testing import ignore_warnings
    from sklearn.exceptions import ConvergenceWarning
from joblib import Parallel,delayed

# from imblearn.under_sampling import RandomUnderSampler
from imblearn.over_sampling import RandomOverSampler
from imblearn.pipeline import make_pipeline # ******for upsampling only******
#==============================================================================
#------------------------------define function 1-------------------------------
from numba import njit, prange
@njit(parallel=True)
def detec_outliers_and_interpolate(data):
    n_timerpoints, n_voxels = data.shape
    for voxel in prange(n_voxels):
        y=data[:,voxel]
        outliers=np.abs(y)>=3
        y[outliers]=np.nan
        
        nans=np.isnan(y)
        if np.all(nans):
            continue
        not_nans=~nans
        x=np.arange(n_timerpoints)
        y[nans]=np.interp(x[nans],x[not_nans],y[not_nans])
        data[:,voxel]=y
    return data


from nilearn import masking
from nilearn.image.resampling import coord_transform
try:
    from nilearn.input_data.nifti_spheres_masker import _apply_mask_and_get_affinity
except:
    from nilearn.maskers.nifti_spheres_masker import _apply_mask_and_get_affinity
def get_sphere_data(mask, BOLD_data, radius,):
    if type(BOLD_data) == str:
        BOLD_data = nilearn.image.load_img(BOLD_data)
    process_mask, process_mask_affine = masking.load_mask_img(mask)
    process_mask_coords = np.where(process_mask !=0)
    process_mask_coords = coord_transform(process_mask_coords[0],
                                          process_mask_coords[1],
                                          process_mask_coords[2],
                                          process_mask_affine,
                                          )
    process_mask_coords = np.asarray(process_mask_coords).T
    
    X, A = _apply_mask_and_get_affinity(seeds=process_mask_coords,
                                        niimg=BOLD_data,
                                        radius=radius,
                                        allow_overlap=True,
                                        mask_img = mask,
                                        )
    return X, A, process_mask, process_mask_affine, process_mask_coords 

# ---------------------------------------seeting1-------------------------------
@ignore_warnings(category=ConvergenceWarning)
def LOO_decode(row, BOLD_signals, labels, cv, chance = False):
    features = BOLD_signals[:, row]
    try:
        if chance:
            #applying VarianceThreshold method (alterative methods:PCA) to extract features: removing low variance voxels
            pipeline = make_pipeline(VarianceThreshold(),
                                     StandardScaler(), # normalize: mean=0, variance=1
                                     DummyClassifier(strategy='uniform',
                                                     random_state=123,
                                                     ))
            
        else:
            svm = LinearSVC(penalty='l1', ##
                            dual=False,
                            class_weight="balanced",
                            random_state=123,
                            )    
            svm = CalibratedClassifierCV(svm, method='sigmoid',
                                         cv=5,)
            # pipeline = make_pipeline(VarianceThreshold(), 
            #                          StandardScaler(),
            #                          svm,
            #                          )
            #up-sampling**********************************************************
            pipeline = make_pipeline(VarianceThreshold(), 
                                     StandardScaler(),
                                     RandomOverSampler(sampling_strategy='all',random_state=123),#******
                                     svm,) 
            
        res = cross_validate(pipeline,
                             features,
                             labels,
                             scoring='roc_auc_ovr', # 'roc_auc',
                             cv=cv,
                             n_jobs=1,
                             verbose=0,
                             return_estimator=True,
                             )
        estimators = res['estimator']    
        idxs_test = [idx_test for _, idx_test in cv.split(features, labels,)]
        scores = [roc_auc_score(labels[idx_test], est.predict_proba(features[idx_test]), 
                                multi_class='ovr', average='macro') for idx_test,est in zip(idxs_test,estimators)]
    except:
        scores = [0.5 for train,test in cv.split(features, labels)]
        print('broken')
    if chance:
        return np.array(scores)
    else:
        outer_coefs = []
        for est in estimators:
            calibrator = est.named_steps['calibratedclassifiercv']
            inner_coefs = np.stack([cc.estimator.coef_
                for cc in calibrator.calibrated_classifiers_
                ],axis=0)
            # mean_inner = np.abs(inner_coefs).mean(axis=0)
            mean_inner = inner_coefs.mean(axis=0)
            outer_coefs.append(mean_inner)
        outer_coefs = np.stack(outer_coefs,axis=0)
        weights_arr = outer_coefs.mean(axis=0).T
        return np.array(scores), weights_arr


#------------------------------define function 2-------------------------------
class MyStruct:
    pass
common_settings = MyStruct()
common_settings.radius = 10 #in mm ######### change*
common_settings.model_name = 'mvpa_Ners_KFolds' 
common_settings.label_map = {0: 0, 1: 1, 2:2} ######### change*
common_settings.concepts_map = {0:'0er', 1:'1er', 2:'2er'} ######### change*

#==============================================================================
# train and test
subs=[1108, 1110, 1114, 1157]
for sub in subs:
    if __name__ == "__main__":
        [os.remove(item) for item in glob('../*/core*')]
        radius = common_settings.radius
        cond_cat  = 'trials'
        condition = 'type_dev_error_1' #(ners) according to experiment conditions to change ################# change*
        feature_type = common_settings.model_name
        sub_num = sub;
        sub = f'sub{sub:03d}'
        os.getcwd()
        working_dir = f'../data/RawData/{sub}' ################## change*
        mask_dir = '../data/barinmask/'         ################## change*
        results_dir = '../results/KFolds_MNI_space_Ners_raws_upsamp10' ################## change*
#   --------------- 
        if not os.path.exists(results_dir):
            os.mkdir(results_dir)
        working_data = os.path.join(working_dir,
                                    'whole_brain_average.nii.gz')
        event_file = os.path.join(working_dir,
                                  f'data_raw_decodingRSA_confidence_{sub_num}.csv')
        wholebrain_mask = os.path.join(mask_dir, 'mask_mu_dif_ner.nii.gz')
        
#define cross validation procedure
        df_data = pd.read_csv(event_file)
        df_data.loc[df_data[condition]>=2,condition]=2 # ******************
        labels = df_data[condition].map(common_settings.label_map).values # label linked to values 
        
        
        masker = NiftiMasker(mask_img = wholebrain_mask,
                             detrend=True,
                             high_pass=1/128,
                             t_r=0.85,
                             standardize=True, #standardize each run separately
                             n_jobs=-1,
                             ).fit() # get the bold data only in brain mask
        if os.path.exists(os.path.join(working_dir,'whole_brain_average.nii.gz')):
            BOLD_signals = nib.load(os.path.join(working_dir,'whole_brain_average.nii.gz'))
        else: 
        # load bold data and csv
            bold_data = []
            rp_files = np.sort(glob(os.path.join('/sharehome/xuchuanyong/Documents/DATA_analysis/MRI_SZU_HR',
                                                 'preprocess_by_spm12','fun_s1_4d0*',f'{sub}','rp_afun_4D.txt')))
            bold_files = np.sort(glob(os.path.join('/sharehome/xuchuanyong/Documents/DATA_analysis/MRI_SZU_HR',
                                                   'preprocess_by_spm12','fun_s1_4d0*',f'{sub}','wrafun_4D.nii')))
            for bold_file, rp_file in zip(bold_files, rp_files):
                nii_bold = nib.load(bold_file)
                #nii_bold = masker.transform(nii_bold)
                #----------------       
                confounds = np.loadtxt(rp_file) #head motions parameters
                nii_bold  = masker.transform(nii_bold, confounds=confounds) # head motion correction
                nii_bold  = detec_outliers_and_interpolate(nii_bold)
                #----------------
                nii_bold = masker.inverse_transform(nii_bold)
                temp_bold = nii_bold.get_fdata()
                bold_data.append(temp_bold)
                del nii_bold
                gc.collect()
            print('Data in mask extracted done')
            bold_data=np.concatenate(bold_data, axis=-1)
            
            del temp_bold
            gc.collect()
            print('Data concatenate done')
        # prepare the whole brain BOLD signal that are averaged from all sessions
            bold_average = []
            for _conf, df_sub in df_data.groupby([cond_cat]):
                for index, tr in df_sub.iterrows():
                    numbers = pd.DataFrame(tr['time_indices'].split('+'), columns=['numbers']).astype(int)
                    temp = bold_data[..., numbers]
                    temp = temp[..., 0]
                    bold_average.append(temp.mean(-1)) # mean for each trial
        
            bold_average = [arr.ravel() for arr in bold_average]
            bold_average = np.vstack(bold_average)
            mask_img_data = masker.mask_img_.get_fdata()
            valid_voxels_mask = mask_img_data>0
            print(np.sum(mask_img_data>0))
            bold_average = bold_average[:,valid_voxels_mask.ravel()] # *******only includ voxels with value > 0   
            BOLD_signals = masker.inverse_transform(bold_average)
            BOLD_image = np.asanyarray(BOLD_signals.dataobj)
            print(BOLD_image.shape)

            nib.save(BOLD_signals, os.path.join(working_dir,'whole_brain_average.nii.gz'))

        # extract the sphere voxels
        print('extracting sphere voxels ...')
        if os.path.exists(os.path.join(working_dir, 'rows_raw_mu10.npy')):
            rows = np.load(os.path.join(working_dir, 'rows_raw_mu10.npy'), allow_pickle=True)
        else:
            (X, A, process_mask, process_mask_affine, process_mask_coords) = get_sphere_data(load_img(wholebrain_mask),
                                                                                             load_img(working_data), 
                                                                                             radius,
                                                                                             )
            rows = A.rows
            np.save(os.path.join(working_dir,'rows_raw_mu10.npy'), rows, allow_pickle=True)
            del X,A,process_mask,process_mask_affine,process_mask_coords

# for MVPA train and test, searchlight
        masker = NiftiMasker(wholebrain_mask,).fit()
        BOLD_signals = masker.transform(working_data)
        #words = df_data[condition].map(common_settings.concepts_map).values # characters of label 
        cv = StratifiedShuffleSplit(n_splits = 100,test_size = .2,random_state = 123) #times of repeat
        row_groups = np.array_split(rows, 1) ##****** parcellate the brain data for faster processing
        
        for ii, row_group in enumerate(row_groups):
            score_file_name = os.path.join(results_dir,
                                           f'{sub}_{condition}_{radius}mm_{ii+1}_score.npy')
            chance_file_name = score_file_name.replace('score','chance')
            weights_file_name = score_file_name.replace('score','weights')
            if not os.path.exists(score_file_name):
                print(f'start {ii+1}')
                # decoding
                gc.collect()
                results = Parallel(n_jobs = -1, verbose = 1)(delayed(LOO_decode)(*[
                    row, BOLD_signals, labels, cv, False]) for row in row_group)
                
                scores, weights=zip(*results)
                scores =  np.stack([r[0] for r in results], axis=0)
                # for idx, arr in enumerate(weights):
                #     print(f"weights[{idx}].shape={arr.shape}")
                weights = np.stack([
                    # np.abs(arr).mean(axis=0)
                    np.mean(arr,axis=0)
                    for arr in weights],axis=0)

                gc.collect()
                # chance level decoding
                chances = Parallel(n_jobs = -1, verbose =1)(delayed(LOO_decode)(*[
                    row, BOLD_signals, labels, cv, True]) for row in row_group)
                gc.collect()
                
                np.save(score_file_name, np.array(scores))
                np.save(chance_file_name, np.array(chances))
                
                # extract the weights of each catogery in SVC prediction
                np.save(weights_file_name, np.array(weights,dtype=float))
            else:
                print(f'{score_file_name} exits')

    #==============================================================================
    # indivudual level: one sample t-test: test - chance, to build a p-value/t-value nii map
    scores_list = []
    scores      = []
    scores_list = np.sort(glob(os.path.join(results_dir, f'{sub}_{condition}_{radius}mm_*_score.npy')))
    scores = np.concatenate([np.load(score, allow_pickle=True) 
                             for score in scores_list ], axis=0)
    chances_list = []
    chances      = []    
    chances_list = np.sort(glob(os.path.join(results_dir, f'{sub}_{condition}_{radius}mm_*_chance.npy')))
    chances = np.concatenate([np.load(chance, allow_pickle=True) 
                              for chance in chances_list ], axis=0)
    
    scores_chances = scores - chances
    scores_chances_brain = np.mean(scores_chances, axis=1)
    scores_chances_brain = masker.inverse_transform(scores_chances_brain.reshape(1,-1))
    scores_chances_brain_name = os.path.join(results_dir,
                                           f'{sub}_{condition}_{radius}mm_scores_chances_brain.nii.gz')
    nib.save(scores_chances_brain, scores_chances_brain_name)
    
    weights_list = []
    weights      = []
    weights_list = np.sort(glob(os.path.join(results_dir, f'{sub}_{condition}_{radius}mm_*_weights.npy')))
    weights = np.concatenate([np.load(weight, allow_pickle=True) 
                             for weight in weights_list ], axis=0)
    weights0 = weights[:,0]
    weights1 = weights[:,1]
    weights2 = weights[:,2]
    weights0 = masker.inverse_transform(weights0.reshape(1,-1))
    weights1 = masker.inverse_transform(weights1.reshape(1,-1))
    weights2 = masker.inverse_transform(weights2.reshape(1,-1))
    weights_brain_name = os.path.join(results_dir,
                                           f'{sub}_{condition}_{radius}mm_weights0_brain.nii.gz')
    nib.save(weights0, weights_brain_name)
    weights_brain_name = os.path.join(results_dir,
                                           f'{sub}_{condition}_{radius}mm_weights1_brain.nii.gz')
    nib.save(weights1, weights_brain_name)      
    weights_brain_name = os.path.join(results_dir,
                                           f'{sub}_{condition}_{radius}mm_weights2_brain.nii.gz')
    nib.save(weights2, weights_brain_name)      
    
    # tvalues, pvalues = [], []
    # tvalues, pvalues = ttest_1samp(scores_chances, 0, axis=1)
    # tvalues_brain = masker.inverse_transform(tvalues.reshape(1, -1))
    # pvalues_brain = masker.inverse_transform(pvalues.reshape(1, -1))
    # t_file_name = os.path.join(results_dir,
    #                                        f'{sub}_{condition}_{radius}mm_tvalues_brain.nii.gz')
    # p_file_name = t_file_name.replace('tvalues_brain','pvalues_brain')
    # nib.save(tvalues_brain,t_file_name)#save coverted nii
    # nib.save(pvalues_brain,p_file_name)#save coverted nii
    
    #plotting------
    # plotting.plot_stat_map(scores_chances_brain,
    #                        wholebrain_mask,
    #                        threshold = 0,
    #                        draw_cross = False,
    #                        cmap = plt.cm.coolwarm,
    #                        vmax = 0.2,
    #                        )
    # plotting.show()

#==============================================================================
# group level: FSL randomise for one sample t-test
