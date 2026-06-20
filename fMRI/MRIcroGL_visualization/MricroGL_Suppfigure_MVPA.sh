############ MVPA decoding precision for Ners
import gl
gl.overlaycloseall
gl.resetdefaults()
gl.view(16)

#open background image
gl.loadimage('/mnt/HR_project_SZU/preprocess_by_spm12/MRIcroGL_visualization/mni_icbm152_t1_tal_nlin_sym_09b_hires_brain.nii.gz')
gl.minmax(0,15, 120)
#open overlay: show positive regions
gl.overlayload('/mnt/HR_project_SZU/preprocess_by_spm12/Decoding_mvpa/results/KFolds_MNI_space_Ners_raws_upsamp10_mu/all_scores_Gmean_p_p095.nii.gz')
gl.overlayloadsmooth(1)

gl.backcolor(255,255,255)
gl.minmax(1,0.5, 0.59)
gl.colorname(1, 'blue2red')
#gl.colorbarposition(1)
gl.colorbarposition(1)
gl.colorbarsize(.1)
gl.colorbarcolor(255,255,255,255)


gl.opacity(1,80)
gl.orthoviewmm(-6,22,44)
#gl.mosaic('S -4')
gl.mosaic('A L+ 44,   C L+ 22,   S L+ -6')


gl.linewidth(0)
#gl.cutout(0.5, 0.5, 0.5, 0, 1,1)
gl.zoomscale(1)

gl.zerointensityinvisible(90,0)
gl.windowposition(200, 100, 1400, 350)
gl.savebmp('/mnt/HR_project_SZU/preprocess_by_spm12/Decoding_mvpa/results/Decoding_precision_ners.png')




############ MVPA decoding precision for Diffs
import gl
gl.overlaycloseall
gl.resetdefaults()
gl.view(16)

#open background image
gl.loadimage('/mnt/HR_project_SZU/preprocess_by_spm12/MRIcroGL_visualization/mni_icbm152_t1_tal_nlin_sym_09b_hires_brain.nii.gz')
gl.minmax(0,15, 120)
#open overlay: show positive regions
gl.overlayload('/mnt/HR_project_SZU/preprocess_by_spm12/Decoding_mvpa/results/KFolds_MNI_space_Diffs_raws10_mu/all_scores_Gmean_p_p095.nii.gz')
gl.overlayloadsmooth(1)

gl.backcolor(255,255,255)
gl.minmax(1,0.5, 0.52)
gl.colorname(1, 'blue2red')
#gl.colorbarposition(1)
gl.colorbarposition(1)
gl.colorbarsize(.1)
gl.colorbarcolor(255,255,255,255)

gl.opacity(1,80)
gl.orthoviewmm(-6,22,44)
#gl.mosaic('S -4')
gl.mosaic('A L+ 44,   C L+ 22,   S L+ -6')


gl.linewidth(0)
#gl.cutout(0.5, 0.5, 0.5, 0, 1,1)
gl.zoomscale(1)

gl.zerointensityinvisible(90,0)
gl.windowposition(200, 100, 1400, 350)
gl.savebmp('/mnt/HR_project_SZU/preprocess_by_spm12/Decoding_mvpa/results/Decoding_precision_diffs.png')


