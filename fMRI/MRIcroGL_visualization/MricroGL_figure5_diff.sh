############################# ners pos.
import gl
gl.resetdefaults()
gl.view(16)

#open background image
gl.loadimage('/mnt/HR_project_SZU/preprocess_by_spm12/MRIcroGL_visualization/mni_icbm152_t1_tal_nlin_sym_09b_hires_brain.nii.gz')
gl.minmax(0,15, 120)
#open overlay: show positive regions
gl.overlayload('/mnt/HR_project_SZU/preprocess_by_spm12/Decoding_mvpa/results/KFolds_MNI_space_Ners_raws_upsamp10_mu/all_onesampT_p_vox_corrp_tstat1.nii.gz')
gl.overlayloadsmooth(0)
gl.overlayload('/mnt/HR_project_SZU/preprocess_by_spm12/Decoding_mvpa/results/KFolds_MNI_space_Diffs_raws10_mu/all_onesampT_p_vox_corrp_tstat1.nii.gz')
gl.overlayloadsmooth(0)

gl.backcolor(255,255,255)
gl.minmax(1,0.95, 1)
gl.minmax(2,0.95, 1)
#gl.minmax(2,0.95, 0.95)
gl.colorname(1, '6warm')
gl.colorname(2, '6bluegrn')
#gl.colorbarposition(1)
gl.colorbarsize(.1)
gl.colorbarcolor(255,255,255,255)
gl.colorbarposition(0)

gl.opacity(1,0)
gl.opacity(2,100)
#gl.orthoviewmm(6,22,46)
#gl.view(2)
gl.mosaic('C 22')
#gl.mosaic('A L+ 50,   C L+ 22,   S L+ -4')


gl.linewidth(0)
#gl.cutout(0.5, 0.5, 0.5, 0, 1,1)
gl.zoomscale(6)

gl.zerointensityinvisible(90,0)
gl.windowposition(200, 100, 1200, 650)
gl.savebmp('/mnt/HR_project_SZU/preprocess_by_spm12/Decoding_mvpa/results/mvpa_in_conf_diff.png')


