############################# all neg.
import gl
gl.overlaycloseall
gl.resetdefaults()
gl.view(16)

#open background image
gl.loadimage('/mnt/HR_project_SZU/preprocess_by_spm12/MRIcroGL_visualization/mni_icbm152_t1_tal_nlin_sym_09b_hires_brain.nii.gz')
gl.minmax(0,15, 120)
#open overlay: overlap only positive regions
#ner pos
gl.overlayload('/mnt/HR_project_SZU/preprocess_by_spm12/GLM2_spm_2nd_nerDif_ner/ner_p_FWEc203.nii')
gl.overlayloadsmooth(1)
#conf pos
gl.overlayload('/mnt/HR_project_SZU/preprocess_by_spm12/GLM2_spm_2nd_muSwitch/mu_p_FWEc227.nii')
gl.overlayloadsmooth(1)
#dif
gl.overlayload('/mnt/HR_project_SZU/preprocess_by_spm12/GLM2_spm_2nd_nerDif_Dif/dif_n_FWEc211.nii')
gl.overlayloadsmooth(1)


gl.backcolor(255,255,255)
gl.minmax(1,3.1, 6.35)
gl.minmax(2,3.1, 5.38)
gl.minmax(3,3.1, 4.84)
gl.colorname(1, '6warm')
gl.colorname(2, '4hot')
gl.colorname(3, '6bluegrn')
#gl.colorfromzero(2,1)

#gl.colorbarposition(1)
gl.colorbarsize(.1)
gl.colorbarcolor(255,255,255,255)
gl.colorbarposition(0)

gl.opacity(1,80)
gl.opacity(2,90)
gl.opacity(3,80)
#gl.orthoviewmm(-4,-14,28)
#gl.mosaic('A L+ 50,   C L+ 22,   S L+ -6')
gl.mosaic('S -8')

gl.linewidth(0)
gl.zoomscale(2)

gl.zerointensityinvisible(90,0)
gl.windowposition(200, 100, 1400, 350)
gl.savebmp('/mnt/HR_project_SZU/preprocess_by_spm12/GLM2_spm_2nd_muSwitch/all_pos_only.png')

