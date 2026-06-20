############################# all neg.
import gl
gl.overlaycloseall
gl.resetdefaults()
gl.view(16)

#open background image
gl.loadimage('/mnt/HR_project_SZU/preprocess_by_spm12/MRIcroGL_visualization/mni_icbm152_t1_tal_nlin_sym_09b_hires_brain.nii.gz')
gl.minmax(0,15, 120)
#open overlay: overlap only positive regions
#dif
gl.overlayload('/mnt/HR_project_SZU/preprocess_by_spm12/GLM2_spm_2nd_nerDif_Dif/dif_n_FWEc211.nii')
gl.overlayloadsmooth(1)


gl.backcolor(255,255,255)
gl.minmax(1,3.1, 4.84)
gl.colorname(1, '6bluegrn')
#gl.colorfromzero(2,1)

gl.colorbarposition(1)
gl.colorbarsize(.05)
gl.colorbarcolor(255,255,255,255)

gl.opacity(1,80)
#gl.orthoviewmm(-4,-14,28)
#gl.mosaic('A L+ 50,   C L+ 22,   S L+ -6')
gl.mosaic('S -4')

gl.linewidth(0)
gl.zoomscale(2)

gl.zerointensityinvisible(90,0)
gl.windowposition(200, 100, 1400, 350)
gl.savebmp('/mnt/HR_project_SZU/preprocess_by_spm12/GLM2_spm_2nd_nerDif_Dif/dif_n_FWEc211.png')

