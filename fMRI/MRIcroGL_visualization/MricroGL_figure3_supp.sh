import gl
gl.resetdefaults()
gl.view(16)

#open background image
gl.loadimage('/mnt/HR_project_SZU/preprocess_by_spm12/MRIcroGL_visualization/mni_icbm152_t1_tal_nlin_sym_09b_hires_brain.nii.gz')
gl.minmax(0,5, 120)
#open overlay: show positive regions
gl.overlayload('/mnt/HR_project_SZU/preprocess_by_spm12/GLM1_spm_2nd/err_infer_instr_FWEc181.nii')
gl.overlayloadsmooth(1)
gl.overlayload('/mnt/HR_project_SZU/preprocess_by_spm12/GLM1_spm_2nd/rew_infer_instr_FWEc219.nii')
gl.overlayloadsmooth(1)

gl.backcolor(255,255,255)
gl.minmax(1,3.16, 5)
gl.minmax(2,3.1, 5.58)
gl.colorname(1, 'HOTIRON')
gl.colorname(2, '5winter')
#gl.colorbarposition(1)
gl.colorbarposition(2)
gl.colorbarsize(0.05)
gl.colorbarcolor(255,255,255,255)

gl.opacity(1,80)
gl.opacity(2,90)
#gl.orthoviewmm(-6,22,50)
#gl.mosaic('S -4')
gl.mosaic('A L+ 40,   C L+ 22,   S L+ -2')


gl.linewidth(0)
#gl.cutout(0.5, 0.5, 0.5, 0, 1,1)
gl.zoomscale(1)

gl.zerointensityinvisible(90,0)
gl.windowposition(200, 100, 1400, 350)
gl.savebmp('/mnt/HR_project_SZU/preprocess_by_spm12/GLM1_spm_2nd/rew_infer_instr_FWEc219.png')

