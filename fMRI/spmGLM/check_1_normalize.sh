cd /mnt/HR_project_SZU/preprocess_by_spm12/fun_s1_4d03 ######
ls -l | grep ^d | awk '{print $9}' > list.txt

for file in  `cat list.txt`
  do
  echo $file
  cd /mnt/HR_project_SZU/preprocess_by_spm12/fun_s1_4d03 ######
  cd $file
  structural_image=/home/hello/fsl/data/standard
  functional_image=/mnt/HR_project_SZU/preprocess_by_spm12/fun_s1_4d03/$file ######
  output_dir=/mnt/HR_project_SZU/preprocess_by_spm12/fun_s1_4d03_check_quality/$file ######

#mkdir
  mkdir -p $output_dir
  cd $output_dir
  fslroi $functional_image/swrafun_4D.nii $output_dir/swrafun_4D_1st 0 1

##flirt for check registration
#  flirt -in $functional_image/swrafun_4D_1st -ref $structural_image/MNI152_T1_2mm.nii.gz -out $output_dir/functional_to_MNI152.nii.gz -omat $output_dir/functional_to_MNI152.mat
#  echo "flirt done"


#photo for check
  slicesdir -p $structural_image/MNI152_T1_2mm.nii.gz $output_dir/swrafun_4D_1st

#plot head motion
  done

  rm /mnt/HR_project_SZU/preprocess_by_spm12/fun_s1_4d03/list.txt ######
#echo "photo has been generated"

