#
cd /mnt/HR_project_SZU/Preprocess/T1Raw
ls -l | grep ^d | awk '{print $9}' > list.txt
#
mkdir T1Img/

for file in `cat list.txt`
do
  echo $file
  mkdir T1Img/$file
  dcm2nii $file/
  cd $file
  mv co*.nii.gz ..
  rm *.nii.gz
  cd ..
  gzip -d co*.nii.gz
  mv co*.nii T1Img/$file
done

#  mv T1Img/ ..
#  cd ..
  cp -r T1Img/ /mnt/HR_project_SZU/preprocess_by_spm12/T1Img01
  cp -r T1Img/ /mnt/HR_project_SZU/preprocess_by_spm12/T1Img02
  cp -r T1Img/ /mnt/HR_project_SZU/preprocess_by_spm12/T1Img03
  cp -r T1Img/ /mnt/HR_project_SZU/preprocess_by_spm12/T1Img04
  cp -r T1Img/ /mnt/HR_project_SZU/preprocess_by_spm12/T1Img05
  mv T1Img/ /mnt/HR_project_SZU/preprocess_by_spm12/T1Img06
