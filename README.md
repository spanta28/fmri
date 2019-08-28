# web_fmri
Runs fmri pipeline using singularity image 

Sample_input_data is provided in this repo.

Optional parameters info: The optional params info is detailed in fmri optional params description.xlsx

json input file is required.
Ex json input to run_script.py:
{"data":[["/home/vagrant/input1","sub1","sess1"],["/home/vagrant/input2","sub2","sess2"]], "options":{    
"reorient_params_x_mm":0, 
"reorient_params_y_mm":0, 
"reorient_params_z_mm":0,
"reorient_params_pitch":0, 
"reorient_params_roll":0, 
"reorient_params_yaw":0,    
"realign_fwhm":8,
    "realign_interp":2,
    "realign_quality":1,
    "realign_register_to_mean":true,
    "realign_separation":4,
    "realign_wrap": [0,0,0],
    "realign_write_interp": 4,
    "realign_write_mask": true,
    "realign_write_which": [2, 1],
    "realign_write_wrap": [0,0,0],
    "slicetime_ref_slice":16,
"normalize_affine_regularization_type":"mni",    
    "normalize_write_bounding_box": [[-78, -112, -70],[78, 76, 85]],
    "normalize_write_interp": 1,
    "normalize_write_voxel_sizes": [3,3,3],
"smoothing_x_mm":10.0,"smoothing_y_mm":10.0,"smoothing_z_mm":10.0,"smoothing_implicit_masking":false
}}



Note: 
trendscenter/fmri_web_image hosted on dockerhub is private. So contact me with your docker ID so that I can give you access to it. 

You can run the code in foll. ways:
1) Build dockerfile file from here 
docker build -t fmri_image_name Dir_containing_Dockerfile


Sometimes you may create a docker container on your local machine and when you move it elsewhere you may not have permissions to run somethings inside that container. So go inside the container and run "chmod -R 755 /" or equivalent to make sure your required software is executable.

To change the cache dir location , incase you dont have permissions to write to the default
in bash:
export SINGULARITY_CACHEDIR="/user/container"

same with writing to default /tmp folder
export SINGULARITY_TMPDIR="/user/container/tmp"

then run singularity build fmri.img trendscenter/fmri_web_image
if you want to create the fmri.img as writable image for editing , then :
singularity build --sandbox fmri.img docker://trendscenter/fmri_web_image
or if you have permissions to use --writable then:
singularity build --writable fmri.img docker://trendscenter/fmri_web_image

This is what I used and it works fine.

Sometimes due to changes in open source code of packages being used in this container, things may not work as expected. 
If you want to use the image that I already built then proceed to step 2) 

2) use trendscenter/fmri_web_image located on dockerhub

To change the cache dir location , incase you dont have permissions to write to the default 
in bash:
export SINGULARITY_CACHEDIR="/export/user/container"

same with writing to default /tmp folder
export SINGULARITY_TMPDIR="/user/container/tmp"

singularity build fmri.img docker://trendscenter/fmri_web_image

if you want to create the fmri.img as writable image for editing , then :
singularity build --sandbox fmri.img docker://trendscenter/fmri_web_image
or if you have permissions to use --writable then:
singularity build --writable fmri.img docker://trendscenter/fmri_web_image

This is what I used and it works fine.

3) Once downloaded to local machine:

3.a) You can use run_script.py and supply required args to have the script call singularity image and run the pipeline or do it manually as in step 3.b)

 python run_script.py -h
usage: run_script.py [-h] -json JSON -o OUTPUT -tmp TMP -image_path IMAGE_PATH
                     [-g SLURM] [-type {nifti,dicoms}]

Run fmri pipeline on mprage T1w dicoms

optional arguments:
  -h, --help            show this help message and exit

required arguments:
  -json JSON, --json JSON
                       
  -o OUTPUT, --output OUTPUT
                        Existing output directory to store vbm outputs
  -tmp TMP              Nipype tmp directory, it can be same as output
                        directory
  -image_path IMAGE_PATH
                        Enter full path to singularity vbm image name

optional arguments:
  -g SLURM, --slurm SLURM
                        Use slurm, True or False, default=False
  -type {nifti,dicoms}, --data_type {nifti,dicoms}
                        input data type. Ex:'nifti','dicoms'

python run_script.py -image_path /container/gsu_fmri.img/ -i /container/input -o /container/output -json /container/json_data/input.json -tmp /output -sub sub1 -sess sess1
args passed are:
Namespace(data_type='nifti', image_path=['/container/gsu_fmri.img/'], input=['/container/input'], output=['/container/output'], sess=['sess1'], slurm=False, smooth=None, sub=['sub1'], tmp=['/container/output'])
Executed command is: singularity exec --contain --workdir /container/output -B /container/input:/data,/container/output:/out /container/gsu_fmri.img/ python3 /computation/run_fmri.py -sub sub1 -sess sess1 -json '/json_data/input.json'

fmri pipeline is running...
fmri output directory: /out/sub1/sess1/anat
fmri preprocessing completed.

3.b) Run singularity
singularity exec --contain --workdir NIPYPE_TMP_DIR -B INPUT_DIR:/data,OUTPUT_DIR:/out fmri.img python3 /computation/run_fmri.py -sub SUBJECT_ID -sess SESSION_ID -json '/json_data/input.json'

    -c/--contain        This option disables the sharing of filesystems on
                        your host (e.g. /dev, $HOME and /tmp).
                        
    --workdir stores tmp files that singulairty image will write , in this case nipype(package that implements fmri pipeline)         writes to this. Nipype will try to write to host /tmp or var_tmp and will have permission issue. Hence usage of this           option.
    
    -B Binds path on host to container
    
    Variables:
    NIPYPE_TMP_DIR : Directory to store nipype tmp files, can be anywhere that the image has write access to. Could be stored     in sessions output dir in the logs as well.
    INPUT_DIR : Full path to directory where input data is located : Data can be a nifti file or dicoms
    The script assumes there is only one nifti file
    OUTPUT_DIR : Directory where outputs will be written to
    
    Required arguments:
    Each input scan is from sepecific subject and sessionthat The foll. arguemnts are used to create output directories using     these args
    SUBJECT_ID : subject id 
    SESSION_ID : session id
    json : full path to accessible json file for the singularity image
  
  Sample output data:
  ls output/fmri_outputs/
outputs_description.txt		quality_control_readme.txt	sub1

ls output/fmri_outputs/sub1/sess1/func/fmri_spm12

QC_Framewise_displacement.txt	meansub2.nii			rsub2.nii			wasub2.nii
asub2.nii			rp_sub2.txt			swasub2.nii

                   
