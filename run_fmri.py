#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This layer includes the interface adapter(IA) for parsing json args to read structural T1w scans (formats:BIDS, nifti files, dicoms)
This layer sends the output to fmri_use_cases_layer with the appropriate inputs to run the pipeine using nipype interface

Sample run examples:
python3 run_fmri.py {"options":{"value":6}, "registration_template":{"value":"/input/local0/simulatorRun/TPM.nii"}, "data":{"value":[["/input/local0/simulatorRun/3D_T1"]]}}
3D_T1 contains T1w dicoms

python3 run_fmri.py python3 run_fmri.py {"options":{"value":6}, "registration_template":{"value":"/input/local0/simulatorRun/TPM.nii"}, "data":{"value":[["/input/local0/simulatorRun/sub1_t1w.nii","/input/local0/simulatorRun/sub1_t1w.nii.gz"]]}}

python3 run_fmri.py {"options":{"value":6}, "registration_template":{"value":"/input/local0/simulatorRun/TPM.nii"}, "data":{"value":[["/input/local0/simulatorRun/BIDS_DIR"]]}}
BIDS_DIR contains bids data

success=True means program finished execution , despite the success or failure of the code
This is to indicate to coinstac that program finished execution
"""
import contextlib


@contextlib.contextmanager
def stdchannel_redirected(stdchannel, dest_filename):
    """
    A context manager to temporarily redirect stdout or stderr

    e.g.:


    with stdchannel_redirected(sys.stderr, os.devnull):
        if compiler.has_function('clock_gettime', libraries=['rt']):
            libraries.append('rt')
    """

    try:
        oldstdchannel = os.dup(stdchannel.fileno())
        dest_file = open(dest_filename, 'w')
        os.dup2(dest_file.fileno(), stdchannel.fileno())

        yield
    finally:
        if oldstdchannel is not None:
            os.dup2(oldstdchannel, stdchannel.fileno())
        if dest_file is not None:
            dest_file.close()


import json, argparse, getopt, re,traceback
import warnings, os, glob, sys
import nibabel as nib
from  math import radians as rad

with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
# Load Nipype spm interface #
from nipype.interfaces import spm
import fmri_use_cases_layer

#Stop printing nipype.workflow info to stdout
from nipype import logging
logging.getLogger('nipype.workflow').setLevel('CRITICAL')

#Create a dictionary to store all paths to softwares,templates & store parameters, names of output files

template_dict = {
    'spm_version':
    '12.7507',
    'matlab_cmd':
    '/opt/spm12/run_spm12.sh /opt/mcr/v95 script',
    'spm_path':
    '/opt/spm12/fsroot',
    'tpm_path':
    '/opt/spm12/fsroot/spm/spm12/toolbox/OldNorm/EPI.nii',
    'transf_mat_path':
    os.path.join('/computation', 'transform.mat'),
    'scan_type':
    'T1w',
    'reorient_params_x_mm': 0,
    'reorient_params_y_mm': 0,
    'reorient_params_z_mm': 0,
    'reorient_params_pitch': 0,
    'reorient_params_roll': 0,
    'reorient_params_yaw': 0,
    'reorient_params_x_scaling': 1,
    'reorient_params_y_scaling': 1,
    'reorient_params_z_scaling': 1,
    'reorient_params_x_affine': 0,
    'reorient_params_y_affine': 0,
    'reorient_params_z_affine': 0,
    'realign_fwhm':8,
    'realign_interp':2,
    'realign_quality':1,
    'realign_register_to_mean':True,
    'realign_separation':4,
    'realign_wrap': [0,0,0],
    'realign_write_interp': 4,
    'realign_write_mask': True,
    'realign_write_which': [2, 1],
    'realign_write_wrap': [0,0,0],
    'slicetime_ref_slice':None,
    'num_slices':None,
    'repetition_time':None,
    'acquisition_order':None,
    'normalize_affine_regularization_type':'mni',
    'normalize_write_bounding_box': [[-78, -112, -70],[78, 76, 85]],
    'normalize_write_interp': 1,
    'normalize_write_voxel_sizes': [3,3,3],
    'FWHM_SMOOTH': [10, 10, 10],
    "smoothing_implicit_masking":False,
    'fmri_output_dirname':
        'fmri_spm12',
    'output_zip_dir':
        'fmri_outputs',
    'display_image_name':
        'wa.png',
    'display_pngimage_name':
        'Normalized Slicetime corrected Image',
    'cut_coords': (0, 0, 0),
    'display_nifti':
        'w*.nii',
    'qc_nifti':
        'wa*nii',
    'fmri_qc_filename':
        'QC_Framewise_displacement.txt',
    'outputs_manual_name':
        'outputs_description.txt',
    'coinstac_display_info':
        'Please read outputs_description.txt for description of pre-processed output files and quality_control_readme.txt for quality control measurement.'
        'These files are placed under the pre-processed data.',
    'bids_outputs_manual_name':
        'outputs_description.txt',
    'nifti_outputs_manual_name':
        'outputs_description.txt',
    'bids_outputs_manual_content':
        "Prefixes descriptions for pre-processed images:"
        "\na-Slicetime corrected\nw-Normalized\ns-Smoothed with fwhm(mm) [6 6 6]\nFor more info. please refer to spm12 manual here: "
        "http://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf and release notes here: http://www.fil.ion.ucl.ac.uk/spm/software/spm12/SPM12_Release_Notes.pdf",
    'nifti_outputs_manual_content':
        "sub-1,sub-2,sub-* denotes each nifti file with respect to the order in the nifti paths given"
        "Prefixes descriptions for segmented images:"
        "\na-Slicetime corrected\nw-Normalized\ns-Smoothed with fwhm(mm) [6 6 6]\nFor more info. please refer to spm12 manual here: "
        "http://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf and release notes here: http://www.fil.ion.ucl.ac.uk/spm/software/spm12/SPM12_Release_Notes.pdf",
    'dicoms_outputs_manual_content':
        "sub-1,sub-2,sub-* denotes each nifti file with respect to the order in the nifti paths given"
        "Prefixes descriptions for segmented images:"
        "\na-Slicetime corrected\nw-Normalized\ns-Smoothed with fwhm(mm) [6 6 6]\nFor more info. please refer to spm12 manual here: "
        "http://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf and release notes here: http://www.fil.ion.ucl.ac.uk/spm/software/spm12/SPM12_Release_Notes.pdf",
    'qc_readme_name':
        'quality_control_readme.txt',
    'qc_readme_content':
        "In each subject's func/fmri_spm12 directory,QC_Framewise_displacement.txt gives the mean of RMS of framewise displacement "
        "\nFramewise Displacement of a time series is defined as the sum of the absolute values of the derivatives of the six realignment parameters "
        "\nRotational displacements are converted from degrees to millimeters by calculating displacement on the surface of a sphere of radius 50 mm "
        "\nFD = 0.15 to 0.2 mm: significant changes begin to be seen"
        "\nFD > 0.5 mm: marked correlation changes observed"
}
"""
More info. on keys in template_dict

spm_path is path to spm software inside docker . 
SPM is Statistical Parametric Mapping toolbox for matlab 
Info. from http://www.fil.ion.ucl.ac.uk/spm/
"Statistical Parametric Mapping refers to the construction and assessment of spatially extended statistical processes used to test hypotheses about functional imaging data. 
These ideas have been instantiated in software that is called SPM.
The SPM software package has been designed for the analysis of brain imaging data sequences. 
The sequences can be a series of images from different cohorts, or time-series from the same subject. 
The current release is designed for the analysis of fMRI, PET, SPECT, EEG and MEG."

tpm_path is the path where the SPM structural template nifti file is stored
This file is used to :
1) Perform segmentation in the fmri pipeline
2) Compute correlation value to smoothed, warped grey matter from output of pipeline, which is stored in the fmri_qc_filename

transf_mat_path is the path to the transformation matrix used in running the reorient step of the pipeline
scan_type is the type of structural scans on which is accepted by this pipeline
FWHM_SMOOTH is the full width half maximum smoothing kernel value in mm in x,y,z directions
fmri_output_dirname is the name of the output directory to which the outputs from this pipeline are written to
fmri_qc_filename is the name of the fmri quality control text file , which is placed in fmri_output_dirname
FWHM_SMOOTH is an optional parameter that can be passed as json in args['input']['opts']

json output description
                    "message"-This string is used by coinstac to display output message to the user on the UI after computation is finished
                    "download_outputs"-Zipped directory where outputs are stored
                    "display"-base64 encoded string of slicetime corrected normalized output nifti
"""
with warnings.catch_warnings():
    warnings.filterwarnings("ignore")


def software_check():
    """This function returns the spm standalone version installed inside the docker
    """
    spm.SPMCommand.set_mlab_paths(
        matlab_cmd=template_dict['matlab_cmd'], use_mcr=True)
    return (spm.SPMCommand().version)


def args_parser():
    parser = argparse.ArgumentParser(
        description='Run vbm pipeline on mprage T1w dicoms')
    required = parser.add_argument_group('required arguments')
    required.add_argument(
        '-sub', required=True, nargs=1, type=str, help='Enter subject_id')
    required.add_argument(
        '-sess', required=True, nargs=1, type=str, help='Enter session')
    required.add_argument(
        '-json',
        '--json',
        nargs=1,
        type=str,
        required=True,
        help=
        'Enter the path to json file. Ex: -json "/data/input.json"'
    )

    return parser.parse_args()

def data_parser(args):


    if 'reorient_params_x_mm' in args['options']:
        template_dict['reorient_params_x_mm'] = float(args['options']['reorient_params_x_mm'])
    if 'reorient_params_y_mm' in args['options']:
        template_dict['reorient_params_y_mm'] = float(args['options']['reorient_params_y_mm'])
    if 'reorient_params_z_mm' in args['options']:
        template_dict['reorient_params_z_mm'] = float(args['options']['reorient_params_z_mm'])
    if 'reorient_params_pitch' in args['options']:
        template_dict['reorient_params_pitch'] = float((args['options']['reorient_params_pitch']))
    if 'reorient_params_roll' in args['options']:
        template_dict['reorient_params_roll'] = float((args['options']['reorient_params_roll']))
    if 'reorient_params_yaw' in args['options']:
        template_dict['reorient_params_yaw'] = float((args['options']['reorient_params_yaw']))
    if 'reorient_params_x_scaling' in args['options']:
        template_dict['reorient_params_x_scaling'] = float(args['options']['reorient_params_x_scaling'])
    if 'reorient_params_y_scaling' in args['options']:
        template_dict['reorient_params_y_scaling'] = float(args['options']['reorient_params_y_scaling'])
    if 'reorient_params_z_scaling' in args['options']:
        template_dict['reorient_params_z_scaling'] = float(args['options']['reorient_params_z_scaling'])
    if 'reorient_params_x_affine' in args['options']:
        template_dict['reorient_params_x_affine'] = float(args['options']['reorient_params_x_affine'])
    if 'reorient_params_y_affine' in args['options']:
        template_dict['reorient_params_y_affine'] = float(args['options']['reorient_params_y_affine'])
    if 'reorient_params_z_affine' in args['options']:
        template_dict['reorient_params_z_affine'] = float(args['options']['reorient_params_z_affine'])


    if 'realign_fwhm' in args['options']:
        template_dict['realign_fwhm']=args['options']['realign_fwhm']
    if 'realign_interp' in args['options']:
        template_dict['realign_interp']=args['options']['realign_interp']
    if 'realign_quality' in args['options']:
        template_dict['realign_quality']=args['options']['realign_quality']
    if 'realign_register_to_mean' in args['options']:
        template_dict['realign_register_to_mean']=args['options']['realign_register_to_mean']
    if 'realign_separation' in args['options']:
        template_dict['realign_separation']=args['options']['realign_separation']
    if 'realign_wrap' in args['options']:
        template_dict['realign_wrap']=args['options']['realign_wrap']
    if 'realign_write_interp' in args['options']:
        template_dict['realign_write_interp']=args['options']['realign_write_interp']
    if 'realign_write_mask' in args['options']:
        template_dict['realign_write_mask']=args['options']['realign_write_mask']
    if 'realign_write_which' in args['options']:
        template_dict['realign_write_which']=args['options']['realign_write_which']
    if 'realign_write_wrap' in args['options']:
        template_dict['realign_write_wrap']=args['options']['realign_write_wrap']

    if 'slicetime_ref_slice' in args['options']:
        template_dict['slicetime_ref_slice']=args['options']['slicetime_ref_slice']
    if 'num_slices' in args['options']:
        template_dict['num_slices']=args['options']['num_slices']
    if 'repetition_time' in args['options']:
        template_dict['repetition_time']=args['options']['repetition_time']
    if 'acquisition_order' in args['options']:
        template_dict['acquisition_order']=args['options']['acquisition_order']


    if 'normalize_affine_regularization_type' in args['options']:
        template_dict['normalize_affine_regularization_type']=args['options']['normalize_affine_regularization_type']
    if 'normalize_write_bounding_box' in args['options']:
        template_dict['normalize_write_bounding_box']=args['options']['normalize_write_bounding_box']
    if 'normalize_write_interp' in args['options']:
        template_dict['normalize_write_interp']=int(args['options']['normalize_write_interp'])
    if 'normalize_write_voxel_sizes' in args['options']:
        template_dict['normalize_write_voxel_sizes']=args['options']['normalize_write_voxel_sizes']


    if 'smoothing_x_mm' in args['options']:
         template_dict['FWHM_SMOOTH'][0]= float(args['options']['smoothing_x_mm'])
    if 'smoothing_y_mm' in args['options']:
         template_dict['FWHM_SMOOTH'][1]= float(args['options']['smoothing_y_mm'])
    if 'smoothing_z_mm' in args['options']:
        template_dict['FWHM_SMOOTH'][2] = float(args['options']['smoothing_z_mm'])

    if 'smoothing_implicit_masking' in args['options']:
        template_dict['implicit_masking']=args['options']['smoothing_implicit_masking']


    if 'registration_template' in args['options']:
        if os.path.isfile(args['options']['registration_template']) and (str(
                ((nib.load(template_dict['tpm_path'])).shape)) == str(
            ((nib.load(args['options']['registration_template'])).shape))):
            template_dict['tpm_path'] = args['options']['registration_template']

def convert_reorientparams_save_to_mat_script():
    from pathlib2 import Path
    import shutil
    shutil.copy('/computation/convert_to_mat_file_template.m', '/computation/convert_to_mat_file.m')
    path = Path('/computation/convert_to_mat_file.m')
    text = path.read_text()
    text = text.replace('x_mm', str(template_dict['reorient_params_x_mm']))
    text = text.replace('y_mm', str(template_dict['reorient_params_y_mm']))
    text = text.replace('z_mm', str(template_dict['reorient_params_z_mm']))
    text = text.replace('pitch', str(template_dict['reorient_params_pitch']))
    text = text.replace('roll', str(template_dict['reorient_params_roll']))
    text = text.replace('yaw', str(template_dict['reorient_params_yaw']))
    text = text.replace('x_scaling', str(template_dict['reorient_params_x_scaling']))
    text = text.replace('y_scaling', str(template_dict['reorient_params_y_scaling']))
    text = text.replace('z_scaling', str(template_dict['reorient_params_z_scaling']))
    text = text.replace('x_affine', str(template_dict['reorient_params_x_affine']))
    text = text.replace('y_affine', str(template_dict['reorient_params_y_affine']))
    text = text.replace('z_affine', str(template_dict['reorient_params_z_affine']))
    path.write_text(text)
    # Run reorient.m script using spm12 standalone and Matlab MCR
    with stdchannel_redirected(sys.stderr, os.devnull):
        spm.SPMCommand.set_mlab_paths(matlab_cmd='/opt/spm12/run_spm12.sh /opt/mcr/v95 script /computation/convert_to_mat_file.m',
                                  use_mcr=True)

if __name__ == '__main__':

    try:
        # Code to check input args
        args = args_parser()
        with open(args.json[0], 'r') as f:
            json_string = json.load(f)
        data_parser(json_string)
    except Exception as e:
        sys.stderr.write('Unable to read input args/data. Error_log:'+str(e)+str(traceback.format_exc()))

    # Code to convert input reorient params options to transform.mat for reorientation, try,except,pass statement is just to supress the error its printing out.
    # Have not been able to come up with a more elegant solution yet. Already tried to supress warnings etc.
    try:
        convert_reorientparams_save_to_mat_script()
    except:
        pass


    try:
        # Check if spm is running
        with stdchannel_redirected(sys.stderr, os.devnull):
            spm_check = software_check()
        if spm_check != template_dict['spm_version']:
            raise EnvironmentError("spm unable to start in vbm docker")

        read_data = '/data'
        output_dir = '/out'

        if args.sub[0] and args.sess[0]:
            template_dict['subject'] = args.sub[0]
            template_dict['session'] = args.sess[0]
            # Check if data has nifti files
            if os.path.isfile(glob.glob(os.path.join(
                    read_data, '*.nii*'))[0]) and os.access(output_dir, os.W_OK):
                nifti_file = glob.glob(os.path.join(read_data, '*.nii*'))[0]
                computation_output = fmri_use_cases_layer.setup_pipeline(
                    data=nifti_file,
                    write_dir=output_dir,
                    data_type='nifti',
                    **template_dict)
                sys.stdout.write(computation_output)
            # Check if inputs are dicoms
            elif os.path.isdir(read_data) and os.listdir(read_data) and os.access(
                    output_dir, os.W_OK):
                dicom_header_info = os.popen('strings' + ' ' + glob.glob(
                    os.path.join(read_data, '*'))[0] + '|grep DICM').read()
                if 'DICM' in dicom_header_info:
                    computation_output = fmri_use_cases_layer.setup_pipeline(
                        data=read_data,
                        write_dir=output_dir,
                        data_type='dicoms',
                        **template_dict)
                    sys.stdout.write(computation_output)
            else:
                sys.stdout.write("No data found")
    except Exception as e:
        sys.stderr.write('Unable to pre-process data. Error_log:'+str(e)+str(traceback.format_exc()))
