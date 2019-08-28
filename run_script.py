#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example:
multiple sessions
python3 run_script.py -o /home/vagrant/output -tmp /home/vagrant/nipype_tmp -image_path /home/vagrant/images/s3_latest  -json '{"smooth":6, "data":[["/home/vagrant/input1","sub1","sess1"],["/home/vagrant/input2","sub2","sess2"]]}'
"""
import json
import argparse
import os
import glob
import shutil
from argparse import ArgumentTypeError as err

# Define parameters for running singluarity
VBM_SCRIPT_NAME = "run_vbm.py"
FMRI_SCRIPT_NAME = "run_fmri.py"
tmp_directory_for_nipype = ''
singularity_image_full_path = ''
input_directory = ''
slurm_cmd = 'sbatch --parsable --partition qBF'
singularity_cmd = 'singularity exec --contain --workdir'
singularity_options_bind_cmd = '-B'
python_cmd_fmt = 'python3 /computation/%s'


class PathType(object):
    def __init__(self, exists=True, type='file', dash_ok=True):
        '''exists:
                True: a path that does exist
                False: a path that does not exist, in a valid parent directory
                None: don't care
           type: file, dir, symlink, None, or a function returning True for valid paths
                None: don't care
           dash_ok: whether to allow "-" as stdin/stdout'''

        assert exists in (True, False, None)
        assert type in ('file', 'dir', 'symlink', None) or hasattr(
            type, '__call__')

        self._exists = exists
        self._type = type
        self._dash_ok = dash_ok

    def __call__(self, string):
        if string == '-':
            # the special argument "-" means sys.std{in,out}
            if self._type == 'dir':
                raise err(
                    'standard input/output (-) not allowed as directory path')
            elif self._type == 'symlink':
                raise err(
                    'standard input/output (-) not allowed as symlink path')
            elif not self._dash_ok:
                raise err('standard input/output (-) not allowed')
        else:
            for sstring in string.split(','):
                e = os.path.exists(sstring)
                if self._exists == True:
                    if not e:
                        raise err("path does not exist: '%s'" % sstring)

                    if self._type is None:
                        pass
                    elif self._type == 'file':
                        if not os.path.isfile(sstring):
                            raise err("path is not a file: '%s'" % sstring)
                    elif self._type == 'symlink':
                        if not os.path.symlink(sstring):
                            raise err("path is not a symlink: '%s'" % sstring)
                    elif self._type == 'dir':
                        if not os.path.isdir(sstring):
                            raise err(
                                "path is not a directory: '%s'" % sstring)
                    elif not self._type(string):
                        raise err("path not valid: '%s'" % sstring)
                else:
                    if self._exists == False and e:
                        raise err("path exists: '%s'" % sstring)

                    p = os.path.dirname(os.path.normpath(sstring)) or '.'
                    if not os.path.isdir(p):
                        raise err("parent path is not a directory: '%s'" % p)
                    elif not os.path.exists(p):
                        raise err("parent directory does not exist: '%s'" % p)

        return string


def args_parser():
    parser = argparse.ArgumentParser(
        description='Run vbm pipeline on mprage T1w dicoms')
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    required.add_argument(
        '-json',
        '--json',
        nargs=1,
        type=str,
        required=True,
        help=
        'Enter the path to json file. Ex: -json "/data/input.json"'
    )
    required.add_argument(
        '-o',
        '--output',
        type=PathType(exists=True, type='dir'),
        required=True,
        nargs=1,
        help='Existing output directory to store vbm outputs')
    required.add_argument(
        '-tmp',
        type=PathType(exists=True, type='dir'),
        required=True,
        nargs=1,
        help='Nipype tmp directory, it can be same as output directory')
    required.add_argument(
        '-image_path',
        required=True,
        nargs=1,
        type=str,
        help='Enter full path to singularity vbm image name')
    optional.add_argument(
        '-g',
        '--slurm',
        type=bool,
        nargs=1,
        default=False,
        help='Use slurm, True or False, default=False')
    optional.add_argument(
        '-type',
        '--data_type',
        type=str,
        nargs=1,
        choices=['nifti', 'dicoms'],
        default='nifti',
        help="input data type. Ex:'nifti','dicoms'")
    return parser.parse_args()


def main(input,
         output,
         singularity_image_full_path,
         tmp_directory_for_nipype,
         sub,
         sess,
         slurm=False,
         logger=None,
         script_name=VBM_SCRIPT_NAME):
    python_cmd = python_cmd_fmt % script_name
    if slurm:
        cmd = slurm_cmd + ' ' + singularity_cmd + ' ' + tmp_directory_for_nipype + ' ' + singularity_options_bind_cmd + ' ' + \
            input + ':/data' + ',' + \
            output + ':/out' + ' ' + singularity_image_full_path + ' ' + python_cmd + ' ' + '-sub' + ' ' + \
            sub + ' ' + '-sess' + ' ' + \
            sess + ' '+ '-json'+' '+json_file
        if logger is not None:
            logger.info('Executed command is: ' + cmd)
        job_id = os.popen(cmd).read()
        if logger is not None:
            logger.info('slurm job_id is:' + job_id)
    else:
        cmd = singularity_cmd + ' ' + tmp_directory_for_nipype + ' ' + singularity_options_bind_cmd + ' ' + \
            input + ':/data' + ',' + \
            output + ':/out' + ' ' + singularity_image_full_path + ' ' + python_cmd + ' ' + '-sub' + ' ' + \
            sub + ' ' + '-sess' + ' ' + \
            sess + ' '+ '-json'+' '+json_file
        if logger is not None:
            logger.info('Executed command is: ' + cmd)
        os.system(cmd)
    shutil.rmtree(os.path.join(tmp_directory_for_nipype, 'tmp'))
    shutil.rmtree(os.path.join(tmp_directory_for_nipype, 'var_tmp'))
    for nii_file in glob.glob(
            os.path.join(tmp_directory_for_nipype, sub, sess, 'anat',
                         '*.nii')):
        if os.path.basename(nii_file) == os.path.basename(input):
            os.remove(nii_file)


if __name__ == '__main__':
    args = args_parser()
    print('args passed are:')
    print(args)
    python_cmd = python_cmd_fmt % VBM_SCRIPT_NAME
    singularity_image_full_path = args.image_path[0]
    tmp_directory_for_nipype = args.tmp[0]

    with open(args.json[0], 'r') as f:
        json_string = json.load(f)
    json_file = args.json[0]

    data = json_string['data']

    for each_subject in data:
        input = each_subject[0]
        sub = each_subject[1]
        sess = each_subject[2]
        if args.slurm == True:
            cmd = slurm_cmd + ' ' + singularity_cmd + ' ' + tmp_directory_for_nipype + ' ' + singularity_options_bind_cmd + ' ' + \
                  input + ':/data' + ',' + \
                  args.output[0] + ':/out' + ' ' + singularity_image_full_path + ' ' + python_cmd + ' ' + '-sub' + ' ' + \
                  sub + ' ' + '-sess' + ' ' + \
                  sess + ' '+ '-json'+' '+json_file
            print('Executed command is: ' + cmd)
            job_id = os.popen(cmd).read()
            print('slurm job_id is:' + job_id)

        else:
            cmd = singularity_cmd + ' ' + tmp_directory_for_nipype + ' ' + singularity_options_bind_cmd + ' ' + \
                  input + ':/data' + ',' + \
                  args.output[0] + ':/out' + ' ' + singularity_image_full_path + ' ' + python_cmd + ' ' + '-sub' + ' ' + \
                  sub + ' ' + '-sess' + ' ' + \
                  sess + ' '+ '-json'+' '+json_file
            print('Executed command is: ' + cmd)
            os.system(cmd)
