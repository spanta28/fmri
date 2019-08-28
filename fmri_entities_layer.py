#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This layer defines the nodes of fmri pre-processing pipeline
"""

import nipype.pipeline.engine as pe
import nipype.interfaces.spm as spm
import nipype.interfaces.spm.utils as spmu
spm.terminal_output = 'file'
from nipype.interfaces.io import DataSink

#Stop printing nipype.workflow info to stdout
from nipype import logging
logging.getLogger('nipype.workflow').setLevel('CRITICAL')

## 1 Reorientation node & settings ##
class Reorient:
    def __init__(self, nifti_file, **template_dict):
        self.node = pe.Node(interface=spmu.ApplyTransform(), name='reorient')
        self.node.inputs.mat = template_dict['transf_mat_path']
        self.node.inputs.paths = template_dict['spm_path']


## 2 Realign node & settings ##
class Realign:
    def __init__(self, **template_dict):
        self.node = pe.Node(interface=spm.Realign(), name='realign')
        self.node.inputs.paths = template_dict['spm_path']
        self.node.inputs.fwhm = template_dict['realign_fwhm']
        self.node.inputs.interp = template_dict['realign_interp']
        self.node.inputs.quality = template_dict['realign_quality']
        self.node.inputs.register_to_mean = template_dict['realign_register_to_mean']
        self.node.inputs.separation = template_dict['realign_separation']
        self.node.inputs.wrap = template_dict['realign_wrap']
        self.node.inputs.write_interp = template_dict['realign_write_interp']
        self.node.inputs.write_mask = template_dict['realign_write_mask']
        self.node.inputs.write_which = template_dict['realign_write_which']
        self.node.inputs.write_wrap = template_dict['realign_write_wrap']

## 3 Slicetiming Node and settings ##
class Slicetiming:
    def __init__(self, **template_dict):
        self.node = pe.Node(interface=spm.SliceTiming(), name='slicetiming')
        self.node.inputs.paths = template_dict['spm_path']

## 4 Normalize Node and settings ##
class Normalize:
    def __init__(self, **template_dict):
        self.node = pe.Node(interface=spm.Normalize12(), name='normalize')
        self.node.inputs.tpm = template_dict['tpm_path']
        self.node.inputs.affine_regularization_type = template_dict['normalize_affine_regularization_type']
        self.node.inputs.write_bounding_box = template_dict['normalize_write_bounding_box']
        self.node.inputs.write_interp = template_dict['normalize_write_interp']
        self.node.inputs.write_voxel_sizes = template_dict['normalize_write_voxel_sizes']


## 5 Smoothing Node & Settings ##
class Smooth:
    def __init__(self, **template_dict):
        self.node = pe.Node(interface=spm.Smooth(), name='smoothing')
        self.node.inputs.paths = template_dict['spm_path']
        self.node.inputs.fwhm = template_dict['FWHM_SMOOTH']
        self.node.inputs.implicit_masking=template_dict['smoothing_implicit_masking']

## 5 Datsink Node that collects segmented, smoothed files and writes to temp_write_dir ##
class Datasink:
    def __init__(self):
        self.node = pe.Node(interface=DataSink(), name='sinker')