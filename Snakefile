#
# Preprocessing of Human Connectome Project data for use in whole-brain modeling
# Viktor Sip, 2021
#


import os
import glob

envvars: "DICERPATH"
DICERPATH = os.environ["DICERPATH"]

SUBJECTS = [os.path.split(d)[1].split('_')[0] for d in glob.glob("data/*_3T_Structural_preproc.zip")]
CONFS = ["all", "ctx"]
NTRACKS = "10M"

scriptdir = os.path.join(workflow.basedir, 'scripts')
etcdir = os.path.join(workflow.basedir, 'etc')


localrules: unzip

# Utils -----------------------------------------------------------------------------------------------------------

rule unzip:
    input: "data/{s}_{package}.zip"
    output: directory("run/{s}/data/{package}")
    shell: """
        unzip -qq {input} -d run/{wildcards.s}/data/{wildcards.package}     && \
        shopt -s dotglob nullglob && mv {output}/{wildcards.s}/* {output}/  && \
        rmdir {output}/{wildcards.s}
    """

# Diffusion processing --------------------------------------------------------------------------------------------

rule dwi_convert:
    input: "run/{s}/data/3T_Diffusion_preproc/"
    output: temp("run/{s}/diff/dwi.mif")
    threads: 4
    shell:
        "mrconvert {input}/T1w/Diffusion/data.nii.gz -fslgrad {input}/T1w/Diffusion/bvecs {input}/T1w/Diffusion/bvals"
        "    -datatype float32 -strides 0,0,0,1 {output[0]} -nthreads {threads}"


rule dwi2response:
    input: "run/{s}/diff/dwi.mif"
    output:
        wm="run/{s}/diff/response_wm.txt",
        gm="run/{s}/diff/response_gm.txt",
        csf="run/{s}/diff/response_csf.txt"
    threads: 4
    shell: "dwi2response dhollander {input} {output.wm} {output.gm} {output.csf} -nthreads {threads}"


rule avg_response:
    input: expand("run/{s}/diff/response_{{tissue}}.txt", s=SUBJECTS)
    output: "run/group/response_{tissue}.txt"
    shell: "responsemean {input} {output}"


rule dwi2fod:
    input:
        dwi="run/{s}/diff/dwi.mif",
        resp_wm="run/group/response_wm.txt",
        resp_gm="run/group/response_gm.txt",
        resp_csf="run/group/response_csf.txt",
        datadir="run/{s}/data/3T_Diffusion_preproc"
    output:
        fod_wm =temp("run/{s}/diff/fod_wm.mif"),
        fod_gm =temp("run/{s}/diff/fod_gm.mif"),
        fod_csf=temp("run/{s}/diff/fod_csf.mif")
    threads: 4
    shell:
        "dwi2fod msmt_csd {input.dwi}"
        "    {input.resp_wm} {output.fod_wm} {input.resp_gm} {output.fod_gm} {input.resp_csf} {output.fod_csf}"
        "    -mask {input.datadir}/T1w/Diffusion/nodif_brain_mask.nii.gz -nthreads {threads}"


rule tckgen:
    input:
        fod_wm="run/{s}/diff/fod_wm.mif",
        datadir="run/{s}/data/3T_Diffusion_preproc"
    output:
        tracks=temp("run/{s}/diff/tracks_all.tck")
    threads: 4
    shell:
        "tckgen {input.fod_wm} {output.tracks} -algorithm iFOD2"
        "    -mask       {input.datadir}/T1w/Diffusion/nodif_brain_mask.nii.gz"
        "    -seed_image {input.datadir}/T1w/Diffusion/nodif_brain_mask.nii.gz"
        "    -select {NTRACKS} -nthreads {threads}"


rule sift:
    input:
        tracks="run/{s}/diff/tracks_all.tck",
        fod_wm="run/{s}/diff/fod_wm.mif"
    output: "run/{s}/diff/tracks_sifted.tck"
    threads: 4
    shell: "tcksift {input.tracks} {input.fod_wm} {output} -nthreads {threads}"


rule labels:
    input:
        datadir="run/{s}/data/3T_Structural_preproc",
        inlut=os.path.join(etcdir, "FreeSurferColorLUT.txt"),
        outlut=os.path.join(etcdir, "lut.dk.{conf}.txt")
    output: "run/{s}/diff/labels.{conf}.nii.gz"
    shell: "labelconvert {input.datadir}/T1w/aparc+aseg.nii.gz {input.inlut} {input.outlut} {output}"


rule conn:
    input:
        tracks="run/{s}/diff/tracks_sifted.tck",
        labels="run/{s}/diff/labels.{conf}.nii.gz",
    output: "run/{s}/res_{conf}/weights.txt",
    shell: "tck2connectome {input.tracks} {input.labels} {output} -symmetric -zero_diagonal"


rule region_list:
    input: os.path.join(etcdir, "lut.dk.{conf}.txt")
    output: "run/{s}/res_{conf}/region_names.txt"
    shell: "cat {input} | grep -Ev '(#.*$)|(^$)' | tr -s ' '    | cut -d' ' -f3 | tail -n +2 > {output}"
    #                   | Remove comments        | Merge delims | Third column  | Skip Unknown


rule sc:
    input: expand(["run/{s}/res_{conf}/weights.txt", "run/{s}/res_{conf}/region_names.txt"], s=SUBJECTS, conf=CONFS)

# fMRI processing using DiCER ------------------------------------------------------------------------------------

rule labels_subcort:
    input:  "run/{s}/data/3T_Structural_preproc"
    output: "run/{s}/rfmri/labels_subcort.dlabel.nii"
    shell:
        "wb_command -cifti-create-label {output} -volume {input}/MNINonLinear/ROIs/Atlas_ROIs.2.nii.gz"
        "    {input}/MNINonLinear/ROIs/Atlas_ROIs.2.nii.gz"


rule labels_cort:
    input:  "run/{s}/data/3T_Structural_preproc"
    output: "run/{s}/rfmri/labels_cort.dlabel.nii"
    shell: "cp {input}/MNINonLinear/fsaverage_LR32k/{wildcards.s}.aparc.32k_fs_LR.dlabel.nii {output}"


rule get_timeseries:
    input: "run/{s}/data/3T_rfMRI_REST_fix/"
    output:
        lr1=temp("run/{s}/rfmri/rfMRI_REST1-LR_orig.dtseries.nii"),
        rl1=temp("run/{s}/rfmri/rfMRI_REST1-RL_orig.dtseries.nii"),
        lr2=temp("run/{s}/rfmri/rfMRI_REST2-LR_orig.dtseries.nii"),
        rl2=temp("run/{s}/rfmri/rfMRI_REST2-RL_orig.dtseries.nii")
    shell:
        "ln -sr {input}/MNINonLinear/Results/rfMRI_REST1_LR/rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii {output.lr1};"
        "ln -sr {input}/MNINonLinear/Results/rfMRI_REST1_RL/rfMRI_REST1_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii {output.rl1};"
        "ln -sr {input}/MNINonLinear/Results/rfMRI_REST2_LR/rfMRI_REST2_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii {output.lr2};"
        "ln -sr {input}/MNINonLinear/Results/rfMRI_REST2_RL/rfMRI_REST2_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii {output.rl2};"


rule get_movement_regressors:
    input: "data/{s}_3T_rfMRI_REST{ses}_preproc.zip"
    output: "run/{s}/rfmri/Movement_Regressors_REST{ses}_{pe}.txt"
    shell: "unzip -p {input} {wildcards.s}/MNINonLinear/Results/rfMRI_REST{ses}_{pe}/Movement_Regressors.txt > {output}"


rule concat_fmri:
    input:
        lr1="run/{s}/rfmri/rfMRI_REST1-LR_orig.dtseries.nii",
        rl1="run/{s}/rfmri/rfMRI_REST1-RL_orig.dtseries.nii",
        lr2="run/{s}/rfmri/rfMRI_REST2-LR_orig.dtseries.nii",
        rl2="run/{s}/rfmri/rfMRI_REST2-RL_orig.dtseries.nii"
    output:
        cifti=temp("run/{s}/rfmri/rfMRI_ALL.dtseries.nii"),
        nifti=temp("run/{s}/rfmri/rfMRI_ALL.nii.gz"),
    shell:
        "wb_command -cifti-merge {output.cifti} -cifti {input.lr1} -cifti {input.rl1} -cifti {input.lr2} -cifti {input.rl2};"
        "wb_command -cifti-convert -to-nifti {output.cifti} {output.nifti}"


rule remove_session_effects:
    input: "run/{s}/rfmri/rfMRI_ALL.nii.gz"
    output: temp("run/{s}/rfmri/rfMRI_ALLFILT.nii.gz")
    resources:
        mem_mb=20000
    shell: "fsl_regfilt -i {input} -d {DICERPATH}/hcp_processing/restStopPoints.tsv -f 1,2,3,4 -o {output}"


rule dicer:
    input:
        rest_nifti="run/{s}/rfmri/rfMRI_ALLFILT.nii.gz"
    output: temp(directory("run/{s}/rfmri/dicer"))
    resources:
        mem_mb=20000
    shell: """
        mkdir {output};                 \
        ln -sr {input} {output};        \
        d=$PWD; cd $DICERPATH;          \
        bash DiCER_lightweight.sh -f -i `basename {input}` -w $d/{output} -s {wildcards.s} -p 5
    """


rule get_dicer_results:
    input:  "run/{s}/rfmri/dicer"
    output: directory("run/{s}/rfmri/dicer-results")
    shell:
        "mkdir {output};"
        "cp {input}/*.{{png,html,css}} {output}/;"
        "cp {input}/{wildcards.s}_dbscan_liberal_regressors.tsv {output}/;"


rule dicer_to_cifti:
    input:
        cifti_template="run/{s}/rfmri/rfMRI_ALL.dtseries.nii",
        dicer_dir="run/{s}/rfmri/dicer/"
    output: temp("run/{s}/rfmri/rfMRI_ALL_dicer.dtseries.nii")
    shell: "wb_command -cifti-convert -from-nifti {input.dicer_dir}/func_temp_correctorder_dbscan.nii.gz {input.cifti_template} {output}"


rule parcellate_fmri:
    input:
        fmri="run/{s}/rfmri/rfMRI_{ses}_{variant}.dtseries.nii",
        labels="run/{s}/rfmri/labels_{type}.dlabel.nii"
    output: "run/{s}/rfmri/rfMRI_{ses}_{variant}.{type}.ptseries.nii"
    shell: "wb_command -cifti-parcellate {input.fmri} {input.labels} COLUMN {output}"


rule join_cort_subcort_ptseries:
    input:
        cort   ="run/{s}/rfmri/rfMRI_{ses}_{variant}.cort.ptseries.nii",
        subcort="run/{s}/rfmri/rfMRI_{ses}_{variant}.subcort.ptseries.nii",
        regions="run/{s}/res_{conf}/region_names.txt"
    output: temp("run/{s}/rfmri/rfMRI_{ses}_{variant}.{conf}.npz")
    shell: "python {scriptdir}/util.py join_cort_subcort {input.cort} {input.subcort} {input.regions} {output}"


rule global_signal_regression:
    input:  "run/{s}/rfmri/rfMRI_{ses}_orig.{conf}.npz"
    output: temp("run/{s}/rfmri/rfMRI_{ses}_gsr.{conf}.npz")
    shell: "python {scriptdir}/util.py gsr {input} {output}"


rule split_fmri:
    input: "run/{s}/rfmri/rfMRI_ALL_{variant}.{conf}.npz"
    output: [temp(f"run/{{s}}/rfmri/rfMRI_{ses}_{{variant}}.{{conf}}.npz")      \
             for ses in ['REST1-LR', 'REST1-RL', 'REST2-LR', 'REST2-RL']]
    shell: "python {scriptdir}/util.py split {input} '{output}'"


rule export:
    input:  "run/{s}/rfmri/rfMRI_{ses}_{var}.{conf}.npz"
    output: "run/{s}/res_{conf}/rfMRI_{ses}_{var}.npz"
    shell: "cp {input} {output}"


rule fmri:
    input:
        expand("run/{s}/rfmri/dicer-results/", s=SUBJECTS),
        expand("run/{s}/res_{conf}/rfMRI_{ses}_{var}.npz", s=SUBJECTS, conf=CONFS, \
               ses=['REST1-LR', 'REST1-RL', 'REST2-LR', 'REST2-RL'], var=['orig', 'gsr', 'dicer']),
        expand("run/{s}/rfmri/Movement_Regressors_REST{ses}_{pe}.txt", s=SUBJECTS, ses=['1', '2'], pe=['LR', 'RL'])
