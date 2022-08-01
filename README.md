# bender_study
Code for analysis of the Bender study phase.

## Analysis procedure

### Installation

Run on Lonestar 6 with Python 3.9.

Steps to run installation (only needs to be run once):

```bash
python3 -m venv ~/software/venv/bender_study
. ~/software/venv/bender_study/bin/activate
pip install -U pip
pip install numpy scipy nibabel scikit-learn
mkdir -p ~/analysis && cd ~/analysis
cd PyMVPA || git clone https://github.com/mortonne/PyMVPA.git && cd PyMVPA
module load swig
python setup.py build_ext
python setup.py install
pip install git+https://github.com/mortonne/pprocess.git
pip install ezlaunch mindstorm
cd bender_study || git clone git@github.com:prestonlab/bender_study.git && cd bender_study
pip install -e .
```

Set up environment to run scripts (run on each login). The `.profile` will need to be edited for your specific setup (e.g., indicating where the data are located and where figures should be saved).

```bash
. ~/analysis/bender_study/.profile
```

### Preprocessing

Basic preprocessing was done using fPrep 1.0.0 (originally called FAT). See the fPrep project for details. After basic preprocessing, functional data were motion-corrected, unwarped, and registered into a common space. Also had anatomical ROIs derived from FreeSurfer 5.3.0.

* Convert DICOM files to NIfTI
  * `bender_init.sh bender_02`
* Reorganize files and merge across days
  * `slaunch -J clean_subj "bender_clean_subj.py {}" $SUBJIDS`
* Preprocess BOLD scans
  * `rlaunch -J prep_bold "prep_bold_run.sh $STUDYDIR/{s}/BOLD/{r}" $SUBJIDS $ALLRUNS` 
* Run cortical reconstruction
  * `slaunch -J run_freesurfer "run_freesurfer.sh {} 32" $SUBJIDS` 
* Convert FreeSurfer output
  * `slaunch -J convert_freesurfer "convert_freesurfer.py {}" $SUBJIDS`
* Register FreeSurfer output to main images
  * `slaunch -J reg_freesurfer "reg_freesurfer.py {}" $SUBJIDS`
* Prepare fieldmap images for unwarping
  * `slaunch -J prep_fieldmap "bender_prep_fieldmap.py {}" $SUBJIDS`
  * May need to adjust the `--dte` option for subsets of participants, as scanning parameters changed about halfway through scanning
* Calculate unwarping of BOLD scans
  * `rlaunch -J epi_reg "bender_epi_reg_run.py {s} {r}" $SUBJIDS $ALLRUNS` 
* Unwarp and co-register BOLD scans
  * `rlaunch -J reg_unwarp "reg_unwarp_bold_run.py {s} {r} study_1" $SUBJIDS $ALLRUNS` 
  
### Smoothing and filtering

* Apply spatial smoothing and temporal filtering to pre-exposure and study data
  * `rlaunch -J smooth_susan "bender_smooth_susan.sh -f 32 -v /work/03206/mortonne/lonestar/bender/{s}/BOLD/antsreg/data/{r} /work/03206/mortonne/lonestar/bender/{s}/BOLD/{r}/fm/brainmask 4.0 /work/03206/mortonne/lonestar/bender/{s}/BOLD/antsreg/data/{r}_hpfsm" $SUBJIDS $PREXRUNS`
  * `rlaunch -J smooth_susan "bender_smooth_susan.sh -f 32 -v /work/03206/mortonne/lonestar/bender/{s}/BOLD/antsreg/data/{r} /work/03206/mortonne/lonestar/bender/{s}/BOLD/{r}/fm/brainmask 4.0 /work/03206/mortonne/lonestar/bender/{s}/BOLD/antsreg/data/{r}_hpfsm" $SUBJIDS $STUDYRUNS`

### Betaseries estimation

* Set up model FSF files using FEAT. This would have involved first creating a template FSF file for one participant, and then modifying that template using string substitution in sed.
* Estimate betaseries images with activation of each item in each run, for both pre-exposure and study data.
  * `slaunch -J estimate_betaseries "bender_estimate_betaseries.py {} prex_stim2 40" $SUBJIDS`
  * `slaunch -J estimate_betaseries "bender_estimate_betaseries.py {} study_stim2 30" $SUBJIDS`

### Searchlight statmaps
* Run a searchlight looking for item reactivation during BC study. 
  * `slaunch -J sl_item_react "bender_sl_item.py {} brainmask cat_react_item2 item_react -o full -c both -d cat -s _stim2 -r 3 -p 100 -n 128" $SUBJIDS`
* Run a searchlight looking for item suppression during BC study. 
  * `slaunch -J sl_item_suppress "bender_sl_item.py {} brainmask item_suppress_gvt item_suppress -o full -c both -d cat -s _stim2 -r 3 -p 100 -n 128" $SUBJIDS`
* Run a searchlight looking for item reactivation that predicts AC accuracy.
  * `slaunch -J sl_item_react_sme "bender_sl_item.py {} brainmask cat_react_item_sme2 item_react_sme -o full -c both -d cat -s _stim2 -r 3 -p 100 -n 128" $SUBJIDS`
* Run a searchlight over study-phase item betaseries images, to find where the A model correlation minus the BC model correlation is greater for correct compared to incorrect trials. 
  * `slaunch -J sl_model "bender_sl_model.py {} brainmask a-bcxy study_wiki_w2v_fix_cont_a_bc_sme -s _stim2 -p 100 -n 128" $SUBJIDS`
* Run a searchlight over study-phase item betaseries images, to find where the AC model correlation minus the B model correlation is greater for correct compared to incorrect trials.
  * `slaunch -J sl_model "bender_sl_model.py {} brainmask ac-bx study_wiki_w2v_fix_cont_ac_bx_sme -s _stim2 -p 100 -n 128" $SUBJIDS`

### Searchlight group analysis

* Run all steps of group voxel threshold analysis.
  * `bender_run_gvt.sh -m /work/03206/mortonne/lonestar/bender/gptemplate/highres_brain_all/gp_template_mni_affine_mask.nii.gz -i BSpline -a 0.01 -p 100 mvpa/cat_react_item2 $SUBJNOS`
* Estimate smoothness for each run within a given mask.
  * `rlaunch -J res_smoothness "bender_res_smoothness.sh {s} {r} study_stim2 $WORK/bender/gptemplate/highres_brain_all/b_gray.nii.gz 128" $SUBJIDS $STUDYRUNS`
* Calculate average smoothness and run 3dClustSim to estimate a null distribution of cluster sizes. 
  * `bender_clustsim_res.sh $SUBJNOS $STUDYRUNS study_stim2 $WORK/bender/gptemplate/highres_brain_all/b_gray.nii.gz`
* Display cluster correction results for the specified ROIs. 
  * `bender_svc.sh mvpa/cat_react_item2 study_stim2 b_phc b_prc b_ifg_insula b_mpfc b_gray`

### Searchlight followup analysis

* Create individual participant masks for significant searchlight clusters.
  * `bender_sl_rois.sh`
* Calculate individual reactivation dissimilarity matrices (pre-exposure to study phase).
  * `bender_indiv_react.py cat_react_item2_lprc_dil1c -s _stim2`
* Calculate individual representational dissimilarity matrices during the study phase.
  * `bender_indiv_rdm.py cat_react_item2_lprc_dil1c -s _stim2`

### Prepare model of semantics

* Convert from MAT-file format to NPZ
  * `bender-convert-model ~/work/bender/batch/models3/allstim/mat_wiki_w2v.mat ~/work/bender/batch/semantics/wiki_w2v.npz` 

### Analysis notebooks

Notebooks are stored in the `jupyter` directory. To run, install Jupyter Lab and the notebook kernel using: 

```bash
pip install jupyterlab
python -m ipykernel install --user --name bender_study
```

You must have your virtual environment installed when you run this. 
Set up necessary environment variables by editing `.profile` for your system; you must indicate where to find data and searchlight results, and where to save figures.
Set up the Bash environment using `source .profile`.
Launch Jupyter Lab using `jupyter lab &`.

* Analyze behavioral test performance and plot.
  * `plot_behav.ipynb`
* Examine reactivation of A item perceptual templates during BC study.
  * `plot_react_stats.ipynb`
* Examine regions with reactivation that predicts subsequent AC inference.
  * `plot_sme_stats.ipynb`
* Examine the relationship between item suppression and AC semantic integration.
  * `plot_react_sem.ipynb`
