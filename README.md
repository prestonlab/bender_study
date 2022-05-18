# bender_study
Code for analysis of the Bender study phase.

## Analysis procedure

### Installation

Run on Lonestar 6 with Python 3.9.

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

* `. ~/analysis/bender_study/.profile`
  * Set environment variables needed to run scripts

### Preprocessing

Basic preprocessing was done using fPrep 1.0.0 (originally called FAT). See the fPrep project for details. After basic preprocessing, functional data were motion-corrected, unwarped, and registered into a common space. Also had anatomical ROIs derived from FreeSurfer 5.3.0.

### Smoothing and filtering

* `rlaunch -J "smooth_susan bender_smooth_susan.sh -f 32 -v /work/03206/mortonne/lonestar/bender/{s}/BOLD/antsreg/data/{r} /work/03206/mortonne/lonestar/bender/{s}/BOLD/{r}/fm/brainmask 4.0 /work/03206/mortonne/lonestar/bender/{s}/BOLD/antsreg/data/{r}_hpfsm" $SUBJIDS $PREXRUNS`
* `rlaunch -J "smooth_susan bender_smooth_susan.sh -f 32 -v /work/03206/mortonne/lonestar/bender/{s}/BOLD/antsreg/data/{r} /work/03206/mortonne/lonestar/bender/{s}/BOLD/{r}/fm/brainmask 4.0 /work/03206/mortonne/lonestar/bender/{s}/BOLD/antsreg/data/{r}_hpfsm" $SUBJIDS $STUDYRUNS`

### Betaseries estimation

* Set up model FSF files using FEAT. This would have involved first creating a template FSF file for one participant, and then modifying that template using string substitution in sed.
* `slaunch -J estimate_betaseries "bender_estimate_betaseries.py {} study_stim2 30" $SUBJIDS`
* `slaunch -J estimate_betaseries "bender_estimate_betaseries.py {} prex_stim2 40" $SUBJIDS`

### Searchlight statmaps
* `slaunch -J sl_rdm_obj "bender_sl_rdm_obj.py {} brainmask a-bcxy study_wiki_w2v_fix_cont_a_bc_sme -m sme -s _stim2 -p 100 -n 48" $SUBJIDS`
  * Run a searchlight over study-phase item betaseries images, to find where the A model correlation minus the BC model correlation is greater for correct compared to incorrect trials. 
* `slaunch -J sl_rdm_obj "bender_sl_rdm_obj.py {} brainmask ac-bx study_wiki_w2v_fix_cont_ac_bx_sme -m sme -s _stim2 -p 100 -n 48" $SUBJIDS`
  * Run a searchlight over study-phase item betaseries images, to find where the AC model correlation minus the B model correlation is greater for correct compared to incorrect trials.

### Searchlight statistical analysis
