# Transcranial Ultrasound Stimulation (TUS) Whole-Brain Modeling

Python code to reproduce the computational model of **transcranial ultrasound stimulation (TUS)** at different stimulation intensities, as proposed in the preprint:

**[Gatica et al., 2025 - "Understanding the high-order network plasticity mechanisms of ultrasound neuromodulation"](https://www.biorxiv.org/content/10.1101/2025.01.11.632528v1)**

---

## 🧠 Overview

This repository contains the simulation pipeline to model TUS-driven changes in brain dynamics using a whole-brain Hopf model, incorporating:

- Homogeneous simulations for control conditions  
- Fitting procedures for global coupling  
- Heterogeneous modeling based on communicability and distance metrics  
- Intensity-dependent stimulation of brain regions  

---

## 📂 Pipeline: How to Run
<p>
Scripts are available in the 
  
<a href="https://github.com/KGatica/Modelling-TUS/tree/f6431cadf0c6e90cf6f1670c8db36f4cbc77828f/scripts" target="_blank">
<code>scripts</code></a> folder.
</p>
Follow the steps below to reproduce the model results.
<ol>
  <li>
    <code>01_main_runHopf_CN_homogeneous</code><br>
    Running the simulations for control conditions and several global coupling parameters (G).
  </li>
  <li>
    <code>02_fitting_CN</code><br>
    Fitting the global coupling parameter (G = 0.16).
  </li>
  <li>
    <code>03main_runHopf_CN_heterogeneous</code><br>
    Adding heterogeneity based on communicability (CMY) or distance vector, testing several bias and scale parameters.
  </li>
  <li>
    <code>04_fitting_stim_het</code><br>
    Fitting the bias and scale parameters for the heterogeneous model.
  </li>
  <li>
    <code>05_main_runHopf_sameparams_alpha</code><br>
    Running the dynamic stimulation model based on communicability or distance, at several stimulation intensities (modulated by alpha).<br>
    The stimulated targets are the thalamus (<code>"thalamusproper"</code>) and the inferior frontal cortex (<code>"parstriangularis"</code>).
  </li>
</ol>

<h2>📄 Original Hopf Model</h2>
<p>
This model is based on and adapted from the Hopf whole-brain simulation framework available at:<br>
<a href="https://github.com/carlosmig/StarCraft-2-Modeling" target="_blank">https://github.com/carlosmig/StarCraft-2-Modeling</a><br>

Coronel-Oliveros C, Medel V, Orellana S, et al. (2024). <em>Gaming expertise induces meso-scale brain plasticity and efficiency mechanisms as revealed by whole-brain modeling.</em> <strong>NeuroImage</strong>, 293:120633.<br>
<a href="https://doi.org/10.1016/j.neuroimage.2024.120633" target="_blank">https://doi.org/10.1016/j.neuroimage.2024.120633</a>
</p>

## 📊 Datasets
<p>
Data is available in the 
<a href="https://github.com/KGatica/Modelling-TUS/tree/f6431cadf0c6e90cf6f1670c8db36f4cbc77828f/datasets" target="_blank">
<code>datasets</code></a> folder.
</p>
The model uses 84 brain regions, parcellated using the Desikan-Killiany atlas, and functional and structural data from 19 subjects.  

The following datasets are provided:

- **Empirical functional data**:
  - Functional connectivity matrices (upper triangular part), shape: `3486 x 19`  
  - Kuramoto Order Parameter (KOP) values, shape: `1 x 19`  
  - **Filenames**: `triu_emp_subset.npy`, `kop_emp_subset.npy`

- **Structural connectomes**:
  - Structural connectivity matrices per subject, shape: `84 x 84 x 19`  
  - **Filename**: `sc_subset.npy`

- **Heterogeneity vectors**:
  - Global communicability: `global_CMY.mat`, shape: `1 x 84`  
  - Euclidean distance vector: `global_distance.mat`, shape: `1 x 84`

- **Stimulated targets (binary mask)**:
  - Binary vectors (1 = stimulated region, 0 = others), shape: `1 x 84`  
  - **Filenames**:  
    - `discrete_local_thalamusproper.mat`  
    - `discrete_local_parstriangularis.mat`
   

<h2>📬 Contact</h2>
<p>
For questions or collaboration inquiries, contact:<br><br>
<strong>Marilyn Gatica</strong><br>
Email: <a href="mailto:marilyn.gatica@nulondon.ac.uk">marilyn.gatica@nulondon.ac.uk</a>
</p>

<h2>📜 License</h2>
<p>
Distributed for academic, non-commercial use. Please cite the relevant papers when using this code.
</p>

<h2>🧾 Citation</h2>
<p>
If you use this code, please cite both of the following works:
</p>

<p>
Gatica M, Atkinson-Clement C, Coronel-Oliveros C, et al. (2025). <em>Differential impact of transcranial ultrasound stimulation intensities on whole-brain dynamics.</em> <strong>bioRxiv</strong>.<br>
<a href="https://www.biorxiv.org/content/10.1101/2025.01.11.632528v1" target="_blank">https://doi.org/10.1101/2025.01.11.632528v1</a>
</p>

<p>
Coronel-Oliveros C, Medel V, Orellana S, et al. (2024). <em>Gaming expertise induces meso-scale brain plasticity and efficiency mechanisms as revealed by whole-brain modeling.</em> <strong>NeuroImage</strong>, 293:120633.<br>
<a href="https://doi.org/10.1016/j.neuroimage.2024.120633" target="_blank">https://doi.org/10.1016/j.neuroimage.2024.120633</a>
</p>

