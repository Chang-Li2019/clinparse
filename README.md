# clinparse
A pipeline to parse ClinVar and perform descriptive analyses on generated PDB files. Specifically, first ClinVar missense pathogenic/benign variants were identified. Second, users can run AlphaFold to generate PDB files for mutant and wild-type protein sequence (gene) of interest. Third, RMSD will be calculated and plDDT scores/contact maps will be visualized for each mutant sequence.

1. Create a conda environment

```
conda create -n clinparse python==3.9
```

2. Activate envrionment and install pip

```
conda activate clinparse
conda install pip
```

3. move to the folder with the setup.py file and install setup

```
cd clinparse
pip3 install .
```

4. Preprocess ClinVar and prepare mutant protein sequences for gene name of interest:

```
python clinparse.py APOE -f data/Clinvar_20220517.vcf.gz 
```

5. (Optional) To run AlphaFold locally, follow the prompt to install localcolabfold: https://github.com/YoshitakaMo/localcolabfold;
To run AlphaFold using Google Colab, please refer to ColabFold https://github.com/sokrypton/ColabFold, or https://github.com/deepmind/alphafold

  *After installation, to run the prediction using cpu:
  
  ```
  colabfold_batch --amber --templates --num-recycle 3 --use-gpu-relax inputfile outputdir/ --cpu
  ```
  
6. Calculate RMSD and visualization using sample files:

```
clinparse % python visualize.py -f ./data/test/APOE_ref_relaxed_rank_1_model_3.pdb -tp ./data/test/TP -tn ./data/test/TN
```

