# clinparse
A pipeline to parse ClinVar and perform descriptive analyses on generated PDB files. Specifically, first ClinVar missense pathogenic/benign variants were identified. Second, users can run AlphaFold to generate PDB files for mutant and wild-type protein sequence (gene) of interest. Third, RMSD will be calculated and plDDT scores/contact maps will be visualized for each mutant sequence.

1. Download package

2. Install requirement packages in a virtual environment

```
cd clinparse
conda env create -f environment.yml
```

3. Activate environment and move to clinparse folder

```
conda activate clinparse
cd clinparse
```

4. Install reference databases for annotation:

```
pyensembl install --release 106 --species human
```

5. Preprocess ClinVar and prepare mutant protein sequences for gene name of interest:

```
python clinparse.py APOE -f data/Clinvar_20220517.vcf.gz 
```

6. (Optional) To run AlphaFold locally, follow the prompt to install localcolabfold: https://github.com/YoshitakaMo/localcolabfold;
To run AlphaFold using Google Colab, please refer to ColabFold https://github.com/sokrypton/ColabFold, or https://github.com/deepmind/alphafold

  *After installation, to run the prediction using cpu:
  
  ```
  colabfold_batch --amber --templates --num-recycle 3 --use-gpu-relax inputfile outputdir/ --cpu
  ```
  
7. Calculate RMSD and visualization using sample files:

```
python visualize.py -f ./data/test/APOE_ref_relaxed_rank_1_model_3.pdb -tp ./data/test/TP -tn ./data/test/TN
```

