# High-Throughput Crystallographic Fragment Screening of Cyclophilin D: A Structural Atlas and Resource for Inhibitor Discovery


Analysis code for the computational parts of *High-Throughput Crystallographic Fragment Screening of Cyclophilin D: A Structural Atlas and Resource for Inhibitor Discovery* (Gebauer et al.).

There are three notebooks. They don't depend on each other and only read from `data/` and `binding_site_conservation/`, so you can run them in any order.

- `CypD Binding Family Analysis Notebook Individual mol_CLEAN.ipynb`: groups the fragment poses into binding-mode families using protein-ligand interaction fingerprints (PLIFs). Figs. 3A, 3B and S3B.
- `CypD Docking Analysis_CLEAN.ipynb`: checks how well GNINA reproduces the crystal poses. Fig. S3A.
- `CypD Binding Site Convservation_CLEAN.ipynb`: scores how conserved the binding-site residues are across the human cyclophilins.

The first two notebooks start with a longer methods write-up. What's below is the short version.

## Setup

```
conda env create -f environment.yml
conda activate cypd_analysis
```

Open the notebooks from the repo root and pick the `cypd_analysis` kernel. Every path is relative to the root, so the file loads will fail if the working directory is somewhere else. The helpers in `utils/` are imported directly. There's nothing to pip install.

Everything was run with Python 3.10.19 on macOS (Apple Silicon).

## Layout

```
data/                        structures, ligand archives, PLIF table
binding_site_conservation/   alignment and residue list for the conservation notebook
utils/                       helper functions shared by the notebooks
processed_data/              tables and per-family SDFs written by the notebooks
final_figures/               figure panels written by the notebooks
environment.yml
```

## Binding-mode families

The starting point is `data/040126_mol_prepared_plifs_w_atom.csv`. It has one row per crystal structure, with the raw PLIF from MOE (0.5 kcal/mol cutoff) and a map from each ligand atom to its residue and chain. A lot of structures have more than one fragment bound across chains A, B and X. The notebook uses that map to split each structure's PLIF into one fingerprint per pose. Poses are named `<crystal ID>-<pose index>`, e.g. `x0022-0`.

Residue numbers are shifted by +42 to match CypD numbering. Poses with fewer than two interactions are dropped.

Poses and interaction bits then go into a bipartite graph, where each pose connects to every bit it has. Bits that show up in fewer than 5 poses are removed. The families are the 7 largest communities from greedy modularity clustering in NetworkX, and anything left over is labelled `Unassigned`. Chemical similarity within each family is Tanimoto on Morgan fingerprints (radius 2, 2048 bits).

The other inputs are two CypD reference structures and two ligand archives:

- `cypd_renumbered.pdb`: used for residue names and numbering
- `cypd_reference_superpose.pdb`: the reference frame for viewing ligands
- `mol_prepared_superpose_lig_all.tar.xz`: one SDF per structure
- `mol_superpose_lig_sep.tar.xz`: one SDF per pose

Outputs:

- `my_molecular_network_colored.html`: the interactive family network. Open it in a browser.
- `processed_data/Neighborhood_Molecules/<family>/`: an SDF for every pose in the family. Next to each folder is `<family>_neighborhood_mols.pdf`, a grid of 2D structures with one entry per crystal ID, sorted by ID.
- `processed_data/cluster_mol_members.xlsx`: SMILES and crystal IDs for each family
- `processed_data/per_cluster_plifs.xlsx`: every interaction in each family
- `final_figures/all_cluster_morgan_sim.png`: within-family similarity heatmaps
- `final_figures/cluster_plifs.png`: most common interactions per family

There's also a py3Dmol cell for looking at one family inside the pocket. Set `CLUST_IDX` to choose which family.

## Docking

The 438 fragments were docked into apo CypD (PDB 8UC5) with GNINA 1.3, using a 27 Å box centred on the protein, exhaustiveness 48, seed 0 and 9 poses per ligand. The docking run itself isn't in this repo. What's here is its output (`102825_gnina_docking_results.tar.xz`), the ligands that went in (`isolated_mols.tar.xz`), and the crystal complexes superposed into the same frame (`103025_superposed_mol.tar.xz`).

GNINA reorders atoms, so each docked pose is mapped back onto the input ligand before computing RMSD. RMSD uses heavy atoms only, and each docked pose is compared against the closest crystal pose in its structure. `x0603` and `x0625` are skipped because of tautomer differences in the heavy-atom count, which leaves 436.

**The raw data containing the inputs and outputs for docking can be found at `data/gnina_docking_raw_input_output.zip`**

Outputs:

- `processed_data/min_docking_rmsd_per_structure.xlsx`: minimum RMSD for each of the 9 ranked poses, per structure
- `final_figures/gnina_min_rmsd_violin_plot.pdf`: RMSD distribution by pose rank
- `final_figures/gnina_rmsd_histogram_by_pose.pdf`: same data as a stacked histogram
- `final_figures/3A_nth_pose_inclusion_102825_gnina_docking.pdf`: number of ligands with a pose under 3 Å when the top N poses are allowed. Change `MIN_RMSD` for a different cutoff. The file name follows it.

## Binding-site conservation

`binding_site_conservation/human_cyclophilins_filtered_aligned_MAFFT.fasta` is a MAFFT alignment of 17 human cyclophilin-domain proteins from UniProt, plus the CypD construct from the crystals (chain X). `residues.csv` lists the binding-site residues, the type of interaction each one makes, and which site (colour) it belongs to.

Each residue is located in the alignment by its position in PPIF_HUMAN (CypD). At that column the notebook reports the average pairwise BLOSUM62 score and the percent identity, counting only sequences that have a residue there rather than a gap. The results are then plotted by site. If the alignment rows aren't all the same length, the notebook stops, and the markdown cell after the loading step explains how to fix it.

Output: `processed_data/conservation_results.xlsx`

## utils

- `plif.py`: builds the bit mapping and binary fingerprints from raw PLIFs
- `mol.py`: loads SDFs and PDBs from tar archives, computes RMSD with atom remapping, Morgan similarity matrices, 2D grid PDFs, SDF export and the py3Dmol viewer
- `metrics.py`: normalized AUC and RMSD helpers
- `graphs.py`: plot styling, colour maps and the heatmap function

## Citation

If you use any of this, please cite the paper above.

## License

Code is released under the MIT License (see `LICENSE`). Data in `data/`, `binding_site_conservation/` and `processed_data/` is released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
