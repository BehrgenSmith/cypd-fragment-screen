import io
import os
import sys
import tarfile
from pathlib import Path
from typing import Dict, List, Optional, Union

import cairosvg
import numpy as np
import prody as pr
import py3Dmol
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, rdFingerprintGenerator
from rdkit.Chem.Draw import rdMolDraw2D
from tqdm import tqdm


def morgan_similarity_matrix(mols, radius=2, nBits=2048, return_order=False):
    """
    Given an iterable of RDKit Mol objects, compute a Tanimoto similarity
    matrix using Morgan fingerprints (new rdFingerprintGenerator API) and
    sort it so that the most similar molecules are in the top-left.

    Returns
    -------
    sims_sorted : (N, N) np.ndarray
        Similarity matrix with rows/cols reordered.
    order : np.ndarray, optional
        Index array giving the new ordering (only if return_order=True).
    """
    # modern generator
    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=nBits)

    # fingerprints
    # fps = [mfpgen.GetFingerprint(m) for m in mols]
    fps = []
    for i, mol in enumerate(mols):
        try:
            fps.append((mfpgen.GetFingerprint(mol)))
        except:
            print(f"Failed with mol {i} : {mol}")
    # fps = np.array(fps)
    n = len(fps)

    sims = np.zeros((n, n), dtype=float)
    for i in range(n):
        sims[i, i] = 1.0
        for j in range(i + 1, n):
            s = DataStructs.TanimotoSimilarity(fps[i], fps[j])
            sims[i, j] = sims[j, i] = s

    # sort by mean similarity so “dense” ones float to top-left
    scores = sims.mean(axis=1)
    order = np.argsort(-scores)
    sims_sorted = sims[order][:, order]

    if return_order:
        return sims_sorted, order
    return sims_sorted


def calculate_matrix_in_memory(complex_name, ref_mol, dock_mol):
    """
    Calculates RMSD matrix for single RDKit Mol objects that contain
    multiple conformers (e.g., loaded via specialized multi-conf functions).
    """
    if not isinstance(ref_mol, Chem.Mol) or not isinstance(dock_mol, Chem.Mol):
        raise ValueError(f"Inputs must be RDKit Mol objects, got {type(ref_mol)}")

    # 1. Get the atomic mapping ONCE
    # The graph match is done on the molecule level, independent of conformers
    match = dock_mol.GetSubstructMatch(ref_mol)
    if not match:
        raise ValueError("Topology mismatch (Isomorphism failed)")

    match_idx = list(match)

    # 2. Extract coordinates directly from the CONFORMERS
    # Iterating over mol.GetConformers() gets all the poses stored inside the Mol object

    # Shape: (N_ref_confs, N_atoms, 3)
    ref_coords = np.array([conf.GetPositions() for conf in ref_mol.GetConformers()])

    # Shape: (N_dock_confs, N_atoms, 3) - Sliced dynamically using the match map
    dock_coords = np.array(
        [conf.GetPositions()[match_idx] for conf in dock_mol.GetConformers()]
    )

    if ref_coords.size == 0 or dock_coords.size == 0:
        raise ValueError("One of the molecules has 0 conformers loaded.")

    # 3. Broad-casted RMSD Math
    # Resulting diff shape: (N_ref_confs, N_dock_confs, N_atoms, 3)
    diff = ref_coords[:, np.newaxis, :, :] - dock_coords[np.newaxis, :, :, :]

    # Square, Sum over XYZ, Mean over Atoms, Square Root
    rmsd_matrix = np.sqrt(np.mean(np.sum(diff**2, axis=-1), axis=-1))

    return rmsd_matrix


def mol_from_sdf_bytes(data: bytes, removeHs: bool = False) -> Optional[Chem.Mol]:
    """
    Converts a multi-record SDF byte-string into a single
    multi-conformer RDKit Mol object.
    """
    # Normalize line endings at the byte level
    clean_data = data.replace(b"\r\n", b"\n").strip()

    # ForwardSDMolSupplier needs a BytesIO object for binary data
    suppl = Chem.ForwardSDMolSupplier(
        io.BytesIO(clean_data), removeHs=removeHs, sanitize=False
    )
    conf_mols = [m for m in suppl if m is not None]

    if not conf_mols:
        return None

    base = Chem.Mol(conf_mols[0])
    base.RemoveAllConformers()

    for i, m in enumerate(conf_mols):
        try:
            conf = m.GetConformer(0)
            conf.SetId(i)
            base.AddConformer(conf, assignId=True)
        except (ValueError, IndexError):
            continue

    return base


def load_sdfs_from_archive(archive_path: str, name_lambda=None) -> Dict[str, Chem.Mol]:
    mol_dict = {}
    archive_path = Path(archive_path)

    with tarfile.open(archive_path, "r:*") as tar:
        members = [m for m in tar.getmembers() if m.name.endswith(".sdf")]

        for member in tqdm(members, desc="Processing Archive"):
            f = tar.extractfile(member)
            if f is None:
                continue

            # Read as raw bytes, do not decode to string
            byte_data = f.read()
            try:
                mol = mol_from_sdf_bytes(byte_data)
            except:
                name = Path(member.name).stem
                print(f"Failed mol read for :{name}")

            if mol:
                name = Path(member.name).stem
                if name_lambda:
                    try:
                        name = name_lambda(name)
                    except Exception as e:
                        print(f"Lambda failed for {name}: {e}")

                mol_dict[name] = mol

    return mol_dict


def pdb_from_bytes(data: bytes) -> Optional[pr.AtomGroup]:
    """
    Parses a PDB from bytes using ProDy's stream parser.
    """
    # Create a stream from the bytes
    pdb_stream = io.BytesIO(data)

    # ProDy's parsePDBStream is the key here
    try:
        ag = pr.parsePDBStream(pdb_stream)
        return ag
    except Exception as e:
        print(f"⚠️ Failed to parse PDB stream: {e}")
        return None


def load_pdbs_from_archive(archive_path: str, name_lambda=None):
    # Silence ProDy noise
    pr.confProDy(verbosity="none")
    pdb_dict = {}
    archive_path = Path(archive_path)

    with tarfile.open(archive_path, "r:*") as tar:
        # Only grab actual files ending in .pdb
        # Grab actual .pdb files and explicitly ignore Mac '._' hidden files
        members = [
            m
            for m in tar.getmembers()
            if m.isfile()
            and m.name.lower().endswith(".pdb")
            and not os.path.basename(m.name).startswith("._")
        ]

        for member in tqdm(members, desc="Loading PDBs"):
            f = tar.extractfile(member)
            if f is None:
                continue

            # Temporary redirect to swallow the <ExFileObject> prints
            with open(os.devnull, "w") as devnull:
                old_stdout = sys.stdout
                sys.stdout = devnull
                try:
                    # CRITICAL FIX: Decode bytes to a regular string!
                    text_data = f.read().decode("utf-8", errors="ignore")
                    pdb_stream = io.StringIO(text_data)

                    ag = pr.parsePDBStream(pdb_stream)
                except Exception:
                    ag = None
                finally:
                    sys.stdout = old_stdout

            # ag will no longer be None!
            if ag is not None:
                # GET JUST THE FILENAME: removes '103025_superposed_mol/' prefix
                raw_filename = os.path.basename(member.name)
                clean_name = Path(raw_filename).stem  # removes '.pdb'

                if name_lambda:
                    try:
                        clean_name = name_lambda(clean_name)
                    except:
                        pass

                pdb_dict[clean_name] = ag

    print(f"✅ Successfully loaded {len(pdb_dict)} PDBs.")
    return pdb_dict


def extract_unl_ligands_heavy_coords(pdb):
    """
    Extracts coordinates of heavy atoms (non-hydrogen) for ligands
    with residue name 'UNL' from a PDB file.

    Parameters
    ----------
    pdb_path : str
        Path to the PDB file.

    Returns
    -------
    List[np.ndarray]
        List of NumPy arrays, one per UNL ligand.
        Each array has shape (N_atoms, 3) for the heavy-atom coordinates.
    """
    structure = pdb
    # Select all UNL residues (ligands)
    unl_sel = structure.select("resname UNL and not element H")
    if unl_sel is None:
        print("No UNL residues found.")
        return []

    lig_coords = []
    for res in structure.iterResidues():
        if res.getResname() == "UNL":
            heavy_atoms = res.select("not element H")
            if heavy_atoms is not None:
                lig_coords.append(heavy_atoms.getCoords())

    return lig_coords


def get_non_h_coords_per_conf(mol: Chem.Mol):
    """
    Returns a list of NumPy arrays containing the 3D coordinates
    of all non-hydrogen atoms for each conformer in the molecule.

    Parameters
    ----------
    mol : Chem.Mol
        RDKit molecule with one or more conformers.

    Returns
    -------
    List[np.ndarray]
        List of (N_atoms, 3) arrays, one per conformer, for non-hydrogen atoms.
    """
    non_h_indices = [
        atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() != "H"
    ]
    coords_list = []

    for conf in mol.GetConformers():
        coords = np.array([list(conf.GetAtomPosition(i)) for i in non_h_indices])
        coords_list.append(coords)

    return coords_list


def mols_to_2d_tile_image(
    mols: List[Chem.Mol],
    legends: Optional[List[str]] = None,
    mols_per_row: int = 4,
    sub_img_size=(300, 300),  # Doubled default size for better detail
    pdf_path: Optional[Union[str, Path]] = None,
):
    # 1. Prepare molecules (2D coords)
    mols_2d = []
    for m in mols:
        if m is None:
            continue
        m2 = Chem.Mol(m)
        AllChem.Compute2DCoords(m2)
        mols_2d.append(m2)

    # 2. Setup the Drawing Canvas (MolDraw2DSVG is the key for Vectors)
    n_mols = len(mols_2d)
    n_rows = (n_mols + mols_per_row - 1) // mols_per_row
    full_width = mols_per_row * sub_img_size[0]
    full_height = n_rows * sub_img_size[1]

    drawer = rdMolDraw2D.MolDraw2DSVG(
        full_width, full_height, sub_img_size[0], sub_img_size[1]
    )

    # Optional: Tweak drawing options for "Publication Quality"
    options = drawer.drawOptions()
    options.addStereoAnnotation = True
    options.prepareMolsBeforeDrawing = True
    options.legendFontSize = 20  # Make legends readable
    options.bondLineWidth = 2  # Thicker bonds look better in PDFs

    # 3. Draw the grid
    drawer.DrawMolecules(mols_2d, legends=legends if legends else None)
    drawer.FinishDrawing()
    svg_text = drawer.GetDrawingText()

    # 4. Save to PDF using cairosvg
    if pdf_path:
        # Ensure we are passing a string/bytes to cairosvg
        cairosvg.svg2pdf(bytestring=svg_text.encode("utf-8"), write_to=str(pdf_path))

    return svg_text


def save_mol_to_sdf(mol, output_dir, pose=None):
    """
    Saves an RDKit molecule to an SDF file using its internal name.
    """
    # 1. Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # 2. Extract the molecule name
    # Fallback to 'unnamed_mol' if the property isn't set
    if mol.HasProp("_Name"):
        mol_name = mol.GetProp("_Name").strip()
    else:
        mol_name = "unnamed_mol"
    if pose is not None:
        mol_name = f"{mol_name}-{pose}"
    # 3. Sanitize filename (remove characters that might break file systems)
    # This replaces spaces with underscores and keeps it simple
    safe_name = "".join(
        [c if c.isalnum() or c in ("-", "_") else "_" for c in mol_name]
    )
    file_path = os.path.join(output_dir, f"{safe_name}.sdf")

    # 4. Write to SDF
    writer = Chem.SDWriter(file_path)
    try:
        writer.write(mol)
        print(f"Successfully saved: {file_path}")
    except Exception as e:
        print(f"Failed to write {mol_name}: {e}")
    finally:
        writer.close()

def visualize_complex(pdb_file, rdkit_dict):
    """
    pdb_file: Path to your reference PDB
    rdkit_dict: {'name': mol_object, ...}
    """
    # 1. Initialize the viewer
    view = py3Dmol.view(width=800, height=600)

    # 2. Load the Reference Protein
    with open(pdb_file, 'r') as f:
        pdb_data = f.read()

    view.addModel(pdb_data, 'pdb')
    view.setStyle({'model': 0}, {"cartoon": {'color': 'spectrum'}})

    # 3. Load the RDKit Ligands
    # We start the index at 1 because the protein is index 0
    for i, (_name, mol) in enumerate(rdkit_dict.items(), start=1):
        # Ensure the molecule has 3D coordinates!
        mb = Chem.MolToMolBlock(mol)

        view.addModel(mb, 'sdf')
        # 'stick' is the py3Dmol equivalent of 'licorice'
        view.setStyle({'model': i}, {'stick': {'colorscheme': 'cyanCarbon'}})

    # 4. Final touches
    view.zoomTo()
    return view.show()