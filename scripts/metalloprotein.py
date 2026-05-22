###################################################################################################################
#   Automated Force Fields for Metals     /$$$$$$$   /$$$$$$  /$$$$$$$  /$$      /$$                              # 
#                                        | $$__  $$ /$$__  $$| $$__  $$| $$$    /$$$                              #
#   /$$$$$$   /$$$$$$   /$$$$$$$ /$$   /$$| $$  \ $$| $$  \ $$| $$  \ $$| $$$$  /$$$$                             #
#  /$$__  $$ |____  $$ /$$_____/| $$  | $$| $$$$$$$/| $$$$$$$$| $$$$$$$/| $$ $$/$$ $$                             #
# | $$$$$$$$  /$$$$$$$|  $$$$$$ | $$  | $$| $$____/ | $$__  $$| $$__  $$| $$  $$$| $$                             #
# | $$_____/ /$$__  $$ \____  $$| $$  | $$| $$      | $$  | $$| $$  \ $$| $$\  $ | $$                             #
# |  $$$$$$$|  $$$$$$$ /$$$$$$$/|  $$$$$$$| $$      | $$  | $$| $$  | $$| $$ \/  | $$                             #
#  \_______/ \_______/|_______/  \____  $$|__/      |__/  |__/|__/  |__/|__/     |__/                             #
#                               /$$  | $$                                                                         #
#                              |  $$$$$$/              Ver. 4.25 - 19 May 2026                                    #
#                               \______/                                                                          #
#                                                                                                                 #
# Developer: Abdelazim M. A. Abdelgawwad.                                                                         #
# Institut de Ciència Molecular (ICMol), Universitat de València, P.O. Box 22085, València 46071, Spain           #
#                                                                                                                 #
#Distributed under the GNU LESSER GENERAL PUBLIC LICENSE Version 2.1, February 1999                               #
#Copyright 2024 Abdelazim M. A. Abdelgawwad, Universitat de València. E-mail: abdelazim.abdelgawwad@uv.es         #
###################################################################################################################

from Bio import PDB
import numpy as np
import re
import sys
import shutil
import os 

def extract_atom_types_from_mol2(file_path):
    """
    Parse a MOL2 file and extract atom type assignments.
    
    Reads the @<TRIPOS>ATOM section of a MOL2 file and builds a dictionary
    mapping atom IDs (int) to their GAFF2/AMBER atom type strings (column 6).
    
    Args:
        file_path (str): Path to the MOL2 file to parse.
        
    Returns:
        dict: {atom_id: atom_type} mapping, e.g. {1: 'ns', 2: 'hn', 133: 'Mn'}.
              Returns empty dict on error.
    """
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
        
        atom_types = {}
        is_atom_section = False
        
        for line in lines:
            if line.startswith("@<TRIPOS>ATOM"):
                is_atom_section = True
                continue
            if line.startswith("@<TRIPOS>BOND"):
                is_atom_section = False
                break
            if is_atom_section and line.strip():
                parts = line.split()
                if len(parts) >= 6:
                    atom_id = int(parts[0])
                    atom_type = parts[5]
                    atom_types[atom_id] = atom_type
        
        return atom_types
    except Exception as e:
        print(f"Error reading {file_path}: {str(e)}")
        return {}

#Read atom type mappings from library file
def read_aminofb15_lib(file_path):
    """
    Read atom type mappings from an AMBER force field library (.lib) file.
    
    Parses the '!entry.<RESNAME>.unit.atoms table' sections to extract
    the mapping from atom names to their force field atom types for each
    standard residue (e.g., HID: {'N': 'N', 'CA': 'XC', 'CB': 'CT', ...}).
    
    Args:
        file_path (str): Path to the .lib file (e.g., ff19SB.lib, fb15.lib).
        
    Returns:
        dict: Nested dict {residue_name: {atom_name: atom_type}}.
              Returns empty dict on error.
    """
    residue_atom_types = {}
    current_residue = None
    
    try:
        with open(file_path, 'r') as file:
            reading_atoms = False
            for line in file:
                if '!entry.' in line and '.unit.atoms table' in line:
                    current_residue = line.split('.')[1]
                    reading_atoms = True
                    residue_atom_types[current_residue] = {}
                elif '!entry.' in line and '.unit.atomspertinfo table' in line:
                    reading_atoms = False
                elif reading_atoms and line.strip() and '"' in line:
                    parts = line.strip().split('"')
                    if len(parts) >= 4:
                        orig_atom_name = parts[1]
                        new_atom_type = parts[3]
                        residue_atom_types[current_residue][orig_atom_name] = new_atom_type
        
        return residue_atom_types
    except Exception as e:
        print(f"Error reading library file: {e}")
        return {}

#Write initial PDB with MOL2 atom types
def write_pdb_with_mol2_types(structure, output_file, atom_types):
    """
    Write a PDB file where the atom name field contains MOL2 atom types.
    
    Iterates through a BioPython Structure object and writes each atom
    in PDB format, but replaces the standard atom name (columns 13-16)
    with the corresponding GAFF2 atom type from the MOL2 file. This creates
    a hybrid PDB used as an intermediate for atom type reassignment.
    
    Args:
        structure: BioPython Structure object containing the extracted residues.
        output_file (str): Path to write the output PDB (e.g., 'QM.pdb').
        atom_types (dict): {sequential_atom_id: mol2_atom_type} from extract_atom_types_from_mol2().
    """
    atom_counter = 1
    with open(output_file, 'w') as f:
        for model in structure:
            for chain in model:
                for residue in chain:
                    hetero, resnum, icode = residue.get_id()
                    for atom in residue:
                        record_type = "HETATM" if residue.get_id()[0] != " " else "ATOM"
                        mol2_type = atom_types.get(atom_counter, atom.name)
                        x, y, z = atom.get_coord()
                        resname = residue.get_resname()
                        chainid = chain.get_id()
                        occupancy = atom.get_occupancy()
                        bfactor = atom.get_bfactor()
                        element = atom.element if hasattr(atom, 'element') else '  '
                        
                        line = f"{record_type:<6}{atom_counter:>5}  {mol2_type:<3} {resname:>3} {chainid}{resnum:>4}    "
                        line += f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{occupancy:>6.2f}{bfactor:>6.2f}          {element:>2}"
                        f.write(line + '\n')
                        atom_counter += 1

#Update atom types in QM.pdb based on REFQM.pdb and library
def update_atom_types_from_library(qm_pdb, ref_pdb, lib_types, terminal_info):
    """
    Replace MOL2 atom types in QM.pdb with protein FF atom types from the library.
    
    For each atom in QM.pdb, looks up its residue name and atom name from REFQM.pdb
    (which has the original PDB atom names), then finds the corresponding protein FF
    atom type (e.g., ff19SB type) from the library. Handles N/C-terminal residues by
    prepending 'N' or 'C' to the residue name for library lookup.
    
    Only standard amino acid atoms that exist in the library get updated; non-standard
    atoms (metals, ligands) retain their original MOL2/GAFF2 types.
    
    Note: Uses zip(qm_lines, ref_lines) which assumes both files have identical line
    counts. REFQM.pdb may contain TER/END records that cause misalignment.
    
    Args:
        qm_pdb (str): Path to QM.pdb (will be overwritten with updated types).
        ref_pdb (str): Path to REFQM.pdb (reference with original atom names).
        lib_types (dict): {residue_name: {atom_name: ff_type}} from read_aminofb15_lib().
        terminal_info (dict): {(chain_id, res_id): 'NTERM'|'CTERM'} for terminal residues.
    """
    with open(ref_pdb, 'r') as ref_file:
        ref_lines = ref_file.readlines()
    
    updated_lines = []
    with open(qm_pdb, 'r') as qm_file:
        qm_lines = qm_file.readlines()
        
        for qm_line, ref_line in zip(qm_lines, ref_lines):
            if qm_line.startswith(('ATOM', 'HETATM')):
                ref_resname = ref_line[17:20].strip()
                ref_atomname = ref_line[12:16].strip()
                ref_chain = ref_line[21:22].strip()
                
                # Parse residue ID
                try:
                    hetero = ' '
                    if ref_line.startswith('HETATM'):
                        hetero = ref_line[17:20].strip()
                    resnum = int(ref_line[22:26].strip())
                    icode = ref_line[26:27].strip() if len(ref_line) > 26 else ' '
                    ref_res_id = (hetero if hetero != ' ' else ' ', resnum, icode)
                    
                    # Check if this residue is a terminal
                    res_key = (ref_chain, ref_res_id)
                    if res_key in terminal_info:
                        terminal_type = terminal_info[res_key]
                        # Add N or C prefix for library lookup
                        lookup_name = ('N' if terminal_type == 'NTERM' else 'C') + ref_resname
                    else:
                        lookup_name = ref_resname
                except:
                    lookup_name = ref_resname
                
                # Look up atom type in library using terminal-aware name
                if lookup_name in lib_types and ref_atomname in lib_types[lookup_name]:
                    new_type = lib_types[lookup_name][ref_atomname]
                    new_line = qm_line[:12] + f"{new_type:<3}" + qm_line[16:]
                    updated_lines.append(new_line)
                else:
                    updated_lines.append(qm_line)
            else:
                updated_lines.append(qm_line)
    
    with open(qm_pdb, 'w') as f:
        f.writelines(updated_lines)

#Update atom types in MOL2 file using atom types from QM.pdb while preserving the original formatting.
def update_mol2_with_qm_types(mol2_file, qm_pdb, output_mol2):

    """
    Update the atom type column in a MOL2 file using types from QM.pdb.
    
    After QM.pdb has been updated with protein FF atom types (via 
    update_atom_types_from_library), this function propagates those types back
    into the MOL2 file. Standard amino acid atoms get ff19SB types while
    non-standard atoms keep their GAFF2 types.
    
    The output MOL2 (COMPLEX_updated.mol2) is a hybrid: protein FF types for
    standard residue atoms, GAFF2 types for the metal center and ligands.
    
    Args:
        mol2_file (str): Path to input MOL2 (NEW_COMPLEX.mol2).
        qm_pdb (str): Path to QM.pdb with updated atom types.
        output_mol2 (str): Path for output MOL2 (COMPLEX_updated.mol2).
    """
    # Read QM.pdb atom types
    qm_atom_types = {}
    with open(qm_pdb, 'r') as qm_file:
        for line in qm_file:
            if line.startswith(('ATOM', 'HETATM')):
                try:
                    atom_number = int(line[6:11].strip())
                    atom_type = line[12:16].strip()
                    qm_atom_types[atom_number] = atom_type
                except ValueError:
                    print(f"Warning: Skipping malformed QM.pdb line: {line.strip()}")

    # Read and update MOL2 file
    updated_mol2_lines = []
    try:
        with open(mol2_file, 'r') as mol2:
            lines = mol2.readlines()
            
            is_atom_section = False
            for line in lines:
                if line.startswith("@<TRIPOS>ATOM"):
                    is_atom_section = True
                    updated_mol2_lines.append(line)
                    continue
                
                if line.startswith("@<TRIPOS>BOND"):
                    is_atom_section = False
                
                if is_atom_section and line.strip():
                    try:
                        # Use regex to parse MOL2 atom line while preserving spacing
                        atom_match = re.match(
                            r"(\s*\d+\s+)(\S+)(\s+-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+\s+)(\S+)(\s+\d+\s+\S+\s+-?\d+\.\d+)", 
                            line
                        )
                        
                        if atom_match:
                            # Extract groups
                            pre_id, atom_name, coords, atom_type, post_type = atom_match.groups()
                            atom_id = int(pre_id.strip())

                            # Update atom type if found in QM.pdb
                            if atom_id in qm_atom_types:
                                atom_type = f"{qm_atom_types[atom_id]:<2}"  # Left-align to match MOL2 format

                            updated_line = f"{pre_id}{atom_name}{coords}{atom_type}{post_type}\n"
                            updated_mol2_lines.append(updated_line)
                        else:
                            updated_mol2_lines.append(line)
                    except Exception as e:
                        print(f"Warning: Skipping malformed MOL2 line: {line.strip()} ({e})")
                        updated_mol2_lines.append(line)  # Keep the original line
                else:
                    updated_mol2_lines.append(line)
        
        # Write updated MOL2 file
        with open(output_mol2, 'w') as output_file:
            output_file.writelines(updated_mol2_lines)
            print("COMPLEX_updated.mol2 has been created with updated atom types.")
            
    except Exception as e:
        print(f"Error updating MOL2 file: {e}")

#Identify N-terminal and C-terminal residues for each chain.
def identify_terminal_residues(structure, standard_residues):
    
    """
    Identify N-terminal and C-terminal residues in each chain.
    
    For each chain in the structure, finds the first and last standard amino
    acid residues and marks them as NTERM/CTERM. This information is used to
    look up terminal-specific atom types from the library (e.g., NHID, CHID).
    
    Args:
        structure: BioPython Structure object.
        standard_residues (set): Set of recognized standard residue names.
        
    Returns:
        dict: {(chain_id, residue_id_tuple): 'NTERM'|'CTERM'} for terminal residues.
    """
    terminal_info = {}

    # Base residue names (without N/C prefix)
    base_residues = {res for res in standard_residues if not res.startswith(('N', 'C'))}

    for model in structure:
        for chain in model:
            residues = [res for res in chain if res.get_resname() in base_residues]

            if len(residues) == 0:
                continue

            # Mark first standard residue as N-terminal
            first_res = residues[0]
            terminal_info[(chain.id, first_res.get_id())] = 'NTERM'

            # Mark last standard residue as C-terminal
            last_res = residues[-1]
            terminal_info[(chain.id, last_res.get_id())] = 'CTERM'

    return terminal_info


def analyze_and_extract_metal_site(input_pdb, mol2_file="NEW_COMPLEX.mol2", lib_file="amber_refernce.lib", 
                                   output_pdb="QM.pdb", distance_cutoff=2.5):
    """
    Core function: analyze metal coordination and extract the metal site from the protein.
    
    This is the main analysis step that:
    1. Parses the protein PDB and identifies all metal atoms
    2. Finds residues within 2.5 Å of each metal (coordination shell)
    3. Finds standard residues covalently bonded (<1.9 Å) to non-standard coordinating residues
    4. Rebuilds a new structure containing only the extracted residues (preserving order)
    5. Writes REFQM.pdb (reference with original atom names) and QM.pdb (with MOL2 types)
    6. Updates QM.pdb atom types using the protein FF library for standard residue atoms
    7. Writes fixed_charges.dat listing backbone atoms (N, CA, C, O) that need fixed charges
    
    Args:
        input_pdb (str): Path to metalloprotein_easyPARM.pdb (cleaned protein structure).
        mol2_file (str): Path to NEW_COMPLEX.mol2 (full QM cluster with GAFF2 types/charges).
        lib_file (str): Basename of the FF library file (looked up in script directory).
        output_pdb (str): Path for QM.pdb output.
        distance_cutoff (float): Metal coordination distance cutoff in Å (default 2.5).
        
    Returns:
        tuple: (metal_coordination dict, terminal_info dict)
            - metal_coordination: {metal_key: {'metal_position': coord, 'coordinating_residues': [...]}}
            - terminal_info: {(chain_id, res_id): 'NTERM'|'CTERM'}
    """
    metals=['MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'W', 'RE', 'OS', 'IR', 'PT', 'AU', 'NA', 'K', 'LI', 'RB', 'CS', 'MG', 'CA', 'SR', 'BA', 'V', 'CR', 'CD', 'HG', 'AL', 'GA', 'IN', 'SN', 'PB', 'BI', 'LA', 'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB', 'LU', 'FE2', 'FE3', 'FE4', 'CU1', 'CU2', 'MN2', 'MN3', 'MN4', 'CO2', 'CO3', 'NI2', 'NI3', 'V2', 'V3', 'V4', 'V5']    
    # Define standard amino acid residues
    standard_residues = {
        # Base residues
        "ALA", "ARG", "ASH", "ASN", "ASP", "CYM", "CYS", "CYX", "GLH", "GLN",
        "GLU", "GLY", "HID", "HIE", "HIP", "HYP", "ILE", "LEU", "LYN", "LYS",
        "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "NHE", "NME",
        "ACE",
        # N-terminal variants
        "NALA", "NARG", "NASH", "NASN", "NASP", "NCYM", "NCYS", "NCYX", "NGLH", "NGLN",
        "NGLU", "NGLY", "NHID", "NHIE", "NHIP", "NHYP", "NILE", "NLEU", "NLYN", "NLYS",
        "NMET", "NPHE", "NPRO", "NSER", "NTHR", "NTRP", "NTYR", "NVAL", "NHIS",
        # C-terminal variants
        "CALA", "CARG", "CASH", "CASN", "CASP", "CCYM", "CCYS", "CCYX", "CGLH", "CGLN",
        "CGLU", "CGLY", "CHID", "CHIE", "CHIP", "CHYP", "CILE", "CLEU", "CLYN", "CLYS",
        "CMET", "CPHE", "CPRO", "CSER", "CTHR", "CTRP", "CTYR", "CVAL", "CHIS"
    }

    # Get the script paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    lib_file_path = os.path.join(script_dir, lib_file)

    # Initialize PDB parser
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('protein', input_pdb)
    
    # Dictionaries and sets to store results
    metal_coordination = {}
    residues_to_extract = set()
    original_order = []  # Maintain the original residue order

    # Step 1: Analyze coordination and track residues
    for model in structure:
        for chain in model:
            for residue in chain:
                residue_key = (chain.id, residue.get_id())
                original_order.append(residue_key)
                
                # Check for metal atoms
                for atom in residue:
                    if atom.element in metals:
                        residues_to_extract.add(residue_key)
                        metal_coord = atom.coord
                        key = f"{atom.element}_{chain.id}_{residue.get_id()}"
                        metal_coordination[key] = {
                            'metal_position': metal_coord,
                            'coordinating_residues': []
                        }

                        # Find nearby residues within cutoff
                        for chain2 in model:
                            for residue2 in chain2:
                                for atom2 in residue2:
                                    distance = np.linalg.norm(metal_coord - atom2.coord)
                                    if distance <= distance_cutoff:
                                        coord_info = {
                                            'chain': chain2.id,
                                            'residue_number': residue2.get_id()[1],
                                            'residue_name': residue2.get_resname(),
                                            'atom_name': atom2.get_name().strip(),
                                            'distance': round(distance, 2),
                                            'atom_coord': atom2.coord,
                                            'element': atom2.element
                                        }
                                        metal_coordination[key]['coordinating_residues'].append(coord_info)
                                        residues_to_extract.add((chain2.id, residue2.get_id()))
    
    # NEW SECTION: Find standard residues linked to non-standard coordinating residues
    # Collect all heavy atoms for neighbor search
    heavy_atoms = []
    non_standard_coordinating = set()
    
    # Identify non-standard coordinating residues
    for model in structure:
        for chain in model:
            for residue in chain:
                res_key = (chain.id, residue.get_id())
                if res_key in residues_to_extract:
                    # Track if this is a non-standard coordinating residue
                    if residue.get_resname() not in standard_residues:
                        non_standard_coordinating.add(res_key)
                
                # Collect all heavy atoms for neighbor search
                for atom in residue:
                    if atom.element != 'H':  # Only include heavy atoms
                        heavy_atoms.append(atom)
    
    # Create neighbor search for all heavy atoms
    ns = PDB.NeighborSearch(heavy_atoms)
    bond_cutoff = 1.9  # Covalent bond distance cutoff in Å
    
    # Check for standard residues bonded to non-standard coordinating residues
    for res_key in non_standard_coordinating:
        chain_id, res_id = res_key
        non_std_residue = structure[0][chain_id][res_id]
        
        # Check each heavy atom in the non-standard residue
        for atom in non_std_residue:
            if atom.element != 'H':
                # Find potential bonding partners
                nearby_atoms = ns.search(atom.coord, bond_cutoff, level='A')
                
                for nearby_atom in nearby_atoms:
                    if nearby_atom == atom:
                        continue  # Skip self
                    
                    nearby_res = nearby_atom.get_parent()
                    nearby_res_key = (nearby_res.get_parent().id, nearby_res.get_id())
                    
                    # If it's a standard residue and not already marked for extraction
                    if (nearby_res.get_resname() in standard_residues and 
                        nearby_res_key != res_key and 
                        nearby_res_key not in residues_to_extract):
                        
                        residues_to_extract.add(nearby_res_key)
                        
                        # If not already in original_order, add it
                        if nearby_res_key not in original_order:
                            original_order.append(nearby_res_key)

    # Step 2: Rebuild structure, explicitly handling hydrogens
    new_structure = PDB.Structure.Structure('metal_site')
    new_model = PDB.Model.Model(0)
    new_structure.add(new_model)
    new_chain = PDB.Chain.Chain('A')
    new_model.add(new_chain)

    # Add residues in original order, including all hydrogens
    for chain_id, res_id in original_order:
        if (chain_id, res_id) in residues_to_extract:
            original_residue = structure[0][chain_id][res_id]
            new_residue = PDB.Residue.Residue(original_residue.id, original_residue.resname, original_residue.segid)

            # Add all atoms, avoiding the skipping of duplicate names
            for atom in original_residue.get_atoms():
                new_atom = PDB.Atom.Atom(
                    atom.get_name(), atom.coord, atom.bfactor, atom.occupancy,
                    atom.altloc, atom.fullname, atom.serial_number, element=atom.element
                )
                new_residue.add(new_atom)
            new_chain.add(new_residue)

    # Step 3: Save the new structure
    ref_output = "REFQM.pdb"
    io = PDB.PDBIO()
    io.set_structure(new_structure)
    io.save(ref_output)
    
    # Identify terminals
    terminal_info = identify_terminal_residues(new_structure, standard_residues)

    # Step 4: Extract specific atoms and write to fixed_charges.dat
    target_atoms = {'N', 'CA', 'C', 'O'}
    # Step 4.1 — read ref_output and store {atom_number: atom_name_ref}
    ref_atoms = {}

    with open(ref_output, 'r') as ref_file:
        for line in ref_file:
            if line.startswith("ATOM"):
                atom_name = line[12:16].strip()
                if atom_name in target_atoms:
                    atom_number = int(line[6:11].strip())
                    ref_atoms[atom_number] = atom_name
    # Step 5: Process MOL2 and library files
    mol2_atom_types = extract_atom_types_from_mol2(mol2_file)
    write_pdb_with_mol2_types(new_structure, output_pdb, mol2_atom_types)

    lib_types = read_aminofb15_lib(lib_file_path)
    update_atom_types_from_library(output_pdb, ref_output, lib_types, terminal_info)

    # Step 5.1 — read output_pdb and match by atom_number
    with open(output_pdb, 'r') as out_file, open("fixed_charges.dat", 'w') as outfile:
        for line in out_file:
            if line.startswith("ATOM"):
                atom_number = int(line[6:11].strip())

                if atom_number in ref_atoms:
                    atom_name_out = line[12:16].strip()
                    residue_name = line[16:20].strip()
                    atom_name_ref = ref_atoms[atom_number]

                    outfile.write(f"{atom_number} {atom_name_ref} {atom_name_out} {residue_name}\n")
    
    return metal_coordination, terminal_info


#Extracts all residue names from the given PDB file.
def get_existing_residues(pdb_file):
    """
    Extract all unique residue names from a PDB file.
    
    Used to check for name collisions when generating unique 3-character
    residue names for metal-coordinating standard residues.
    
    Args:
        pdb_file (str): Path to a PDB file (typically metalloprotein_easyPARM.pdb).
        
    Returns:
        set: Set of residue name strings found in the file.
    """
    existing_residues = set()
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                resname = line[17:20].strip()
                existing_residues.add(resname)
    return existing_residues

#Generates a unique residue name based on base_name. Increments a counter until a unique name is found.

def generate_unique_residue_name(base_name, existing_residues, residue_type_count):
    """
    Generate a unique 3-character residue name that doesn't conflict with existing names.
    
    Uses the first 2 characters of the original residue name as a base (e.g., 'HI' for HID,
    'CY' for CYS) and appends/prepends digits in three patterns:
      - Pattern 1: HI1, HI2, ..., HI9
      - Pattern 2: 1HI, 2HI, ..., 9HI  
      - Pattern 3: H1I, H2I, ..., H9I
    Cycles through patterns until a unique name is found.
    
    Args:
        base_name (str): 2-character base (e.g., 'HI', 'CY', 'AS', 'GL').
        existing_residues (set): Names already in use (modified in-place when new name added).
        residue_type_count (dict): Counter for each base_name (modified in-place).
        
    Returns:
        str: Unique 3-character residue name (e.g., 'HI1', '2CY', 'A1S').
    """
    # Initialize counter if not already done
    if base_name not in residue_type_count:
        residue_type_count[base_name] = 1
    else:
        residue_type_count[base_name] += 1

    while True:
        count = residue_type_count[base_name]

        # Determine which pattern and digit to use
        if count <= 9:
            # Pattern 1: CY1-CY9 (digit at position 2)
            candidate = f"{base_name}{count}"
        elif count <= 18:
            # Pattern 2: 1CY-9CY (digit at position 0)
            digit = count - 9
            candidate = f"{digit}{base_name}"
        elif count <= 27:
            # Pattern 3: C1Y-C9Y (digit at position 1)
            digit = count - 18
            candidate = f"{base_name[0]}{digit}{base_name[1]}"
        else:
            # Cycle back: repeat the patterns
            cycle_position = ((count - 1) % 27) + 1
            if cycle_position <= 9:
                candidate = f"{base_name}{cycle_position}"
            elif cycle_position <= 18:
                digit = cycle_position - 9
                candidate = f"{digit}{base_name}"
            else:
                digit = cycle_position - 18
                candidate = f"{base_name[0]}{digit}{base_name[1]}"

        if candidate not in existing_residues:
            # Add candidate to the set to avoid future conflicts
            existing_residues.add(candidate)
            return candidate
        residue_type_count[base_name] += 1

#Extract non-standard residues (including those standard residues near metals),
#standard residues linked to non-standard residues, and generate QM PDB and XYZ files.
def extract_non_standard_residues_from_ref(ref_pdb, output_pdb="part_QM.pdb", output_xyz="part_QM.xyz", 
                                           standard_residues=None, all_residues_pdb="metalloprotein_easyPARM.pdb",
                                           terminal_info=None):
    
    """
    Split the extracted metal site into non-standard (metal+ligand) and standard (amino acid) parts.
    
    Reads REFQM.pdb and separates residues into:
    - Non-standard residues (metal centers, non-AA ligands) → part_QM.pdb + part_QM.xyz
    - Standard amino acid residues near metals → individual PDB/XYZ files with unique names
      (e.g., HI1.pdb, HI1.xyz, CY1.pdb, CY1.xyz)
    
    Standard residues are renamed to unique 3-char names to avoid conflicts with the protein
    FF definitions (since they will get custom charges from the QM calculation).
    
    Args:
        ref_pdb (str): Path to REFQM.pdb (extracted metal site with original names).
        output_pdb (str): Path for non-standard residues PDB output (part_QM.pdb).
        output_xyz (str): Path for non-standard residues XYZ output (part_QM.xyz).
        standard_residues (set): Optional set of standard residue names.
        all_residues_pdb (str): Full protein PDB for checking name collisions.
        terminal_info (dict): Terminal residue information.
        
    Returns:
        dict: {(resname, chain, resnum): new_3char_name} mapping for renamed standard residues.
    """
    metals=['MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 'MO', 'TC', 'RU', 'RH', 'PD', 'AG', 'W', 'RE', 'OS', 'IR', 'PT', 'AU', 'NA', 'K', 'LI', 'RB', 'CS', 'MG', 'CA', 'SR', 'BA', 'V', 'CR', 'CD', 'HG', 'AL', 'GA', 'IN', 'SN', 'PB', 'BI', 'LA', 'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB', 'LU', 'FE2', 'FE3', 'FE4', 'CU1', 'CU2', 'MN2', 'MN3', 'MN4', 'CO2', 'CO3', 'NI2', 'NI3', 'V2', 'V3', 'V4', 'V5']

    element_converter = {
        'ZN': 'Zn', 'FE': 'Fe', 'CO': 'Co', 'RU': 'Ru', 'IR': 'Ir', 'PT': 'Pt', 'AU': 'Au', 'AG': 'Ag', 
        'CU': 'Cu', 'MG': 'Mg', 'MN': 'Mn', 'NI': 'Ni', 'PD': 'Pd', 'CD': 'Cd', 'HG': 'Hg', 'BR': 'Br',
        'CL': 'Cl', 'AL': 'Al', 'GA': 'Ga', 'IN': 'In', 'SB': 'Sb', 'TL': 'Tl', 'PB': 'Pb', 'BI': 'Bi',
        'AS': 'As', 'SE': 'Se', 'SR': 'Sr', 'MO': 'Mo', 'TC': 'Tc', 'RE': 'Re', 'OS': 'Os', 'RH': 'Rh', 
        'NA': 'Na', 'BA': 'Ba', 'SI': 'Si'
    }

    if standard_residues is None:
        standard_residues = {
            # Base residues
            "ALA", "ARG", "ASH", "ASN", "ASP", "CYM", "CYS", "CYX", "GLH", "GLN",
            "GLU", "GLY", "HID", "HIE", "HIP", "HYP", "ILE", "LEU", "LYN", "LYS",
            "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "NHE", "NME",
            "ACE",
            # N-terminal variants
            "NALA", "NARG", "NASH", "NASN", "NASP", "NCYM", "NCYS", "NCYX", "NGLH", "NGLN",
            "NGLU", "NGLY", "NHID", "NHIE", "NHIP", "NHYP", "NILE", "NLEU", "NLYN", "NLYS",
            "NMET", "NPHE", "NPRO", "NSER", "NTHR", "NTRP", "NTYR", "NVAL", 
            # C-terminal variants
            "CALA", "CARG", "CASH", "CASN", "CASP", "CCYM", "CCYS", "CCYX", "CGLH", "CGLN",
            "CGLU", "CGLY", "CHID", "CHIE", "CHIP", "CHYP", "CILE", "CLEU", "CLYN", "CLYS",
            "CMET", "CPHE", "CPRO", "CSER", "CTHR", "CTRP", "CTYR", "CVAL"
        }

    if terminal_info is None:
        terminal_info = {}
    # Get complete set of residue names from the input PDB file
    existing_residues = get_existing_residues(all_residues_pdb)

    # Parse structure for structural analysis
    parser = PDB.PDBParser(QUIET=True)
    try:
        structure = parser.get_structure('protein', ref_pdb)
    except Exception as e:
        print(f"Error parsing PDB file {ref_pdb}: {e}")
        raise

    # Store metal coordinates
    metal_coords = []
    # For standard residues that will be renamed
    standard_res_atoms = {}
    # For non-standard residues (or standard ones not near metals)
    nonstandard_res_atoms = []
    # Mapping from original residue identifier (tuple) to new name
    residue_name_mapping = {}
    # Counter for each base name (e.g., "CY" or "HI")
    residue_type_count = {}
    # Track residues to extract
    residues_to_extract = set()
    # Track non-standard residues for later bond checking
    non_standard_res_keys = set()

    # First pass: collect metal coordinates and identify coordinating residues
    with open(ref_pdb, 'r') as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                element = line.split()[-1]
                resname = line[17:20].strip()
                chain = line[21:22].strip()
                resnum = int(line[22:26])
                res_key = (resname, chain, resnum)
                
                # Track all non-standard residues
                if resname not in standard_residues:
                    non_standard_res_keys.add(res_key)
                
                if element.upper() in metals:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    metal_coords.append(np.array([x, y, z]))
                    # Get residue info for this metal
                    residues_to_extract.add(res_key)

    # Build dictionary of heavy atoms for neighbor search
    heavy_atoms = []
    residue_atoms = {}  # Maps residue key to its atoms
    atom_residue_map = {}  # Maps atom object to its residue key
    
    for model in structure:
        for chain in model:
            for residue in chain:
                res_id = residue.get_id()
                resname = residue.get_resname()
                chain_id = chain.id
                resnum = res_id[1]
                res_key = (resname, chain_id, resnum)
                
                # Initialize residue atoms list if not already there
                if res_key not in residue_atoms:
                    residue_atoms[res_key] = []
                
                # Check if this residue is near a metal
                for atom in residue:
                    if atom.element != 'H':  # Only use heavy atoms
                        # Add to heavy atoms list
                        heavy_atoms.append(atom)
                        # Store atom in residue atoms list
                        residue_atoms[res_key].append(atom)
                        # Map this atom to its residue
                        atom_residue_map[atom.get_full_id()] = res_key
                        
                        # Check proximity to metals
                        coord = atom.get_coord()
                        for metal_coord in metal_coords:
                            if np.linalg.norm(metal_coord - coord) <= 2.5:
                                residues_to_extract.add(res_key)
                                break

    # Create neighbor search for all heavy atoms
    ns = PDB.NeighborSearch(heavy_atoms)
    bond_cutoff = 1.9  # Covalent bond distance cutoff

    # Check for standard residues bonded to non-standard residues (not just those coordinating metals)
    for ns_res_key in non_standard_res_keys:
        if ns_res_key in residue_atoms:  # Make sure we have atoms for this residue
            # Process each heavy atom in the non-standard residue
            for atom in residue_atoms[ns_res_key]:
                # Find potential bonding partners
                nearby_atoms = ns.search(atom.get_coord(), bond_cutoff, level='A')
                
                for nearby_atom in nearby_atoms:
                    if nearby_atom != atom:  # Skip self
                        nearby_atom_id = nearby_atom.get_full_id()
                        if nearby_atom_id in atom_residue_map:
                            nearby_res_key = atom_residue_map[nearby_atom_id]
                            
                            # If it's a standard residue and not already marked for extraction
                            if (nearby_res_key[0] in standard_residues and 
                                nearby_res_key != ns_res_key and 
                                nearby_res_key not in residues_to_extract):
                                
                                # Extract standard residues linked to non-standard residues
                                residues_to_extract.add(nearby_res_key)

    # Second pass: process and rename residues
    with open(ref_pdb, 'r') as f:
        current_res = None
        current_res_atoms = []
        
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                resname = line[17:20].strip()
                chain = line[21:22].strip()
                resnum = int(line[22:26])
                res_key = (resname, chain, resnum)
                
                # New residue encountered
                if res_key != current_res:
                    # Process previous residue
                    if current_res is not None and current_res_atoms:
                        if current_res in residues_to_extract:
                            if current_res[0] in standard_residues:
                                # Use the first two letters of original name as base
                                base_name = current_res[0][:2].upper()
                                new_name = generate_unique_residue_name(base_name, existing_residues, residue_type_count)
                                residue_name_mapping[current_res] = new_name
                                updated_atoms = [atom_line[:17] + f"{new_name:<3}" + atom_line[20:] for atom_line in current_res_atoms]
                                standard_res_atoms[current_res] = updated_atoms
                            else:
                                nonstandard_res_atoms.extend(current_res_atoms)
                    
                    # Reset for new residue
                    current_res = res_key
                    current_res_atoms = []
                current_res_atoms.append(line)
        
        # Process last residue in file
        if current_res is not None and current_res_atoms:
            if current_res in residues_to_extract:
                if current_res[0] in standard_residues:
                    base_name = current_res[0][:2].upper()
                    new_name = generate_unique_residue_name(base_name, existing_residues, residue_type_count)
                    residue_name_mapping[current_res] = new_name
                    updated_atoms = [atom_line[:17] + f"{new_name:<3}" + atom_line[20:] for atom_line in current_res_atoms]
                    standard_res_atoms[current_res] = updated_atoms
                else:
                    nonstandard_res_atoms.extend(current_res_atoms)

    # Write non-standard residues to output_pdb and output_xyz
    if nonstandard_res_atoms:
        with open(output_pdb, 'w') as f:
            for line in nonstandard_res_atoms:
                f.write(line)
        
        xyz_atoms = []
        for line in nonstandard_res_atoms:
            atom_name = line[12:16].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            element = line.split()[-1]
            if not element:
                alpha_chars = ''.join(c for c in atom_name if c.isalpha())
                element = (alpha_chars[:2].upper() if len(alpha_chars) >= 2 else alpha_chars[0].upper()) if alpha_chars else '' 
            element = element_converter.get(element.upper(), element)
            xyz_atoms.append((element, x, y, z))
        
        with open(output_xyz, 'w') as f:
            f.write(f"{len(xyz_atoms)}\n")
            f.write(f"Generated from {ref_pdb} - Non-standard residues\n")
            for element, x, y, z in xyz_atoms:
                f.write(f"{element:<2} {x:>10.6f} {y:>10.6f} {z:>10.6f}\n")

    # Write each renamed standard residue to its own PDB and XYZ file
    for res_key, atoms in standard_res_atoms.items():
        new_name = residue_name_mapping[res_key]
        pdb_filename = f"{new_name}.pdb"
        with open(pdb_filename, 'w') as f:
            for line in atoms:
                f.write(line)
        
        xyz_atoms = []
        for line in atoms:
            atom_name = line[12:16].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            element = line.split()[-1]
            if not element:
                alpha_chars = ''.join(c for c in atom_name if c.isalpha())
                element = (alpha_chars[:2].upper() if len(alpha_chars) >= 2 else alpha_chars[0].upper()) if alpha_chars else '' 
            element = element_converter.get(element.upper(), element)
            xyz_atoms.append((element, x, y, z))
        
        xyz_filename = f"{new_name}.xyz"
        with open(xyz_filename, 'w') as f:
            f.write(f"{len(xyz_atoms)}\n")
            f.write(f"Generated from {ref_pdb} - {new_name} {res_key[1]}:{res_key[2]}\n")
            for element, x, y, z in xyz_atoms:
                f.write(f"{element:<2} {x:>10.6f} {y:>10.6f} {z:>10.6f}\n")
    
    return residue_name_mapping

#Generate a PDB file containing all residues except those in part_QM.pdb,
#with updated residue names for metal-coordinating residues.    
def generate_nonstand_pdb(input_pdb, part_qm_pdb, output_pdb="nonstand.pdb"):
    """
    Generate the main protein PDB for tleap with metal region removed and residues renamed.
    
    Creates nonstand.pdb by:
    1. Excluding all residues that appear in part_QM.pdb (the metal+non-standard part)
    2. Renaming metal-coordinating standard residues to their unique names (HI1, CY1, etc.)
       by scanning for 3-char-named PDB files in the working directory
    
    The output is used as input for tleap to build the full protein topology, where
    the renamed residues are loaded as custom residues with QM-derived charges.
    
    Args:
        input_pdb (str): Path to metalloprotein_easyPARM.pdb (full protein).
        part_qm_pdb (str): Path to part_QM.pdb (residues to exclude).
        output_pdb (str): Output path (default: nonstand.pdb).
    """
    try:
        # First get the mapping from HI*.pdb, CY*.pdb files to know which residues to rename
        rename_mapping = {}
        for filename in os.listdir():
            if (filename.endswith('.pdb') and 
                len(filename) == 7):  # e.g., "HI1.pdb"
                    if filename[2].isdigit():  # CY1
                        new_name = filename[:3]
                    elif filename[0].isdigit():  # 1CY
                        new_name = filename[:3]
                    elif filename[1].isdigit():  # C1Y
                        new_name = filename[:3]
                    else:
                        # If none of the patterns match, skip this file
                        continue

                # Read the first ATOM/HETATM line to get residue info
                    try:
                        with open(filename, 'r') as f:
                            for line in f:
                                if line.startswith(("ATOM", "HETATM")):
                                    resnum = int(line[22:26])
                                    # Store mapping by resnum
                                    rename_mapping[resnum] = new_name
                                    break
                    except Exception as e:
                        print(f"Warning: Could not read {filename}: {e}")
                        continue 
        # Read residues to exclude from part_QM.pdb
        exclude_residues = set()
        with open(part_qm_pdb, 'r') as qm_file:
            for line in qm_file:
                if line.startswith(("ATOM", "HETATM")):
                    resnum = int(line[22:26])
                    exclude_residues.add(resnum)
        
        # Process input PDB and write output
        current_resnum = None
        with open(input_pdb, 'r') as infile, open(output_pdb, 'w') as outfile:
            for line in infile:
                if line.startswith("END"):
                    continue
                    
                if line.startswith(("ATOM", "HETATM")):
                    resnum = int(line[22:26])
                    
                    # Skip if residue should be excluded
                    if resnum in exclude_residues:
                        continue
                    
                    # Check if this residue needs to be renamed
                    if resnum in rename_mapping:
                        new_name = rename_mapping[resnum]
                        # Replace residue name while keeping the rest of the line intact
                        line = line[:17] + f"{new_name:<3}" + line[20:]
                    
                    outfile.write(line)
                
                elif not line.startswith(("ATOM", "HETATM")):
                    outfile.write(line)
            
            outfile.write("END\n")
            
    except Exception as e:
        print(f"Error generating nonstand.pdb: {e}")
        raise

#Extract charges and atom types for atoms in PDB files from MOL2 file.    
def extract_qm_charges(ref_pdb, mol2_file, residue_name_mapping=None):
    """
    Extract charges and atom types from the MOL2 file for each split-out residue.
    
    For each PDB file (part_QM.pdb and each renamed standard residue like HI1.pdb),
    looks up the corresponding atom IDs in NEW_COMPLEX.mol2 and extracts their charges
    and atom types. Writes:
    - charge_qm.dat: charges for the metal+non-standard part
    - charge_HI1.dat, charge_CY1.dat, etc.: charges for each renamed standard residue
    - charges_all.dat: index mapping mol2 files to charge files
    - qm.xyz: coordinates of the metal+non-standard part (for Seminario method)
    
    Args:
        ref_pdb (str): Path to REFQM.pdb (not directly used, but defines atom ordering context).
        mol2_file (str): Path to NEW_COMPLEX.mol2 (source of charges and coordinates).
        residue_name_mapping (dict): {(resname, chain, resnum): new_name} from extract_non_standard_residues_from_ref().
    """
    if residue_name_mapping is None:
        residue_name_mapping = {}
    
    try:
        # First, identify all PDB files to process
        pdb_files = []
         
        # Add part_QM.pdb first (for non-standard residues) - with XYZ file
        if os.path.exists("part_QM.pdb"):
            pdb_files.append(("part_QM.pdb", "charge_qm.dat", "qm.xyz", True))
        
        # Add all renamed standard residue PDB files - WITHOUT XYZ files
        for key, new_name in residue_name_mapping.items():
            pdb_filename = f"{new_name}.pdb"
            if os.path.exists(pdb_filename):
                charge_filename = f"charge_{new_name}.dat"
                pdb_files.append((pdb_filename, charge_filename, None, False))
        
        # Open charges_all.dat for writing the mapping
        with open("charges_all.dat", 'w') as charges_all_file:
            # First, handle part_QM.pdb for non-standard residues
            if os.path.exists("part_QM.pdb"):
                charges_all_file.write("QM.mol2 charge_qm.dat\n")
            
            # Write mappings for renamed standard residues
            for key, new_name in residue_name_mapping.items():
                pdb_filename = f"{new_name}.pdb"
                if os.path.exists(pdb_filename):
                    charges_all_file.write(f"{new_name}.mol2 charge_{new_name}.dat\n")
        
        # Process each PDB file
        for pdb_info in pdb_files:
            pdb_file, charge_output, xyz_output, create_xyz = pdb_info
            
            # Read atom IDs, atom names, and residue information from PDB
            qm_atom_ids = []
            atom_residue_mapping = {}
            atom_pdb_names = {}
            
            with open(pdb_file, 'r') as pdb:
                for line in pdb:
                    if line.startswith(("ATOM", "HETATM")):
                        try:
                            atom_id = int(line[6:11].strip())
                            pdb_atom_name = line[12:16].strip()
                            resname = line[17:20].strip()
                            chain = line[21:22].strip()
                            resnum = int(line[22:26])
                            
                            # Check if this residue was renamed
                            original_key = (resname, chain, resnum)
                            new_resname = residue_name_mapping.get(original_key, resname)
                            
                            qm_atom_ids.append(atom_id)
                            atom_residue_mapping[atom_id] = new_resname
                            atom_pdb_names[atom_id] = pdb_atom_name
                        except ValueError as e:
                            print(f"Warning: Could not parse atom ID from PDB line: {line.strip()}")
                            continue
            
            # Read charges, atom types, atom names, and coordinates from MOL2 file
            charges_and_types = {}
            atom_data = {}
            is_atom_section = False
            
            with open(mol2_file, 'r') as mol2:
                for line in mol2:
                    if "@<TRIPOS>ATOM" in line:
                        is_atom_section = True
                        continue
                    elif "@<TRIPOS>" in line and "ATOM" not in line:
                        is_atom_section = False
                        continue
                    
                    if is_atom_section and line.strip():
                        try:
                            parts = line.split()
                            atom_id = int(parts[0])
                            if atom_id in qm_atom_ids:
                                atom_name = parts[1]
                                x = float(parts[2])
                                y = float(parts[3])
                                z = float(parts[4])
                                atom_type = parts[5]
                                charge = float(parts[-1])
                                
                                charges_and_types[atom_id] = {
                                    'charge': charge, 
                                    'atom_type': atom_type
                                }
                                
                                # Extract element by removing numbers from atom name
                                element = ''.join([c for c in atom_name if not c.isdigit()])
                                
                                atom_data[atom_id] = {
                                    'element': element,
                                    'x': x,
                                    'y': y,
                                    'z': z
                                }
                        except (ValueError, IndexError) as e:
                            print(f"Warning: Could not parse MOL2 line: {line.strip()}")
                            continue
            
            # Write charges to output file (charge, atom_type, PDB_atom_name)
            with open(charge_output, 'w') as charge_file:
                for atom_id in qm_atom_ids:
                    if atom_id in charges_and_types:
                        data = charges_and_types[atom_id]
                        pdb_name = atom_pdb_names.get(atom_id, "???")
                        charge_file.write(f"{data['charge']:.6f} {data['atom_type']} {pdb_name}\n")
                    else:
                        print(f"Warning: No charge found for atom ID {atom_id}")
            
            # Write XYZ coordinate file only for QM residue
            if create_xyz and xyz_output:
                with open(xyz_output, 'w') as xyz_file:
                    # Line 1: number of atoms
                    xyz_file.write(f"{len(qm_atom_ids)}\n")
                    # Line 2: blank line
                    xyz_file.write("\n")
                    # Lines 3+: element name (without numbers) and xyz coordinates
                    for atom_id in qm_atom_ids:
                        if atom_id in atom_data:
                            data = atom_data[atom_id]
                            xyz_file.write(f"{data['element']:2s} {data['x']:12.6f} {data['y']:12.6f} {data['z']:12.6f}\n")
                        else:
                            print(f"Warning: No coordinate data found for atom ID {atom_id}")
                 
    except Exception as e:
        print(f"Error extracting charges: {e}")
        raise

#Generate easyPARM_residues.dat file containing information about renamed standard residues.
#Format: standard_residue.pdb standard_residue.mol2 standard_residue_name
#Only includes residues that were originally standard residues (like HID, CYS, etc.)
def generate_easyparm_residues(residue_name_mapping, output_file="easyPARM_residues.dat", output_file2="ALL_RESIDUE_tleap.input"):

    """
    Generate configuration files listing renamed standard residues for downstream processing.
    
    Creates two files:
    1. easyPARM_residues.dat — lists each renamed residue in format:
       "HI1.pdb HI1.mol2 HI1" (used by the bash script to run antechamber on each)
    2. ALL_RESIDUE_tleap.input — tleap script to load all renamed residue mol2 files
       and save them to COMPLEX.lib
    
    Only includes residues that were originally standard amino acids (HID, CYS, ASP, etc.).
    
    Args:
        residue_name_mapping (dict): {(resname, chain, resnum): new_name} mapping.
        output_file (str): Path for easyPARM_residues.dat output.
        output_file2 (str): Path for ALL_RESIDUE_tleap.input output.
    """
    try:
        # Define standard residues
        standard_residues = {
            "ALA", "ARG", "ASH", "ASN", "ASP", "CYM", "CYS", "CYX", "GLH", "GLN",
            "GLU", "GLY", "HID", "HIE", "HIP", "HYP", "ILE", "LEU", "LYN", "LYS",
            "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", 'HIS'
        }

        # List to hold the special residues
        special_residues = []

        # Go through the residue name mapping
        for original_res, new_name in residue_name_mapping.items():
            resname = original_res[0]

            # Only include if the original residue is standard
            if resname in standard_residues:
                pdb_file = f"{new_name}.pdb"
                mol2_file = f"{new_name}.mol2"
                residue_name = new_name

                special_residues.append((pdb_file, mol2_file, residue_name))

        # Sort the list for consistent output
        special_residues.sort()

        # Write to output file
        with open(output_file, 'w') as f:
            for pdb, mol2, name in special_residues:
                f.write(f"{pdb} {mol2} {name}\n")
        
        with open(output_file2, 'w') as f:
            f.write(f"source leaprc.gaff\n")
            for pdb, mol2, name in special_residues:
                f.write(f"loadamberparams COMPLEX.frcmod\n")
                f.write(f"{name} = loadmol2 {mol2}\n")
                f.write(f"saveoff {name} COMPLEX.lib\n")
            f.write(f"quit\n")


    except Exception as e:
        print(f"Error generating {output_file}: {e}")
        raise

#Generate a PDB file containing all residues except those in part_QM.pdb
def generate_easyPARM_nonstand_pdb(input_pdb, part_qm_pdb, output_pdb="easynonstands.pdb"):
    """
    Generate an alternative nonstand PDB using coordinate-based exclusion.
    
    Similar to generate_nonstand_pdb() but excludes atoms by matching their
    XYZ coordinates (rounded to 3 decimals) against part_QM.pdb, rather than
    by residue number. This is more robust when atom numbering differs between
    the protein PDB and the QM extraction.
    
    Output (easynonstands.pdb) is used later when assembling the final system
    with pdb4amber.
    
    Args:
        input_pdb (str): Path to metalloprotein_easyPARM.pdb (full protein).
        part_qm_pdb (str): Path to part_QM.pdb (atoms to exclude by coordinates).
        output_pdb (str): Output path (default: easynonstands.pdb).
    """
    try:
        # Build exclude set from coordinates in part_qm_pdb
        # Coordinates are the only reliable common field between the two file formats
        exclude_coords = set()
        with open(part_qm_pdb, 'r') as qm_file:
            for line in qm_file:
                if line.startswith(("ATOM", "HETATM")):
                    try:
                        x = round(float(line[30:38].strip()), 3)
                        y = round(float(line[38:46].strip()), 3)
                        z = round(float(line[46:54].strip()), 3)
                        exclude_coords.add((x, y, z))
                    except ValueError:
                        continue  # skip malformed lines

        with open(input_pdb, 'r') as infile, open(output_pdb, 'w') as outfile:
            for line in infile:
                if line.startswith("END"):
                    continue

                if line.startswith(("ATOM", "HETATM")):
                    try:
                        x = round(float(line[30:38].strip()), 3)
                        y = round(float(line[38:46].strip()), 3)
                        z = round(float(line[46:54].strip()), 3)
                        coords = (x, y, z)
                    except ValueError:
                        outfile.write(line)  # keep malformed lines as-is
                        continue

                    if coords not in exclude_coords:
                        outfile.write(line)

                else:
                    outfile.write(line)

    except Exception as e:
        print(f"Error generating nonstand.pdb: {e}")
        raise

def main():
    """
    Main entry point for metalloprotein.py.
    
    Orchestrates the full metalloprotein processing workflow:
    1. analyze_and_extract_metal_site() — identify metals, extract coordination shell
    2. update_mol2_with_qm_types() — create COMPLEX_updated.mol2 with hybrid atom types
    3. extract_non_standard_residues_from_ref() — split into metal part + standard residues
    4. generate_nonstand_pdb() — create protein PDB with metal region removed
    5. extract_qm_charges() — extract per-residue charges from MOL2
    6. generate_easyparm_residues() — write config files for bash script
    7. generate_easyPARM_nonstand_pdb() — coordinate-based exclusion PDB
    
    Command-line usage:
        python3 metalloprotein.py <path_to_ff_library.lib>
    """
    try:
        # Check if library file argument is provided
        if len(sys.argv) < 2:
            print("Example: python3 code.py amber_fb15.lib")
            sys.exit(1)

        # Get library file from command line argument
        lib_file = sys.argv[1]

        # Check if library file exists
        if not os.path.exists(lib_file):
            print(f"Error: Library file '{lib_file}' not found")
            sys.exit(1)
        input_pdb = "metalloprotein_easyPARM.pdb"
        mol2_file = "NEW_COMPLEX.mol2"
        output_mol2 = "COMPLEX_updated.mol2"
        part_qm_pdb = "part_QM.pdb"
        part_qm_xyz = "part_QM.xyz"
        qm_pdb = "QM.pdb"
        nonstand_pdb = "nonstand.pdb"
        charge_file = "charge_qm.dat"
        easynonstand_pdb = "easynonstands.pdb"

        # Analyze metal site and generate files
        metal_coordination = analyze_and_extract_metal_site(input_pdb, mol2_file, lib_file, qm_pdb)

        # Update MOL2 file with QM.pdb atom types
        update_mol2_with_qm_types(mol2_file, qm_pdb, output_mol2)
        
        # Extract non-standard residues
        residue_name_mapping = extract_non_standard_residues_from_ref("REFQM.pdb", part_qm_pdb, part_qm_xyz)
                                                                      
        # Generate nonstand.pdb excluding the metal and ligand residues
        generate_nonstand_pdb(input_pdb, part_qm_pdb, nonstand_pdb)
 
        # Extract charges for QM region
        extract_qm_charges("REFQM.pdb", "NEW_COMPLEX.mol2", residue_name_mapping) 
        # After running extract_non_standard_residues_from_ref and generate_nonstand_pdb
        generate_easyparm_residues(residue_name_mapping)
        generate_easyPARM_nonstand_pdb(input_pdb, part_qm_pdb, easynonstand_pdb)
    except FileNotFoundError as e:
        print(f"Error: File not found - {e}")
    except Exception as e:
        print(f"Error processing files: {e}")
        raise

if __name__ == "__main__":
    main()
