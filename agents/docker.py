import os
import subprocess
from vina import Vina
from meeko import MoleculePreparation
from rdkit import Chem

# --- CONFIG ---
PDB_ID = "6B95"
LIGAND_SDF = "data/3d_structures/benzofuran_series.sdf"
PROTEIN_PDB = f"data/protein/{PDB_ID}.pdb"
OUTPUT_DIR = "data/docking_results"

# Coordinates from Targeter (Using your exact values)
CENTER_X = 135.717
CENTER_Y = 11.911
CENTER_Z = -6.150
SIZE_X = 24.0
SIZE_Y = 48.0
SIZE_Z = 20.0

class DockerAgent:
    def __init__(self):
        print("⚔️ [DOCKER] Initializing Physics Simulation...")
        if not os.path.exists(OUTPUT_DIR):
            os.makedirs(OUTPUT_DIR)

    def prepare_receptor(self):
        print("   ...Preparing Protein Receptor (adding charges)...")
        output_pdbqt = f"{OUTPUT_DIR}/receptor.pdbqt"
        
        # Use OpenBabel to convert PDB -> PDBQT
        # -xr: output rigid molecule
        cmd = [
            "obabel", PROTEIN_PDB, "-O", output_pdbqt, 
            "--partialcharge", "gasteiger", "-xr"
        ]
        
        try:
            subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            print("   ✅ Receptor Ready.")
            return output_pdbqt
        except Exception as e:
            print(f"❌ OpenBabel Failed. Error: {e}")
            return None

    def prepare_ligand(self):
        print("   ...Preparing Benzofuran Ligand...")
        suppl = Chem.SDMolSupplier(LIGAND_SDF, removeHs=False)
        mol = next(iter(suppl)) # Take the first molecule (best sibling)
        
        # Use Meeko to prep ligand
        preparator = MoleculePreparation()
        preparator.prepare(mol)
        pdbqt_string = preparator.write_pdbqt_string()
        
        output_pdbqt = f"{OUTPUT_DIR}/ligand.pdbqt"
        with open(output_pdbqt, "w") as f:
            f.write(pdbqt_string)
            
        print("   ✅ Ligand Ready.")
        return output_pdbqt

    def run_simulation(self, receptor_path, ligand_path):
        print(f"   ...Starting Vina Docking @ ({CENTER_X}, {CENTER_Y}, {CENTER_Z})...")
        
        v = Vina(sf_name='vina')
        v.set_receptor(receptor_path)
        v.set_ligand_from_file(ligand_path)
        
        # Set the Box
        v.compute_vina_maps(center=[CENTER_X, CENTER_Y, CENTER_Z], 
                            box_size=[SIZE_X, SIZE_Y, SIZE_Z])
        
        # Dock!
        v.dock(exhaustiveness=8, n_poses=5)
        
        # Save Results
        out_pdbqt = f"{OUTPUT_DIR}/docked_poses.pdbqt"
        v.write_poses(out_pdbqt, n_poses=5, overwrite=True)
        
        print("\n🏆 [DOCKING COMPLETE]")
        print(f"   Saved poses to: {out_pdbqt}")
        print("-" * 30)
        print("   BINDING AFFINITIES (kcal/mol):")
        
        energies = v.energies(n_poses=5)
        for i, energy in enumerate(energies):
            print(f"   Pose {i+1}: {energy[0]:.2f} kcal/mol")

if __name__ == "__main__":
    bot = DockerAgent()
    rec = bot.prepare_receptor()
    lig = bot.prepare_ligand()
    
    if rec and lig:
        bot.run_simulation(rec, lig)