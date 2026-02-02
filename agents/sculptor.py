import pandas as pd
import os
from rdkit import Chem
from rdkit.Chem import AllChem

# --- CONFIG ---
INPUT_FILE = "data/processed/ptp1b_allosteric_family.csv"
OUTPUT_DIR = "data/3d_structures"

class SculptorAgent:
    def __init__(self):
        print("🗿 [SCULPTOR] Initializing 3D Conformer Generator...")
        if not os.path.exists(OUTPUT_DIR):
            os.makedirs(OUTPUT_DIR)

    def generate_3d(self):
        print(f"   ...Reading family from {INPUT_FILE}")
        try:
            df = pd.read_csv(INPUT_FILE)
        except:
            print("❌ File not found.")
            return

        print(f"   ...Sculpting {len(df)} candidates into 3D...")
        
        success_count = 0
        writer = Chem.SDWriter(f"{OUTPUT_DIR}/benzofuran_series.sdf")

        for index, row in df.iterrows():
            cid = row['chembl_id']
            smiles = row['smiles']
            
            # 1. Create Mol from SMILES
            mol = Chem.MolFromSmiles(smiles)
            mol = Chem.AddHs(mol) # Add Hydrogens (Crucial for 3D)
            
            # 2. Embed in 3D Space (Inflate the balloon)
            params = AllChem.ETKDGv3()
            res = AllChem.EmbedMolecule(mol, params)
            
            if res == 0:
                # 3. Energy Minimize (Relax the molecule)
                AllChem.MMFFOptimizeMolecule(mol)
                
                # Set Name
                mol.SetProp("_Name", cid)
                mol.SetProp("IC50_nM", str(row['ic50_nm']))
                
                # Save
                writer.write(mol)
                success_count += 1
                print(f"   ✅ Sculpted {cid} -> 3D Ready")
            else:
                print(f"   ⚠️ Failed to generate 3D for {cid}")

        writer.close()
        print(f"\n🏆 [SUCCESS] Saved {success_count} 3D structures to {OUTPUT_DIR}/benzofuran_series.sdf")
        print("   -> Next Step: Load this .sdf file into PyMOL or AutoDock Vina.")

if __name__ == "__main__":
    bot = SculptorAgent()
    bot.generate_3d()