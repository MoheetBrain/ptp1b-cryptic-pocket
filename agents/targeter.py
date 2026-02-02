import os
import urllib.request
import numpy as np

# --- CONFIG ---
PDB_ID = "6B95"  # The Crystal Structure of PTP1B with Allosteric Inhibitor
OUTPUT_DIR = "data/protein"
PDB_URL = f"https://files.rcsb.org/download/{PDB_ID}.pdb"

class TargeterAgent:
    def __init__(self):
        print(f"🎯 [TARGETER] Initializing. Target System: {PDB_ID}")
        if not os.path.exists(OUTPUT_DIR):
            os.makedirs(OUTPUT_DIR)

    def fetch_protein(self):
        # 1. Download the PDB file
        pdb_path = f"{OUTPUT_DIR}/{PDB_ID}.pdb"
        if not os.path.exists(pdb_path):
            print(f"   ...Downloading {PDB_ID} from RCSB...")
            urllib.request.urlretrieve(PDB_URL, pdb_path)
        else:
            print(f"   ...Found {PDB_ID} locally.")
        return pdb_path

    def locate_pocket(self, pdb_path):
        print("   ...Scanning protein for Allosteric Ligand (BB3)...")
        
        ligand_coords = []
        
        # 2. Parse the PDB to find the ligand atoms (HETATM)
        # In 6B95, the allosteric ligand is often labeled 'BB3' or similar.
        # We will scan for HETATM lines that are NOT water (HOH).
        
        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith("HETATM"):
                    res_name = line[17:20].strip()
                    # Exclude Water and Ions
                    if res_name not in ["HOH", "EDO", "MG", "NA"]:
                        # Extract X, Y, Z
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        ligand_coords.append([x, y, z])

        if not ligand_coords:
            print("❌ Could not auto-detect ligand. Manual coordinates required.")
            return

        # 3. Calculate Center of Mass (The Bullseye)
        coords = np.array(ligand_coords)
        center = np.mean(coords, axis=0)
        
        min_c = np.min(coords, axis=0)
        max_c = np.max(coords, axis=0)
        size = max_c - min_c + 10.0 # Add 10 Angstrom buffer

        print("\n✅ [TARGET ACQUIRED] Cryptic Pocket Coordinates Found:")
        print(f"   Center X: {center[0]:.3f}")
        print(f"   Center Y: {center[1]:.3f}")
        print(f"   Center Z: {center[2]:.3f}")
        print("-" * 30)
        print("📦 [GRID BOX CONFIGURATION]")
        print(f"   size_x = {size[0]:.1f}")
        print(f"   size_y = {size[1]:.1f}")
        print(f"   size_z = {size[2]:.1f}")
        print(f"   center_x = {center[0]:.3f}")
        print(f"   center_y = {center[1]:.3f}")
        print(f"   center_z = {center[2]:.3f}")
        
        # Save these config params for later
        with open(f"{OUTPUT_DIR}/docking_config.txt", "w") as f:
            f.write(f"center_x = {center[0]:.3f}\n")
            f.write(f"center_y = {center[1]:.3f}\n")
            f.write(f"center_z = {center[2]:.3f}\n")
            f.write(f"size_x = {size[0]:.1f}\n")
            f.write(f"size_y = {size[1]:.1f}\n")
            f.write(f"size_z = {size[2]:.1f}\n")

if __name__ == "__main__":
    bot = TargeterAgent()
    pdb = bot.fetch_protein()
    bot.locate_pocket(pdb)