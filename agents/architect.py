import pandas as pd
import os
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold   # <--- THE CORRECT IMPORT
from rdkit.Chem import Descriptors, Lipinski
from rdkit.Chem import FilterCatalog 

# --- CONFIG ---
INPUT_FILE = "data/processed/ptp1b_high_confidence.csv"
OUTPUT_FILE = "data/processed/ptp1b_architect_designs.csv"

class ArchitectAgent:
    def __init__(self):
        print("🏛️ [ARCHITECT] Initializing Generative Design Logic...")
        # Initialize PAINS filter
        params = FilterCatalog.FilterCatalogParams()
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
        self.pains_filter = FilterCatalog.FilterCatalog(params)

    def analyze_structure(self, smile):
        mol = Chem.MolFromSmiles(smile)
        if not mol: return None, False
        if self.pains_filter.HasMatch(mol):
            return mol, False 
        return mol, True

    def get_scaffold(self, mol):
        # <--- YOUR FIX IS HERE
        core = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(core)

    def execute_design_phase(self):
        print(f"📂 [ARCHITECT] Loading mining data from {INPUT_FILE}...")
        
        if not os.path.exists(INPUT_FILE):
            print("❌ [ERROR] Input file missing.")
            return

        df = pd.read_csv(INPUT_FILE)
        
        # Check if chembl_id exists in input
        if 'chembl_id' not in df.columns:
            print("❌ [ERROR] Input CSV is missing 'chembl_id' column.")
            return

        unique_scaffolds = {}

        print("📐 [ARCHITECT] Extracting scaffolds and filtering PAINS...")
        
        for index, row in df.iterrows():
            smile = row['smiles']
            ic50 = row['ic50_nm']
            cid = row['chembl_id']
            
            mol, is_clean = self.analyze_structure(smile)
            if not is_clean: continue 

            scaffold = self.get_scaffold(mol)
            
            # Logic: Keep the BEST binder for each scaffold
            if scaffold not in unique_scaffolds:
                unique_scaffolds[scaffold] = {
                    'chembl_id': cid,
                    'smiles': smile,
                    'scaffold': scaffold,
                    'ic50_nm': ic50,
                    'mw': Descriptors.MolWt(mol),
                    'logp': Descriptors.MolLogP(mol),
                    'hbd': Lipinski.NumHDonors(mol)
                }
            else:
                if ic50 < unique_scaffolds[scaffold]['ic50_nm']:
                     unique_scaffolds[scaffold].update({
                        'chembl_id': cid,
                        'smiles': smile,
                        'ic50_nm': ic50,
                        'mw': Descriptors.MolWt(mol),
                        'logp': Descriptors.MolLogP(mol)
                     })

        architect_picks = list(unique_scaffolds.values())
        
        # Filter for Cryptic Pocket (Hydrophobic + Small)
        final_designs = []
        for pick in architect_picks:
            if 2.0 <= pick['logp'] <= 5.0 and pick['mw'] < 550:
                final_designs.append(pick)

        results_df = pd.DataFrame(final_designs)
        results_df.sort_values(by='ic50_nm', inplace=True)
        
        results_df.to_csv(OUTPUT_FILE, index=False)
        print(f"\n✅ [ARCHITECT] Fixed Data Saved to: {OUTPUT_FILE}")
        print(f"   ...Count: {len(results_df)}")

if __name__ == "__main__":
    bot = ArchitectAgent()
    bot.execute_design_phase() 