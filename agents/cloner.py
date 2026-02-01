import pandas as pd
from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
from rdkit.Chem.Scaffolds import MurckoScaffold

# --- CONFIG ---
TARGET_ID = "CHEMBL598832"    # The "Unicorn"
SIMILARITY_THRESHOLD = 70     # 70% match

class ClonerAgent:
    def __init__(self):
        print(f"🧬 [CLONER] Initializing. Target: {TARGET_ID}")
        self.similarity = new_client.similarity
        self.molecule = new_client.molecule
        self.activity = new_client.activity

    def get_lead_structure(self):
        try:
            mol = self.molecule.get(TARGET_ID)
            return mol['molecule_structures']['canonical_smiles']
        except Exception as e:
            return None

    def expand_series(self):
        print(f"   ...Fetching structure for {TARGET_ID}...")
        lead_smiles = self.get_lead_structure()
        
        if not lead_smiles: return

        print(f"   ...Scanning ChEMBL for siblings >{SIMILARITY_THRESHOLD}% similar...")
        try:
            res = self.similarity.filter(smiles=lead_smiles, similarity=SIMILARITY_THRESHOLD)
            print(f"   ...Found {len(res)} structural siblings.")
        except: return

        candidates = []
        print("   ...Filtering for PTP1B activity...")
        
        for match in res:
            cid = match['molecule_chembl_id']
            similarity = match['similarity']
            
            acts = self.activity.filter(molecule_chembl_id=cid, target_chembl_id="CHEMBL335", standard_type="IC50")
            if not acts: continue

            best_ic50 = float('inf')
            for a in acts:
                try:
                    val = float(a['standard_value'])
                    if val < best_ic50: best_ic50 = val
                except: continue
            if best_ic50 == float('inf'): continue

            # Properties & Scaffold (CRITICAL FOR DASHBOARD)
            smiles = match['molecule_structures']['canonical_smiles']
            mol = Chem.MolFromSmiles(smiles)
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            hbd = Lipinski.NumHDonors(mol)
            
            try:
                core = MurckoScaffold.GetScaffoldForMol(mol)
                scaffold = Chem.MolToSmiles(core)
            except: scaffold = smiles

            candidates.append({
                'chembl_id': cid,
                'similarity': similarity,
                'ic50_nm': best_ic50,
                'mw': round(mw, 2),
                'logp': round(logp, 2),
                'hbd': hbd,
                'scaffold': scaffold,
                'smiles': smiles
            })
            print(f"   ✅ Sibling: {cid} ({similarity}%) -> {best_ic50} nM")

        if candidates:
            df = pd.DataFrame(candidates)
            df.sort_values(by='ic50_nm', inplace=True)
            
            # Save to the file the Dashboard reads
            output_file = "data/processed/ptp1b_architect_designs.csv" 
            df.to_csv(output_file, index=False)
            
            print(f"\n🏆 [SUCCESS] Family Tree Saved for Dashboard.")
            print(f"   ...Best Sibling: {df.iloc[0]['chembl_id']} ({df.iloc[0]['ic50_nm']} nM)")
        else:
            print("❌ No active siblings found.")

if __name__ == "__main__":
    bot = ClonerAgent()
    bot.expand_series()