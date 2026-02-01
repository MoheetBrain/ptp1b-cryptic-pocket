import pandas as pd
from chembl_webresource_client.new_client import new_client
from rdkit import Chem
from rdkit.Chem import Descriptors

class DeepDiveAgent:
    def __init__(self):
        print("🤿 [DEEP DIVE] Analyzing the True Allosteric hits...")
        self.molecule = new_client.molecule
        self.activity = new_client.activity

    def analyze_candidate(self, chembl_id):
        # 1. Get Structure
        mol_data = self.molecule.get(chembl_id)
        smiles = mol_data['molecule_structures']['canonical_smiles']
        
        # 2. Get Potency (IC50)
        acts = self.activity.filter(
            molecule_chembl_id=chembl_id, 
            target_chembl_id="CHEMBL335",
            standard_type="IC50"
        )
        
        # Find the best IC50
        best_ic50 = float('inf')
        for a in acts:
            try:
                val = float(a['standard_value'])
                if val < best_ic50:
                    best_ic50 = val
            except:
                continue

        # 3. Calculate Properties
        mol = Chem.MolFromSmiles(smiles)
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)

        return {
            'chembl_id': chembl_id,
            'ic50_nm': best_ic50,
            'mw': round(mw, 2),
            'logp': round(logp, 2),
            'smiles': smiles
        }

    def execute(self):
        # The Two Unicorns found by the Pivot
        targets = ['CHEMBL598832', 'CHEMBL4128743']
        
        results = []
        for t in targets:
            print(f"   ...Processing {t}")
            data = self.analyze_candidate(t)
            results.append(data)
            
        # Report
        print("\n🏆 --- THE TRUE ALLOSTERIC SCAFFOLDS ---")
        df = pd.DataFrame(results)
        print(df[['chembl_id', 'ic50_nm', 'mw', 'logp']].to_string(index=False))
        
        # Save for visualization
        df.to_csv("data/processed/ptp1b_true_gold.csv", index=False)

if __name__ == "__main__":
    bot = DeepDiveAgent()
    bot.execute()