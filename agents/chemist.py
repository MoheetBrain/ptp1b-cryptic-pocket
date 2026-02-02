import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, BRICS, Lipinski

# --- CONFIG ---
INPUT_FILE = "data/processed/ptp1b_allosteric_family.csv" 
OUTPUT_FILE = "data/processed/ptp1b_next_gen_leads.csv"

class ChemistAgent:
    def __init__(self):
        print("⚗️ [CHEMIST] Initializing Sledgehammer Protocol (Fixed)...")

    def fragment_monster(self, mol, parent_id):
        print(f"   🪓 Chopping oversized molecule {parent_id}...")
        try:
            frags = BRICS.BreakBRICSBonds(mol)
            frag_mols = Chem.GetMolFrags(frags, asMols=True)
            valid_cores = []
            for f in frag_mols:
                mw = Descriptors.MolWt(f)
                if 200 < mw < 500:
                    try:
                        Chem.SanitizeMol(f)
                        valid_cores.append({
                            'parent': parent_id,
                            'modification': 'BRICS_Core_Cut',
                            'smiles': Chem.MolToSmiles(f),
                            'mw': round(mw, 2),
                            'logp': round(Descriptors.MolLogP(f), 2)
                        })
                    except: continue
            return valid_cores
        except: return []

    def optimize_lead(self):
        try:
            df = pd.read_csv(INPUT_FILE)
            df.sort_values(by='ic50_nm', inplace=True)
            print(f"   ...Loaded {len(df)} candidates for optimization.")
        except:
            print("❌ Input file missing.")
            return

        # GREEDY REACTIONS: Target ANY aromatic carbon [c:1]
        reactions = [
            (AllChem.ReactionFromSmarts('[c:1]>>[c:1](F)'), "+Fluorine"),
            (AllChem.ReactionFromSmarts('[c:1]>>[c:1](Cl)'), "+Chlorine"),
            (AllChem.ReactionFromSmarts('[c:1]>>[c:1](C)'), "+Methyl"),
            (AllChem.ReactionFromSmarts('[c:1]>>[c:1](O)'), "+Hydroxyl"),
            (AllChem.ReactionFromSmarts('[c:1]>>[c:1](C#N)'), "+Nitrile"),
            (AllChem.ReactionFromSmarts('[c:1]>>[c:1](C(F)(F)F)'), "+CF3")
        ]

        analogs = []  # <--- FIXED LINE 60

        for index, row in df.iterrows():
            parent_id = row['chembl_id']
            parent_mw = row['mw']
            mol = Chem.MolFromSmiles(row['smiles'])
            if not mol: continue

            if parent_mw > 600:
                analogs.extend(self.fragment_monster(mol, parent_id))
                continue 

            print(f"   ...Reacting {parent_id}...")
            
            for rxn, tag in reactions:
                try:
                    # Run on EVERY aromatic atom
                    ps = rxn.RunReactants((mol,))
                    for p in ps:
                        try:
                            new_mol = p[0]
                            Chem.SanitizeMol(new_mol)
                            mw = Descriptors.MolWt(new_mol)
                            if mw < 600:
                                analogs.append({
                                    'parent': parent_id,
                                    'modification': tag,
                                    'smiles': Chem.MolToSmiles(new_mol),
                                    'mw': round(mw, 2),
                                    'logp': round(Descriptors.MolLogP(new_mol), 2)
                                })
                        except: continue
                except: continue

        if analogs:
            res_df = pd.DataFrame(analogs).drop_duplicates(subset=['smiles'])
            res_df.to_csv(OUTPUT_FILE, index=False)
            print(f"\n🏆 [OPTIMIZATION COMPLETE]")
            print(f"   Generated {len(res_df)} Next-Gen Analogs.")
            print(f"   Saved to: {OUTPUT_FILE}")
            print("\n   --- Modification Distribution ---")
            print(res_df['modification'].value_counts().to_string())
        else:
            print("   ❌ No valid analogs generated.")

if __name__ == "__main__":
    bot = ChemistAgent()
    bot.optimize_lead()