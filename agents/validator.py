import pandas as pd
from chembl_webresource_client.new_client import new_client

# --- CONFIG ---
INPUT_FILE = "data/processed/ptp1b_architect_designs.csv"

class ValidatorAgent:
    def __init__(self):
        print("⚖️ [VALIDATOR] Initializing Mechanism Check...")
        self.mechanism = new_client.mechanism

    def check_mechanism(self, chembl_id):
        """Asks ChEMBL: 'How does this kill the enzyme?'"""
        try:
            # Fetch mechanism records
            mechs = self.mechanism.filter(molecule_chembl_id=chembl_id)
            
            if not mechs:
                return "❓ Unknown / Not Annotated"
            
            # Look for keywords in the mechanism description
            details = []
            for m in mechs:
                desc = m.get('mechanism_of_action', 'N/A')
                site = m.get('binding_site_comment', 'N/A')
                details.append(f"{desc} | Site: {site}")
                
            return " | ".join(details)

        except Exception as e:
            return f"⚠️ Error: {e}"

    def execute_validation(self):
        print(f"📂 [VALIDATOR] Reading candidates from {INPUT_FILE}...")
        
        try:
            df = pd.read_csv(INPUT_FILE)
        except FileNotFoundError:
            print("❌ File not found.")
            return

        # Let's check the Top 5
        print("\n🔍 --- MECHANISM AUDIT ---")
        print(f"{'CHEMBL ID':<20} {'IC50 (nM)':<15} {'MECHANISM VERDICT'}")
        print("-" * 80)
        
        for index, row in df.head(5).iterrows():
            cid = row['chembl_id']
            ic50 = row['ic50_nm']
            
            # The Verdict
            verdict = self.check_mechanism(cid)
            print(f"{cid:<20} {ic50:<15} {verdict}")

if __name__ == "__main__":
    bot = ValidatorAgent()
    bot.execute_validation()