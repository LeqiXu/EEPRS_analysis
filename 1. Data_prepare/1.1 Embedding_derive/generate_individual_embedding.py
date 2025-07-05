import pandas as pd
import numpy as np
import pickle
import os

def load_eids(eid_path, eid_col=0, sep='\t'):
    """
    Load eids from a text/csv file.
    """
    eids = pd.read_csv(eid_path, sep=sep, header=None)
    return eids.iloc[:, eid_col].tolist()

def load_embeddings(file_path, from_pkl=True):
    """
    Load ICD-10 code embeddings from a .pkl or .txt/.csv file.
    Returns dict: {code: vector}
    """
    if from_pkl:
        with open(file_path, 'rb') as file:
            return pickle.load(file)
    else:
        df = pd.read_csv(file_path, sep=None, engine='python', index_col=0)
        return {index: row.values for index, row in df.iterrows()}

def get_ind_embed(eid, icd_embed, hesin_df, icd_col='diag_icd10'):
    icds = hesin_df[hesin_df['eid'] == eid][icd_col]
    vectors = [icd_embed[icd] for icd in icds if icd in icd_embed]
    if vectors:
        return np.mean(np.vstack(vectors), axis=0)
    else:
        # Return zero vector of correct dimension
        return np.zeros(len(next(iter(icd_embed.values()))))

def get_group_embed(eid_list, hesin_df, icd_embed, icd_col='diag_icd10'):
    group_embed = {}
    for i, eid in enumerate(eid_list):
        group_embed[eid] = get_ind_embed(eid, icd_embed, hesin_df, icd_col)
        if (i+1) % 1000 == 0:
            print(f"{i+1} individuals embedded.")
    return group_embed

def main(
    hesin_path,
    eid_path,
    icd_embedding_path,
    output_prefix,
    eid_col=0,
    sep_eid='\t',
    icd_embedding_from_pkl=True,
    icd_col='diag_icd10'
):
    print("Loading hesin diagnosis file...")
    hesin = pd.read_csv(hesin_path, sep='\t')
    hesin = hesin[['eid', icd_col]].drop_duplicates()

    print(f"Loading EIDs from {eid_path} ...")
    eids = load_eids(eid_path, eid_col=eid_col, sep=sep_eid)

    print(f"Filtering hesin for selected EIDs (N={len(eids)}) ...")
    hesin_subset = hesin[hesin['eid'].isin(eids)]

    print("Loading ICD-10 code embeddings ...")
    icd_embed = load_embeddings(icd_embedding_path, from_pkl=icd_embedding_from_pkl)

    print("Generating individual embeddings ...")
    group_embed = get_group_embed(eids, hesin_subset, icd_embed, icd_col=icd_col)

    save_path = f"{output_prefix}_embed.pkl"
    with open(save_path, "wb") as f:
        pickle.dump(group_embed, f)
    print(f"Saved: {save_path}")

# =========================
# Example usage as a script
# =========================
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Generate individual phenotype embeddings from hesin + ICD-10 embeddings.")
    parser.add_argument('--hesin', required=True, help='Path to hesin diagnosis file (.txt or .csv)')
    parser.add_argument('--eid', required=True, help='Path to EID list (one column)')
    parser.add_argument('--icd_embedding', required=True, help='Path to ICD-10 code embedding file (.pkl or .csv/.txt)')
    parser.add_argument('--output_prefix', required=True, help='Prefix for output .pkl file')
    parser.add_argument('--eid_col', type=int, default=0, help='Column index for EIDs in eid file (default: 0)')
    parser.add_argument('--sep_eid', default='\t', help='Separator for EID file (default: tab)')
    parser.add_argument('--icd_embedding_from_pkl', action='store_true', help='Set if ICD embeddings are in .pkl format')
    parser.add_argument('--icd_col', default='diag_icd10', help='Column name for ICD-10 codes in hesin file (default: diag_icd10)')
    args = parser.parse_args()

    main(
        hesin_path=args.hesin,
        eid_path=args.eid,
        icd_embedding_path=args.icd_embedding,
        output_prefix=args.output_prefix,
        eid_col=args.eid_col,
        sep_eid=args.sep_eid,
        icd_embedding_from_pkl=args.icd_embedding_from_pkl,
        icd_col=args.icd_col
    )