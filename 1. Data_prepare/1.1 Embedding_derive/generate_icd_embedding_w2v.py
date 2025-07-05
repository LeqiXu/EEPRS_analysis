import pandas as pd
import numpy as np
from gensim.models import Word2Vec
from nltk.tokenize import word_tokenize
from nltk.corpus import stopwords
import string
import argparse
import os

# If running for the first time, uncomment the following lines to download NLTK resources
#import nltk
#nltk.download('punkt')
#nltk.download('stopwords')

def clean_tokenize(description):
    description = description.lower()
    description = description.translate(str.maketrans('', '', string.punctuation))
    words = word_tokenize(description)
    words = [word for word in words if word not in stopwords.words('english')]
    return words

def get_vector(word, model):
    try:
        return model.wv[word]
    except KeyError:
        return np.zeros(model.vector_size)

def average_embedding(tokens, model):
    vectors = [get_vector(token, model) for token in tokens]
    vectors = [vec for vec in vectors if np.any(vec)]
    if vectors:
        return np.mean(vectors, axis=0)
    else:
        return np.zeros(model.vector_size)

def generate_icd_embeddings(
    hesin_path,
    icd_mapping_path,
    output_path,
    eid_list_path=None,
    eid_column=None,
    sep='\t'
):
    print(f"Loading hesin file: {hesin_path}")
    hesin_df = pd.read_table(hesin_path)
    # If using a subset of individuals, filter here
    if eid_list_path:
        print(f"Filtering with eid list: {eid_list_path}")
        eid_list = pd.read_csv(eid_list_path, sep=sep)
        if eid_column is not None:
            eids = eid_list[eid_column]
        else:
            eids = eid_list.iloc[:, 0]
        hesin_df = hesin_df[hesin_df['eid'].isin(eids)]
    
    print(f"Aggregating descriptions...")
    df_aggregated = hesin_df.groupby('eid')['meaning'].apply(' '.join).reset_index()
    df_aggregated['tokens'] = df_aggregated['meaning'].apply(clean_tokenize)
    sentences = df_aggregated['tokens'].tolist()

    print(f"Training Word2Vec on {len(sentences)} sentences...")
    model = Word2Vec(sentences, vector_size=100, window=5, min_count=1, workers=4)

    print(f"Loading ICD-10 mapping: {icd_mapping_path}")
    icd_description = pd.read_csv(icd_mapping_path)
    icd_description['tokens'] = icd_description['meaning'].apply(clean_tokenize)
    icd_description['average_embedding'] = icd_description['tokens'].apply(
        lambda tokens: average_embedding(tokens, model)
    )

    embeddings_df = pd.DataFrame(
        icd_description['average_embedding'].tolist(),
        index=icd_description['coding']
    )
    embeddings_df.columns = [f'Dimension_{i+1}' for i in range(embeddings_df.shape[1])]

    print(f"Saving embeddings to: {output_path}")
    embeddings_df.to_csv(output_path, sep='\t', header=True, index=True)

    print(f"Done. Embeddings shape: {embeddings_df.shape}")
    return embeddings_df

# ----------------- CLI for easy use -----------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate ICD-10 code embeddings from hesin descriptions.")
    parser.add_argument("--hesin", type=str, required=True, help="Path to hesin description file (TSV).")
    parser.add_argument("--icd_mapping", type=str, required=True, help="Path to ICD-10 mapping CSV.")
    parser.add_argument("--output", type=str, required=True, help="Path to save the embeddings.")
    parser.add_argument("--eid_list", type=str, default=None, help="Optional: Path to file containing list of eids.")
    parser.add_argument("--eid_column", type=int, default=None, help="Optional: Index of eid column in eid_list file (default: 0).")
    parser.add_argument("--eid_sep", type=str, default='\t', help="Optional: Separator for eid_list file (default: tab).")

    args = parser.parse_args()

    generate_icd_embeddings(
        hesin_path=args.hesin,
        icd_mapping_path=args.icd_mapping,
        output_path=args.output,
        eid_list_path=args.eid_list,
        eid_column=args.eid_column,
        sep=args.eid_sep
    )