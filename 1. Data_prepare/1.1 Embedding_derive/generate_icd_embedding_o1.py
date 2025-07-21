import requests
from bs4 import BeautifulSoup
import html2text
import mygene
import json
import pickle
mg = mygene.MyGeneInfo()

import pandas as pd
import scanpy as sc
import pickle


import openai
import time
delay_sec = 5
# remember to set your open AI API key!
openai.api_key = 'set api'

import numpy as np

df = pd.read_csv("./hesin_icd10_descrip_embed.txt", sep='\t', index_col=0)

EMBED_DIM = 3072 # embedding dim from GPT-3.5

def get_gpt_embedding(text, model="text-embedding-3-large"):
    text = text.replace("\n", " ")
    return np.array(openai.Embedding.create(input = [text], model=model)['data'][0]['embedding'])

gene_name_to_GPT_response = {}
gene_name_getembedding = {}

df_setdisease = sorted(set(df['meaning'].values))
gene_completion_test = df_setdisease

gene_completion_test = sorted(set(df_setdisease))

# gene_completion_test = list(GPT_3_5_gene_embeddings.keys())
for gene in gene_completion_test:
    print(gene)
    try:
        input_data = f"Please provide a concise, one-paragraph summary of disease {gene}, incorporating established medical facts about its causes, clinical features, diagnosis, and treatment. If any aspect is unclear or not well-documented, simply state, 'I do not know.'"
        completion = openai.ChatCompletion.create(model="o1-preview", 
                    messages=[{"role": "user", 
                               "content": f"[Instruction]:{input_data}\n[Agent Description]:",
                               }], seed = 0)
        gene_name_to_GPT_response[gene] = completion.choices[0].message.content
        gene_name_getembedding[gene] = get_gpt_embedding(gene_name_to_GPT_response[gene])
        time.sleep(1)
    except (openai.APIError, 
                    openai.error.APIError, 
                    openai.error.APIConnectionError, 
                    openai.error.RateLimitError, 
                    openai.error.ServiceUnavailableError, 
                    openai.error.Timeout) as e:
        #Handle API error here, e.g. retry or log
        print(f"OpenAI API returned an API Error: {e}")
        pass


with open('ensem_describe_diseaseiccode.pickle', 'wb') as handle:
    pickle.dump(gene_name_to_GPT_response, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
with open('ensem_emb_diseaseiccode.pickle', 'wb') as handle:
    pickle.dump(gene_name_getembedding, handle, protocol=pickle.HIGHEST_PROTOCOL)
