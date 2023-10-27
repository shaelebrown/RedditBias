
# this scipt was based on code I recieved from Colleen
# used to generate data for computing persistence diagrams

# import modules
import pandas as pd
import numpy as np
from sentence_transformers import SentenceTransformer

# get data and fill missing values of text
mydata = pd.read_csv('data/sampled_bias.csv',encoding='latin')
mydata['Reddit_Comments'] = mydata['Reddit_Comments'].fillna(value=".")

# strip to text for input into BERT model
text_list = list(mydata.Reddit_Comments)

# get BERT--384 vectors
sbert_model = SentenceTransformer('all-MiniLM-L6-v2')

# encode data with BERT
encoded_text = sbert_model.encode(text_list)

# wrangle to array
BERT_array = np.array([x for x in encoded_text])

# convert to dataframe
BERT_df = pd.DataFrame(BERT_array)

# split by bias type
race_df = BERT_df[mydata['Bias_Type'] == 'race']
gender_df = BERT_df[mydata['Bias_Type'] == 'gender']
orientation_df = BERT_df[mydata['Bias_Type'] == 'orientation']
religion_df = BERT_df[mydata['Bias_Type'] == 'religion']

# save to data directory
race_df.to_csv('data/race_df.csv')
gender_df.to_csv('data/gender_df.csv')
orientation_df.to_csv('data/orientation_df.csv')
religion_df.to_csv('data/religion_df.csv')
