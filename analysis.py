
# this scipt was based on code I recieved from Colleen

# import modules
import pandas as pd
import numpy as np
from sentence_transformers import SentenceTransformer
from ripser import ripser
from persim import plot_diagrams

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

# compute persistence diagrams
race_diag = ripser(race_df, maxdim = 2)['dgms']
gender_diag = ripser(gender_df, maxdim = 2)['dgms']
orientation_diag = ripser(orientation_df, maxdim = 2)['dgms']
religion_diag = ripser(religion_df, maxdim = 2)['dgms']

# visualize
plot_diagrams(race_diag, show=True, title = 'Race')
plot_diagrams(gender_diag, show=True, title = 'Gender')
plot_diagrams(orientation_diag, show=True, title = 'Orientation')
plot_diagrams(religion_diag, show=True, title = 'Religion')