import streamlit as st
import pandas as pd
from PIL import Image
import subprocess
import os
import base64
import pickle
import csv

st.set_page_config(
   page_title="Acetylcholine",
   page_icon=Image.open('icons.ico'),
   layout="wide",
#    initial_sidebar_state="expanded",
)
# Molecular descriptor calculator
def desc_calc():
    # Performs the descriptor calculation
    bashCommand = "java -Xms2G -Xmx2G -Djava.awt.headless=true -jar ./PaDEL-Descriptor/PaDEL-Descriptor.jar -removesalt -standardizenitro -fingerprints -descriptortypes ./PaDEL-Descriptor/PubchemFingerprinter.xml -dir ./ -file descriptors_output.csv"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    os.remove('molecule.smi')

# File download
def filedownload(df):
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()  # strings <-> bytes conversions
    href = f'<a href="data:file/csv;base64,{b64}" download="prediction.csv">Download Predictions</a>'
    return href

# Model building
def acetylcholinesterase_build_model(input_data):
    # Reads in saved regression model
    load_model = pickle.load(open('acetylcholinesterase_model.pkl', 'rb'))
    # Apply model to make predictions
    prediction = load_model.predict(input_data)
    st.header('**Prediction output**')
    prediction_output = pd.Series(prediction, name='pIC50')
    canonical_smile =  load_data.iloc[:,[0]]
    df = pd.concat([canonical_smile, prediction_output], axis=1)
    df['IC50'] = (10**(9-(df['pIC50'])))
    bioactivity_class = []
    for i in df.IC50:
        if float(i) >= 10000:
            bioactivity_class.append("inactive")
        elif float(i) <= 1000:
            bioactivity_class.append("active")
        else:
            bioactivity_class.append("intermediate")
    bioactivity_class = pd.Series(bioactivity_class, name='Bioactivity_class')
    df = pd.concat([df, bioactivity_class], axis=1)

    st.write(df)
    st.markdown(filedownload(df), unsafe_allow_html=True)


def acetylcholinesterase_build_Ki_model(input_data):
    # Reads in saved regression model
    load_model = pickle.load(open('acetylcholinesterase_Ki_model.pkl', 'rb'))
    # Apply model to make predictions
    prediction = load_model.predict(input_data)
    st.header('**Prediction output**')
    prediction_output = pd.Series(prediction, name='pKi')
    canonical_smile =  load_data.iloc[:,[0]]
    df = pd.concat([canonical_smile, prediction_output], axis=1)
    df['Binding Affinity'] = round(float(10**(9-(df['pKi']))),2)
    st.write(df)
    st.markdown(filedownload(df), unsafe_allow_html=True)


# Logo image
image = Image.open('logo.png')

st.image(image, use_column_width=True)

# Page title
st.markdown("""
# Compounds Bioactivity Prediction App for Drug Target in Metabolic Syndromes

This app allows you to predict the bioactivity towards inhibiting several protein ligands and enzymes.

**Credits**
- Built using `Python` and `Streamlit` by Nutritional and industrial Laboratory, Department of Biochemistry, University of Ibadan
- Descriptors were calculated using [PaDEL-Descriptor](http://www.yapcwsoft.com/dd/padeldescriptor/) [[Reference Paper]](https://doi.org/10.1002/jcc.21707).
- Chanin Nantasenamat (Data Professor)
---
""")
smile = st.text_input('Enter canonical smile')
smiles = smile.split()
dsf = pd.DataFrame(smiles, columns=['Canonical smile'])
dsf['Name'] ='Compound'
#st.write(dsf)
receptor_mode = st.selectbox('Choose a Target Receptor Activity',['Acetylcholinesterase IC50', 'Acetylcholinesterase Binding Affinity (Ki)'])
if st.button('Predict'):
    load_data = dsf
    load_data.to_csv('molecule.smi', sep = '\t', header = False, index = False)
    with st.spinner("Calculating descriptors in 4 seconds..."):
        desc_calc()
    desc = pd.read_csv('descriptors_output.csv')
    desc = desc.loc[desc['Name'] == 'Compound']
    desc = desc.reset_index()
    # Read descriptor list used in previously built model
    if receptor_mode=='Acetylcholinesterase IC50':
        # st.header('**Subset of descriptors from previously built Acetylcholinesterase models**')
        Xlist = list(pd.read_csv('descriptor_list.csv').columns)
        desc_subset = desc[Xlist]
        acetylcholinesterase_build_model(desc_subset)
    elif receptor_mode=='Acetylcholinesterase Binding Affinity (Ki)':
        Xlist = list(pd.read_csv('acetylcholinesterase_Ki_descriptor_list.csv').columns)
        desc_subset = desc[Xlist]
        acetylcholinesterase_build_Ki_model(desc_subset)
    
    else:
        st.info('Reload Page and choose a target Receptor!')

    
else:
    st.info('Input canonical smile to start!')
