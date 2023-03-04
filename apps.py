import streamlit as st
from streamlit_option_menu import option_menu
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
############################################################################################################################################################################################
def acetylcholinesterase_single_build_IC50_model(input_data):
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


def acetylcholinesterase_single_build_Ki_model(input_data):
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

#######################################################################################################################################################################
def acetylcholinesterase_multiple_build_IC50_model(input_data):
    # Reads in saved regression model
    load_model = pickle.load(open('acetylcholinesterase_model.pkl', 'rb'))
    # Apply model to make predictions
    prediction = load_model.predict(input_data)
    st.header('**Prediction output**')
    prediction_output = pd.Series(prediction, name='pIC50')
    molecule_name = dfss.iloc[:,[0]]
    canonical_smile = dfss.iloc[:,[1]]
    df = pd.concat([molecule_name,canonical_smile, prediction_output], axis=1)
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

def acetylcholinesterase_multiple_build_Ki_model(input_data):
    # Reads in saved regression model
    load_model = pickle.load(open('acetylcholinesterase_Ki_model.pkl', 'rb'))
    # Apply model to make predictions
    prediction = load_model.predict(input_data)
    st.header('**Prediction output**')
    prediction_output = pd.Series(prediction, name='pKi')
    molecule_name = dfss.iloc[:,[0]]
    canonical_smile = dfss.iloc[:,[1]]
    df = pd.concat([molecule_name,canonical_smile, prediction_output], axis=1)
    df['Ki'] = (10**(9-(df['pKi'])))
    st.write(df)
    st.markdown(filedownload(df), unsafe_allow_html=True)

#######################################################################################################################################################################

# import streamlit as st

selected = option_menu(
            menu_title=None,  # required
            options=["Home", "Multiple Prediction", "Contact", 'Credits'],  # required
            icons=["house", "book", "envelope","hand-thumbs-up"],  # optional
            menu_icon="cast",  # optional
            default_index=0,  # optional
            orientation="horizontal",
        )
########################################################################################################################################################################################
if selected =='Multiple Prediction':
    st.header('Upload your CSV data')
    uploaded_file = st.file_uploader("Upload your input file", type=['csv'])
    st.write('Input CSV file format is shown below: Compound name or ID on first column, Canonical smile on second column:')
    st.markdown("""
    [Example input file](https://raw.githubusercontent.com/Macaulay001/Machine_Learning/main/example_file.csv)
    """)

    receptor_mode = st.selectbox('Choose a Target Receptor Activity',['Acetylcholinesterase IC50', 'Acetylcholinesterase Binding Affinity (Ki)'])

    if st.button('Predict'):
        dfs = pd.read_csv(uploaded_file,  sep=',')
        load_data = dfs.iloc[:,[1]]
        load_data.to_csv('molecule.smi', sep = '\t', header = False, index = False)
        lr = load_data.shape
        percentag = str(int(lr[0]) * 0.0067)
        percentage = str(round(float(percentag), 2))
        st.header('**Original input data**')
        st.write(load_data)

        with st.spinner("Calculating descriptors in " + percentage+"mins..."):
            desc_calc()

        # Read in calculated descriptors and display the dataframe
        # st.header('**Calculated molecular descriptors**')
        desc = pd.read_csv('descriptors_output.csv')
        input_column = dfs.iloc[:, 0:2]
        descs = pd.concat([desc, input_column], axis =1)
        descs = descs.dropna(how='any')
        dfss = descs.iloc[:, 882:]
        desc = descs.iloc[:, 0:882]
        dfss.reset_index(drop=True, inplace=True)
        desc.reset_index(drop=True, inplace=True)
        # st.write(desc)
        # st.write(desc.shape)

        # Read descriptor list used in previously built model
        if receptor_mode=='Acetylcholinesterase IC50':
            # st.header('**Subset of descriptors from previously built Acetylcholinesterase models**')
            Xlist = list(pd.read_csv('descriptor_list.csv').columns)
            desc_subset = desc[Xlist]
            # st.write(desc_subset)
            # st.write(desc_subset.shape)
            # Apply trained model to make prediction on query compounds
            acetylcholinesterase_multiple_build_IC50_model(desc_subset)
        elif receptor_mode=='Acetylcholinesterase Binding Affinity (Ki)':
            # st.header('**Subset of descriptors from previously built Acetylcholinesterase models**')
            Xlist = list(pd.read_csv('acetylcholinesterase_Ki_descriptor_list.csv').columns)
            desc_subset = desc[Xlist]
            # st.write(desc_subset)
            # st.write(desc_subset.shape)
            # Apply trained model to make prediction on query compounds
            acetylcholinesterase_multiple_build_Ki_model(desc_subset)
        else:
            st.info('Reload Page and choose a target Receptor!')

        
    else:
        st.info('Upload input data to start!')


####################################################################################################################################################################################################################################################################
if selected =='Credits':
    st.markdown("""
    # **Credits**
    - Built using `Python` and `Streamlit` by Nutritional and industrial Laboratory, Department of Biochemistry, University of Ibadan
    - Descriptors were calculated using [PaDEL-Descriptor](http://www.yapcwsoft.com/dd/padeldescriptor/) [[Reference Paper]](https://doi.org/10.1002/jcc.21707).
    - Chanin Nantasenamat (Data Professor)
    ---
    """)
####################################################################################################################################################################################################################################
if selected =='Home':
    
    # # Logo image
    # image = Image.open('logo.png')

    # st.image(image, use_column_width=True)

    # Page title
    st.markdown("""
    # Compounds Bioactivity Prediction App for Drug Target in Metabolic Syndromes

    This app allows you to predict the bioactivity of compounds towards acetylcholinesterase enzymes.

    **Credits**
    - Built using `Python` and `Streamlit` by Nutritional and industrial Laboratory, Department of Biochemistry, University of Ibadan
    """)
    smile = st.text_input('Enter canonical smile')
    smiles = smile.split()
    dsf = pd.DataFrame(smiles, columns=['Canonical smile'])
    dsf['Name'] ='Compound'
    # st.write(dsf)
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
            acetylcholinesterase_single_build_IC50_model(desc_subset)
        elif receptor_mode=='Acetylcholinesterase Binding Affinity (Ki)':
            Xlist = list(pd.read_csv('acetylcholinesterase_Ki_descriptor_list.csv').columns)
            desc_subset = desc[Xlist]
            acetylcholinesterase_single_build_Ki_model(desc_subset)
        
        else:
            st.info('Reload Page and choose a target Receptor!')

        
    else:
        st.info('Upload input data to start!')
####################################################################################################################################################################################################################################
if selected =='Contact':
    st.write('Contacts will be available soon')