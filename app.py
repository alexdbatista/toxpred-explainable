"""
ToxPred-Explainable: Interpretable Molecular Toxicity Screening
Main Streamlit Application with Explainability Engine

Author: Alex Domingues Batista
"""

import io
import os
import pickle
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st
from PIL import Image
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.Draw import SimilarityMaps, rdMolDraw2D

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from config import ASSAY_INFO, MODEL_PATH
from utils import (
    calculate_lipinski,
    explain_molecule,
    get_morgan_fingerprint,
    name_to_smiles,
    validate_smiles,
)

# Page config
st.set_page_config(
    page_title="ToxPred-Explainable", 
    page_icon="🧪", 
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS - Modern Design
st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;600;700&display=swap');
    
    .main {
        font-family: 'Inter', sans-serif;
        max-width: 1400px;
        padding: 1rem;
    }
    
    /* Sidebar styling for better readability */
    [data-testid="stSidebar"] {
        background: linear-gradient(180deg, #1a1a2e 0%, #16213e 100%);
        padding: 1rem 0.5rem;
    }
    
    [data-testid="stSidebar"] * {
        color: white !important;
    }
    
    [data-testid="stSidebar"] h2 {
        font-size: 1.5rem !important;
    }
    
    [data-testid="stSidebar"] h3 {
        color: white !important;
        font-weight: 600;
        font-size: 1rem !important;
    }
    
    [data-testid="stSidebar"] .stRadio > label {
        color: rgba(255,255,255,0.95) !important;
        font-size: 0.9rem !important;
    }
    
    [data-testid="stSidebar"] [data-testid="stMarkdownContainer"] p {
        color: rgba(255,255,255,0.9) !important;
        font-size: 0.9rem !important;
    }
    
    .hero-header {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        padding: 1.5rem 2rem;
        border-radius: 16px;
        text-align: center;
        margin-bottom: 1.5rem;
        box-shadow: 0 8px 24px rgba(0,0,0,0.15);
    }
    
    .hero-title {
        font-size: 2.2rem;
        font-weight: 700;
        color: white;
        margin: 0;
        text-shadow: 2px 2px 4px rgba(0,0,0,0.2);
    }
    
    .hero-subtitle {
        font-size: 1rem;
        color: rgba(255,255,255,0.95);
        margin-top: 0.5rem;
    }
    
    .section-header {
        font-size: 1.5rem;
        font-weight: 700;
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        margin: 1rem 0 0.5rem 0;
    }
    
    .feature-card {
        background: white;
        padding: 1.2rem;
        border-radius: 12px;
        box-shadow: 0 2px 12px rgba(0,0,0,0.08);
        border: 1px solid #e9ecef;
        margin-bottom: 0.8rem;
    }
    
    .feature-card h3 {
        font-size: 1.2rem !important;
        margin-bottom: 1rem;
    }
    
    .feature-card h4 {
        font-size: 1rem !important;
        margin: 1rem 0 0.5rem 0;
    }
    
    .toxic-box {
        background: linear-gradient(135deg, #ff6b6b 0%, #ee5a6f 100%);
        color: white;
        padding: 1rem 1.5rem;
        border-radius: 12px;
        box-shadow: 0 4px 16px rgba(255,107,107,0.3);
        font-weight: 600;
        font-size: 1rem;
        text-align: center;
        margin: 0.8rem 0;
    }
    
    .safe-box {
        background: linear-gradient(135deg, #51cf66 0%, #37b24d 100%);
        color: white;
        padding: 1rem 1.5rem;
        border-radius: 12px;
        box-shadow: 0 4px 16px rgba(55,178,77,0.3);
        font-weight: 600;
        font-size: 1rem;
        text-align: center;
        margin: 0.8rem 0;
    }
    
    .info-card {
        background: linear-gradient(135deg, #4facfe 0%, #00f2fe 100%);
        color: white;
        padding: 1.2rem;
        border-radius: 12px;
        box-shadow: 0 6px 20px rgba(79,172,254,0.25);
        margin: 0.8rem 0;
    }
    
    .info-card h4 {
        font-size: 1rem !important;
        margin: 0 0 0.5rem 0;
    }
    
    .info-card p {
        font-size: 0.95rem;
        margin: 0;
    }
    
    .stats-box {
        background: rgba(255,255,255,0.95);
        padding: 1rem;
        border-radius: 12px;
        box-shadow: 0 3px 12px rgba(0,0,0,0.1);
        margin: 0.8rem 0;
    }
    
    .badge {
        display: inline-block;
        padding: 0.4rem 1rem;
        border-radius: 16px;
        font-size: 0.8rem;
        font-weight: 600;
        margin: 0.3rem;
    }
    
    .badge-success {
        background: linear-gradient(135deg, #51cf66, #37b24d);
        color: white;
    }
    
    .badge-info {
        background: linear-gradient(135deg, #4facfe, #00f2fe);
        color: white;
    }
    
    .badge-warning {
        background: linear-gradient(135deg, #ffd93d, #ffc93d);
        color: #333;
    }
    
    .color-legend {
        display: flex;
        align-items: center;
        margin: 10px 0;
        padding: 10px;
        background: rgba(255,255,255,0.1);
        border-radius: 10px;
    }
    
    .color-box {
        width: 30px;
        height: 30px;
        border-radius: 8px;
        margin-right: 15px;
        box-shadow: 0 2px 8px rgba(0,0,0,0.2);
    }
</style>
""", unsafe_allow_html=True)

# Load model
@st.cache_resource
def load_model():
    """Load model, training it first if it doesn't exist"""
    if not os.path.exists(MODEL_PATH):
        st.warning("⏳ Model not found. Training now... This will take ~30 seconds.")
        try:
            # Import and run training
            from src.train_model import train_and_save_model
            train_and_save_model()
            st.success("✅ Model trained successfully!")
        except Exception as e:
            st.error(f"❌ Error training model: {str(e)}")
            return None
    
    try:
        with open(MODEL_PATH, 'rb') as f:
            return pickle.load(f)
    except Exception as e:
        st.error(f"❌ Error loading model: {str(e)}")
        return None

model = load_model()

# Sidebar
with st.sidebar:
    st.markdown("""
    <div style='text-align: center; padding: 20px 0;'>
        <div style='font-size: 4rem; margin-bottom: 10px;'>&#x1F9EA;</div>
        <h2 style='color: white; margin: 0; font-weight: 700;'>ToxPred</h2>
        <p style='color: rgba(255,255,255,0.8); font-size: 0.9rem; margin: 5px 0;'>Explainable</p>
        <div style='background: rgba(255,255,255,0.3); height: 2px; width: 60px; margin: 15px auto; border-radius: 2px;'></div>
        <p style='color: rgba(255,255,255,0.9); font-size: 0.95rem;'><strong>Interpretable AI for Drug Safety</strong></p>
    </div>
    """, unsafe_allow_html=True)
    
    st.markdown("<br>", unsafe_allow_html=True)
    
    # Model Info Card
    st.markdown("""
    <div class='stats-box'>
        <h3 style='color: #667eea; margin-top: 0; font-size: 1.2rem;'>&#x1F4CA; Model Stats</h3>
        <div style='margin: 15px 0;'>
            <div style='display: flex; justify-content: space-between; margin: 8px 0;'>
                <span style='color: #666;'><strong>Assay:</strong></span>
                <span style='color: #333; font-size: 0.85rem;'>SR-ARE</span>
            </div>
            <div style='display: flex; justify-content: space-between; margin: 8px 0;'>
                <span style='color: #666;'><strong>Training Data:</strong></span>
                <span style='color: #333;'><strong>5,832</strong> molecules</span>
            </div>
            <div style='display: flex; justify-content: space-between; margin: 8px 0;'>
                <span style='color: #666;'><strong>Test Accuracy:</strong></span>
                <span style='color: #51cf66; font-weight: 700;'>86.6%</span>
            </div>
            <div style='display: flex; justify-content: space-between; margin: 8px 0;'>
                <span style='color: #666;'><strong>ROC-AUC:</strong></span>
                <span style='color: #4facfe; font-weight: 700;'>0.822</span>
            </div>
        </div>
        <div style='background: linear-gradient(135deg, #667eea, #764ba2); padding: 10px; border-radius: 8px; text-align: center; color: white; font-weight: 600; margin-top: 15px;'>
            &#x2728; Atom-Level Explainability
        </div>
    </div>
    """, unsafe_allow_html=True)
    
    st.markdown("<br>", unsafe_allow_html=True)
    
    # Model explanation expander
    with st.expander("ℹ️ About SR-ARE Assay", expanded=False):
        st.markdown("""
        **Stress Response - Antioxidant Response Element**
        
        The SR-ARE assay detects compounds that activate cellular stress response pathways. 
        Activation indicates potential toxicity through oxidative stress.
        
        **Why it matters:**
        - Early toxicity warning signal
        - Predicts liver and kidney damage
        - Essential for drug safety screening
        """)
    
    # Navigation
    st.markdown("<h3 style='color: white; font-size: 1.1rem;'>&#x1F3AF; Navigation</h3>", unsafe_allow_html=True)
    page = st.radio(
        "Select Function:",
        ["Single Prediction", "Batch Analysis", "About"],
        label_visibility="collapsed"
    )
    
    st.markdown("<br><br>", unsafe_allow_html=True)
    
    # Quick Tips
    st.markdown("""
    <div style='background: rgba(255,255,255,0.1); padding: 15px; border-radius: 10px; color: white;'>
        <h4 style='margin-top: 0; font-size: 1rem;'>&#x1F4A1; Quick Tips</h4>
        <ul style='font-size: 0.85rem; line-height: 1.6; padding-left: 20px;'>
            <li>Red atoms = toxic contributors</li>
            <li>Green atoms = safe contributors</li>
            <li>Try example molecules first</li>
            <li>Download heatmaps for reports</li>
        </ul>
    </div>
    """, unsafe_allow_html=True)

# Hero Header
st.markdown("""
<div class='hero-header'>
    <h1 class='hero-title'>&#x1F9EA; ToxPred-Explainable</h1>
    <p class='hero-subtitle'>Interpretable Molecular Toxicity Screening with Atom-Level Attribution</p>
    <div style='margin-top: 1.5rem;'>
        <span class='badge badge-success'>&#x2713; Production Ready</span>
        <span class='badge badge-info'>&#x1F3AF; 86.6% Accuracy</span>
        <span class='badge badge-warning'>&#x1F525; Explainable AI</span>
    </div>
</div>
""", unsafe_allow_html=True)

# Page routing
if page == "Single Prediction":
    st.markdown("<h2 class='section-header'>&#x1F52C; Single Molecule Analysis</h2>", unsafe_allow_html=True)
    
    col1, col2 = st.columns([1.2, 1])
    
    with col1:
        st.markdown("""
        <div class='feature-card'>
            <h3 style='color: #667eea; margin-top: 0;'>&#x1F4DD; Input Molecule</h3>
        """, unsafe_allow_html=True)
        
        # Input method selection
        input_method = st.radio(
            "Search by:",
            ["Chemical Name", "SMILES String"],
            horizontal=True,
            help="Choose to search by common name (e.g., 'aspirin') or SMILES notation"
        )
        
        # Initialize session state for inputs
        if 'example_value' not in st.session_state:
            st.session_state.example_value = ""
        
        if input_method == "Chemical Name":
            chemical_name = st.text_input(
                "Enter Chemical Name:",
                value=st.session_state.example_value if st.session_state.get('example_type') == 'name' else "",
                placeholder="e.g., aspirin, caffeine, ibuprofen",
                help="Common name or IUPAC name"
            )
            smiles_input = None
        else:
            smiles_input = st.text_input(
                "Enter SMILES:",
                value=st.session_state.example_value if st.session_state.get('example_type') == 'smiles' else "",
                placeholder="e.g., CCOc1ccc2nc(S(N)(=O)=O)sc2c1",
                help="Simplified Molecular Input Line Entry System"
            )
            chemical_name = None
        
        # Example molecules
        st.markdown("<h4 style='color: #667eea; margin-top: 20px;'>&#x1F48A; Example Molecules</h4>", unsafe_allow_html=True)
        
        if input_method == "Chemical Name":
            examples = {
                "Aspirin": "aspirin",
                "Caffeine": "caffeine",
                "Ibuprofen": "ibuprofen",
                "Paracetamol": "paracetamol"
            }
        else:
            examples = {
                "Aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
                "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
                "Benzoic Acid": "c1ccccc1C(=O)O",
                "Toxic Example": "CCOc1ccc2nc(S(N)(=O)=O)sc2c1"
            }
        
        example_cols = st.columns(4)
        for i, (name, value) in enumerate(examples.items()):
            with example_cols[i]:
                if st.button(name, key=f"ex_{i}_{input_method}"):
                    st.session_state.example_value = value
                    st.session_state.example_type = 'name' if input_method == "Chemical Name" else 'smiles'
                    st.rerun()
        
        st.markdown("</div>", unsafe_allow_html=True)
        
        if st.button("Analyze Molecule", type="primary", use_container_width=True):
            # Convert chemical name to SMILES if needed
            if input_method == "Chemical Name":
                if not chemical_name:
                    st.error("Please enter a chemical name!")
                    st.stop()
                
                with st.spinner(f"Looking up '{chemical_name}' in PubChem database..."):
                    smiles_input, status = name_to_smiles(chemical_name)
                    
                    if status != "success":
                        st.error(f"&#x274C; {status}")
                        st.info("&#x1F4A1; Try using SMILES input method instead, or check the spelling")
                        st.stop()
                    else:
                        st.success(f"&#x2705; Found: `{smiles_input}`")
            
            if not smiles_input:
                st.error("Please enter a SMILES string!")
            elif not validate_smiles(smiles_input):
                st.error("Invalid SMILES string! Please check your input.")
            elif model is None:
                st.error("Model not loaded!")
            else:
                with st.spinner("Analyzing molecule..."):
                    st.session_state.smiles = smiles_input
                    st.session_state.analyzed = True
    
    with col2:
        if 'analyzed' in st.session_state and st.session_state.analyzed:
            smiles = st.session_state.smiles
            mol = Chem.MolFromSmiles(smiles)
            
            st.markdown("""
            <div class='feature-card'>
                <h3 style='color: #667eea; margin-top: 0;'>&#x1F52C; Molecular Structure</h3>
            """, unsafe_allow_html=True)
            
            # Generate high-quality molecular image
            drawer = rdMolDraw2D.MolDraw2DCairo(600, 600)
            drawer.drawOptions().addAtomIndices = False
            drawer.drawOptions().bondLineWidth = 2
            drawer.DrawMolecule(mol)
            drawer.FinishDrawing()
            img_data = drawer.GetDrawingText()
            
            # Convert to PIL Image for display
            img = Image.open(io.BytesIO(img_data))
            st.image(img, width=350)
            st.markdown("</div>", unsafe_allow_html=True)
            
            # Lipinski's Rule
            lipinski = calculate_lipinski(mol)
            st.markdown("""
            <div class='feature-card'>
                <h3 style='color: #667eea; margin-top: 0;'>&#x1F48A; Drug-Likeness (Lipinski's Rule)</h3>
            """, unsafe_allow_html=True)
            
            lipinski_cols = st.columns(4)
            metrics = [
                ("MW", lipinski['MW'], "< 500"),
                ("LogP", lipinski['LogP'], "< 5"),
                ("HBD", lipinski['HBD'], "< 5"),
                ("HBA", lipinski['HBA'], "< 10")
            ]
            
            for col, (name, value, threshold) in zip(lipinski_cols, metrics):
                with col:
                    st.metric(name, f"{value:.1f}", threshold)
            
            passes = lipinski['Passes']
            if passes:
                st.markdown('<div class="safe-box">&#x2705; Passes Lipinski\'s Rule - Drug-Like Molecule!</div>', 
                           unsafe_allow_html=True)
            else:
                st.markdown('<div class="toxic-box">&#x26A0;&#xFE0F; Violates Lipinski\'s Rule</div>', 
                           unsafe_allow_html=True)
            
            st.markdown("</div>", unsafe_allow_html=True)
    
    # Prediction and Explainability
    if 'analyzed' in st.session_state and st.session_state.analyzed:
        st.markdown("<h2 class='section-header'>&#x1F3AF; Toxicity Prediction</h2>", unsafe_allow_html=True)
        
        smiles = st.session_state.smiles
        mol = Chem.MolFromSmiles(smiles)
        
        # Get prediction
        fp = get_morgan_fingerprint(mol)
        prediction_proba = model.predict_proba(fp.reshape(1, -1))[0][1]
        prediction_class = "TOXIC" if prediction_proba > 0.5 else "SAFE"
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.metric("Classification", prediction_class)
        with col2:
            st.metric("Toxicity Probability", f"{prediction_proba:.1%}")
        with col3:
            confidence = max(prediction_proba, 1 - prediction_proba)
            st.metric("Confidence", f"{confidence:.1%}")
        
        if prediction_proba > 0.5:
            st.markdown(f'<div class="toxic-box">&#x2620;&#xFE0F; <b>TOXIC COMPOUND DETECTED</b><br><small>{ASSAY_INFO["name"]} positive response</small></div>', 
                       unsafe_allow_html=True)
            
            with st.expander("🔍 What does this mean?"):
                st.markdown("""
                **Positive SR-ARE Response Detected**
                
                This compound is predicted to activate cellular stress response pathways, which may indicate:
                
                - ⚠️ Potential oxidative stress induction
                - ⚠️ Risk of cellular damage
                - ⚠️ May require further safety evaluation
                
                **Confidence Level:** {:.1f}%
                
                **Recommendation:** This prediction suggests caution. Consider:
                - Further *in vitro* testing
                - Structure modification to reduce toxicity
                - Alternative compounds if available
                """.format(confidence * 100))
        else:
            st.markdown(f'<div class="safe-box">&#x2705; <b>SAFE COMPOUND</b><br><small>No {ASSAY_INFO["name"]} response detected</small></div>', 
                       unsafe_allow_html=True)
            
            with st.expander("✅ What does this mean?"):
                st.markdown("""
                **Negative SR-ARE Response Detected**
                
                This compound is predicted to NOT activate stress response pathways, suggesting:
                
                - ✅ Lower risk of oxidative stress
                - ✅ Potentially safer toxicity profile
                - ✅ Good candidate for further development
                
                **Confidence Level:** {:.1f}%
                
                **Note:** This prediction indicates lower risk in this specific assay. Complete safety assessment requires multiple toxicity tests.
                """.format(confidence * 100))
        
        # Explainability
        st.markdown("<h2 class='section-header'>&#x1F525; Explainability: Atom-Level Attribution</h2>", unsafe_allow_html=True)
        st.markdown("""
        <div class='info-card'>
            <h4 style='margin-top: 0; font-size: 1.1rem;'>&#x1F3A8; Visual Attribution Map</h4>
            <p style='margin: 0;'>See which molecular substructures drive the toxicity prediction</p>
        </div>
        """, unsafe_allow_html=True)
        
        with st.spinner("Generating explainability heatmap..."):
            heatmap_img, max_weight = explain_molecule(mol, model)
        
        col1, col2 = st.columns([2, 1])
        
        with col1:
            st.markdown("<div class='feature-card'>", unsafe_allow_html=True)
            st.image(heatmap_img, caption="Toxicity Attribution Heatmap", width=600)
            st.markdown("</div>", unsafe_allow_html=True)
        
        with col2:
            # Interpretation guide with native Streamlit
            st.markdown("### 🔍 Interpretation Guide")
            
            # Color legends
            col_a, col_b = st.columns([1, 3])
            with col_a:
                st.markdown('<div style="width:40px;height:40px;background:linear-gradient(135deg,#ff6b6b,#ee5a6f);border-radius:8px;"></div>', unsafe_allow_html=True)
            with col_b:
                st.markdown("**Red/Orange**  \nToxic contribution")
            
            col_c, col_d = st.columns([1, 3])
            with col_c:
                st.markdown('<div style="width:40px;height:40px;background:linear-gradient(135deg,#51cf66,#37b24d);border-radius:8px;"></div>', unsafe_allow_html=True)
            with col_d:
                st.markdown("**Green/Blue**  \nSafe contribution")
            
            col_e, col_f = st.columns([1, 3])
            with col_e:
                st.markdown('<div style="width:40px;height:40px;background:#ffffff;border:2px solid #ccc;border-radius:8px;"></div>', unsafe_allow_html=True)
            with col_f:
                st.markdown("**White**  \nNeutral atoms")
            
            st.markdown("---")
            st.markdown("#### Max Attribution Score")
            st.metric("", f"{max_weight:.4f}", help="Higher = stronger influence on prediction")
            st.caption("Higher weight = stronger influence on prediction")
            
            # Download button
            buf = io.BytesIO()
            heatmap_img.save(buf, format='PNG')
            buf.seek(0)
            
            st.download_button(
                label="Download Heatmap",
                data=buf,
                file_name=f"toxpred_heatmap_{smiles[:20]}.png",
                mime="image/png",
                use_container_width=True
            )

elif page == "Batch Analysis":
    st.markdown("<h2 class='section-header'>&#x1F4CA; Batch Molecule Analysis</h2>", unsafe_allow_html=True)
    
    st.markdown("""
    Upload a CSV file with SMILES strings to analyze multiple molecules at once.
    
    **Required format:** CSV with a column named `smiles`
    """)
    
    uploaded_file = st.file_uploader("Upload CSV file", type=['csv'])
    
    if uploaded_file is not None:
        df = pd.read_csv(uploaded_file)
        
        if 'smiles' not in df.columns:
            st.error("CSV must contain a 'smiles' column!")
        else:
            st.success(f"Loaded {len(df)} molecules")
            st.dataframe(df.head())
            
            if st.button("Analyze All Molecules", type="primary"):
                with st.spinner(f"Analyzing {len(df)} molecules..."):
                    results = []
                    progress_bar = st.progress(0)
                    
                    for idx, row in df.iterrows():
                        smiles = row['smiles']
                        mol = Chem.MolFromSmiles(smiles)
                        
                        if mol is not None:
                            fp = get_morgan_fingerprint(mol)
                            proba = model.predict_proba(fp.reshape(1, -1))[0][1]
                            pred = "TOXIC" if proba > 0.5 else "SAFE"
                            
                            results.append({
                                'SMILES': smiles,
                                'Prediction': pred,
                                'Toxicity_Probability': proba,
                                'Confidence': max(proba, 1 - proba)
                            })
                        else:
                            results.append({
                                'SMILES': smiles,
                                'Prediction': 'INVALID',
                                'Toxicity_Probability': np.nan,
                                'Confidence': np.nan
                            })
                        
                        progress_bar.progress((idx + 1) / len(df))
                    
                    results_df = pd.DataFrame(results)
                    
                    st.markdown("### Results Summary")
                    
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        toxic_count = (results_df['Prediction'] == 'TOXIC').sum()
                        st.metric("Toxic Molecules", toxic_count)
                    with col2:
                        safe_count = (results_df['Prediction'] == 'SAFE').sum()
                        st.metric("Safe Molecules", safe_count)
                    with col3:
                        invalid_count = (results_df['Prediction'] == 'INVALID').sum()
                        st.metric("Invalid SMILES", invalid_count)
                    
                    st.markdown("### Detailed Results")
                    st.dataframe(results_df, use_container_width=True)
                    
                    # Download results
                    csv = results_df.to_csv(index=False)
                    st.download_button(
                        label="Download Results CSV",
                        data=csv,
                        file_name="toxpred_batch_results.csv",
                        mime="text/csv",
                        use_container_width=True
                    )

elif page == "About":
    st.markdown("<h2 class='section-header'>&#x2139;&#xFE0F; About ToxPred-Explainable</h2>", unsafe_allow_html=True)
    
    # Overview
    st.markdown("""
    <div class='info-card'>
        <h3 style='margin-top: 0;'>🎯 Mission</h3>
        <p style='font-size: 1.05rem; line-height: 1.6;'>
        Transform molecular toxicity prediction from a black-box process into an <strong>interpretable, 
        trustworthy system</strong> where scientists can understand exactly which molecular features 
        drive toxicity predictions.
        </p>
    </div>
    """, unsafe_allow_html=True)
    
    # Technical Details
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("### 🧠 Model Architecture")
        st.markdown("""
        **Algorithm:** Random Forest Classifier
        - 100 decision trees (ensemble learning)
        - Class-balanced weights (handles imbalanced data)
        - Out-of-bag validation for robust estimates
        
        **Feature Engineering:**
        - Morgan Fingerprints (ECFP4)
        - 2048-bit circular fingerprints
        - Radius = 2 (captures substructures up to 4 bonds)
        - Encodes molecular connectivity and topology
        
        **Why Random Forest?**
        - Handles high-dimensional data (2048 features)
        - Naturally provides feature importance
        - Resistant to overfitting
        - Fast inference for real-time predictions
        """)
        
        st.markdown("### 📊 Dataset: Tox21 Challenge")
        st.markdown("""
        **Source:** EPA/FDA Toxicology in the 21st Century
        
        **SR-ARE Assay Specifics:**
        - Total molecules: 7,831
        - Training samples: 5,832 (after QC)
        - Positive (toxic): 945 compounds (16.2%)
        - Negative (safe): 4,887 compounds (83.8%)
        
        **Data Quality:**
        - Experimental validation from HTS screening
        - Multiple concentration tests
        - Quality control filters applied
        - Reproducible measurements
        """)
    
    with col2:
        st.markdown("### 🎯 Performance Metrics")
        
        # Performance table
        perf_data = {
            "Metric": ["Training Accuracy", "Test Accuracy", "ROC-AUC", "Precision", "Recall", "F1-Score"],
            "Value": ["99.85%", "86.63%", "0.822", "High", "Balanced", "Strong"],
            "Interpretation": [
                "Model learns patterns well",
                "Good generalization",
                "Strong discrimination ability",
                "Few false positives",
                "Catches most toxic compounds",
                "Overall robust performance"
            ]
        }
        st.dataframe(perf_data, hide_index=True, use_container_width=True)
        
        st.markdown("""
        **What makes this model production-ready?**
        
        ✅ **Accuracy:** 86.6% correct predictions on unseen data  
        ✅ **ROC-AUC 0.822:** Strong ability to distinguish toxic from safe  
        ✅ **Balanced:** Handles imbalanced data effectively  
        ✅ **Explainable:** Atom-level attribution for transparency  
        """)
        
        st.markdown("### 🔬 SR-ARE Assay Explained")
        st.markdown("""
        **Biological Mechanism:**
        
        The Stress Response - Antioxidant Response Element (SR-ARE) pathway:
        
        1. **Normal State:** Nrf2 protein stays in cytoplasm
        2. **Stress Activation:** Toxic compounds cause oxidative stress
        3. **Response:** Nrf2 moves to nucleus, binds to ARE
        4. **Result:** Antioxidant genes activate (defense mechanism)
        
        **Clinical Relevance:**
        - Predicts liver toxicity (hepatotoxicity)
        - Indicates oxidative damage potential
        - Correlates with drug-induced organ damage
        
        **In Drug Development:**
        - 🔴 High SR-ARE = early warning signal
        - 🟡 Medium SR-ARE = requires investigation  
        - 🟢 Low SR-ARE = safer toxicity profile
        """)
    
    # Use Cases
    st.markdown("### 🚀 Real-World Applications")
    
    use_case_cols = st.columns(3)
    
    with use_case_cols[0]:
        st.markdown("""
        **💊 Drug Discovery**
        - Virtual screening of compound libraries
        - Prioritize safe candidates early
        - Reduce late-stage failures
        - Save millions in development costs
        """)
    
    with use_case_cols[1]:
        st.markdown("""
        **🔬 Medicinal Chemistry**
        - Identify toxic substructures
        - Guide structure optimization
        - Balance potency vs. safety
        - Rational drug design decisions
        """)
    
    with use_case_cols[2]:
        st.markdown("""
        **📚 Research & Education**
        - Understand structure-toxicity relationships
        - Teach explainable AI concepts
        - Publish reproducible science
        - Regulatory documentation
        """)
    
    # Technical Stack
    st.markdown("---")
    st.markdown("### 🛠️ Technology Stack")
    
    tech_cols = st.columns(4)
    with tech_cols[0]:
        st.markdown("**RDKit 2025.9**  \nCheminformatics")
    with tech_cols[1]:
        st.markdown("**Scikit-Learn**  \nMachine Learning")
    with tech_cols[2]:
        st.markdown("**Streamlit 1.41**  \nWeb Interface")
    with tech_cols[3]:
        st.markdown("**Python 3.13**  \nCore Language")
    
    st.markdown("---")
    st.markdown("""
    <div style='text-align: center; color: #666; padding: 20px;'>
        <p><strong>ToxPred-Explainable</strong> | Built with ❤️ for the scientific community</p>
        <p style='font-size: 0.9rem;'>⚠️ For research purposes only. Not for clinical or regulatory decisions without proper validation.</p>
    </div>
    """, unsafe_allow_html=True)

# Footer
st.markdown("---")
st.markdown("""
<div style='text-align: center; color: #666; padding: 20px;'>
    <p><strong>ToxPred-Explainable</strong> | Production-Ready Interpretable ML for Drug Safety</p>
    <p>&#x26A0;&#xFE0F; For research purposes only. Not for clinical or regulatory decisions.</p>
</div>
""", unsafe_allow_html=True)
