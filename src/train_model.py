"""
Model Training Script for ToxPred-Explainable
Train Random Forest on Tox21 SR-ARE data
"""

import os
import pickle
import urllib.request

import numpy as np
import pandas as pd
from config import ASSAY_INFO, DATA_PATH, FP_PARAMS, MODEL_PARAMS, MODEL_PATH, BBB_MODEL_PATH
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, classification_report, roc_auc_score
from sklearn.model_selection import train_test_split


def download_data():
    """Download Tox21 dataset if not present."""
    if not os.path.exists(DATA_PATH):
        print("📥 Downloading Tox21 dataset...")
        os.makedirs(os.path.dirname(DATA_PATH), exist_ok=True)
        
        url = "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/tox21.csv.gz"
        urllib.request.urlretrieve(url, DATA_PATH)
        print("✅ Dataset downloaded!")
    else:
        print("✅ Dataset already exists!")


def get_fingerprint_arr(mol):
    """Convert molecule to Morgan fingerprint array."""
    if mol is None:
        return np.zeros(FP_PARAMS['nBits'])
    return np.array(
        AllChem.GetMorganFingerprintAsBitVect(
            mol, 
            FP_PARAMS['radius'], 
            nBits=FP_PARAMS['nBits']
        )
    )


def load_and_clean_data():
    """Load and clean Tox21 data."""
    print(f"\n🔄 Loading Tox21 dataset...")
    
    # Load data
    tox21_df = pd.read_csv(DATA_PATH, compression='gzip')
    print(f"   - Loaded {len(tox21_df)} molecules")
    
    # Focus on target assay
    target_col = ASSAY_INFO['target_col']
    clean_df = tox21_df.dropna(subset=[target_col]).copy()
    
    print(f"\n✅ Data cleaned for {target_col}!")
    print(f"   - Valid samples: {len(clean_df)}")
    print(f"   - Toxic (1): {int(clean_df[target_col].sum())}")
    print(f"   - Safe (0): {int(len(clean_df) - clean_df[target_col].sum())}")
    
    return clean_df, target_col


def featurize_data(clean_df, target_col):
    """Convert SMILES to fingerprints."""
    print("\n🔄 Generating fingerprints for all molecules...")
    print("   (This may take 1-2 minutes)")
    
    X_smiles = clean_df['smiles'].values
    y = clean_df[target_col].values
    
    # Convert SMILES to molecules
    mols = [Chem.MolFromSmiles(s) for s in X_smiles]
    
    # Convert to fingerprints
    X = np.array([get_fingerprint_arr(m) for m in mols])
    
    print(f"✅ Featurization complete!")
    print(f"   - Feature matrix shape: {X.shape}")
    
    return X, y


def train_model(X, y):
    """Train Random Forest model."""
    print("\n🚀 Training Random Forest model...")
    
    # Split data
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=y
    )
    
    print(f"   - Training samples: {len(y_train)}")
    print(f"   - Test samples: {len(y_test)}")
    
    # Train model
    model = RandomForestClassifier(**MODEL_PARAMS)
    model.fit(X_train, y_train)
    
    # Evaluate
    train_acc = accuracy_score(y_train, model.predict(X_train))
    test_acc = accuracy_score(y_test, model.predict(X_test))
    
    y_pred_proba = model.predict_proba(X_test)[:, 1]
    test_roc_auc = roc_auc_score(y_test, y_pred_proba)
    
    print("\n✅ Model Training Complete!")
    print(f"\n📊 Performance Metrics:")
    print(f"   - Training Accuracy: {train_acc:.2%}")
    print(f"   - Test Accuracy: {test_acc:.2%}")
    print(f"   - Test ROC-AUC: {test_roc_auc:.3f}")
    
    # Classification report
    print(f"\n📋 Classification Report:")
    print(classification_report(y_test, model.predict(X_test), 
                                target_names=['Safe', 'Toxic']))
    
    return model


def save_model(model):
    """Save trained model."""
    os.makedirs(os.path.dirname(MODEL_PATH), exist_ok=True)
    
    with open(MODEL_PATH, 'wb') as f:
        pickle.dump(model, f)
    
    print(f"\n💾 Model saved to: {MODEL_PATH}")


def train_and_save_model():
    """Complete training pipeline - callable from app.py"""
    print("=" * 70)
    print("🧪 ToxPred-Explainable Model Training")
    print("=" * 70)
    
    # Download data
    download_data()
    
    # Load and clean
    clean_df, target_col = load_and_clean_data()
    
    # Featurize
    X, y = featurize_data(clean_df, target_col)
    
    # Train
    model = train_model(X, y)
    
    # Save
    save_model(model)
    
    print("\n" + "=" * 70)
    print("🎉 Training Complete!")
    print("=" * 70)
    return model


def download_bbb_data():
    """Download BBBP dataset if not present."""
    bbb_data_path = os.path.join(os.path.dirname(DATA_PATH), 'bbbp_data.csv')
    if not os.path.exists(bbb_data_path):
        print("📥 Downloading BBBP dataset...")
        os.makedirs(os.path.dirname(bbb_data_path), exist_ok=True)
        
        url = "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/BBBP.csv"
        urllib.request.urlretrieve(url, bbb_data_path)
        print("✅ BBBP dataset downloaded!")
    else:
        print("✅ BBBP dataset already exists!")
    return bbb_data_path


def train_and_save_bbb_model():
    """Train and save BBB permeability model."""
    print("=" * 70)
    print("🧠 Training Blood-Brain Barrier Permeability Model")
    print("=" * 70)
    
    # Download data
    bbb_data_path = download_bbb_data()
    
    # Load and clean data
    print("\n🔄 Loading BBBP dataset...")
    bbb_df = pd.read_csv(bbb_data_path)
    print(f"   - Loaded {len(bbb_df)} molecules")
    
    # Clean data - keep only valid SMILES
    bbb_df = bbb_df.dropna(subset=['smiles', 'p_np']).copy()
    print(f"   - Valid samples: {len(bbb_df)}")
    print(f"   - BBB+ (permeable): {int(bbb_df['p_np'].sum())}")
    print(f"   - BBB- (non-permeable): {int(len(bbb_df) - bbb_df['p_np'].sum())}")
    
    # Featurize
    print("\n🔄 Generating fingerprints...")
    X_smiles = bbb_df['smiles'].values
    y = bbb_df['p_np'].values
    
    mols = [Chem.MolFromSmiles(s) for s in X_smiles]
    X = np.array([get_fingerprint_arr(m) for m in mols])
    
    print(f"✅ Featurization complete! Shape: {X.shape}")
    
    # Train
    print("\n🚀 Training Random Forest model...")
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=y
    )
    
    model = RandomForestClassifier(**MODEL_PARAMS)
    model.fit(X_train, y_train)
    
    # Evaluate
    test_acc = accuracy_score(y_test, model.predict(X_test))
    y_pred_proba = model.predict_proba(X_test)[:, 1]
    test_roc_auc = roc_auc_score(y_test, y_pred_proba)
    
    print(f"\n✅ BBB Model Training Complete!")
    print(f"   - Test Accuracy: {test_acc:.2%}")
    print(f"   - Test ROC-AUC: {test_roc_auc:.3f}")
    
    # Save
    os.makedirs(os.path.dirname(BBB_MODEL_PATH), exist_ok=True)
    with open(BBB_MODEL_PATH, 'wb') as f:
        pickle.dump(model, f)
    
    print(f"\n💾 BBB Model saved to: {BBB_MODEL_PATH}")
    print("=" * 70)
    
    return model


def main():
    """Main training pipeline."""
    train_and_save_model()
    
    print("\n🚀 Next steps:")
    print("   1. Run the Streamlit app: streamlit run app.py")
    print("   2. Test explainability on molecules")
    print("   3. Export heatmaps for your portfolio")
    print("\n")


if __name__ == "__main__":
    main()
