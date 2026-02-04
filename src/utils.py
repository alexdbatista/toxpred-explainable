"""
Utility Functions for ToxPred-Explainable
Molecular featurization and explainability functions
"""

import io

import numpy as np
import requests
from PIL import Image
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.Draw import SimilarityMaps, rdMolDraw2D


def validate_smiles(smiles: str) -> bool:
    """
    Validate if a SMILES string is valid.
    
    Args:
        smiles: SMILES string to validate
        
    Returns:
        True if valid, False otherwise
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None
    except:
        return False


def name_to_smiles(chemical_name: str) -> tuple[str, str]:
    """
    Convert chemical name to SMILES using PubChem API.
    
    Args:
        chemical_name: Common or IUPAC name of chemical
        
    Returns:
        Tuple of (SMILES string, status message)
    """
    try:
        # Clean the input
        chemical_name = chemical_name.strip()
        
        # PubChem REST API endpoint - try IsomericSMILES first, fallback to CanonicalSMILES
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{chemical_name}/property/IsomericSMILES,CanonicalSMILES/JSON"
        
        response = requests.get(url, timeout=10)
        
        if response.status_code == 200:
            data = response.json()
            
            # Check if the response has the expected structure
            if 'PropertyTable' in data and 'Properties' in data['PropertyTable']:
                properties = data['PropertyTable']['Properties']
                if len(properties) > 0:
                    # Try IsomericSMILES first (more specific), then CanonicalSMILES, then any SMILES field
                    smiles = (properties[0].get('IsomericSMILES') or 
                             properties[0].get('CanonicalSMILES') or
                             properties[0].get('ConnectivitySMILES'))
                    
                    if smiles:
                        return smiles, "success"
                    else:
                        return None, "SMILES not found in PubChem response"
                else:
                    return None, "No properties returned from PubChem"
            else:
                return None, "Unexpected PubChem response format"
                
        elif response.status_code == 404:
            return None, "Chemical name not found in PubChem database"
        else:
            return None, f"PubChem API error: {response.status_code}"
            
    except requests.exceptions.Timeout:
        return None, "Request timeout - please try again"
    except requests.exceptions.ConnectionError:
        return None, "Connection error - check your internet connection"
    except KeyError as e:
        return None, f"Data parsing error: missing key {str(e)}"
    except Exception as e:
        return None, f"Error: {str(e)}"


def get_morgan_fingerprint(mol, radius=2, nBits=2048):
    """
    Generate Morgan fingerprint for a molecule.
    
    Args:
        mol: RDKit molecule object
        radius: Fingerprint radius (default: 2)
        nBits: Number of bits (default: 2048)
        
    Returns:
        numpy array of fingerprint
    """
    if mol is None:
        return np.zeros(nBits)
    
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits)
    return np.array(fp)


def calculate_lipinski(mol):
    """
    Calculate Lipinski's Rule of Five properties.
    
    Args:
        mol: RDKit molecule object
        
    Returns:
        dict with MW, LogP, HBD, HBA, and Passes flag
    """
    if mol is None:
        return {
            'MW': 0,
            'LogP': 0,
            'HBD': 0,
            'HBA': 0,
            'Passes': False
        }
    
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    
    # Lipinski's Rule: MW<500, LogP<5, HBD<5, HBA<10
    passes = (mw < 500) and (logp < 5) and (hbd < 5) and (hba < 10)
    
    return {
        'MW': mw,
        'LogP': logp,
        'HBD': hbd,
        'HBA': hba,
        'Passes': passes
    }


def explain_molecule(mol, model, radius=2, nBits=2048):
    """
    Generate explainability heatmap for a molecule.
    
    Args:
        mol: RDKit molecule object
        model: Trained sklearn model
        radius: Fingerprint radius
        nBits: Number of bits in fingerprint
        
    Returns:
        tuple: (PIL Image of heatmap, max_weight)
    """
    if mol is None:
        raise ValueError("Invalid molecule")
    
    # Define fingerprint function for explainability
    def fp_func(mol, atomId=-1):
        return SimilarityMaps.GetMorganFingerprint(mol, atomId, radius=radius, nBits=nBits)
    
    # Create high-quality drawing object
    draw2d = rdMolDraw2D.MolDraw2DCairo(800, 800)
    draw2d.drawOptions().bondLineWidth = 2
    
    # Generate heatmap
    fig, maxweight = SimilarityMaps.GetSimilarityMapForModel(
        mol,
        fp_func,
        lambda fp: model.predict_proba([fp])[0][1],  # Probability of toxic class
        draw2d=draw2d
    )
    
    # Convert to PIL Image
    draw2d.FinishDrawing()
    img_data = draw2d.GetDrawingText()
    img = Image.open(io.BytesIO(img_data))
    
    return img, maxweight


def smiles_to_image(smiles, size=(300, 300)):
    """
    Convert SMILES to molecule image.
    
    Args:
        smiles: SMILES string
        size: Image size tuple
        
    Returns:
        PIL Image
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    from rdkit.Chem import Draw
    img = Draw.MolToImage(mol, size=size)
    return img


def batch_predict(smiles_list, model):
    """
    Predict toxicity for a batch of SMILES strings.
    
    Args:
        smiles_list: List of SMILES strings
        model: Trained sklearn model
        
    Returns:
        list of dicts with predictions
    """
    results = []
    
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        
        if mol is not None:
            fp = get_morgan_fingerprint(mol)
            proba = model.predict_proba(fp.reshape(1, -1))[0][1]
            pred = "TOXIC" if proba > 0.5 else "SAFE"
            
            results.append({
                'smiles': smiles,
                'prediction': pred,
                'probability': proba,
                'confidence': max(proba, 1 - proba)
            })
        else:
            results.append({
                'smiles': smiles,
                'prediction': 'INVALID',
                'probability': np.nan,
                'confidence': np.nan
            })
    
    return results


def get_toxic_substructures(mol, model, threshold=0.3):
    """
    Identify toxic substructures in a molecule.
    
    Args:
        mol: RDKit molecule object
        model: Trained model
        threshold: Attribution threshold for toxicity
        
    Returns:
        list of atom indices with high toxic attribution
    """
    if mol is None:
        return []
    
    # Get atom-level attributions
    def fp_func(mol, atomId=-1):
        return SimilarityMaps.GetMorganFingerprint(mol, atomId, radius=2, nBits=2048)
    
    # Calculate per-atom contributions
    weights = []
    for atomId in range(mol.GetNumAtoms()):
        fp = fp_func(mol, atomId)
        contrib = model.predict_proba([fp])[0][1]
        weights.append(contrib)
    
    # Get atoms above threshold
    toxic_atoms = [i for i, w in enumerate(weights) if w > threshold]
    
    return toxic_atoms
