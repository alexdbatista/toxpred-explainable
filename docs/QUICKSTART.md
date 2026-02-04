# Quick Start Guide

## Installation (5 minutes)

```bash
# 1. Clone repo
git clone https://github.com/yourusername/toxpred-explainable.git
cd toxpred-explainable

# 2. Create environment
conda create -n toxpred python=3.13
conda activate toxpred

# 3. Install dependencies
conda install -c conda-forge rdkit
pip install -r requirements.txt

# 4. Train model
cd src && python train_model.py && cd ..

# 5. Run app
streamlit run app.py
```

## Usage Examples

### Example 1: Aspirin (Safe)
```
SMILES: CC(=O)Oc1ccccc1C(=O)O
Expected: SAFE (low probability)
Key features: Ester, carboxylic acid groups
```

### Example 2: Benzene (Toxic)
```
SMILES: c1ccccc1
Expected: TOXIC (higher probability)
Key features: Aromatic ring may show some toxic attribution
```

### Example 3: Ethanol (Safe)
```
SMILES: CCO
Expected: SAFE (very low probability)
Key features: Simple alcohol, drug-like
```

## Troubleshooting

**Issue**: RDKit import error
**Fix**: Install via conda: `conda install -c conda-forge rdkit`

**Issue**: Model file not found
**Fix**: Run training script: `cd src && python train_model.py`

**Issue**: Port already in use
**Fix**: Specify different port: `streamlit run app.py --server.port 8502`
