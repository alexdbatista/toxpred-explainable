"""
Configuration file for ToxPred-Explainable
"""

import os

# Paths
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODEL_PATH = os.path.join(BASE_DIR, 'models', 'rf_toxpred_sr_are.pkl')
DATA_PATH = os.path.join(BASE_DIR, 'data', 'tox21_data.csv.gz')

# Model parameters
MODEL_PARAMS = {
    'n_estimators': 100,
    'class_weight': 'balanced',
    'random_state': 42,
    'n_jobs': -1
}

# Fingerprint parameters
FP_PARAMS = {
    'radius': 2,
    'nBits': 2048
}

# Assay information
ASSAY_INFO = {
    'name': 'SR-ARE (Stress Response - Antioxidant Response Element)',
    'description': 'Measures cellular stress response, a key indicator of toxicity',
    'target_col': 'SR-ARE'
}

# App settings
APP_SETTINGS = {
    'title': 'ToxPred-Explainable',
    'icon': '🧪',
    'layout': 'wide',
    'theme': 'light'
}

# Lipinski's Rule thresholds
LIPINSKI_THRESHOLDS = {
    'MW': 500,
    'LogP': 5,
    'HBD': 5,
    'HBA': 10
}
