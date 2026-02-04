# 🧪 ToxPred-Explainable

**Explainable AI for Molecular Toxicity Prediction**

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://your-app-url-here.streamlit.app)
[![Python 3.13+](https://img.shields.io/badge/python-3.13+-blue.svg)](https://www.python.org/downloads/)
[![Streamlit](https://img.shields.io/badge/streamlit-1.41+-FF4B4B.svg)](https://streamlit.io)
[![RDKit](https://img.shields.io/badge/rdkit-2023.9+-brightgreen.svg)](https://www.rdkit.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A production-ready web application for molecular toxicity prediction with atom-level explainability. Search by chemical name or SMILES, get instant predictions, and visualize which atoms contribute to toxicity.

**🎥 [Live Demo](https://your-app-url-here.streamlit.app)** | **📊 [Model Card](docs/model_card.md)** | **📖 [Technical Report](docs/technical_report.md)**

---

## ✨ Key Features

- 🔍 **Dual Search**: Chemical name (PubChem API) OR SMILES structure input
- 🎯 **Instant Predictions**: Random Forest classifier with 86.6% test accuracy
- 🗺️ **Explainability Heatmaps**: Atom-level attribution showing toxic substructures
- 📊 **Drug-Likeness**: Lipinski's Rule of Five validation (MW, LogP, HBD, HBA)
- 📁 **Batch Processing**: Upload CSV files for bulk predictions
- 💾 **Export Ready**: Download heatmaps and predictions for reports
- 📚 **Educational Content**: Comprehensive explanations of SR-ARE assay and predictions

---

## 🚀 Quick Start

### Try the Live App
Visit **[your-app-url-here.streamlit.app](https://your-app-url-here.streamlit.app)** to try it instantly!

### Run Locally

```bash
# 1. Clone the repository
git clone https://github.com/yourusername/toxpred-explainable.git
cd toxpred-explainable

# 2. Create conda environment (RDKit requires conda)
conda create -n toxpred python=3.13
conda activate toxpred
conda install -c conda-forge rdkit

# 3. Install dependencies
pip install -r requirements.txt

# 4. Launch app
streamlit run app.py
```

The app will open at `http://localhost:8501`

---

## 📊 Model Performance

**Dataset**: EPA/FDA Tox21 Challenge - SR-ARE Assay (Stress Response)
- **Training**: 4,665 molecules (16.2% toxic, 83.8% safe)
- **Testing**: 1,167 molecules

**Metrics**:
| Metric | Training | Testing |
|--------|----------|---------|
| Accuracy | 99.85% | 86.63% |
| Precision | 99.86% | 55.68% |
| Recall | 99.07% | 42.86% |
| F1-Score | 99.46% | 48.44% |
| ROC-AUC | 99.99% | 0.822 |

**Architecture**: Random Forest (100 trees, class-balanced)  
**Features**: Morgan Fingerprints (2048-bit ECFP4, radius=2)

---

## 🛠️ Tech Stack

| Component | Technology | Purpose |
|-----------|------------|---------|
| **ML Framework** | scikit-learn 1.3+ | Random Forest classifier |
| **Cheminformatics** | RDKit 2023.9+ | Molecular processing, fingerprints, visualization |
| **Web Framework** | Streamlit 1.28+ | Interactive web application |
| **API Integration** | PubChem REST API | Chemical name → SMILES conversion |
| **Data Processing** | pandas 2.0+, NumPy 1.24+ | Data manipulation |
| **Visualization** | Matplotlib 3.7+, PIL 10.0+ | Heatmaps and molecular rendering |

---

## 🎮 How to Use

### 1. Search by Chemical Name
```
Input: "aspirin" → System fetches SMILES from PubChem → Prediction
```

### 2. Search by SMILES
```
Input: "CC(=O)Oc1ccccc1C(=O)O" → Direct prediction
```

### 3. Try Examples
Click example molecules (Aspirin, Caffeine, etc.) to auto-populate input

### 4. Interpret Results
- **Prediction**: Safe/Toxic with probability
- **Heatmap**: Red atoms = high toxicity contribution
- **Lipinski**: Drug-likeness validation
- **Download**: Export heatmap as PNG

### 5. Batch Analysis
Upload CSV with `smiles` column → Download predictions

**CSV Format:**
```csv
smiles
CCO
CC(=O)Oc1ccccc1C(=O)O
CN1C=NC2=C1C(=O)N(C(=O)N2C)C
```

---

## 📁 Project Structure

```
toxpred-explainable/
│
├── app.py                        # Main Streamlit application
├── requirements.txt              # Python dependencies (flexible constraints)
├── README.md                     # Project documentation
├── .gitignore                    # Git ignore patterns
│
├── .streamlit/
│   └── config.toml              # Streamlit theme configuration
│
├── src/
│   ├── __init__.py
│   ├── config.py                # Paths and constants
│   ├── utils.py                 # Molecular processing, PubChem API, explainability
│   └── train_model.py           # Model training pipeline
│
├── models/
│   └── rf_toxpred_sr_are.pkl    # Trained Random Forest (86.6% test accuracy)
│
├── data/
│   └── tox21_data.csv.gz        # Tox21 SR-ARE dataset (auto-downloaded)
│
├── notebooks/
│   └── toxpred_explainability_prototype.ipynb  # Development notebook
│
└── docs/
    ├── model_card.md            # Model documentation
    └── technical_report.md      # Detailed analysis
```

---

## 🧬 Dataset Details

**Source**: EPA/FDA Tox21 Challenge  
**Assay**: SR-ARE (Stress Response - Antioxidant Response Element)  
**Biological Target**: Nuclear factor erythroid 2-related factor 2 (Nrf2) pathway

**What it measures**: Cellular stress response activation. Positive results indicate compounds that trigger oxidative stress defense mechanisms, which can signal:
- Toxicity concerns (excessive stress response)
- Potential therapeutic effects (moderate activation)

**Data Quality**:
- 7,831 total molecules screened
- 5,832 with valid SR-ARE measurements (1,999 missing data excluded)
- High-throughput screening (HTS) with quality control
- Binary classification: Active (toxic) vs Inactive (safe)

---

## 🤖 Model Architecture

### Why Random Forest?

1. **Handles Imbalance**: Class weighting (16.2% toxic) without SMOTE
2. **Feature Importance**: Built-in explainability via Gini importance
3. **No Feature Scaling**: Works directly with binary fingerprints
4. **Production-Ready**: Fast predictions (~10ms per molecule)
5. **Robust**: Less prone to overfitting than deep learning on small datasets

### Why Morgan Fingerprints?

- **Industry Standard**: Used in drug discovery pipelines
- **Structural Encoding**: Captures circular substructures (radius=2 = atoms + neighbors)
- **Fixed Size**: 2048 bits for consistent model input
- **Collision Resistance**: Low chance of different structures → same fingerprint

### Training Details

```python
RandomForestClassifier(
    n_estimators=100,
    class_weight='balanced',  # Handles 16.2% vs 83.8% imbalance
    random_state=42,
    n_jobs=-1
)
```

**Cross-validation**: 5-fold CV during development  
**Final model**: Trained on full 4,665-molecule training set

---

## 🎯 Use Cases

### 1. Drug Discovery
- **Early-stage screening**: Filter out toxic candidates before synthesis
- **Lead optimization**: Identify problematic substructures to modify
- **Cost savings**: Reduce expensive in vitro/in vivo testing

### 2. Medicinal Chemistry
- **Structure-activity relationships**: Understand toxicity drivers
- **Molecular design**: Avoid toxic motifs in new compounds
- **Patent analysis**: Assess safety profiles of competitive compounds

### 3. Regulatory Compliance
- **Pre-submission screening**: Check compounds before regulatory filing
- **Safety assessment**: Support ICH M7 mutagenicity evaluations
- **Risk prioritization**: Rank compounds for deeper toxicology studies

### 4. Education & Research
- **Teaching tool**: Demonstrate explainable AI in chemistry
- **Benchmarking**: Compare new models against established baseline
- **Hypothesis generation**: Explore oxidative stress mechanisms

---

## ⚠️ Limitations & Disclaimer

**This tool is for research purposes only. NOT for clinical or regulatory decisions.**

### Model Limitations
- ✋ **Single assay**: Only predicts SR-ARE toxicity (oxidative stress)
- ✋ **In vitro only**: Cell-based assay doesn't capture full organism effects
- ✋ **Limited chemical space**: Trained on Tox21 library (may not generalize to all drugs)
- ✋ **Imbalanced performance**: 55.7% precision, 42.9% recall on toxic class

### When NOT to use
- ❌ Drug approval decisions
- ❌ Patient treatment planning
- ❌ Regulatory submissions (without additional validation)
- ❌ Replacing experimental toxicology

### Recommended Next Steps
1. Validate high-risk predictions with in vitro assays
2. Test similar compounds (series) for consistency
3. Consult toxicology experts for interpretation
4. Consider multi-assay panels (hERG, CYP450, AMES, etc.)

---

## 🚀 Deployment

### Streamlit Cloud

1. **Push to GitHub**:
```bash
git init
git add .
git commit -m "Initial commit: ToxPred-Explainable"
git remote add origin https://github.com/yourusername/toxpred-explainable.git
git push -u origin main
```

2. **Deploy**:
- Visit [share.streamlit.io](https://share.streamlit.io)
- Connect GitHub repository
- Set main file: `app.py`
- Click "Deploy"

3. **Update README badge** with your Streamlit Cloud URL

### Local Deployment

```bash
# Production mode with gunicorn (optional)
pip install streamlit
streamlit run app.py --server.port 8501 --server.headless true
```

---

## 🤝 Contributing

Contributions welcome! Areas for improvement:
- [ ] Multi-assay predictions (hERG, AMES, CYP450)
- [ ] Uncertainty quantification (conformal prediction)
- [ ] Attention-based explainability (GNNExplainer)
- [ ] Active learning for data-efficient retraining
- [ ] API endpoint for programmatic access

---

## 📄 License

MIT License - see [LICENSE](LICENSE) file

---

## 📧 Contact

**Author**: Alex Domingues Batista  
**Portfolio**: [github.com/yourusername](https://github.com/yourusername)  
**Email**: your.email@example.com

---

## 🙏 Acknowledgments

- **Tox21 Consortium**: EPA, NIH, FDA for open dataset
- **RDKit Community**: Cheminformatics tools
- **Streamlit Team**: Web framework
- **PubChem**: Chemical structure database

---

## 📚 References

1. Tox21 Challenge: https://tripod.nih.gov/tox21/challenge/
2. Morgan Fingerprints: Rogers & Hahn (2010) J. Chem. Inf. Model.
3. Explainable AI in Drug Discovery: Jiménez-Luna et al. (2020) Nat. Mach. Intell.
4. Lipinski's Rule of Five: Lipinski et al. (1997) Adv. Drug Deliv. Rev.

---

**⭐ If you found this project helpful, please give it a star!**
│
├── notebooks/
│   └── explainability_engine_prototype.ipynb  # Development notebook
│
└── docs/
    └── demo_screenshots/      # App screenshots for README
```

---

## 🧠 How It Works

### 1. Featurization
- Input: SMILES string (e.g., `CCO`)
- Convert to RDKit molecule object
- Generate **Morgan Fingerprint** (ECFP4): 2048-bit circular fingerprint encoding molecular structure

### 2. Prediction
- Random Forest classifier (100 trees)
- Class-weighted to handle imbalanced data (16% toxic)
- Output: Probability of toxicity (0-1)

### 3. Explainability
- For each atom: perturb fingerprint → measure change in prediction
- High contribution → red color (toxic)
- Low contribution → blue color (safe)
- Uses RDKit's `SimilarityMaps.GetSimilarityMapForModel`

### 4. Lipinski's Rule
Drug-likeness filter used in pharma:
- ✅ Molecular Weight < 500 Da
- ✅ LogP < 5
- ✅ H-bond donors < 5
- ✅ H-bond acceptors < 10

---

## 📈 Example Results

### Toxic Molecule: TOX3021
**SMILES**: `CCOc1ccc2nc(S(N)(=O)=O)sc2c1`

| Property | Value |
|----------|-------|
| Prediction | TOXIC |
| Probability | 0.38 (38%) |
| Lipinski Pass | ✅ Yes |
| Toxic Substructure | Sulfonamide group (S(N)(=O)=O) |

**Heatmap**: Red highlighting on sulfonamide group indicates this functional group drives the toxicity prediction.

### Safe Molecule: Ethanol
**SMILES**: `CCO`

| Property | Value |
|----------|-------|
| Prediction | SAFE |
| Probability | 0.08 (8%) |
| Lipinski Pass | ✅ Yes |

---

## 🎓 Educational Value

This project demonstrates:

1. **Explainable AI**: Moving beyond black-box models to interpretable predictions
2. **Domain Knowledge**: Understanding cheminformatics, drug-likeness rules
3. **Software Engineering**: Modular code, configuration management, professional structure
4. **Data Science Pipeline**: Download → Clean → Featurize → Train → Deploy
5. **Production Readiness**: Web app, batch processing, exports, error handling

---

## 🔮 Future Enhancements

- [ ] Multi-assay prediction (all 12 Tox21 assays)
- [ ] Deep learning models (Graph Neural Networks)
- [ ] Molecular property optimization (suggest modifications)
- [ ] Integration with PubChem/ChEMBL databases
- [ ] Docker containerization for easy deployment
- [ ] REST API for programmatic access

---

## 📚 References

- **Tox21 Dataset**: [EPA/FDA Tox21 Challenge](https://tripod.nih.gov/tox21/challenge/)
- **RDKit**: [Open-source cheminformatics toolkit](https://www.rdkit.org/)
- **Morgan Fingerprints**: Rogers & Hahn, JCIM 2010
- **SimilarityMaps**: Riniker & Landrum, JCIM 2013

---

## 🤝 Contributing

This is a portfolio project, but suggestions welcome! Feel free to:
- Open issues for bugs/improvements
- Submit PRs for enhancements
- Use as a template for your own projects

---

## 📜 License

MIT License - feel free to use for learning/portfolio purposes.

---

## 👤 Author

**Alex Domingues Batista**
- LinkedIn: [Your LinkedIn](https://linkedin.com/in/yourprofile)
- Portfolio: [Your Portfolio Site](https://yourportfolio.com)
- GitHub: [@yourusername](https://github.com/yourusername)

---

## 🙏 Acknowledgments

- EPA/FDA for the Tox21 dataset
- RDKit development team
- Streamlit for the awesome web framework

---

*Built with ❤️ for drug discovery and explainable AI*
