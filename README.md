# Deep-AmPEP30: Deep Learning-Based Antimicrobial Peptide Prediction

![R](https://img.shields.io/badge/R-276DC3.svg?style=flat&logo=r&logoColor=white)
![Keras](https://img.shields.io/badge/Keras-FF0000.svg?style=flat&logo=keras&logoColor=white)
![TensorFlow](https://img.shields.io/badge/TensorFlow-FF6F00.svg?style=flat&logo=tensorflow&logoColor=white)

**Deep-AmPEP30** is a state-of-the-art machine learning system for predicting antimicrobial peptides (AMPs) from protein sequences. This system serves as the core prediction engine for the AxPEP web server, providing accurate and efficient AMP prediction for peptides ranging from 5 to 30 amino acids in length.

## 🌟 Key Features

- **Dual Model Architecture**: Implements both Convolutional Neural Network (CNN) and Random Forest (RF) models for robust predictions
- **Advanced Feature Engineering**: Utilizes PseKRAAC (Pseudo K-tuple Reduced Amino Acid Composition) and multiple amino acid classification schemes
- **High Performance**: Optimized for large-scale sequence processing with sub-second per sequence prediction times
- **Production Ready**: Designed for real-world applications with comprehensive error handling and validation
- **Batch Processing**: Supports FASTA format input for processing multiple sequences simultaneously

## 🛠️ Technical Architecture

### Machine Learning Models

#### CNN Model (Convolutional Neural Network)
- **Framework**: Keras/TensorFlow
- **Architecture**: 128+128 filters, 3×3 convolution kernels, 10 dense units
- **Training Time**: ~3 minutes (3,297 training sequences)
- **Prediction Speed**: ~22 seconds (10,000 sequences)

#### Random Forest Model
- **Configuration**: 1,200 decision trees
- **Training Time**: ~27 seconds
- **Prediction Speed**: ~23 seconds (10,000 sequences)

### Feature Extraction
- **Amino Acid Composition (AAC)**: Frequency analysis of 20 standard amino acids
- **PseKRAAC**: Pseudo K-tuple Reduced Amino Acid Composition
- **Multiple Classification Systems**: Type 8-17, Type 3A-19, Type 12-18, Type 7-15, Type 12-17
- **Optimized Feature Selection**: Best 5 feature combination strategy

## 📊 Performance Metrics

### Model Accuracy
- **UniProt Non-AMP Dataset**: 75.4% accuracy (24.6% error rate)
- **Positive Sample Error Rate**: 5.9%
- **Large-scale Prediction**: 38% sequences predicted as AMPs

### Computational Performance
- **Feature Generation**: 0.52±0.02 seconds (10,000 sequences)
- **CNN Training**: 197.5±0.8 seconds
- **RF Training**: 28.1±0.3 seconds
- **Prediction Speed**: ~1.7 seconds (10,000 sequences) for both models

## 📋 Requirements

### R Dependencies
```r
# Core packages
library(seqinr)      # Sequence analysis
library(keras)       # Deep learning
library(kerasR)      # Keras interface
library(caret)       # Machine learning
library(randomForest) # Random forest
library(ROCR)        # Performance evaluation
library(protr)       # Protein analysis
```

### System Requirements
- **R**: Version ≥ 3.6.0
- **Python**: Version ≥ 3.6 (for TensorFlow backend)
- **Memory**: ≥ 4GB RAM recommended
- **Storage**: ≥ 1GB free space for models and datasets

## 🚀 Installation

1. **Clone the repository**:
```bash
git clone https://github.com/your-username/Deep-AmPEP30.git
cd Deep-AmPEP30
```

2. **Install R dependencies**:
```r
# Install required packages
install.packages(c("seqinr", "caret", "randomForest", "ROCR", "protr"))

# Install Keras and TensorFlow
install.packages("keras")
library(keras)
install_keras()
```

3. **Initialize Keras**:
```r
library(kerasR)
keras_init()
```

## 📖 Usage

### Basic Prediction

#### Using CNN Model
```r
source("Deep-AmPEP30.R")

# Predict antimicrobial peptides using CNN
results <- cnn_yan_predict("your_sequences.fasta")
```

#### Using Random Forest Model
```r
# Predict using Random Forest
results <- rf_yan_predict("your_sequences.fasta")
```

### Training Custom Models

#### Train CNN Model
```r
# Train CNN model with default parameters
develop_cnn_mdl_yan(rf = FALSE)
```

#### Train Random Forest Model
```r
# Train Random Forest model
develop_cnn_mdl_yan(rf = TRUE)
```

### Feature Generation
```r
# Generate features for custom sequences
features <- totally_psekraac_best5(
    test_path = "your_sequences.fasta",
    test = TRUE,
    without_fasta_name = FALSE,
    check_len = TRUE
)
```

## 📁 Dataset Structure

```
dataset/
├── model_po.fasta          # Positive training samples (3,303 AMPs)
├── model_ne.fasta          # Negative training samples (3,293 non-AMPs)
├── train_po.fasta          # Additional positive training data
├── train_ne.fasta          # Additional negative training data
├── test_po.fasta           # Positive test samples
└── test_ne.fasta           # Negative test samples
```

### Data Format
Sequences should be provided in FASTA format:
```
>sequence_id_1
ACSAG
>sequence_id_2
FRWWHR
>sequence_id_3
RKKWFW
```

## 🔬 Model Validation

The system uses 10-fold cross-validation for model evaluation:

```r
# Cross-validation results are automatically generated during training
# Performance metrics include:
# - Accuracy
# - Sensitivity (True Positive Rate)
# - Specificity (True Negative Rate)
# - AUC (Area Under Curve)
```

## 📈 Performance Benchmarks

All performance benchmarks were conducted on sequences of length 5-30 amino acids:

| Operation | Time (seconds) | Sequences | Model |
|-----------|----------------|-----------|--------|
| Feature Generation | 0.52±0.02 | 10,000 | - |
| CNN Training | 197.5±0.8 | 3,297 | CNN |
| RF Training | 28.1±0.3 | 3,297 | RF |
| CNN Prediction | 1.73±0.06 | 10,000 | CNN |
| RF Prediction | 1.70±0.04 | 10,000 | RF |

## 🔧 Configuration

### Amino Acid Classification Systems
The system uses multiple classification schemes optimized for AMP prediction:

```r
# Available classification types
type8.cluster17 <- c('AT','C','DE','F','G','H','IV','K','L','M','N','P','Q','R','S','V','W')
type3a.cluster19 <- c('FA','P','G','S','T','D','E','Q','N','K','R','H','W','Y','M','L','I','V','C')
type12.cluster18 <- c('TVLI','M','F','W','Y','C','A','H','G','N','Q','P','R','K','S','T','D','E')
type7.cluster15 <- c('C','K','R','W','Y','A','FILV','M','D','E','Q','H','TP','GS','N')
type12.cluster17 <- c('TVLI','M','F','W','Y','C','A','H','G','N','Q','P','R','K','S','T','DE')
```

## 🧪 Applications

This system is valuable for various research and industrial applications:

- **Drug Discovery**: Screening novel antimicrobial peptide candidates
- **Proteomics**: Large-scale protein function annotation
- **Antibiotic Alternatives**: Identifying new antimicrobial therapies
- **Biotechnology**: Designing peptides with specific antimicrobial activities

## 📚 File Structure

```
Deep-AmPEP30/
├── Deep-AmPEP30.R                 # Main prediction script
├── dataset/                       # Training and test datasets
├── note                          # Development notes and observations
├── run_time_performance.txt      # Detailed performance benchmarks
├── rf_cnn_mdl_use_time.note     # Model timing comparisons
└── C_glabarta_protein/          # Example application data
```

## 🤝 Contributing

We welcome contributions to improve Deep-AmPEP30:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/improvement`)
3. Commit your changes (`git commit -am 'Add new feature'`)
4. Push to the branch (`git push origin feature/improvement`)
5. Create a Pull Request

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- **Protein Sequence Analysis**: Built using the `seqinr` and `protr` packages
- **Machine Learning**: Powered by `keras`, `caret`, and `randomForest`
- **Performance Optimization**: Extensive benchmarking and optimization for production use

## 📞 Contact

For questions, suggestions, or collaboration opportunities, please contact:
- **Project Maintainer**: [Your Name]
- **Email**: [your.email@domain.com]
- **Issue Tracker**: [GitHub Issues](https://github.com/your-username/Deep-AmPEP30/issues)

## 🔗 Related Projects

- **AxPEP Web Server**: [Link to web interface]
- **AMP Database**: [Link to database]
- **Publication**: [Link to associated research paper]

---

*Deep-AmPEP30: Advancing antimicrobial peptide prediction through deep learning and comprehensive feature engineering.*