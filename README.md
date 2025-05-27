# OpenCRAVAT Novel Variants Discovery Pipeline

[![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://python.org)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Documentation](https://img.shields.io/badge/docs-latest-brightgreen.svg)](docs/)

A comprehensive computational pipeline for discovering potentially pathogenic novel variants from OpenCRAVAT-annotated whole exome sequencing data, with focus on population-specific analysis and avoidance of database bias.

## 🔬 Project Overview

This pipeline analyzes large-scale OpenCRAVAT-annotated WES datasets to identify novel pathogenic variants using population genetics principles and computational prediction tools. The pipeline specifically avoids reliance on existing clinical databases (ClinVar) to discover truly novel variants that may not yet be documented in clinical literature.

### Key Features
- ✅ **Multi-dimensional scoring system** (0-100 points)
- ✅ **Population-specific variant discovery**
- ✅ **Ultra-rare variant identification** (frequency <0.001%)
- ✅ **Loss-of-function variant characterization**
- ✅ **Cross-population frequency analysis**
- ✅ **Memory-efficient chunk processing**

## 📊 Pipeline Capabilities

- **Large dataset processing**: Handles millions of variants efficiently
- **Comprehensive annotations**: Integrates 100+ annotation columns
- **Population focus**: Customizable for different ancestry groups
- **Scalable analysis**: Chunk-based processing for memory optimization

## 🧬 Methodology

### Novel Variant Scoring Formula

```
Novel_Score = CADD_Component + REVEL_Component + Frequency_Component + 
              LoF_Bonus + Consensus_Score + Population_Specificity_Bonus
```

**Component Details:**
- **CADD Component**: (CADD_score / 40.0) × 35 points (max 35)
- **REVEL Component**: REVEL_score × 25 points (max 25)
- **Frequency Component**: (-log₁₀(max_population_freq + 1e-8) / 8) × 20 points (max 20)
- **LoF Bonus**: 10 points for stop_gained, frameshift, splice variants
- **Consensus Score**: 3 points (PolyPhen2) + 5 points (REVEL > 0.5)
- **Population Specificity**: log₁₀(other_pop_freq / target_pop_freq) (max 10)

### Population Data Sources

**European Populations:**
- gnomAD Non-Finnish European (NFE)
- gnomAD Finnish (FIN)
- 1000 Genomes: British (GBR), Spanish (IBS), Italian (TSI), Utah/CEPH (CEU)

**Cross-Population Validation:**
- gnomAD: African (AFR), East Asian (EAS), South Asian (SAS), Latino/American (AMR)

## 🏆 Example Output

### Sample High-Scoring Variants

| Gene | Position | Type | Score | CADD | REVEL | Status |
|------|----------|------|-------|------|-------|--------|
| **GENE1** | chr1:12345678 | missense | 75.2 | 32.0 | 0.95 | Ultra-rare |
| **GENE2** | chr2:23456789 | missense | 72.1 | 30.5 | 0.92 | Ultra-rare |
| **GENE3** | chr3:34567890 | missense | 70.8 | 29.8 | 0.88 | Very rare |

*Note: These are example outputs. Actual results depend on your input data.*

## 📁 Repository Structure

```
├── scripts/                    # Analysis pipeline scripts
│   ├── european_novel_variants_discovery.py    # Main discovery pipeline
│   ├── chromosome_wise_analysis.py             # Chromosome-based analysis
│   ├── high_impact_variants_analysis.py        # High-impact variant detection
│   └── run_analysis.py                         # Pipeline runner
├── reports/                    # Example analysis reports
│   ├── sample_analysis_report.md               # Example findings report
│   └── methodology_summary.md                  # Analysis approach summary
├── docs/                      # Documentation and methodology
├── data/                      # Sample data files and format specs
├── results/                   # Analysis output files
└── visualizations/           # Generated plots and figures
```

## 🚀 Quick Start

### Prerequisites

```bash
pip install -r requirements.txt
```

### Basic Usage

```python
# Run the complete novel variants discovery pipeline
python scripts/european_novel_variants_discovery.py

# Analyze high-impact variants specifically
python scripts/high_impact_variants_analysis.py

# Run chromosome-wise analysis
python scripts/chromosome_wise_analysis.py
```

### Input Data Format

The pipeline expects OpenCRAVAT-annotated CSV files with the following key columns:
- `cadd.phred`: CADD deleteriousness scores
- `revel.score`: REVEL pathogenicity scores
- `polyphen2.hdiv_pred`: PolyPhen2 predictions
- `gnomad.*.af`: Population frequency data
- `thousandgenomes.*.af`: 1000 Genomes frequencies

## 📈 Analysis Features

### Core Capabilities
- ✅ **Population-specific variant discovery**
- ✅ **Multi-predictor consensus scoring**
- ✅ **Loss-of-function variant identification**
- ✅ **Cross-population frequency analysis**
- ✅ **Ultra-rare variant prioritization**
- ✅ **Chunk-based processing for large datasets**

### Computational Optimizations
- Memory-efficient chunk processing (configurable chunk size)
- Parallel processing support
- Progress tracking and logging
- Error handling and validation

## 🔍 Expected Outputs

### Variant Categories
- **Ultra-rare variants**: frequency <0.001%
- **Population-specific variants**: differential frequency patterns
- **Loss-of-function variants**: stop, frameshift, splice
- **High-scoring variants**: CADD >30, REVEL >0.8
- **Consensus variants**: multiple predictor agreement

### Clinical Prioritization Tiers
- **Tier 1**: Highest scoring variants requiring immediate investigation
- **Tier 2**: High-priority candidates for follow-up
- **Tier 3**: Research-interest variants for population studies

## 📊 Quality Metrics

The pipeline includes quality control measures:
- Cross-database consistency validation
- Multi-predictor agreement assessment
- Population stratification verification
- Processing efficiency monitoring

## 🔬 Biological Applications

### Target Gene Categories
- **Cardiac/Muscle function**: ion channels, structural proteins
- **Metabolic pathways**: enzymes, transporters
- **Neurological function**: neurotransmitters, development
- **Immune system**: signaling, receptors
- **Cell cycle/DNA repair**: tumor suppressors, repair genes

### Population Genetics Applications
- Founder effect identification
- Population-specific selective pressure analysis
- Pharmacogenomic variant discovery
- Ancestry-specific risk assessment

## 📚 Documentation

Detailed documentation available in the [`docs/`](docs/) directory:
- [Methodology Guide](docs/methodology.md)
- [Installation Instructions](docs/installation.md)
- [Tutorial Examples](docs/tutorials/)
- [Troubleshooting Guide](docs/troubleshooting.md)

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## ⚠️ Important Notes

### Data Privacy
- **Never upload actual patient data**
- **Remove all personally identifiable information**
- **Follow institutional data sharing policies**
- **Comply with relevant regulations (GDPR, HIPAA, etc.)**

### Usage Guidelines
- This pipeline is for research purposes
- Computational predictions require experimental validation
- Results should be interpreted by qualified professionals
- Consider population-specific factors in interpretation

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📞 Contact

**Project Lead**: Sinan SAGLAM 
**Email**: sinansaglam6@gmail.com  
**Institution**: University of Pisa

## 🙏 Acknowledgments

- OpenCRAVAT team for annotation pipeline
- gnomAD consortium for population data
- 1000 Genomes Project for reference populations
- CADD and REVEL prediction tools

## 📈 Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{opencravat_novel_variants_2024,
  author = {Sinan Saglam]},
  title = {Novel Variants Discovery Pipeline},
  year = {2025},
  url = {https://github.com/[ssaglam2]/novel-variants},
  
}
```

---

**Last Updated**: May 2025 
**Version**: 1.0.0  
**Status**: Template/Framework for Novel Variant Discovery 
