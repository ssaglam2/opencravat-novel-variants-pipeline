# Methodology Summary Report
## OpenCRAVAT Novel Variants Discovery Pipeline

### Overview

This document provides a concise summary of the methodology used in the OpenCRAVAT Novel Variants Discovery Pipeline for identifying potentially pathogenic novel variants from whole exome sequencing data.

---

## 🎯 Objective

**Primary Goal**: Identify novel pathogenic variants using computational prediction tools while avoiding bias from existing clinical databases.

**Key Innovation**: Population-specific analysis that avoids ClinVar dependence to discover truly novel variants.

---

## 🧬 Methodological Approach

### 1. Data Input Requirements

**Essential Annotations:**
- CADD deleteriousness scores
- REVEL pathogenicity predictions
- PolyPhen2 functional predictions
- Population frequency data (gnomAD, 1000 Genomes)
- Variant functional annotations

### 2. Scoring Algorithm

**Multi-Component Scoring System (0-100 points):**

```
Novel_Score = CADD_Component (35) + REVEL_Component (25) + 
              Frequency_Component (20) + LoF_Bonus (10) + 
              Consensus_Score (8) + Population_Specificity (10)
```

### 3. Component Calculations

#### CADD Component (0-35 points)
- Formula: `(CADD_score / 40.0) × 35`
- Rationale: CADD >30 considered highly deleterious
- Normalization: Scaled to 40 for extreme scores

#### REVEL Component (0-25 points)
- Formula: `REVEL_score × 25`
- Rationale: REVEL >0.5 indicates pathogenic
- Range: Direct multiplication of 0-1 score

#### Frequency Component (0-20 points)
- Formula: `(-log₁₀(max_freq + 1e-8) / 8) × 20`
- Rationale: Rare variants more likely pathogenic
- Scaling: Logarithmic to handle wide frequency range

#### Loss-of-Function Bonus (0-10 points)
- Criteria: stop_gained, frameshift_variant, splice variants
- Rationale: Clear functional impact mechanism
- Implementation: Binary bonus for qualifying variants

#### Consensus Score (0-8 points)
- PolyPhen2 damaging: +3 points
- REVEL >0.5: +5 points
- Rationale: Multiple predictor agreement

#### Population Specificity (0-10 points)
- Formula: `log₁₀(other_pop_freq / target_pop_freq)`
- Rationale: Identifies population-specific variants
- Applications: Pharmacogenomics, founder effects

---

## 🌍 Population Analysis Framework

### Target Population Strategy

**European Ancestry Focus:**
- Primary: gnomAD NFE, gnomAD FIN
- Secondary: 1000G GBR, IBS, TSI, CEU
- Validation: Cross-population comparison

**Frequency Calculations:**
- European frequency: Minimum across European populations
- Non-European frequency: Maximum across other populations
- Specificity ratio: Other/European frequency

### Population Categories

1. **Ultra-rare**: <0.001% frequency
2. **Very rare**: 0.001-0.01% frequency  
3. **Rare**: 0.01-1% frequency
4. **Population-specific**: Differential frequency patterns

---

## 🔬 Quality Control Measures

### Data Validation
- Column presence verification
- Score range validation
- Missing data handling
- Duplicate variant identification

### Statistical Validation
- Cross-database consistency checks
- Multi-predictor correlation analysis
- Population stratification verification

### Computational Validation
- Memory-efficient chunk processing
- Performance monitoring
- Error handling and logging

---

## 📊 Output Classification

### Functional Categories
- **Loss-of-function**: Clear impact variants
- **Missense damaging**: High prediction scores
- **Consensus variants**: Multiple predictor agreement
- **Population-specific**: Frequency differentials

### Clinical Prioritization Tiers
- **Tier 1**: Score >70, immediate investigation
- **Tier 2**: Score 60-70, high priority
- **Tier 3**: Score 50-60, research interest

---

## 🎯 Validation Strategy

### Computational Validation
- Score distribution analysis
- Predictor agreement assessment
- Population frequency consistency

### Biological Validation
- Gene function enrichment
- Pathway analysis
- Literature correlation

### Clinical Validation
- Phenotype association
- Family segregation
- Functional assays

---

## ⚙️ Technical Implementation

### Processing Optimization
- **Chunk size**: Configurable (default: 100K variants)
- **Memory usage**: <8GB RAM for large datasets
- **Processing time**: Approximately 1 minute per 100K variants
- **Output formats**: CSV, summary statistics, visualizations

### Software Dependencies
- Python 3.8+
- pandas, numpy for data processing
- matplotlib, seaborn for visualization
- scipy for statistical analysis

---

## 📈 Performance Metrics

### Expected Performance
- **Computational efficiency**: Linear scaling with dataset size
- **Memory footprint**: Optimized for large datasets
- **Accuracy metrics**: Cross-database correlation >0.9
- **Consensus rate**: >75% multi-predictor agreement

### Quality Indicators
- Score distribution normality
- Population frequency consistency
- Predictor correlation strength
- Processing completion rate

---

## 🚨 Limitations and Considerations

### Technical Limitations
- Dependence on annotation quality
- Population representation bias
- Predictor tool limitations
- Missing annotation handling

### Biological Limitations
- Computational predictions need validation
- Clinical correlation not guaranteed
- Population coverage constraints
- Temporal database changes

### Ethical Considerations
- Privacy protection requirements
- Institutional policy compliance
- Consent and data sharing guidelines
- Clinical interpretation responsibilities

---

## 🔄 Future Enhancements

### Algorithmic Improvements
- Machine learning integration
- Additional predictor tools
- Tissue-specific annotations
- Pathway-based scoring

### Data Integration
- Expanded population coverage
- Functional genomics data
- Clinical phenotype data
- Pharmacogenomic annotations

---

## 📚 References and Resources

### Key Publications
- CADD: Kircher et al., Nat Genet 2014
- REVEL: Ioannidis et al., Am J Hum Genet 2016
- gnomAD: Karczewski et al., Nature 2020
- 1000 Genomes: 1000 Genomes Consortium, Nature 2015

### Software Resources
- OpenCRAVAT: https://opencravat.org/
- gnomAD: https://gnomad.broadinstitute.org/
- Pipeline repository: [Your GitHub URL]

---

## 📞 Support and Contact

For methodology questions or technical support:

**Pipeline Development**: [Your Name/Team]  
**Institution**: [Your Institution]  
**Documentation**: See `docs/methodology.md` for complete details  
**Issues**: GitHub repository issue tracker  

---

**Document Version**: 1.0  
**Last Updated**: [Current Date]  
**Methodology Version**: OpenCRAVAT Novel Variants Pipeline v1.0 