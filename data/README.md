# Data Directory

This directory contains input data files for the OpenCRAVAT Novel Variants Discovery Pipeline.

## Required Data Format

The pipeline expects OpenCRAVAT-annotated CSV files with the following structure:

### Essential Columns

#### Variant Information
- `chrom`: Chromosome (e.g., chr1, chr2, ..., chrX, chrY)
- `pos`: Genomic position
- `ref`: Reference allele
- `alt`: Alternative allele
- `hugo`: Gene symbol
- `so`: Sequence Ontology term (variant type)

#### Deleteriousness Scores
- `cadd.phred`: CADD Phred-scaled scores (higher = more deleterious)
- `revel.score`: REVEL pathogenicity scores (0-1, higher = more pathogenic)

#### Prediction Tools
- `polyphen2.hdiv_pred`: PolyPhen2 HumDiv predictions
- `polyphen2.hdiv_score`: PolyPhen2 HumDiv scores

#### Population Frequencies

**gnomAD v3.1:**
- `gnomad.af`: Overall allele frequency
- `gnomad.af_nfe`: Non-Finnish European frequency
- `gnomad.af_fin`: Finnish frequency
- `gnomad.af_afr`: African frequency
- `gnomad.af_eas`: East Asian frequency
- `gnomad.af_sas`: South Asian frequency
- `gnomad.af_amr`: Latino/American frequency

**1000 Genomes Project:**
- `thousandgenomes.af`: Overall 1000G frequency
- `thousandgenomes.eur_af`: European frequency
- `thousandgenomes.gbr_af`: British frequency
- `thousandgenomes.ibs_af`: Spanish frequency
- `thousandgenomes.tsi_af`: Italian frequency
- `thousandgenomes.ceu_af`: Utah/CEPH frequency
- `thousandgenomes.afr_af`: African frequency
- `thousandgenomes.eas_af`: East Asian frequency
- `thousandgenomes.sas_af`: South Asian frequency

#### Clinical Annotations (Optional)
- `clinvar.sig`: ClinVar clinical significance
- `clinvar.review`: ClinVar review status

## Data Size Considerations

- **Large files**: Use Git LFS for files >100MB
- **CSV format**: Recommended for analysis
- **Compression**: .gz files supported
- **Memory**: Pipeline uses chunk processing for large datasets

## Sample Data

A small sample dataset is provided for testing:
- `sample_variants.csv`: 1000 representative variants
- Column headers and format example

## Data Sources

### Primary Source
- **OpenCRAVAT**: https://opencravat.org/
- Comprehensive variant annotation platform
- Integrates multiple databases and predictors

### Population Databases
- **gnomAD**: https://gnomad.broadinstitute.org/
- **1000 Genomes**: https://www.internationalgenome.org/

### Prediction Tools
- **CADD**: https://cadd.gs.washington.edu/
- **REVEL**: https://sites.google.com/site/revelgenomics/
- **PolyPhen2**: http://genetics.bwh.harvard.edu/pph2/

## Usage Notes

1. Place your OpenCRAVAT-annotated CSV file in this directory
2. Update the file path in the analysis scripts
3. Ensure all required columns are present
4. Large files (>1GB) may require chunk processing

## Privacy and Ethics

- Remove personally identifiable information
- Follow institutional data sharing policies
- Comply with relevant regulations (GDPR, HIPAA, etc.)
- Consider data use agreements

For questions about data format or requirements, please refer to the main documentation. 