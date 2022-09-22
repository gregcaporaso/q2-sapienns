# q2-sapienns

q2-sapienns is a set of tools that can be used for preparing [BioBakery3](https://doi.org/10.7554/eLife.65088) data for use in [QIIME 2](https://qiime2.org). As QIIME 2 expands support for metagenomics data analysis, this will provide a framework for working with processed BioBakery3 data, and for comparing other methods to BioBakery3.

**q2-sapienns is comprehensively unit tested, but hasn't been field-tested much yet. If you notice any issues, please post to [the issue tracker](https://github.com/gregcaporaso/q2-sapienns/issues).**

Some basic usage examples are provided below.

## HUMAnN 3

### Pathway abundance tables

Import a HUMANn 3 _Pathway abundance file_. See the [HUMAnN 3 User Manual and Tutorial](https://huttenhower.sph.harvard.edu/humann) for details on this file and how to create it. There can be one or more samples in this file. If using default reference data with HUMANn 3, the pathway annotations will refer to [MetaCyc pathways](https://metacyc.org/).

```
qiime tools import --input-path humann-pathabundance-2.tsv --output-path humann-pathabundance-2.qza --type HumannPathAbundanceTable
```

Create `FeatureTable[Frequency]` and `FeatureData[Taxonomy]` artifacts from the imported table.

```
qiime sapienns humann-pathway --i-pathway-table humann-pathabundance-2.qza --o-table table.qza --o-taxonomy feature-data.qza
```

Create `FeatureTable[Frequency]` and `FeatureData[Taxonomy]` artifacts from the imported table, dropping pathway annotations with taxonomic information (destratified).

```
qiime sapienns humann-pathway --i-pathway-table humann-pathabundance-2.qza --o-table table-destratified.qza --o-taxonomy feature-data-destratified.qza --p-destratify
```

Summarize created artifacts for viewing.

```
qiime feature-table summarize --i-table table-destratified.qza --o-visualization table-destratified.qzv --m-sample-metadata-file sample-metadata.tsv
qiime metadata tabulate --m-input-file feature-data-destratified.qza --o-visualization feature-data-destratified.qzv
```

### Gene family tables

Import a HUMANn 3 _Gene families file_ and create `FeatureTable[Frequency]` and `FeatureData[Taxonomy]` artifacts from the imported table. See the [HUMAnN 3 User Manual and Tutorial](https://huttenhower.sph.harvard.edu/humann) for details on this file and how to create it. There can be one or more samples in this file. If using the default reference with HUMANn 3, the gene family annotations will refer to [UniRef50](https://www.uniprot.org/).

```
qiime tools import --input-path humann-genefamilies-2.tsv --output-path humann-genefamilies-2.tsv --type HumannGeneFamilyTable
```

```
qiime sapienns humann-genefamily --i-genefamily-table humann-genefamilies-2.qza --o-table table.qza --o-taxonomy feature-data.qza
qiime sapienns humann-genefamily --i-genefamily-table humann-genefamilies-2.qza --o-table table-destratified.qza --o-taxonomy feature-data-destratified.qza --p-destratify
```

## MetaPhlAn 3

Import a MetaPhlAn 3 taxonomy file and create `FeatureTable[RelativeFrequency]` and `FeatureData[Taxonomy]` artifacts from the imported table.See the [MetaPhlAn 3 documentation](https://huttenhower.sph.harvard.edu/metaphlan) for details on this file and how to create it. There can be one or more samples in this file. If using the default reference with MetaPhlAn 3, the taxnomic ids will refer to the [NCBI taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy).

```
qiime tools import --input-path metaphlan-merged-abundance-1.tsv --output-path metaphlan-merged-abundance-1.qza --type MetaphlanMergedAbundanceTable
```

```
qiime sapienns metaphlan-taxon --i-stratified-table metaphlan-merged-abundance-1.qza --p-level 7 --o-table species-table.qza --o-taxonomy taxonomy.qza
```
