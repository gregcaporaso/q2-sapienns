# q2-sapienns

q2-sapienns is the working title for a set of tools that can be used for preparing [BioBakery3](https://doi.org/10.7554/eLife.65088) data for use in [QIIME 2](https://qiime2.org). As QIIME 2 expands support for metagenomics data analysis, this will provide a framework for working with processed BioBakery3 data, and for comparing other methods to BioBakery3.

**At present this is under-tested, and not ready for use in real analyses. Consider this "alpha" software.** There are test data files in this repository in the `q2-sapienns/tests/data` directory, but no detailed unit tests yet.

Some basic usage examples are provided below.

## HUMAnN 3

### Pathway abundance tables

Import a HUMANn 3 _Pathway abundance file_. See the [HUMAnN 3 User Manual and Tutorial](https://huttenhower.sph.harvard.edu/humann) for details on this file and how to create it. There can be one or more samples in this file.

```
qiime tools import --input-path humann-pathabundance-2.tsv --output-path humann-pathabundance-2.qza --type HumannPathAbundanceTable
```

Create `FeatureTable[Frequency]` and `FeatureData[Taxonomy]` artifacts from the imported table. **The semantic type of the `FeatureData` artifact is likely to change.**

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

Import a HUMANn 3 _Gene families file_ and create `FeatureTable[Frequency]` and `FeatureData[Taxonomy]` artifacts from the imported table. **The semantic type of the `FeatureData` artifact is likely to change.**. See the [HUMAnN 3 User Manual and Tutorial](https://huttenhower.sph.harvard.edu/humann) for details on this file and how to create it. There can be one or more samples in this file.

```
qiime tools import --input-path humann-genefamilies-2.tsv --output-path humann-genefamilies-2.tsv --type HumannGeneFamilyTable
```

```
qiime sapienns humann-genefamily --i-genefamily-table humann-genefamilies-2.qza --o-table table.qza --o-taxonomy feature-data.qza
qiime sapienns humann-genefamily --i-genefamily-table humann-genefamilies-2.qza --o-table table-destratified.qza --o-taxonomy feature-data-destratified.qza --p-destratify
```

## MetaPhlAn 3

Import a MetaPhlAn 3 taxonomy file and create `FeatureTable[RelativeFrequency]` and `FeatureData[Taxonomy]` artifacts from the imported table.See the [MetaPhlAn 3 documentation](https://huttenhower.sph.harvard.edu/metaphlan) for details on this file and how to create it. There can be one or more samples in this file.

```
qiime tools import --input-path metaphlan-merged-abundance-1.tsv --output-path metaphlan-merged-abundance-1.qza --type MetaphlanMergedAbundanceTable
```

```
qiime sapienns metaphlan-taxon --i-stratified-table metaphlan-merged-abundance-1.qza --p-level 7 --o-table species-table.qza --o-taxonomy taxonomy.qza
```
