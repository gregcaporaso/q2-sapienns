# q2-sapienns

q2-sapienns is a set of tools that can be used for preparing [BioBakery3](https://doi.org/10.7554/eLife.65088) data for use in [QIIME 2](https://qiime2.org). As QIIME 2 expands support for metagenomics data analysis, this will provide a framework for working with processed BioBakery3 data, and for comparing other methods to BioBakery3.

**q2-sapienns is comprehensively unit tested, but hasn't been field-tested much yet. If you notice any issues, please post to [the issue tracker](https://github.com/gregcaporaso/q2-sapienns/issues).** Basic usage examples are provided below.

## Installation

First, create and/or activate a QIIME 2 environment by following the [QIIME 2 install instructions](https://docs.qiime2.org/2022.8/install/native/). q2-sapienns has most recently been tested with QIIME 2 2022.8.

Then, install q2-sapienns using `pip` as follows...

```bash
pip install git+https://github.com/gregcaporaso/q2-sapienns.git
```

... refresh your QIIME 2 environment...

```bash
qiime dev refresh-cache
```

... and you should now see `sapienns` in your list of available QIIME 2 plugins:

```bash
$ qiime --help
Usage: qiime [OPTIONS] COMMAND [ARGS]...
...
  sample-classifier   Plugin for machine learning prediction of sample metadata.
  sapienns            Plugin for interacting with biobakery data.
  taxa                Plugin for working with feature taxonomy annotations.
...
```

## Usage: HUMAnN 3

There seem to have not been changes to the HUMAnN file formats used here between version 2-3.5 (and likely version 4), so these tools should work with all of those versions [(source)](https://forum.biobakery.org/t/human-and-metaphlan-file-formats/4024/3?u=gregcaporaso). If you notice any issues, please let me know!

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

## Usage: MetaPhlAn 3

There may be relevant changes to the file formats used here between versions of MetaPhlAn, though those changes may not be relevant to the _Merged Abundance Table_ [(source)](https://forum.biobakery.org/t/human-and-metaphlan-file-formats/4024/3?u=gregcaporaso). This functionality was developed for the MetaPhlAn format that contains exactly two columns (`clade_name` and `NCBI_tax_id`) before the sample abundance columns. I recommend looking at the column headers for the first three columns in your input file before attempting to use this code. The file should look something like:

```
$ head -5 metaphlan-merged-abundance.tsv
#mpa_v30_CHOCOPhlAn_201901
clade_name	NCBI_tax_id	sample1	sample_2
k__Archaea	2157	9.75907	0.02352
k__Archaea|p__Euryarchaeota	2157|28890	9.75907	0.02352
k__Archaea|p__Euryarchaeota|c__Methanobacteria	2157|28890|183925	9.75907	0.02352
```

q2-sapienns _should_ fail if you try to import data in a format different than the one it's expecting, but I can't be sure that format validation will work in all cases. It won't hurt to look at your data before using it with q2-sapienns.

### Merged abundance table

Import a MetaPhlAn 3 taxonomy file and create `FeatureTable[RelativeFrequency]` and `FeatureData[Taxonomy]` artifacts from the imported table.See the [MetaPhlAn 3 documentation](https://huttenhower.sph.harvard.edu/metaphlan) for details on this file and how to create it. There can be one or more samples in this file. If using the default reference with MetaPhlAn 3, the taxnomic ids will refer to the [NCBI taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy).

```
qiime tools import --input-path metaphlan-merged-abundance-1.tsv --output-path metaphlan-merged-abundance-1.qza --type MetaphlanMergedAbundanceTable
```

```
qiime sapienns metaphlan-taxon --i-stratified-table metaphlan-merged-abundance-1.qza --p-level 7 --o-table species-table.qza --o-taxonomy taxonomy.qza
```
