![BuildStatus](https://github.com/h-s-miller/phagepy/workflows/test.yml/badge.svg?event=push)

PhagePy is a library that streamlines the processing of Phage Immunoprecipitation Sequencing (PhIP-Seq) data. I am working on more documentation, so stay stuned :) 


# Downloading PhagePy
```
git clone https://github.com/h-s-miller/phagepy.git
cd phagepy
pip install .
```

# Using PhagePy
### data structure
Phagepy uses the [anndata](https://anndata.readthedocs.io/en/latest/) data structure from [scanpy](https://scanpy.readthedocs.io/en/stable/) to store metadata, sequencing counts, and gene annotation information in one object. Reading this documentation will be very helpful to understanding how to manipulate adata objects for this analysis! 
### how to load the data
Phagepy has a command, `create_anndata()` to load metadata and sequencing counts into an anndata object. `create_anndata()` has 3 parameters:
1. `counts_file`: a csv file with counts for each sample in rpk. 
2. `metadata_file`: a metadata file with annotations for each sample
    1. **note**: the sample ids in the first column of the metadata file **must** match the sample ids in the counts file
3.`transpose`(default=True): a boolean parameter.
    1. if your input counts file is (peptides)x(samples), then `transpose` should be set to True. Note: this is how data downloaded from PhageDB is formatted. 
    2. if your input counts file is (samples)x(peptides), then `transpose` should be set to True.
```
counts='sample.csv'
meta='sample_meta.csv'
adata=pp.create_anndata(counts, meta)
```

# More tutorial coming...
