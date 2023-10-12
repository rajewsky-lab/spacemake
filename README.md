<a href="https://pypi.org/project/spacemake/">
   <img src="https://img.shields.io/pypi/v/spacemake.svg" / ></a>
   
<a href="https://spacemake.readthedocs.io/">
   <img src="https://readthedocs.org/projects/spacemake/badge/?version=latest" / ></a>
   
 <a href="https://pepy.tech/project/spacemake">
   <img src="https://pepy.tech/badge/spacemake" / ></a>

# Spacemake: processing and analysis of large-scale spatial transcriptomics data

<img src="https://raw.githubusercontent.com/rajewsky-lab/spacemake/master/docs/graphical_abstract_twitter.png" width="400" />

Spacemake is a modular, robust, and scalable spatial transcriptomics pipeline built in `Snakemake` and `Python`. Spacemake is designed to handle all major spatial transcriptomics datasets and can be readily configured for other technologies. It can process and analyze several samples in parallel, even if they stem from different experimental methods. Spacemake's unified framework enables reproducible data processing from raw sequencing data to automatically generated downstream analysis reports. Spacemake is built with a modular design and offers additional functionality such as sample merging, saturation analysis, and analysis of long reads as separate modules.

If you find Spacemake useful in your work, consider citing it: 

```
Spacemake: processing and analysis of large-scale spatial transcriptomics data
Tamas Ryszard Sztanka-Toth, Marvin Jens, Nikos Karaiskos, Nikolaus Rajewsky
GigaScience, Volume 11, 2022, giac064
```

Documentation can be found [here](https://spacemake.readthedocs.io/en/latest/).

## Unit testing

We are committed to achieving a high code coverage with unit tests. The master branch utilizes the `unittest` module to run spacemake with small test data sets. On the current development branches, we have switched to `pytest` and cover a much broader range of the code. This work is ongoing.

To run the currently implemented tests on master, run `python spacemake/unittests.py`. This will create a directory `spacemake/_tests/` inside which a minimal spacemake directory structure will be created using `spacemake init` and subsequently some of the core functionality (adding genomes/species, samples, changing configuration, etc.) will be executed. All output will be logged to `spacemake/_tests/run_spacemake.out.log`. If you encounter any weird behavior, please make sure to include the content of this file in your ticket on the issue tracker. Thank you!
...
