
<!-- README.md is generated from README.Rmd. Please edit that file -->

## Research Article

This repository contains the code and data used to carry out analyses
for this research article:

**Doré et al., 2021 - Anthropogenic pressures coincide with Neotropical
biodiversity hotspots in a flagship butterfly group. Diversity and
Distributions.**

[![DOI](https://zenodo.org/badge/DOI/10.1111/ddi.13455.svg)](https://doi.org/10.1111/ddi.13455)

[![DOI](https://zenodo.org/badge/DOI/10.1111/ddi.13455.svg)](https://doi.org/10.1111/ddi.13455).

In this study, we provide new maps of the distribution of the different
facets of Ithomiini diversity in the entire Neotropics in order to
identify areas of both evolutionary and ecological importance for
conservation. Our study build on a large and previously unpublished
dataset of ca. 30 000 geolocalized occurrences, and a recently published
phylogeny of the entire tribe to present diversity analyses and a
conservation assessment to an unprecedent taxonomic and geographic scale
for a non-Vertebrate taxon in the Neotropics. Specifically, we map
mimicry ring and species richness, geographic rarity, as well as
phylogenetic diversity. We evaluate spatial congruence among those
diversity indices and assess current anthropogenic threats to diversity
hotspots.

All content is available
[here](https://github.com/MaelDore/ithomiini_diversity).

## Contents

  - [:file\_folder: **extra\_scripts**](extra_scripts/) directory
<<<<<<< HEAD
    contains additional scripts to conduct supplementary analyses
    provided in the Supplementary Information.

  - [:file\_folder: **for\_article**](for_article/) directory contains
    final Figures produced for the research article.
=======
    contains supplementary scripts for complementary analyses not
    included in the article.

  - [:file\_folder:**for\_article**](for_article/) directory contains
    the final figures from the article.
>>>>>>> 82754e16074c0b3851c4d15b9420290112bdc65c

  - [:file\_folder: **functions**](functions/) directory contains
    homemade functions used in the analyses.

  - [:file\_folder: **graphs**](graphs/) directory contains all the
    graphs generated by the scripts during the analyses.

  - [:file\_folder: **input\_data**](input_data/) directory contains the
    raw and transformed data used in the analyses. Sub-folders include
    environmental data from [MERRAclim
    v.2.0](https://doi.org/10.5061/dryad.s2v81) (Vega et al., 2017),
    anthropogenic pressures data from the [Human Footprint
    map](https://doi.org/10.1038/ncomms12558) (Venter et al., 2016),
    elevation data from the [SRTM v.4.1](http://srtm.csi.cgiar.org)
    (Farr et al., 2007), and forest cover data from the [Landsat Tree
    Cover Continuous Fields
    dataset](https://developers.google.com/earth-engine/datasets/catalog/NASA_MEASURES_GFCC_TC_v3)
    (Sexton et al., 2013).

  - [:file\_folder: **models**](models/) directory contains output of
    species distribution models (SDMs) once generated.

  - [:file\_folder: **phylogenies**](phylogenies/) directory contains
    all phylogeny graphs generated by the scripts during the analyses.

  - [:file\_folder: **outputs**](outputs/) directory contains all files
    generated by the scripts that are not maps, phylogenies or SDM
    outputs.

  - [:file\_folder: **renv**](renv/) directory contains the full library
    of packages used for this analyses to ensure reproductibility.

  - [:file\_folder: **scripts**](scripts/) directory contains the
    scripts used to run the analyses

  - [:file\_folder: **supplementaries**](supplementaries/) directory
    contains the outputs used to generate the Figures for the
    Supplementary Materials

## How to run it

This research has been developed using the statistical programming
language R. To run the analyses, you will need installed on your
computer the [R software](https://cloud.r-project.org/) itself and
optionally [RStudio
Desktop](https://rstudio.com/products/rstudio/download/).

You can download the entire project as a `.zip` from [this
URL](/archive/master.zip). After unzipping:

  - Open the `ithomiini_diversity.Rproj` file, found at the root of the
    project, in RStudio

  - Run sequentially the scripts found in the [:file\_folder:
    **scripts**](scripts/) folder. It will rebuild the model outputs,
    maps, graphs and figures, including the final ones presented in the
    main text and Supplementaries of the article.

## How to cite

Please cite this research article as:

> Doré, M., Willmott, K., Leroy, B., Chazot, N., Mallet, J., Freitas, A.
> V. L., Hall, J. P. W., Lamas, G., Dasmahapatra, K. K., Fontaine, C., &
> Elias, M. (2021). Anthropogenic pressures coincide with Neotropical
> biodiversity hotspots in a flagship butterfly group. Diversity and
> Distributions, 00, 1–19. <https://doi.org/10.1111/ddi.13455>

## Associated online archives

<<<<<<< HEAD
The dataset used to run the analyses is made publicy available at 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4696055.svg)](https://doi.org/10.5281/zenodo.4696055)

The mimicry classification applied in this study can be found at 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5497876.svg)](https://doi.org/10.5281/zenodo.5497876)
=======
The occurrence dataset used to run the analyses is made publicly
available at
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4696055.svg)](https://doi.org/10.5281/zenodo.4696055).
>>>>>>> 82754e16074c0b3851c4d15b9420290112bdc65c

The mimicry classification of ithomiine subspecies defined for this
study is made publicly available at
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5497876.svg)](https://doi.org/10.5281/zenodo.5497876).

The distribution maps generated throughout the analyses are made
publicly available at 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4673446.svg)](https://doi.org/10.5281/zenodo.4673446)
