# Lab journal
Samuel Morabito
smorabit@uci.edu


---
### 03/07/19

I am starting this lab journal markdown document in order to have a better record
of what I do / what I hope to get done in lab each day.

At yesterday's lab meeting we discussed the paper [Integrative transcriptome analyses of the aging brain implicate altered splicing in Alzheimer's disease susceptibility](https://www.nature.com/articles/s41588-018-0238-1?utm_source=feedburner&utm_medium=feed&utm_campaign=Feed%3A+ng%2Frss%2Fcurrent+%28Nature+Genetics+-+Issue%29#Bib1). In this paper intron clusters of interest were computed from their data using a package called [Leafcutter](http://davidaknowles.github.io/leafcutter/). Figure 2b in the aforementioned paper showed the variance explained within each of these intron clusters based on 4 criteria: AD diagnosis, amyloid burden, tangle burden, and neuritic plaques. We thought that it would be interesting to use the same package [variancePartition](https://www.bioconductor.org/packages/devel/bioc/vignettes/variancePartition/inst/doc/variancePartition.pdf) to demonstrate the variance explained by each of these criteria in our modules instead of intron cluters.

I began working on a script to use variancePartition on the Mayo dataset MEs
