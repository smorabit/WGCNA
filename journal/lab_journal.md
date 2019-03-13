# Bioinformatics lab journal
Samuel Morabito
smorabit@uci.edu


---
### 03/07/19

I am starting this lab journal markdown document in order to have a better record
of what I do / what I hope to get done in lab each day.

At yesterday's lab meeting we discussed the paper [Integrative transcriptome analyses of the aging brain implicate altered splicing in Alzheimer's disease susceptibility](https://www.nature.com/articles/s41588-018-0238-1?utm_source=feedburner&utm_medium=feed&utm_campaign=Feed%3A+ng%2Frss%2Fcurrent+%28Nature+Genetics+-+Issue%29#Bib1). In this paper intron clusters of interest were computed from their data using a package called [Leafcutter](http://davidaknowles.github.io/leafcutter/). Figure 2b in the aforementioned paper showed the variance explained within each of these intron clusters based on 4 criteria: AD diagnosis, amyloid burden, tangle burden, and neuritic plaques. We thought that it would be interesting to use the same package [variancePartition](https://www.bioconductor.org/packages/devel/bioc/vignettes/variancePartition/inst/doc/variancePartition.pdf) to demonstrate the variance explained by each of these criteria in our modules instead of intron cluters.

I began working on a script to use variancePartition on the Mayo dataset MEs, producing the following figures:

![image1](figures/mayo_variance_fraction.png?raw=true =475x) ![image2](figures/mayo_variance_partition.png?raw=true =475x)

The next step is to project the MEs from the Mayo dataset onto the ROSMAP data, since there is metadata pertaining to immunohistochemistry for things like amyloid burden, tangle burden, neuritic plaques etc.

---
### 03/12/19

Today I worked on projecting MEs from the Mayo dataset onto ROSMAP data, and then examining the variance within the data attributed to several features describing the ROSMAP samples. The features we are particularly interested are burden of neuritic plaques and tangles.

I was able to project the modules onto the ROSMAP data, and then run the variancePartition package on these MEs, but the output did not look so good. Every feature tested was hovering around zero variance with residuals at nearly 100%.

Next I tried to run variancePartition just using the ROSMAP gene expression data, but I couldn't get it to run. It isn't providing any error messages which is frustrating. I think next I am going to try running module preservation on the the ROSMAP/MAYO MEs. Maybe the modules are just not preserved and don't really mean anything in the context of ROSMAP...? Although I doubt it.

In the time that I started writing this entry I just began running module preservation so we'll se how that goes tomorrow morning.
