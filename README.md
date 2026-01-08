This is the data and code associated with "Population  Dynamics of Native Benthic Macroinvertebrates Utilizing the Zebra Mussel (Dreissena polymorpha) as Habitat in the Mississippi River Headwater Region" by Amanda K.
Schuermann, Andrew W. Hafs, Richard W. Koch, and Debbie L. Guelda.

The code was written by Amanda S. and Andrew H. in 2024, and updated in late 2025 for the most recent version of R packages. It contains not only the code used to visualize and understand the data, but also
the code used to create the figures in the final manuscript.

The data was collected by Amanda S. in the summer of 2024. Field site selection and methodologies can be found in the publication.

Data sheets explanation: Due to the collaborative nature of this project, and the fact that several analyses were done at different times, there are several data sheets associated with this project.

Data_key_taxa_NMDS.csv: This data sheet contains the counts of the key taxa explained in the publication in a format that made the NMDS analysis easier.

diversity_corr_updated_1.csv: This data sheet contains the extrapolated values for the zebra mussel infested site using the surface area measurements of the druses.

diversity_PERMANOVA.upt.csv: This data sheet contains the counts of the key taxa explained in the publication in a format that made the PERMANOVA analysis easier.

data_AFDM_csv.csv: Data from ash free dry mass conducted after field season.

Column name explanation
These data sheets all contain similar column names:

Sample ID: arbitrary number used for laboratory ordering purposes assigned to each individual sample.

Subseason/Site: The site and subseason number assocated with the sample.

Taxa names: the count of each taxa found within specific sample.

Taxa names _corr: The correct number of individuals found in that sample using the surface area of the zebra mussel druse as outlined in the manuscript.

Average length mm: average length of zebra mussel in sample.

Maco count: total number of macroinvertebrates in the sample.

All questions about the code should be directed to the corresponding author, Amanda S, schue621@umn.edu

Macro per zeb: number of macroinvertebrates/number of zebra mussels

SA druse  mm: Surface area of druse in mm based on 2D approximation using foil method.

zeb_presence: 1/0 denoting Y/N for zebra mussel presence.
