# Load the FamAgg package
install.packages("FamAgg")
library("tinytex")

# Load the genotype data

genotype_data <- read.table("obesitycohort.bed")

# Load the phenotype data

phenotype_data <- read.table("obesitycohort.fam")

# Perform the family-based rare variant association test

result <- FamAgg(genotype_data, phenotype_data, method="TDT")

# Print the results

print(result)
