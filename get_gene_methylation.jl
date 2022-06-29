#read metdense file to get average methylation for each gene and cell
using CSV
using DataFrames

allGenes = CSV.File("data/all_genes.txt") |> DataFrame

sum((allGenes.tss .!= allGenes.start) .& (allGenes.tss .!= allGenes.fin))