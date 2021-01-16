# Genetic algorithms applied to 2020 SARS-CoV2 pandemic
## Theory
This kind of optimisation algorithms tries to emulate evolution by generating many different candidate solutions and then trying to improve them separately through applying a selective-pressure-like force on loci (parameters). The effect of this selective process is that beneficial alleles survive along time, while deletereous ones are renewed with new candidates.

New alleles are introduced in a population through the following processes:
- Mutation
- Immigration (from a different population)

Existing alleles are removed from population through:
- Selection
- Genetic drift
- (Emigration is negligible)

New combinations of alleles can arise by crossing-over the existing ones.
Remember that natural selection do not act on alleles themselves, but on their phenotypic effect.

## Implementation

Search space must be discretised (split) into 2^n parts (where n is a positive integer). Thus, a \"phenotypic effect\" (a value belonging this discretised search space) can be encoded as a \"gene\" by rewritting it as a binary number of n digits corresponding to this value.

Each \"organism\" has a bunch of binary sequences associated that act as its \"chromosomes\".










