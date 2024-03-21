# GDD-agreement

This repository contains the Python scripts to compute the GDD-agreement.


## What's the GDD-agreement?

The *GDD-agreement* is the similarity between two networks defined in the following articles:

- N. Pr&zcaron;ulj (2007). Biological network comparison using graphlet degree distribution. *Bioinform.*, **23**, 2, 177&ndash;183.<br><https://doi.org/10.1093/bioinformatics/btl301>


## Usage

```python
import networkx as nx

# Generate two networks
G = nx.fast_gnp_random_graph(20, 0.3, seed=1)
H = nx.fast_gnp_random_graph(15, 0.7, seed=2)

from GDDA import GDDA
gdda = GDDA()

# Compute the GDD-agreement.
a = gdda.agreement(G, H, method="arith")

print(a)# 0.63057...
```
