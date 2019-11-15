# Porch
![porch](img/porch.jpg "Photo (c) Sonja Lovas, https://www.flickr.com/photos/sonjalovas/4038233322 ")  

A method to explore your personal pathways.  

## Background

Many experiments in biology and medicine aim at detecting and quantifying differences in concentrations of biomolecules between experimental conditions. However, such concentration differences may stem from common biological factors and are often correlational. Consequently, one frequently uses so-called pathway analysis to estimate enrichments of structural, functional or locational relations. This is dependent on precomputed relations, thus there is a wide range of publicly available pathway databases.

## Porch
Porch is a method for pathway analysis that estimates a single quantitative pathway expression value for each pathway and sample.  Our method does not only spectacularly outperform current state-of-the-art methods, like GSEA, by orders of magnitude, but we also show how our method can be combined with general linear models to test more complex relations than regular case-control scenarios.


```eval_rst
.. toctree::
   :maxdepth: 2
   :caption: Contents:
```

## Porch API

```eval_rst
.. automodule:: porch.porch
   :members:
```
## The `qvalue` submodule

```eval_rst
.. automodule:: porch.qvalue
   :members:
```
