# Expected output

The numeric output of the scripts, as produced on our systems, so that a rerun can be checked
against it without having to compare figures by eye. 72 KB in total: only the tables the figure
panels rest on, not the intermediate files and not the figures themselves, which the scripts
regenerate.

```bash
diff expected/fig1E_correlations.tsv results/fig1E_correlations.tsv
```

One file will not match exactly. `figS2E_gsea.tsv` carries permutation p-values, reproducible
only with the same seed, the same single process and the same gene set library, so check the
numbers rather than the bytes.
