# Bioinformatics Portfolio

[![Tests](https://github.com/yazanus94-dot/Bioinformatics-Statistical-Genomics-Portfolio/actions/workflows/ci.yml/badge.svg)](https://github.com/yazanus94-dot/Bioinformatics-Statistical-Genomics-Portfolio/actions/workflows/ci.yml)
![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue)
![License](https://img.shields.io/badge/license-MIT-green)

Two self-contained computational biology projects, each with tested code, reproducible results and a worked example problem with a verified answer.

| Project | Question it answers | Methods |
|---------|---------------------|---------|
| [Case-control variant association](project1-case-control-variant-association/) | How much does genotyping-error QC change which variants are significant? | exact Hardy-Weinberg test, allelic Fisher test, odds ratios, Bonferroni/BH, genomic inflation |
| [RNA-seq differential expression](project2-rnaseq-differential-expression/) | Which genes and pathways respond to a treatment? | median-of-ratios normalisation, PCA, Welch DE test, hypergeometric enrichment |

Both projects use only NumPy, SciPy, pandas and Matplotlib, run in seconds, and are checked by automated tests (`pytest`) on every push via GitHub Actions.

## Run everything

```bash
for p in project1-case-control-variant-association project2-rnaseq-differential-expression; do
  (cd "$p" && pip install -r requirements.txt && python -m pytest && python run_analysis.py && python verify_task.py)
done
```

## License

MIT, see [LICENSE](LICENSE).
