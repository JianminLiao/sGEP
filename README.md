# A Successive Two-stage Method for Sparse Generalized Eigenvalue Problems
### The source code is provided the purpose of providing reproducibility information for the paper.
## Authors: 
#### Qia Li, School of Computer Science and Engineering, Guangdong Province Key Laboratory of Computational Science, Sun Yat-sen University, Guangzhou, Guangdong Province 510275, P.R.China (liqia@mail.sysu.edu.cn);
#### Jianmin Liao, Department of Mathematics, Syracuse University, Syracuse, NY 13244, USA (jliao21@syr.edu);
#### Lixin Shen, Department of Mathematics, Syracuse University, Syracuse, NY 13244, USA (lshen03@syr.edu);
#### Na Zhang, Department of Applied Mathematics, College of Mathematics and Informatics, South China Agricultural University, Guangzhou, Guangdong Province 510642, P.R.China (nzhsysu@gmail.com).

## Quick Demos 

This project includes quick demos for sparse Principal Component Analysis (sPCA) and sparse Canonical Correlation Analysis (sCCA). To run these demos, you can use the following scripts:

- `demo_sPCA.m`: For sparse Principal Component Analysis
- `demo_sCCA.m`: For sparse Canonical Correlation Analysis

## Proposed Algorithms

This project introduces two new algorithms, which are implemented in the following scripts:

- `SA_TPM.m`
- `SA_PGSA_ML.m`

### Associated Files

The following associated files are required to run the algorithms:

- `TPower_modified.m`
- `SGEP_PGSA_M.m`
- `sa_partial.m`
- `sa_partial_SGEP.m`

## To reproduce all numerical results presented in the paper

Run the corresponding .m files to reproduce the numerical results in the ``Numerical Experiments in the Paper'' folder. The numerical experiments in the paper are repeated many times to provide meaningful results, which implies that the experiments typically take much longer time than quick demos.

