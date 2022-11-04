# RRA-WGCNA

* **Environment**

```
conda create -n geo-py python=3.8
conda activate geo-py
pip install -r env/requirements.txt --user
```

* **Function**

  ```
  Step 1 : Gene Expression Omnibus(GEO) preprocessing - GEO Database & GEO Meta Database
  Step 2 : Differential expression analysis & Robust rank aggregation - run_DEA_RRA.R
  Step 3 : Weighted gene correlation network analysis(WGCNA) - run_WCGNA_TCGA.R
  Step 4 : Machine learning validation - run_ML_validation.R
  ```

* **Workflow**

  ![workflow](https://user-images.githubusercontent.com/42958809/167774216-553390d6-0e07-470f-996e-420ac8b1afc3.png)
