# Load module
# module
import numpy as np
from numpy import interp
import pandas as pd
import os
import datetime
import json
import codecs
from itertools import cycle
from requests import get
from pathlib import Path
from functools import reduce
from retry import retry
from sqlalchemy import create_engine
import subprocess
import pymysql

from sklearn.preprocessing import StandardScaler, label_binarize
from sklearn.pipeline import Pipeline, make_pipeline
from sklearn.feature_selection import RFECV
from sklearn.model_selection import train_test_split, GridSearchCV, RepeatedStratifiedKFold, StratifiedKFold, RandomizedSearchCV, KFold, cross_val_score
from sklearn.linear_model import Lasso
from sklearn.impute import KNNImputer
from sklearn.svm import SVC,SVR
from sklearn.multiclass import OneVsRestClassifier
from sklearn.metrics import roc_auc_score, roc_curve, auc
import matplotlib.pyplot as plt
from imblearn.over_sampling import SMOTE
import requests
from requests.adapters import HTTPAdapter, Retry

# roc curve
def roc_auc_function(y_test, y_score, num_class):
    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    for i in range(len(num_class)): 
        fpr[i], tpr[i], _ = roc_curve(y_test[:, i], y_score[:, i]) 
        roc_auc[i] = auc(fpr[i], tpr[i])

    # Compute micro-average ROC curve and ROC area 
    fpr["micro"], tpr["micro"], _ = roc_curve(y_test.ravel(), y_score.ravel()) 
    roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

    # First aggregate all false positive rates
    all_fpr = np.unique(np.concatenate([fpr[i] for i in range(len(num_class))]))

    # Then interpolate all ROC curves at this points
    mean_tpr = np.zeros_like(all_fpr)
    for i in range(len(num_class)):
        mean_tpr += interp(all_fpr, fpr[i], tpr[i])
    # Finally average it and compute AUC
    mean_tpr /= len(num_class)
    fpr["macro"] = all_fpr
    tpr["macro"] = mean_tpr
    roc_auc["macro"] = auc(fpr["macro"], tpr["macro"])

    return fpr, tpr, roc_auc
def roc_acu_calculator(DF, feature_name, log_save, over_sampling, module_name):
  X = DF.iloc[:, 1:]
  y = DF.iloc[:, 0]
  
  # SMOTE oversampling for minority
  if over_sampling:
    sm = SMOTE("minority", random_state=331)
    X, y = sm.fit_resample(X,y)
    
  # multi-class detection
  num_class = set(y)
  if len(num_class) > 2:
    y = label_binarize(y, classes=list(num_class))
  
  # Learn to predict each class against the other 
  X_train, X_test, y_train, y_test = train_test_split(X, y, test_size = .3, random_state=331)
  classifier = OneVsRestClassifier(SVC(kernel='linear', probability=True, random_state=331)) 
  y_score = classifier.fit(X_train, y_train).decision_function(X_test)
  
  if len(num_class) > 2:    
      
    # Compute ROC curve and ROC area for each class
    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    for i in range(len(num_class)): 
        fpr[i], tpr[i], _ = roc_curve(y_test[:, i], y_score[:, i]) 
        roc_auc[i] = auc(fpr[i], tpr[i])

    # Compute micro-average ROC curve and ROC area 
    fpr["micro"], tpr["micro"], _ = roc_curve(y_test.ravel(), y_score.ravel()) 
    roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])

    # First aggregate all false positive rates
    all_fpr = np.unique(np.concatenate([fpr[i] for i in range(len(num_class))]))

    # Then interpolate all ROC curves at this points
    mean_tpr = np.zeros_like(all_fpr)
    for i in range(len(num_class)):
        mean_tpr += interp(all_fpr, fpr[i], tpr[i])
    # Finally average it and compute AUC
    mean_tpr /= len(num_class)
    fpr["macro"] = all_fpr
    tpr["macro"] = mean_tpr
    roc_auc["macro"] = auc(fpr["macro"], tpr["macro"])

    # Plot all ROC curves
    lw = 2
    plt.figure(figsize=(10,10))
    plt.plot(fpr["micro"], tpr["micro"],
             label='micro-average ROC curve (area = {0:0.4f})'
                   ''.format(roc_auc["micro"]),
             color='deeppink', linestyle=':', linewidth=4)

    plt.plot(fpr["macro"], tpr["macro"],
             label='macro-average ROC curve (area = {0:0.4f})'
                   ''.format(roc_auc["macro"]),
             color='navy', linestyle=':', linewidth=4)

    colors = cycle(['aqua', 'darkorange', 'darkgreen', 'violet', 'peru', 'gold'])
    for i, color in zip(range(len(num_class)), colors):
        plt.plot(fpr[i], tpr[i], color=color, lw=lw,
                 label='ROC curve of class {0} (AUC={1:0.4f})'
                 ''.format(i, roc_auc[i]))

    plt.plot([0, 1], [0, 1], 'k--', lw=lw)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate (FPR)')
    plt.ylabel('True Positive Rate (TPR)')
    plt.title(module_name + "_" + feature_name + " (Multi-class)")
    plt.legend(loc="lower right")
    plt.savefig(log_save + "/" + module_name + "_" + feature_name + '_ROC_AUC.png')
    plt.clf()
    
    # calculation
    y_prob = classifier.predict_proba(X_test)
    macro_roc_auc_ovr = roc_auc_score(y_test, y_prob, multi_class="ovr", average="macro")
    micro_roc_auc_ovr = roc_auc_score(
        y_test, y_prob, multi_class="ovr", average="micro")
    return {'macro' : macro_roc_auc_ovr, 'micro': micro_roc_auc_ovr}
    return macro_roc_auc_ovr
  else :
    fpr, tpr, thres = roc_curve(y_test, y_score)
    roc_auc = roc_auc_score(y_test, y_score)

    # roc curve
    lw = 2
    plt.plot(fpr, tpr, color='darkorange',
             lw=lw, label='ROC curve (area = %0.4f)' % roc_auc)
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlabel('False Positive Rate (FPR)')
    plt.ylabel('True Positive Rate (TPR)')
    plt.title(module_name + "_" + feature_name + " (Binary-class")
    plt.legend(loc="best")
    plt.savefig(log_save + "/" + module_name + "_" + feature_name + '_ROC_AUC.png')
    plt.clf()
    
    return roc_auc
  
def roc_acu_calculator_cv(DF, feature_name, log_save, over_sampling, module_name):
    X = DF.iloc[:, 1:]
    y = DF.iloc[:, 0]

    # SMOTE oversampling for minority
    if over_sampling:
        sm = SMOTE(sampling_strategy="minority", random_state=331)
        X, y = sm.fit_resample(X,y)

    # multi-class detection
    num_class = set(y)
    if len(num_class) > 2:
        y = label_binarize(y, classes=list(num_class))

    fold_cnt = 5
    kfold = KFold(n_splits=fold_cnt)
    scores = []
    pipe_line = make_pipeline(StandardScaler(),
                             OneVsRestClassifier(SVC(kernel='linear', probability=True, random_state=331)))

    roc_auc_fold = dict()
    fpr_fold = dict()
    tpr_fold = dict()

    for k, (train, test) in enumerate(kfold.split(DF)):
        y_score = pipe_line.fit(X.iloc[train, :], y[train]).decision_function(X.iloc[test, :])

        if len(num_class) > 2:
            fpr, tpr, roc_auc = roc_auc_function(y_test=y[test], y_score=y_score, num_class=num_class)
            roc_auc_fold[k] = roc_auc['micro']
            fpr_fold[k] = fpr['micro']
            tpr_fold[k] = tpr['micro']

        else :
            fpr, tpr, thres = roc_curve(y[test], y_score, drop_intermediate=False)
            roc_auc_fold[k] = roc_auc_score(y[test], y_score)
            fpr_fold[k] = fpr
            tpr_fold[k] = tpr


    # First aggregate all false positive rates
    all_fpr = np.unique(np.concatenate([fpr_fold[i] for i in range(len(num_class))]))

    # Then interpolate all ROC curves at this points
    mean_tpr = np.zeros_like(all_fpr)
    # for i in range(len(num_class)):
    for i in range(fold_cnt):
        mean_tpr += interp(all_fpr, fpr_fold[i], tpr_fold[i])

    # Finally average it and compute AUC
    # mean_tpr /= len(num_class)
    mean_tpr /= fold_cnt
    fpr_fold["macro"] = all_fpr
    tpr_fold["macro"] = mean_tpr
    roc_auc_fold["macro"] = auc(fpr_fold["macro"], tpr_fold["macro"])

    # Plot all ROC curves
    lw = 2
    plt.figure(figsize=(10,10))
    plt.plot(fpr_fold["macro"], tpr_fold["macro"],
             label='macro-average ROC curve (AUC = {0:0.3f})'
                   ''.format(roc_auc_fold["macro"]),
             color='red', linestyle=':', linewidth=4)

    colors = cycle(['aqua', 'darkorange', 'darkgreen', 'violet', 'peru', 'gold'])
    for i, color in zip(range(fold_cnt), colors):
        plt.plot(fpr_fold[i], tpr_fold[i], color=color, lw=lw, alpha=0.2,
                 label='ROC curve fold-{0} (AUC={1:0.3f})'
                 ''.format(i + 1, roc_auc_fold[i]))

    plt.plot([0, 1], [0, 1], 'k--', lw=lw)
    # plt.xlim([0.0, 1.0])
    # plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate (FPR)')
    plt.ylabel('True Positive Rate (TPR)')
    plt.title(module_name + "_" + feature_name + "_Multi-class-CV")
    plt.legend(loc="lower right")
    plt.savefig(log_save + "/Step4_" + module_name + "_" + feature_name + '_ROC_AUC_CV.png')
    plt.close()

    return roc_auc_fold["macro"]
  
# gene selection
def feature_selection_svm_rfecv(X, Y, over_sampling):
  # SMOTE oversampling for minority
  if over_sampling:
    sm = SMOTE("minority",random_state=331)
    X, y = sm.fit_resample(X,y)  

  # model conf.
  k_fold = StratifiedKFold(n_splits=5, shuffle=True, random_state=331)
  svc = SVC(kernel = 'linear', probability = True)
  rfecv = RFECV(estimator=svc, cv=k_fold, scoring='roc_auc')
  
  ## search parameter
  param_grid = {
      # SVC - linear kernel의 경우 'C'만 parameter로 확인하면 됨
      'estimator__C': np.arange(0.01,100,0.01)
  }
  
  # CV conf.
  CV_rfc = RandomizedSearchCV(estimator=rfecv, 
                      param_distributions=param_grid, 
                      n_iter=50,
                      cv= k_fold, 
                      scoring = 'roc_auc', 
                      verbose=1,
                      random_state=331,
                      n_jobs=10)
  CV_rfc.fit(X, Y.values.ravel())
  
  return CV_rfc.best_estimator_.support_

def feature_selection_LASSO(X, y, over_sampling=None):
    num_class = set(y)
    
    # SMOTE oversampling for minority
    if over_sampling:
      sm = SMOTE("minority",random_state=331)
      X, y = sm.fit_resample(X,y)
    
    if len(num_class) > 2:
        y = label_binarize(y, classes=list(num_class))

        # ML pipeline
        pipeline = Pipeline([
                     ('scaler',StandardScaler()),
                     ('model', OneVsRestClassifier(Lasso(max_iter=99999)))])

        # grid search
        search = GridSearchCV(pipeline,
                      {'model__estimator__alpha':np.arange(0.01,3,0.001)},
                      cv = 10, scoring="neg_mean_squared_error",verbose=1
                      )

        search.fit(X,y)

        print(search.best_params_)
        coefficients = pd.DataFrame(search.best_estimator_.named_steps.model.coef_)
        return coefficients.abs().sum().to_list()

    else :
        pipeline = Pipeline([
                     ('scaler',StandardScaler()),
                     ('model',Lasso(max_iter=99999))])

        # grid search
        search = GridSearchCV(pipeline,
                          {'model__alpha':np.arange(0.01,3,0.001)},
                          cv = 10, scoring="neg_mean_squared_error",verbose=1
                          )

        search.fit(X,y.values.ravel())

        print(search.best_params_)
        coefficients = search.best_estimator_.named_steps['model'].coef_
        return coefficients
      
def non_zero_column(DF):
  sample_cnt = int(len(DF.columns) * 0.2)
  zero_row = dict(DF.isin([0]).sum(axis=1))
  non_remove_feature = list()

  for key, value in zero_row.items():
      if value < sample_cnt:
          non_remove_feature.append(key)
  
  return non_remove_feature
def cancer_select(cols, cancer_type, raw_path):
    # phenotype
    phe1 = pd.read_csv(raw_path + "GDC-PANCAN.basic_phenotype.tsv", sep="\t")
    phe1 = phe1.loc[phe1.program == "TCGA", :].loc[:, ['sample', 'sample_type', 'project_id']].drop_duplicates(['sample'])
    phe1['sample'] =  phe1.apply(lambda x : x['sample'][:-1], axis=1)
    phe2 = pd.read_csv(raw_path + "TCGA_phenotype_denseDataOnlyDownload.tsv", sep="\t")
    ph_join = pd.merge(left = phe2 , right = phe1, how = "left", on = "sample").dropna(subset=['project_id'])
    
    if cancer_type == "PAN" or cancer_type == "PANCAN":
        filterd = ph_join.loc[ph_join['sample_type_y'] == "Primary Tumor", :]
        sample_barcode = filterd["sample"].tolist()
    else:
        filterd = ph_join.loc[((ph_join['sample_type_y'] == "Primary Tumor")|(ph_join['sample_type_y'] == "Solid Tissue Normal")) & (ph_join['project_id'] == "TCGA-" + cancer_type) , :]
        sample_barcode = filterd["sample"].tolist()
        
    intersect_ = list(set(cols).intersection(sample_barcode))
    
    return intersect_
def load_tcga_dataset(pkl_path, raw_path, cancer_type):   
  # subfunction
  # main
  if os.path.isfile(pkl_path + cancer_type + "_rna.pkl"):
    # sep
    rna = pd.read_pickle(pkl_path  + cancer_type + "_rna.pkl")
  else :
      # create dir
      Path(pkl_path).mkdir(parents=True, exist_ok=True)
      Path(raw_path).mkdir(parents=True, exist_ok=True)
      
      # file name
      mrna_f = 'tcga_RSEM_gene_fpkm.gz'
      
      # file url
      mrna_url = "https://toil-xena-hub.s3.us-east-1.amazonaws.com/download/tcga_RSEM_gene_fpkm.gz"
       # to make GET request
      def download(url, file_name):
          with open(file_name, "wb") as file:   # open in binary mode
              response = get(url)               # get request
              file.write(response.content)      # write to file
  
      # mrna
      if os.path.isfile(raw_path + mrna_f) == False:
          download(mrna_url, raw_path + mrna_f)
      
      # RNA gene expression
      if os.path.isfile(pkl_path + cancer_type + "_rna.pkl") == False:
          col = pd.read_csv(raw_path + mrna_f, sep = "\t", index_col=0, nrows=0).columns.to_list()
          use_col = ['sample'] + cancer_select(cols=col, cancer_type=cancer_type, raw_path=raw_path)
          df_chunk = pd.read_csv(raw_path + mrna_f,
                       sep = "\t", index_col=0, iterator=True, chunksize=50000, usecols=use_col)
          rna = pd.concat([chunk for chunk in df_chunk])
          rna = rna[rna.index.isin(non_zero_column(rna))]
  
          rna.to_pickle(pkl_path + cancer_type + "_rna.pkl")
      else : 
          rna = pd.read_pickle(pkl_path  + cancer_type + "_rna.pkl")
      
      # set same column for merge
  
      # pickle save
      rna.to_pickle(pkl_path + "/" + cancer_type + "_rna.pkl")
  
  # set index
  rna_index = rna.index.to_list()
  
  # missing impute
  imputer = KNNImputer(n_neighbors=10)
  rna_impute = imputer.fit_transform(rna)
  
  omics = pd.DataFrame(rna_impute, columns=rna.columns)
  rna.index = rna_index
  
  return rna

# Uniprot, PDB ID Search
POLLING_INTERVAL = 3
API_URL = "https://rest.uniprot.org"
tsv_loader = "?compressed=true&fields=accession%2Creviewed%2Cid%2Cprotein_name%2Corganism_name&format=tsv"

retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))


def check_response(response):
    try:
        response.raise_for_status()
    except requests.HTTPError:
        print(response.json())
        raise


def submit_id_mapping(from_db, to_db, taxid, ids):
    request = requests.post(
        f"{API_URL}/idmapping/run",
        data={"from": from_db, "to": to_db, "taxId": taxid, "ids": ",".join(ids)},
    )
    check_response(request)
    return request.json()["jobId"]


def get_next_link(headers):
    re_next_link = re.compile(r'<(.+)>; rel="next"')
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)


def check_id_mapping_results_ready(job_id):
    while True:
        request = session.get(f"{API_URL}/idmapping/status/{job_id}")
        check_response(request)
        j = request.json()
        if "jobStatus" in j:
            if j["jobStatus"] == "RUNNING":
                print(f"Retrying in {POLLING_INTERVAL}s")
                time.sleep(POLLING_INTERVAL)
            else:
                raise Exception(request["jobStatus"])
        else:
            return bool(j["results"] or j["failedIds"])


def get_batch(batch_response, file_format, compressed):
    batch_url = get_next_link(batch_response.headers)
    while batch_url:
        batch_response = session.get(batch_url)
        batch_response.raise_for_status()
        yield decode_results(batch_response, file_format, compressed)
        batch_url = get_next_link(batch_response.headers)


def combine_batches(all_results, batch_results, file_format):
    if file_format == "json":
        for key in ("results", "failedIds"):
            if key in batch_results and batch_results[key]:
                all_results[key] += batch_results[key]
    elif file_format == "tsv":
        return all_results + batch_results[1:]
    else:
        return all_results + batch_results
    return all_results


def get_id_mapping_results_link(job_id):
    url = f"{API_URL}/idmapping/details/{job_id}"
    request = session.get(url)
    check_response(request)
    return request.json()["redirectURL"]


def decode_results(response, file_format, compressed):
    if compressed:
        decompressed = zlib.decompress(response.content, 16 + zlib.MAX_WBITS)
        if file_format == "json":
            j = json.loads(decompressed.decode("utf-8"))
            return j
        elif file_format == "tsv":
            return [line for line in decompressed.decode("utf-8").split("\n") if line]
        elif file_format == "xlsx":
            return [decompressed]
        elif file_format == "xml":
            return [decompressed.decode("utf-8")]
        else:
            return decompressed.decode("utf-8")
    elif file_format == "json":
        return response.json()
    elif file_format == "tsv":
        return [line for line in response.text.split("\n") if line]
    elif file_format == "xlsx":
        return [response.content]
    elif file_format == "xml":
        return [response.text]
    return response.text


def get_xml_namespace(element):
    m = re.match(r"\{(.*)\}", element.tag)
    return m.groups()[0] if m else ""


def merge_xml_results(xml_results):
    merged_root = ElementTree.fromstring(xml_results[0])
    for result in xml_results[1:]:
        root = ElementTree.fromstring(result)
        for child in root.findall("{http://uniprot.org/uniprot}entry"):
            merged_root.insert(-1, child)
    ElementTree.register_namespace("", get_xml_namespace(merged_root[0]))
    return ElementTree.tostring(merged_root, encoding="utf-8", xml_declaration=True)


def print_progress_batches(batch_index, size, total):
    n_fetched = min((batch_index + 1) * size, total)
    print(f"Fetched: {n_fetched} / {total}")

def get_id_mapping_results_stream(url):
    if "/stream/" not in url:
        url = url.replace("/results/", "/results/stream/")
    request = session.get(url)
    check_response(request)
    parsed = urlparse(url)
    query = parse_qs(parsed.query)
    file_format = query["format"][0] if "format" in query else "json"
    compressed = (
        query["compressed"][0].lower() == "true" if "compressed" in query else False
    )
    return decode_results(request, file_format, compressed)

def mergeDictionary(dict_):
    dict_3 = defaultdict(list)
    for d in dict_: # you can list as many input dicts as you want here
        for key, value in d.items():
            dict_3[key].append(value)
            # list(itertools.chain(*list_of_lists))
    return dict_3


def get_pdb_structure_id(gene_list):
    # Gene to Uniprot
    job_id = submit_id_mapping(
        from_db="Gene_Name", to_db="UniProtKB", taxid=9606, ids=gene_list
    )
    
    print(job_id)
    
    link = get_id_mapping_results_link(job_id)
    mapping_result = get_id_mapping_results_stream(link + tsv_loader)
    mapping_result = list(map(lambda x : x.split("\t"), mapping_result))
    mapping_result_df = pd.DataFrame(mapping_result[1:], columns=mapping_result[0])
    m_v = mapping_result_df[mapping_result_df['Reviewed'] == "reviewed"]

    # Uniprot to PDB
    m_v['pdb'] = m_v.apply(lambda x : Query(x['Entry']).search(), axis=1)
    # m_v['pdbIDCount'] = m_v.apply(lambda x : 0 if x['pdb'] is None else int(len(x['pdb'])), axis=1)
    m_v['pdbID'] = m_v.apply(lambda x : None if x['pdb'] is None else ';'.join(x['pdb']), axis=1)
    m_v = m_v[['From', 'Entry', 'Reviewed', 'pdbID']]
    m_v.columns = ['GENE_NAME', 'UniprotKB', 'Reviewed', 'pdbID']
    
    # post processing
    pdbID_collapse = m_v.groupby('GENE_NAME').apply(lambda x : ';'.join(filter(None, x['pdbID'])))
    pdbID_collapse.name = "pdbID"

    UniprotKB_collapse = m_v.groupby('GENE_NAME').apply(lambda x : ';'.join(filter(None, x['UniprotKB'])))
    UniprotKB_collapse.name = "UniprotKB"

    pdb_collapse = pd.merge(UniprotKB_collapse, pdbID_collapse, right_index = True,
                   left_index = True)

    pdb_collapse['pdbID'] = pdb_collapse.apply(lambda x : x['pdbID'].split(';'), axis=1)

    pdb_collapse['pdbCount'] = pdb_collapse.apply(lambda x : len(set(x['pdbID'])) - 1, axis=1)

    pdb_collapse['pdbID'] = pdb_collapse.apply(lambda x : ';'.join(set(x['pdbID'])),axis = 1)
    
    return pdb_collapse

# textmining
aws_mariadb_url = 'mysql+pymysql://root:sempre813!@192.168.0.91:3306/Textmining'
engine_mariadb = create_engine(aws_mariadb_url)  

def db_query(x):   
    q1 = pd.read_sql_table(table_name=x, con=engine_mariadb)
    q1.columns = ['gene', x + '_TYPE', x + '_SUPPORT', x + '_CONFIDENCE', x + '_LIFT', x + '_COUNT']
    return q1

def symbol2pdb(gene_list):
    headers = {'content-type': 'application/x-www-form-urlencoded'}
    params = 'q='+ ','.join(gene_list)+' &scopes=symbol&fields=taxid,pdb'
    res = requests.post('http://mygene.info/v3/query', data=params, headers=headers)
    out = pd.DataFrame(json.loads(codecs.decode(bytes(res.text, 'utf-8'), 'utf-8-sig')))
    out = out[out.taxid == 9606]
    out = out[out.pdb.notna()]
    out = out.loc[:, ['query', 'pdb']]

    out['pdbCount'] = out.apply(lambda x : len(x['pdb']) if isinstance(x['pdb'], list) else 1, axis=1)
    out['pdb'] = out.apply(lambda x : ';'.join(x['pdb']) if isinstance(x['pdb'], list) else x['pdb'],axis = 1)
    
    return out

# ONCOKB
def query_mariadb(db, query):
    # Connect to MariaDB (mariadb.connect 대신 pymysql.connect 로 해도 된다)
    dbconn = pymysql.connect(
        user="root",
        password="sempre813!",
        host="192.168.0.91",
        port=3306,
        database=db
    )
 
    # dbconn = mydb.cursor()  # 이 명령어는 불필요.
    # mariaDB Query to Pandas DataFrame
    query_result= pd.read_sql(query,dbconn)
    dbconn.close()
 
    return query_result

def oncokb_allcuratedGenes():

    ONCOKB_TOKEN = query_mariadb(db="TOKEN", query="SELECT * FROM ONCOKB").TOKEN.to_list()[0]

    proc = subprocess.run(["curl",  "-X", "GET",  
                           'https://www.oncokb.org/api/v1/utils/allCuratedGenes.txt?includeEvidence=true',
                           '-H', 'accept: application/json',
                           '-H', 'Authorization: Bearer '+ ONCOKB_TOKEN[0]

                          ],
                       stdout=subprocess.PIPE, encoding='utf-8')

    oncokb_curated = proc.stdout
    oncokb_curated_df = pd.DataFrame(list(map(lambda x: x.split("\t"), oncokb_curated.split('\n'))))
    oncokb_curated_df.columns = oncokb_curated_df.iloc[0]
    oncokb_curated_df.drop(oncokb_curated_df.index[0], inplace=True)
    oncokb_curated_df = oncokb_curated_df.loc[:, ['Hugo Symbol', 'Is Oncogene', 'Is Tumor Suppressor Gene', 'Highest Level of Evidence(sensitivity)',
           'Highest Level of Evidence(resistance)', 'Background']]
    oncokb_curated_df.columns = ['gene'] + ["OncoKB_" + value for value in oncokb_curated_df.columns[1:]]
    
    return oncokb_curated_df

