import requests
import pandas as pd


# Endpoint base
BASE_URL = "https://www.cbioportal.org/api"
HEADERS = {"accept": "application/json"}


def buscar_nombres_datasets_con_keyword(keywords: list):
    """
    Busca estudios en cBioPortal usando una keywords específicas.
    Devuelve un DataFrame con los nombres de los estudios y sus descripciones.
    """
    # Lista para guardar los resultados
    results = []

    # Iterar sobre cada keyword
    for kw in keywords:
        print(f"Buscando estudios para la keyword: {kw}")
        params = {
            "keyword": kw,
            "projection": "DETAILED",
            "pageSize": 10000000,
            "pageNumber": 0,
            "direction": "ASC"
        }
        HEADERS = {"accept": "application/json"}
        
        response = requests.get(BASE_URL+"/studies", params=params, headers=HEADERS)
        response.raise_for_status()
        data = response.json()
        
        # Cada keyword puede devolver varios estudios
        for study in data:
            results.append({
                "name": study.get("name"),
                "description": study.get("description"),
                "studyId": study.get("studyId"),
                "cancerType": study.get("cancerType", {}).get("name")
            })

    # Convertir a DataFrame
    df = pd.DataFrame(results)

    # Mostrar tabla
    print(df)

    # Guardar en CSV
    df.to_csv("cbioportal_studies.csv", index=False)
    return df


def get_studies():
    """Obtiene todos los estudios disponibles en cBioPortal"""
    url = f"{BASE_URL}/studies?projection=SUMMARY&pageSize=10000"
    response = requests.get(url, headers=HEADERS)
    response.raise_for_status()
    return response.json()


def get_molecular_profiles(study_id):
    """Obtiene los perfiles moleculares de un estudio"""
    url = f"{BASE_URL}/studies/{study_id}/molecular-profiles"
    response = requests.get(url, headers=HEADERS)
    response.raise_for_status()
    return response.json()


def buscar_estudios_con_rna_seq():
    studies = get_studies()
    rna_seq_studies = []

    for study in studies:
        study_id = study["studyId"]
        profiles = get_molecular_profiles(study_id)
        
        # Filtrar por perfiles de tipo RNA-seq
        for profile in profiles:
            if "rna_seq" in profile["molecularAlterationType"].lower() \
               or "rna" in profile["molecularAlterationType"].lower():
                rna_seq_studies.append(study["name"])
                break  # Si ya tiene uno, no hace falta seguir buscando

    print("Estudios con datos de RNA-seq:")
    for s in rna_seq_studies:
        print(f"- {s}")
        


keys = ["Cancer Cell Line Encyclopedia (Broad, 2019)", "Breast Invasive Carcinoma (TCGA, Firehose Legacy)", "Breast Invasive Carcinoma (TCGA, PanCancer Atlas)", "Breast Invasive Carcinoma (TCGA, Cell 2015)", "Chronic Lymphocytic Leukemia (Broad, Nature Genetics 2022)", "Acute Myeloid Leukemia (OHSU, Cancer Cell 2022)", "Colorectal Adenocarcinoma (TCGA, PanCancer Atlas)", "Pan-cancer Analysis of Advanced and Metastatic Tumors (BCGSC, Nature Cancer 2020)", "Ovarian Cancer (Gray Foundation, Cancer Discov 2024)", "Kidney Renal Clear Cell Carcinoma (TCGA, Firehose Legacy)", "Brain Lower Grade Glioma (TCGA, Firehose Legacy)", "Uterine Corpus Endometrial Carcinoma (TCGA, PanCancer Atlas)", "Head and Neck Squamous Cell Carcinoma (TCGA, Firehose Legacy)", "Lung Adenocarcinoma (TCGA, Firehose Legacy)", "Head and Neck Squamous Cell Carcinoma (TCGA, PanCancer Atlas)", "Brain Lower Grade Glioma (TCGA, PanCancer Atlas)", "Kidney Renal Clear Cell Carcinoma (TCGA, PanCancer Atlas)", "Lung Adenocarcinoma (TCGA, PanCancer Atlas)", "Thyroid Carcinoma (TCGA, Firehose Legacy)", "Lung Squamous Cell Carcinoma (TCGA, Firehose Legacy)", "Thyroid Carcinoma (TCGA, PanCancer Atlas)", "Prostate Adenocarcinoma (TCGA, Firehose Legacy)", "Prostate Adenocarcinoma (TCGA, PanCancer Atlas)", "Lung Squamous Cell Carcinoma (TCGA, PanCancer Atlas)", "Papillary Thyroid Carcinoma (TCGA, Cell 2014)", "Skin Cutaneous Melanoma (TCGA, Firehose Legacy)", "Acute Myeloid Leukemia (OHSU, Nature 2018)", "Skin Cutaneous Melanoma (TCGA, PanCancer Atlas)", "Kidney Renal Clear Cell Carcinoma (TCGA, Nature 2013)", "Stomach Adenocarcinoma (TCGA, Firehose Legacy)", "Stomach Adenocarcinoma (TCGA, PanCancer Atlas)", "Bladder Urothelial Carcinoma (TCGA, Firehose Legacy)", "Bladder Cancer (TCGA, Cell 2017)", "Bladder Urothelial Carcinoma (TCGA, PanCancer Atlas)", "Colorectal Adenocarcinoma (TCGA, Firehose Legacy)", "Liver Hepatocellular Carcinoma (TCGA, Firehose Legacy)", "Liver Hepatocellular Carcinoma (TCGA, PanCancer Atlas)", "Diffuse Glioma (GLASS Consortium)", "Uterine Corpus Endometrial Carcinoma (TCGA, Nature 2013)", "Ovarian Serous Cystadenocarcinoma (TCGA, Firehose Legacy)", "Cervical Squamous Cell Carcinoma and Endocervical Adenocarcinoma (TCGA, Firehose Legacy)", "Ovarian Serous Cystadenocarcinoma (TCGA, PanCancer Atlas)", "Mature B-Cell Neoplasms (Simon Fraser University, Blood 2023)", "Bladder Cancer (MSK/TCGA, 2020)", "Cervical Squamous Cell Carcinoma (TCGA, PanCancer Atlas)", "Kidney Renal Papillary Cell Carcinoma (TCGA, Firehose Legacy)", "Prostate Adenocarcinoma (TCGA, Cell 2015)", "Kidney Renal Papillary Cell Carcinoma (TCGA, PanCancer Atlas)", "Head and Neck Squamous Cell Carcinoma (TCGA, Nature 2015)", "Stomach Adenocarcinoma (TCGA, Nature 2014)", "Sarcoma (TCGA, Firehose Legacy)", "Sarcoma (TCGA, PanCancer Atlas)", "Colorectal Adenocarcinoma (TCGA, Nature 2012)", "Pediatric Preclinical Testing Consortium (CHOP, Cell Rep 2019)", "Lung Adenocarcinoma (TCGA, Nature 2014)", "Adult Soft Tissue Sarcomas (TCGA, Cell 2017)", "Pediatric Acute Lymphoid Leukemia - Phase II (TARGET, 2018)", "Normal Melanocytes (UCSF, Nature 2020)", "Pediatric Brain Cancer (CPTAC/CHOP, Cell 2020)", "Esophageal Carcinoma (TCGA, Firehose Legacy)", "Pheochromocytoma and Paraganglioma (TCGA, Firehose Legacy)", "Esophageal Adenocarcinoma (TCGA, PanCancer Atlas)", "Lung Adenocarcinoma (OncoSG, Nat Genet 2020)", "Urothelial Carcinoma (BCAN/HCRN, Nat Commun 2022)", "Pancreatic Adenocarcinoma (TCGA, Firehose Legacy)", "Pheochromocytoma and Paraganglioma (TCGA, PanCancer Atlas)", "Lung Squamous Cell Carcinoma (TCGA, Nature 2012)", "Pancreatic Adenocarcinoma (TCGA, PanCancer Atlas)", "Uterine Corpus Endometrial Carcinoma (TCGA, Firehose Legacy)", "Acute Myeloid Leukemia (TCGA, Firehose Legacy)", "Acute Myeloid Leukemia (TCGA, NEJM 2013)", "Acute Myeloid Leukemia (TCGA, PanCancer Atlas)", "Glioblastoma Multiforme (TCGA, Firehose Legacy)", "Glioblastoma Multiforme (TCGA, PanCancer Atlas)", "The Metastatic Breast Cancer Project (Provisional, December 2021)", "Testicular Germ Cell Cancer (TCGA, Firehose Legacy)", "Glioblastoma (TCGA, Cell 2013)", "Testicular Germ Cell Tumors (TCGA, PanCancer Atlas)", "The Metastatic Breast Cancer Project (Archived, 2020)", "Pediatric Neuroblastoma (TARGET, 2018)", "Pancreatic Ductal Adenocarcinoma (CPTAC, Cell 2021)", "Pediatric Wilms Tumor (TARGET, 2018)", "Bladder Urothelial Carcinoma (TCGA, Nature 2014)", "Metastatic Breast Cancer (AURORA US Network, Nat Cancer 2023)", "Proteogenomic landscape of breast cancer (CPTAC, Cell 2020)", "Thymoma (TCGA, Firehose Legacy)", "Thymoma (TCGA, PanCancer Atlas)", "Metastatic Prostate Cancer (SU2C/PCF Dream Team, Cell 2015)", "Prostate Cancer (DKFZ, Cancer Cell 2018)", "Normal Keratinocytes (UCSF, BioRxiv 2024)", "Lung Adenocarcinoma (CPTAC, Cell 2020)", "Lung Squamous Cell Carcinoma (CPTAC, Cell 2021)", "Colon Cancer (CPTAC-2 Prospective, Cell 2019)", "Glioblastoma (CPTAC, Cell 2021)", "Pancreatic Adenocarcinoma (QCMG, Nature 2016)", "Endometrial Carcinoma (CPTAC, Cell 2020)", "Prostate Cancer MDA PCa PDX (MD Anderson, Clin Cancer Res 2024)", "Mesothelioma (TCGA, Firehose Legacy)", "Mesothelioma (TCGA, PanCancer Atlas)", "Uveal Melanoma (TCGA, PanCancer Atlas)", "Uveal Melanoma (TCGA, Firehose Legacy)", "Diffuse Glioma (GLASS Consortium, Nature 2019)", "Adrenocortical Carcinoma (TCGA, Firehose Legacy)", "Adrenocortical Carcinoma (TCGA, PanCancer Atlas)", "Brain Tumor PDXs (Mayo Clinic, Clin Cancer Res 2020)", "Kidney Chromophobe (TCGA, Firehose Legacy)", "Kidney Chromophobe (TCGA, Cancer Cell 2014)", "Kidney Chromophobe (TCGA, PanCancer Atlas)", "Uterine Carcinosarcoma (TCGA, PanCancer Atlas)", "Uterine Carcinosarcoma (TCGA, Firehose Legacy)", "Pre-cancer Colorectal Polyps (HTAN Vanderbilt, Cell 2021)", "Diffuse Large B-Cell Lymphoma (TCGA, PanCancer Atlas)", "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma (TCGA, Firehose Legacy)", "Pediatric Acute Myeloid Leukemia (TARGET, 2018)", "Pediatric Rhabdoid Tumor (TARGET, 2018)", "Metastatic Melanoma (DFCI, Science 2015)", "Cholangiocarcinoma (TCGA, Firehose Legacy)", "Cholangiocarcinoma (TCGA, PanCancer Atlas)", "Acral Melanoma (TGEN, Genome Res 2017)", "Upper Tract Urothelial Carcinoma (Cornell/Baylor/MDACC, Nat Commun 2019)", "Atypical Small Cell Lung Cancer (MSK, Cancer Discov 2024)", "Metastatic Melanoma (UCLA, Cell 2016)", "Normal Fibroblasts (UCSF, BioRxiv 2024)", "Schwannoma (Children s Tumor Foundation, Acta Neuropathologica 2020)", "Lung Cancer (SMC, Cancer Research 2016)", "Melanoma (MSK, NEJM 2014)", "Prostate Adenocarcinoma Organoids (MSK, Cell 2014)", "Acute Lymphoblastic Leukemia (St Jude, Nat Genet 2015)"]
# keys = ["(TCGA, GDC)"]
buscar_nombres_datasets_con_keyword(keywords=keys)
# buscar_estudios_con_rna_seq()
