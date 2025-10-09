# coding: utf-8
import pandas as pd
import os
from GOTool.GeneOntology import GeneOntology
import pickle

DATA_DIRECTORY = "/media/marcelo_baez/HD_Disc1/data/bacteria_selection/"
GOA_DIRECTORY = DATA_DIRECTORY + "selection_goa/"
if os.path.exists(DATA_DIRECTORY + "proteomes_df.pkl"):
    proteomes_df = pd.read_pickle(DATA_DIRECTORY + "proteomes_df.pkl")
else:
    url = "https://www.ebi.ac.uk/GOA/proteomes"
    url = "https://www.ebi.ac.uk/inc/drupal/goa/proteomes_release.html"
    proteomes_df = pd.read_html(url)[0]
    # proteomes_df.columns = proteomes_df.iloc[0]
    proteomes_df = proteomes_df.reindex(
        proteomes_df.index.drop(0),
    )
    proteomes_df["Tax ID"] = proteomes_df["Tax ID"].astype("str")
    proteomes_df.to_pickle(DATA_DIRECTORY + "proteomes_df.pkl")

if os.path.exists(DATA_DIRECTORY + "taxonomy_df.pkl"):
    taxonomy_df = pd.read_pickle(DATA_DIRECTORY + "taxonomy_df.pkl")
else:
    taxonomy_df = pd.read_csv(DATA_DIRECTORY + "proteomes-all.tab", sep="\t")
    # this is a tab separated file that
    # has to be downloaded from http://www.uniprot.org/proteomes/
    taxonomy_df.drop(columns=["Proteome ID", "Organism"], inplace=True)
    taxonomy_df.drop_duplicates(inplace=True)
    taxonomy_df["Organism ID"] = taxonomy_df["Organism ID"].astype("str")
    taxonomy_df.to_pickle(DATA_DIRECTORY + "taxonomy_df.pkl")

found_taxons = taxonomy_df[
    taxonomy_df["Organism ID"].isin(proteomes_df["Tax ID"].astype("str"))
]
proteomes_df = proteomes_df.merge(
    found_taxons, left_on="Tax ID", right_on="Organism ID"
)
proteomes_df = proteomes_df[
    proteomes_df["Total entries"] == proteomes_df["Protein count"]
]
proteomes_df["Superkingdom"] = proteomes_df["Taxonomic lineage"].str.partition(
    ","
)[  # noqa
    0
]  # noqa
bacteria_df = proteomes_df[proteomes_df["Superkingdom"] == "Bacteria"]

print(f"Downloading {bacteria_df.shape[0]} files...")
# ## Download GOAS

# create the selection_goa directory if it doesn't exist
if not os.path.exists(GOA_DIRECTORY):
    os.makedirs(GOA_DIRECTORY)
found = 0
not_found = []
for f in proteomes_df[proteomes_df["Superkingdom"] == "Bacteria"]["File"]:
    if not os.path.isfile(GOA_DIRECTORY + f):
        url = (
            "ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/{file}".format(  # noqa
                file=f
            )
        )  # noqa
        command = "wget -P " + GOA_DIRECTORY + ' "' + url + '"'
        not_found.append(command)
    else:
        found += 1
print(found, "proteomes found")
print(len(not_found), "proteomes not found")

# try to download all not found, but probably by hand is better
for command in not_found:
    os.system(command)


# ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/4116045.A_yellow_vein_Taiwan_virus-[Taiwan].goa
for f in proteomes_df[
    (proteomes_df["File"].isin([a.split("/")[-1] for a in not_found]))
    & (proteomes_df["Superkingdom"] == "Bacteria")
]["File"]:
    print('wget "ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/' + f + '"')  # noqa


bacteria_df.head()


if os.path.exists(DATA_DIRECTORY + "selection_count_df.pkl"):
    selection_count_df = pd.read_pickle(DATA_DIRECTORY + "selection_count_df.pkl")
else:
    go = GeneOntology(DATA_DIRECTORY + "go.obo")
    go.build_structure()
    data = []
    counter = 0
    for i, organism in bacteria_df.iterrows():
        d = {}
        go.load_annotation_file(
            GOA_DIRECTORY + organism["File"],
            organism["Organism"],
            GeneOntology.ALL_EVIDENCE_CODES,
        )
        go.load_annotation_file(
            GOA_DIRECTORY + organism["File"],
            organism["Organism"] + " exp",
            GeneOntology.EXPERIMENTAL_EVIDENCE_CODES,
        )
        go.up_propagate_annotations(organism["Organism"] + " exp")

        # any evidence code
        annotated_terms = [
            t for t in go.terms.values() if organism["Organism"] in t.annotations
        ]
        annotated_genes = set()
        for t in annotated_terms:
            annotated_genes |= set(t.annotations[organism["Organism"]].keys())

        d["Tax ID"] = organism["Tax ID"]
        d["all genes"] = len(annotated_genes)

        # experimental annotations
        org = organism["Organism"] + " exp"
        annotated_terms = [t for t in go.terms.values() if org in t.annotations]
        annotated_genes = set()
        annotations_by_domain: dict[str, set[str]] = {
            "biological_process": set(),
            "cellular_component": set(),
            "molecular_function": set(),
        }
        terms_by_domain: dict[str, set[str]] = {
            "biological_process": set(),
            "cellular_component": set(),
            "molecular_function": set(),
        }
        for t in annotated_terms:
            annotated_genes |= set(t.annotations[org].keys())
            annotations_by_domain[t.domain] |= set(t.annotations[org].keys())
            terms_by_domain[t.domain].add(t)

        # terms annotated with 3 genes or more
        popular_terms = [
            t for t in annotated_terms if len(t.annotations[org]) >= 3
        ]  # noqa
        # genes annotated to those terms
        popular_genes = set()
        popular_by_domain: dict[str, set[str]] = {
            "biological_process": set(),
            "cellular_component": set(),
            "molecular_function": set(),
        }
        popular_terms_by_domain: dict[str, set[str]] = {
            "biological_process": set(),
            "cellular_component": set(),
            "molecular_function": set(),
        }
        for t in popular_terms:
            popular_genes |= set(t.annotations[org].keys())
            popular_by_domain[t.domain] |= set(t.annotations[org].keys())
            popular_terms_by_domain[t.domain].add(t)

        d["exp annotations"] = len(annotated_genes)
        d["popular genes"] = len(popular_genes)
        d["annotated terms"] = len(annotated_terms)
        d["popular terms"] = len(popular_terms)
        d["terms bp"] = len(terms_by_domain["biological_process"])
        d["terms mf"] = len(terms_by_domain["molecular_function"])
        d["terms cc"] = len(terms_by_domain["cellular_component"])
        d["pop bp"] = len(popular_terms_by_domain["biological_process"])
        d["pop mf"] = len(popular_terms_by_domain["molecular_function"])
        d["pop cc"] = len(popular_terms_by_domain["cellular_component"])

        data.append(d)

        counter += 1
        if counter % 1000 == 0:
            print((float(counter) / len(bacteria_df["Tax ID"])) * 100.0, "%")
    selection_count_df = pd.DataFrame(data)
    selection_count_df.to_pickle(DATA_DIRECTORY + "selection_count_df.pkl")
    pickle.dump(data, open(os.path.join(DATA_DIRECTORY, "data.pkl"), "wb"))
selection_df = bacteria_df.merge(
    selection_count_df, left_on="Tax ID", right_on="Tax ID"
)
selection_df.to_pickle(DATA_DIRECTORY + "selection_df.pkl")

ontology_tau = 8
annotations_tau = 10
res = selection_df[
    (selection_df["exp annotations"] >= annotations_tau)
    & (selection_df["pop bp"] >= ontology_tau)
    & (selection_df["pop mf"] >= ontology_tau)
    & (selection_df["pop cc"] >= ontology_tau)
    #    (results_is_a_part_of['popular genes'] > 50) &
    #    (results_is_a_part_of['popular terms'] > 50)
][
    [
        "Tax ID",
        "Organism",
        "File",
        "all genes",
        "annotated terms",
        "exp annotations",
        "pop bp",
        "pop cc",
        "pop mf",
        "popular genes",
        "popular terms",
        "terms bp",
        "terms cc",
        "terms mf",
    ]
]

res_untrimmed = selection_df[
    (selection_df["exp annotations"] >= annotations_tau)
    & (selection_df["pop bp"] >= ontology_tau)
    & (selection_df["pop mf"] >= ontology_tau)
    & (selection_df["pop cc"] >= ontology_tau)
    #    (results_is_a_part_of['popular genes'] > 50) &
    #    (results_is_a_part_of['popular terms'] > 50)
]

res.to_pickle(os.path.join(DATA_DIRECTORY, "selected_bacteria_df.pkl"))
