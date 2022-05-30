import requests
import gzip
import io

from BCBio.GFF import GFFParser
from Bio.SeqFeature import SeqFeature
from Bio.SeqUtils import seq1
from mutalyzer_retriever.parsers.gff3 import _create_record_model, parse

MAIN = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases"

GRCH_37_ANNOTATIONS = {"105.20190906", "105.20220307"}
GRCH_37_P = "GCF_000001405.25_GRCh37.p13"
GRCH_38_ANNOTATIONS = {
    "109.20190607",
    "109.20190905",
    "109.20191205",
    # "109.20200228",
    # "109.20200522",
    # "109.20200815",
    # "109.20201120",
    # "109.20210226",
    # "109.20210514",
    # "109.20211119",
}
GRCH_38_P = "GCF_000001405.39_GRCh38.p13"
ENDING = "_genomic.gff.gz"

SOURCES = {
    "GCF_000001405.39_GRCh38": {".p13"}
}


def retrieve_gz(annotations, patch):
    for annotation in annotations:
        url = f"{MAIN}/{annotation}/{patch}/{patch}{ENDING}"
        print(url)
        r = requests.get(url)
        open(f"{patch}_{annotation}{ENDING}", "wb").write(r.content)


def extract_grch_37_records():
    gff_parser = GFFParser()
    records = {}
    for annotation_id in GRCH_37_ANNOTATIONS:
        with gzip.open(f"{GRCH_37_P}_{annotation_id}{ENDING}", "rb") as f:
            content = f.read().decode()
            gff = gff_parser.parse(io.StringIO(content))
            print(f"gff for: {annotation_id}")
            for record in gff:
                if record.id == "NC_000001.10":
                    parsed_record = _create_record_model(record)
                    if parsed_record["id"] not in records:
                        records[parsed_record["id"]] = {}
                    if parsed_record.get("features"):
                        for gene in parsed_record["features"]:
                            if gene.get("features"):
                                for feature in gene["features"]:
                                    if (
                                        feature["id"]
                                        not in records[parsed_record["id"]]
                                    ):
                                        records[parsed_record["id"]][feature["id"]] = {}
                                    records[parsed_record["id"]][feature["id"]][
                                        annotation_id
                                    ] = feature
                else:
                    print(f"  skipping: {record.id}")
    # print(len(records))
    # for r_id in records:
    #     print(f"{r_id}: {len(records[r_id])}")
    return records


def _location(annotations):
    locations = [
        (
            annotations[i]["location"]["start"]["position"],
            annotations[i]["location"]["start"]["position"],
        )
        for i in annotations
    ]
    return locations[:-1] == locations[1:]


def _exons(annotations):
    exons = []
    for annotation in annotations:
        exon = []
    locations = [
        (
            annotations[i]["location"]["start"]["position"],
            annotations[i]["location"]["start"]["position"],
        )
        for i in annotations
    ]
    return locations[:-1] == locations[1:]


def get_stats(records):
    not_equal = {}
    for r_id in records:
        for f_id in records[r_id]:
            if len(records[r_id][f_id]) == len(GRCH_37_ANNOTATIONS):
                if _location(records[r_id][f_id]):
                    records[r_id][f_id] = {"equal": True}
                else:
                    print("  - different:", r_id, f_id)
                    if r_id not in not_equal:
                        not_equal[r_id] = []
                    if f_id not in not_equal[r_id]:
                        not_equal[r_id].append(f_id)

    for r_id in records:
        to_print = f"{r_id}: {len(records[r_id])}"
        if r_id in not_equal:
            to_print += f"; not_equal: {len(not_equal[r_id])}"
        print(to_print)


def split_gff(annotations, patch):
    for annotation_id in annotations:
        with gzip.open(f"{patch}_{annotation_id}{ENDING}", "rb") as f:
            print(annotation_id)
            current_id = ""
            current_content = ""
            extras = ""
            for line in f:
                s_line = line.decode()
                if s_line.startswith("#!"):
                    extras += s_line
                elif s_line.startswith("##sequence-region"):
                    if current_id:
                        open(f"split/{current_id}_{annotation_id}.gff3", "w").write(
                            current_content
                        )
                    current_id = s_line.split(" ")[1]
                    current_content = f"##gff-version 3\n{extras}{s_line}"
                elif s_line.startswith("##species") or s_line.startswith(current_id):
                    current_content += s_line


def _hgnc_check(hgnc, features):
    for gene in features.get("features"):
        if gene.get("qualifiers") and gene["qualifiers"].get("HGNC") and gene["qualifiers"]["HGNC"] == hgnc:
            return gene["id"]


def _synonym_check(synonyms, features):
    print(synonyms)
    for synonym in synonyms:
        for i, gene in enumerate(features.get("features")):
            if synonym == gene["id"]:
                return gene["id"], i
            if gene.get("qualifiers") and gene["qualifiers"].get("synonym") and synonym in gene["qualifiers"]["synonym"]:
                return gene["id"], i
    return None, None


def merge_models(new, old):
    print(len(new.get("features")))
    print(len(old.get("features")))

    new_genes = {}
    old_genes = {}
    for i, gene in enumerate(new.get("features")):
        new_genes[gene["id"]] = i
    for i, gene in enumerate(old.get("features")):
        old_genes[gene["id"]] = i
    not_in = set()
    in_different = set()
    in_same = set()
    for i, gene in enumerate(old.get("features")):
        if gene["id"] not in new_genes:
            not_in.add((gene['id'], i))
            # print(gene["qualifiers"])
            # if gene.get("features"):
            #     for f in gene.get("features"):
            #         print(f["id"])
            # else:
            #     print("no features")
        else:
            if gene != new["features"][new_genes[gene["id"]]]:
                in_different.add(gene['id'])
            else:
                in_same.add(gene['id'])
    hgnc_found = set()
    for i in not_in:
        gene = old["features"][i[1]]
        if gene["qualifiers"].get("HGNC"):
            new_gene_id = _hgnc_check(gene["qualifiers"]["HGNC"], new)
            if new_gene_id:
                print(f"not in, but found through hgnc: {gene['id']}, {gene['qualifiers']['HGNC']}, {new_gene_id}")
                hgnc_found.add(i)

    not_in = not_in - hgnc_found
    synonym_found = set()
    for i in not_in:
        gene = old["features"][i[1]]
        synonyms = [gene["id"]]
        if gene["qualifiers"].get("synonym"):
            synonyms += gene["qualifiers"]["synonym"]
        new_gene_id, n_i = _synonym_check(synonyms, new)
        if new_gene_id:
            print(f"not in, but found through synonyms: {gene['id']}, {synonyms}, {new_gene_id}, {new['features'][n_i]}")
            synonym_found.add(i)
    not_in = not_in - synonym_found
    print(f" - not in: {len(not_in)}")
    print(f" - HGNC found: {len(hgnc_found)}")
    print(f" - synonym found: {len(synonym_found)}")
    print(f" - in_different: {len(in_different)}")
    print(f" - in_same: {len(in_same)}")
    print(set(new_genes) - set(old_genes))


def sandbox():
    out = {}
    files = [
        "GCF_000001405.25_GRCh37.p13_105.20190906_genomic.gff.gz",
        "GCF_000001405.25_GRCh37.p13_105.20220307_genomic.gff.gz",
        "GCF_000001405.39_GRCh38.p13_109.20190607_genomic.gff.gz",
        "GCF_000001405.39_GRCh38.p13_109.20190905_genomic.gff.gz"
    ]
    for file in files:
        assembly_id = file.split(".p")[0]
        patch_id = file.split(".")[2].split("_")[0]
        annotation_id = file.split(".")[3].split("_")[0]
        print(assembly_id, patch_id, annotation_id)

        with gzip.open(file, "rb") as f:
            current_id = ""
            current_content = ""
            extras = ""
            for line in f:
                s_line = line.decode()
                if s_line.startswith("#!"):
                    extras += s_line
                elif s_line.startswith("##sequence-region"):
                    if current_id and current_id in ["NC_000001.10", "NC_000001.11"]:
                        current_model = parse(current_content)
                        print(current_id)
                        print(current_model["qualifiers"])
                        if current_id not in out:
                            out[current_id] = current_model
                        else:
                            merge_models(current_model, out[current_id])
                    current_id = s_line.split(" ")[1]
                    current_content = f"##gff-version 3\n{extras}{s_line}"
                elif s_line.startswith("##species") or s_line.startswith(current_id):
                    current_content += s_line


if __name__ == "__main__":
    # retrieve_gz(GRCH_37_ANNOTATIONS, GRCH_37_P)
    # retrieve_gz(GRCH_38_ANNOTATIONS, GRCH_38_P)
    # records = extract_grch_37_records()
    # get_stats(records)
    # split_gff(GRCH_37_ANNOTATIONS, GRCH_37_P)
    # split_gff(GRCH_38_ANNOTATIONS, GRCH_38_P)
    sandbox()
