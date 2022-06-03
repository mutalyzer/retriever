import gzip
import io
import json
import xml.etree.ElementTree as ET
from copy import deepcopy
from ftplib import FTP, error_perm

import requests
from BCBio.GFF import GFFParser

from mutalyzer_retriever.parsers.gff3 import _create_record_model, parse

MAIN = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases"

GRCH_37_ANNOTATIONS = {"105.20190906", "105.20220307"}
GRCH_37_P = "GCF_000001405.25_GRCh37.p13"
GRCH_38_ANNOTATIONS = {
    "109.20190607",
    "109.20190905",
    "109.20191205",
    "109.20200228",
    "109.20200522",
    "109.20200815",
    "109.20201120",
    "109.20210226",
    "109.20210514",
    "109.20211119",
}
GRCH_38_P = "GCF_000001405.39_GRCh38.p13"
ENDING = "_genomic.gff.gz"

SOURCES = {"GCF_000001405.39_GRCh38": {".p13"}}


# def retrieve_gz(annotations, patch):
#     for annotation in annotations:
#         url = f"{MAIN}/{annotation}/{patch}/{patch}{ENDING}"
#         print(url)
#         r = requests.get(url)
#         open(f"{patch}_{annotation}{ENDING}", "wb").write(r.content)


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
        if (
            gene.get("qualifiers")
            and gene["qualifiers"].get("HGNC")
            and gene["qualifiers"]["HGNC"] == hgnc
        ):
            return gene["id"]


def _synonym_check(synonyms, features):
    for synonym in synonyms:
        for gene in features.get("features"):
            if synonym == gene["id"]:
                return gene["id"]
            if (
                gene.get("qualifiers")
                and gene["qualifiers"].get("synonym")
                and synonym in gene["qualifiers"]["synonym"]
            ):
                return gene["id"]
    return None


def _get_gene_transcript_ids(gene):
    transcripts = []
    if gene.get("features"):
        for feature in gene["features"]:
            transcripts.append(feature["id"])
    return transcripts


def _get_transcript_ids(model):
    transcripts = []
    if model.get("features"):
        for feature in model["features"]:
            transcripts += _get_gene_transcript_ids(feature)
    return transcripts


def _in(item_id, features):
    for feature in features:
        if feature.get("id") and feature["id"] == item_id:
            return True
    return False


def _get_genes_for_transcript(t_id, model):
    genes = []
    if model.get("features"):
        for feature in model["features"]:
            if feature.get("features"):
                for transcript in feature["features"]:
                    if transcript["id"] == t_id:
                        genes.append(feature)
    return genes


def _merge_genes(gene_new, gene_old, already_new_transcripts=[]):
    if gene_old.get("features") and gene_new.get("features"):
        old_transcripts = {t["id"]: i for i, t in enumerate(gene_old["features"])}
        new_transcripts = {t["id"]: i for i, t in enumerate(gene_new["features"])}
        for t_id in set(old_transcripts) - set(new_transcripts):
            if already_new_transcripts and t_id in already_new_transcripts:
                continue
            gene_new["features"].append(gene_old["features"][old_transcripts[t_id]])
        return set(old_transcripts) - set(new_transcripts)
    return set()


def merge(new, old):
    def _genes_summary():
        for idx, gene in enumerate(new.get("features")):
            genes_new[gene["id"]] = idx
        for idx, gene in enumerate(old.get("features")):
            genes_old[gene["id"]] = idx
        for gene in old.get("features"):
            if gene["id"] not in genes_new:
                genes_old_not_in.add(gene["id"])
            else:
                genes_old_in.add(gene["id"])

    def _merge_same_hgnc():
        for gene_old_id in genes_old_not_in:
            gene_old = old["features"][genes_old[gene_old_id]]
            if gene_old["qualifiers"].get("HGNC"):
                gene_new_id = _hgnc_check(gene_old["qualifiers"]["HGNC"], new)
                if gene_new_id:
                    hgnc_found.add(gene_old_id)
                    genes_old_in_new_id[gene_old_id] = gene_new_id
                    transcripts_old_added.update(
                        _merge_genes(
                            new["features"][genes_new[gene_new_id]],
                            gene_old,
                            _get_transcript_ids(new),
                        )
                    )
                    # print("\n")
                    # print(json.dumps(old["features"][genes_old[gene_old_id]],
                    #                  indent=2))
                    # print("=====")
                    # print(json.dumps(new["features"][genes_new[gene_new_id]],
                    #                  indent=2))
                    # print("=====")

    def _synonym_lookup():
        for gene_old_id in genes_old_not_in:
            gene = old["features"][genes_old[gene_old_id]]
            synonyms = [gene["id"]]
            if gene["qualifiers"].get("synonym"):
                synonyms += gene["qualifiers"]["synonym"]
            gene_new_id = _synonym_check(synonyms, new)
            if gene_new_id:
                synonym_found.add(gene_old_id)
                # print("\n")
                # print(json.dumps(old["features"][genes_old[gene_old_id]], indent=2))
                # print("=====")
                # print(json.dumps(new["features"][genes_new[gene_new_id]], indent=2))
                # print("=====")

    def _transcript_lookup():
        for gene_old_id in genes_old_not_in:
            gene_old = old["features"][genes_old[gene_old_id]]
            transcripts_old = _get_gene_transcript_ids(gene_old)
            for gene_new in new.get("features"):
                transcripts_new = _get_gene_transcript_ids(gene_new)
                transcripts_intersection = set(transcripts_old).intersection(
                    set(transcripts_new)
                )
                if transcripts_intersection:
                    if (
                        len(transcripts_intersection) == len(transcripts_new)
                        or set(transcripts_old) == transcripts_intersection
                    ):
                        transcript_found.add(gene_old_id)
                    # else:
                    #     print(
                    #         f"  - old gene {gene_old['id']} - new gene {gene_new['id']}:\n    {transcripts_old}\n    {transcripts_new}"
                    #     )

                    # print(f"  - old gene {gene_old['id']} - new gene {gene_new['id']} transcripts overlap: {transcripts_intersection}")
                    # print("---")
                    # print(json.dumps(gene_old, indent=2))
                    # print("---")
                    # print(json.dumps(gene_new, indent=2))

    def _merge():
        for gene_id in genes_old_in:
            gene_old = old["features"][genes_old[gene_id]]
            gene_new = new["features"][genes_new[gene_id]]
            if gene_new == gene_old:
                genes_equal.add(gene_id)
            else:
                genes_different.add(gene_id)
                transcripts_old_added.update(
                    _merge_genes(gene_new, gene_old, _get_transcript_ids(new))
                )

    def _add_not_in():
        already_new_transcripts = _get_transcript_ids(new)
        for gene_old_id in genes_old_not_in:
            # print("add", gene_old_id, genes_old[gene_old_id])
            gene_not_in = deepcopy(old["features"][genes_old[gene_old_id]])
            if gene_not_in.get("features"):
                to_del = []
                for i, feature in enumerate(gene_not_in["features"]):
                    if feature["id"] in already_new_transcripts:
                        to_del.append(i)
                for i in reversed(to_del):
                    del gene_not_in["features"][i]

            new["features"].append(gene_not_in)

    print("------")
    # print(f"new genes: {len(new.get('features'))}")
    # print(f"old genes: {len(old.get('features'))}")

    genes_new = {}
    genes_old = {}

    genes_old_not_in = set()
    genes_old_in = set()
    genes_old_in_new_id = {}
    hgnc_found = set()
    synonym_found = set()
    transcript_found = set()
    genes_equal = set()
    genes_different = set()
    transcripts_old_added = set()

    _genes_summary()

    transcripts = _get_transcript_ids(new)
    transcripts_start = f"{len(transcripts)} vs {len(set(transcripts))}"

    _merge_same_hgnc()
    genes_old_not_in -= hgnc_found

    _synonym_lookup()
    # genes_old_not_in = genes_old_not_in - synonym_found

    _transcript_lookup()
    genes_old_not_in -= transcript_found

    transcripts = _get_transcript_ids(new)
    transcripts_before_merge = f"{len(transcripts)} vs {len(set(transcripts))}"

    _merge()

    transcripts = _get_transcript_ids(new)
    transcripts_before_not_in = f"{len(transcripts)} vs {len(set(transcripts))}"

    # seen = set()
    # dupes = [x for x in transcripts if x in seen or seen.add(x)]
    # print(f"\ndupes: {dupes}\n")

    # for dupe in dupes:
    #     print(f" - {dupe}")
    # for gene in _get_genes_for_transcript(dupe, new):
    #     print(json.dumps(gene, indent=2))
    #     print("--=-=-")
    # print(_get_genes_for_transcript(dupe, new))
    # print("----")
    # print(_get_genes_for_transcript(dupe, old))

    _add_not_in()

    transcripts = _get_transcript_ids(new)
    transcripts_after_not_in = f"{len(transcripts)} vs {len(set(transcripts))}"

    print("----")
    print(f"- not in genes: {len(genes_old_not_in)}")
    print(f"- HGNC found genes: {len(hgnc_found)}")
    print(f"- synonym found genes: {len(synonym_found)}")
    print(f"- transcript found genes: {len(transcript_found)}")
    print(f"- genes in in total: {len(genes_old_in)}")
    print(f"  - equal: {len(genes_equal)}")
    print(f"  - different: {len(genes_different)}")
    print(f"    - old transcripts added: {len(transcripts_old_added)}")
    print(f" - transcripts start: {transcripts_start}")
    print(f" - transcripts before merge: {transcripts_before_merge}")
    print(f" - transcripts before not in: {transcripts_before_not_in}")
    print(f" - transcripts after adding not in: {transcripts_after_not_in}")
    print("------")


def get_models():
    out = {}
    files = [
        # "GCF_000001405.25_GRCh37.p13_105.20190906_genomic.gff.gz",
        # "GCF_000001405.25_GRCh37.p13_105.20220307_genomic.gff.gz",
        "GCF_000001405.38_GRCh38.p12_109_genomic.gff.gz",
        "GCF_000001405.39_GRCh38.p13_109.20190607_genomic.gff.gz",
        "GCF_000001405.39_GRCh38.p13_109.20190905_genomic.gff.gz",
        "GCF_000001405.39_GRCh38.p13_109.20191205_genomic.gff.gz",
        "GCF_000001405.39_GRCh38.p13_109.20200228_genomic.gff.gz",
        "GCF_000001405.39_GRCh38.p13_109.20200522_genomic.gff.gz",
        "GCF_000001405.39_GRCh38.p13_109.20200815_genomic.gff.gz",
        "GCF_000001405.39_GRCh38.p13_109.20201120_genomic.gff.gz",
        "GCF_000001405.39_GRCh38.p13_109.20210226_genomic.gff.gz",
        "GCF_000001405.39_GRCh38.p13_109.20210514_genomic.gff.gz",
        "GCF_000001405.39_GRCh38.p13_109.20211119_genomic.gff.gz",
        "GCF_000001405.40_GRCh38.p14_110_genomic.gff.gz",
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
                    # if current_id and current_id in ["NC_000001.10", "NC_000001.11", "NC_000024.10"]:
                    # if current_id and current_id in ["NC_000001.11", "NC_000002.12"]:
                    if current_id and current_id in ["NC_000001.11"]:
                        # if current_id and current_id in ["NC_000024.10"]:
                        # if current_id and current_id in ["NC_000003.12"]:
                        current_model = parse(current_content)
                        print(current_id)
                        # print(current_model["qualifiers"])
                        if current_id not in out:
                            out[current_id] = current_model
                        else:
                            merge(current_model, out[current_id])
                            out[current_id] = current_model

                    current_id = s_line.split(" ")[1]
                    current_content = f"##gff-version 3\n{extras}{s_line}"
                elif s_line.startswith("##species") or s_line.startswith(current_id):
                    current_content += s_line
    print("\n\n\n\n")
    for r_id in out:
        print(r_id)
        transcripts = _get_transcript_ids(out[r_id])
        print(f"{len(transcripts)} vs {len(set(transcripts))}")
        # print(len(_get_transcript_ids(out[r_id])))


def get_ftp_locations():
    annotations = []
    main_url = "ftp.ncbi.nlm.nih.gov"
    main_dir = "genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases"
    with FTP(main_url) as ftp:
        ftp.login()
        ftp.cwd(main_dir)
        for d_a in ftp.nlst():
            try:
                ftp.cwd(d_a)
            except error_perm:
                continue
            annotation = {"id": d_a}
            # print(f" - {d_a}")
            for d_d in ftp.nlst():
                # print(f"   - {d_d}")
                if d_d.endswith("annotation_report.xml"):
                    annotation["annotation_report"] = d_d
                if d_d.startswith("GCF_") and "GRCh" in d_d:
                    annotation["dir"] = d_d
                    try:
                        ftp.cwd(d_d)
                    except error_perm:
                        continue
                    for d_f in ftp.nlst():
                        # print(f"     - {d_f}")
                        if d_f.endswith("_genomic.gff.gz"):
                            annotation["file"] = d_f
                    ftp.cwd("..")
            ftp.cwd("..")
            annotations.append(annotation)
    return {"url": main_url, "dir": main_dir, "annotations": annotations}


def _gff_file_name(annotation):
    return (
        annotation["file"].split("_genomic")[0]
        + "_"
        + annotation["id"]
        + "_genomic.gff.gz"
    )


def _report_file_name(annotation):
    return annotation["annotation_report"]


def retrieve_files(annotations):
    main_url = "https://" + annotations["url"] + "/" + annotations["dir"]
    for annotation in annotations["annotations"]:
        url = f"{main_url}/{annotation['id']}/{annotation['dir']}/{annotation['file']}"
        print(url)
        r = requests.get(url)
        open(_gff_file_name(annotation), "wb").write(r.content)

        url = f"{main_url}/{annotation['id']}/{annotation['annotation_report']}"
        print(url)
        r = requests.get(url)
        open(_report_file_name(annotation), "wb").write(r.content)


def _report(report_file):
    tree = ET.parse(report_file)
    root = tree.getroot()
    return {
        "freeze_date_id": root.find("./BuildInfo/FreezeDateId").text,
        "assembly_name": root.find("./AssembliesReport/FullAssembly/Name").text,
        "assembly_accession": root.find(
            "./AssembliesReport/FullAssembly/Accession"
        ).text,
    }


def report_updates(annotations):
    for annotation in annotations["annotations"]:
        annotation.update(_report(_report_file_name(annotation)))


def group_by_accession(annotations):
    groups = {}
    sorted_annotations = sorted(
        annotations["annotations"], key=lambda d: d["freeze_date_id"]
    )
    for annotation in sorted_annotations:
        assembly = annotation["assembly_name"].split(".")[0]
        if assembly not in groups:
            groups[assembly] = []
        groups[assembly].append(annotation)
    return groups


def main():
    # retrieve_gz(GRCH_37_ANNOTATIONS, GRCH_37_P)
    # retrieve_gz(GRCH_38_ANNOTATIONS, GRCH_38_P)
    # records = extract_grch_37_records()
    # get_stats(records)
    # split_gff(GRCH_37_ANNOTATIONS, GRCH_37_P)
    # split_gff(GRCH_38_ANNOTATIONS, GRCH_38_P)

    # annotations = get_ftp_locations()
    # open("annotations.json", "w").write(json.dumps(annotations, indent=2))
    # print(json.dumps(annotations, indent=2))
    # print(json.dumps(annotations, indent=2))
    # retrieve_files(annotations)
    # get_models()

    annotations = json.loads(open("annotations.json", "r").read())
    report_updates(annotations)
    # print(json.dumps(annotations, indent=2))
    print(json.dumps(group_by_accession(annotations), indent=2))
    # _report("Homo_sapiens_AR105.20220307_annotation_report.xml")


if __name__ == "__main__":
    main()
