import requests
import gzip
import io

from BCBio.GFF import GFFParser
from Bio.SeqFeature import SeqFeature
from Bio.SeqUtils import seq1
from mutalyzer_retriever.parsers.gff3 import _create_record_model

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


if __name__ == "__main__":
    # retrieve_gz(GRCH_37_ANNOTATIONS, GRCH_37_P)
    retrieve_gz(GRCH_38_ANNOTATIONS, GRCH_38_P)
    # records = extract_grch_37_records()
    # get_stats(records)
    # split_gff(GRCH_37_ANNOTATIONS, GRCH_37_P)
    split_gff(GRCH_38_ANNOTATIONS, GRCH_38_P)
