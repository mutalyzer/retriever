import gzip
import json
import xml.etree.ElementTree as ET
from copy import deepcopy
from pathlib import Path

import requests

from mutalyzer_retriever.parser import parse

from ..parser import parse
from .ncbi import fetch


def _common_url():
    return "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/"


def _annotations_urls():
    common = _common_url()
    annotations_urls = {
        "GRCh37": [
            [
                common + "105.20190906/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gff.gz",
                common + "105.20190906/Homo_sapiens_AR105_annotation_report.xml",
            ],
            [
                common + "105.20220307/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gff.gz",
                common + "105.20220307/Homo_sapiens_AR105.20220307_annotation_report.xml",
            ],
            [
                common + "GCF_000001405.25-RS_2024_09/GCF_000001405.25_GRCh37.p13_genomic.gff.gz",
                common + "GCF_000001405.25-RS_2024_09/GCF_000001405.25-RS_2024_09_annotation_report.xml",
            ],
        ],
        "GRCh38": [
            [
                common + "109/GCF_000001405.38_GRCh38.p12/GCF_000001405.38_GRCh38.p12_genomic.gff.gz",
                common + "109/GCF_000001405.38_GRCh38.p12/Homo_sapiens_AR109_annotation_report.xml"
            ],
            [
                common + "109.20190607/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz",
                common + "109.20190607/GCF_000001405.39_GRCh38.p13/Homo_sapiens_AR109.20190607_annotation_report.xml"
            ],
            [
                common + "109.20190905/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz",
                common + "109.20190905/GCF_000001405.39_GRCh38.p13/Homo_sapiens_AR109.20190905_annotation_report.xml"
            ],
            [
                common + "109.20191205/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz",
                common + "109.20191205/GCF_000001405.39_GRCh38.p13/Homo_sapiens_AR109.20191205_annotation_report.xml"
            ],
            [
                common + "109.20200228/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz",
                common + "109.20200228/GCF_000001405.39_GRCh38.p13/Homo_sapiens_AR109.20200228_annotation_report.xml"
            ],
            [
                common + "109.20200522/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz",
                common + "109.20200522/GCF_000001405.39_GRCh38.p13/Homo_sapiens_AR109.20200522_annotation_report.xml"
            ],
            [
                common + "109.20200815/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz",
                common + "109.20200815/GCF_000001405.39_GRCh38.p13/Homo_sapiens_AR109.20200815_annotation_report.xml"
            ],
            [
                common + "109.20201120/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz",
                common + "109.20201120/GCF_000001405.39_GRCh38.p13/Homo_sapiens_AR109.20201120_annotation_report.xml"
            ],
            [
                common + "109.20210226/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz",
                common + "109.20210226/GCF_000001405.39_GRCh38.p13/Homo_sapiens_AR109.20210226_annotation_report.xml"
            ],
            [
                common + "109.20210514/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz",
                common + "109.20210514/GCF_000001405.39_GRCh38.p13/Homo_sapiens_AR109.20210514_annotation_report.xml"
            ],
            [
                common + "109.20211119/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz",
                common + "109.20211119/GCF_000001405.39_GRCh38.p13/Homo_sapiens_AR109.20211119_annotation_report.xml"
            ],
            [
                common + "110/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz",
                common + "110/Homo_sapiens_AR110_annotation_report.xml"
            ],
            [
                common + "GCF_000001405.40-RS_2023_03/GCF_000001405.40_GRCh38.p14_genomic.gff.gz",
                common + "GCF_000001405.40-RS_2023_03/GCF_000001405.40-RS_2023_03_annotation_report.xml"
            ],
            [
                common + "GCF_000001405.40-RS_2023_10/GCF_000001405.40_GRCh38.p14_genomic.gff.gz",
                common + "GCF_000001405.40-RS_2023_10/GCF_000001405.40-RS_2023_10_annotation_report.xml"
            ],
            [
                common + "GCF_000001405.40-RS_2024_08/GCF_000001405.40_GRCh38.p14_genomic.gff.gz",
                common + "GCF_000001405.40-RS_2024_08/GCF_000001405.40-RS_2024_08_annotation_report.xml"
            ]
        ],
        "T2T-CHM13v2": [
            [
                common + "GCF_009914755.1-RS_2023_03/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz",
                common + "GCF_009914755.1-RS_2023_03/GCF_009914755.1-RS_2023_03_annotation_report.xml"
            ],
            [
                common + "GCF_009914755.1-RS_2023_10/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz",
                common + "GCF_009914755.1-RS_2023_10/GCF_009914755.1-RS_2023_10_annotation_report.xml"
            ],
            [
                common + "GCF_009914755.1-RS_2024_08/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz",
                common + "GCF_009914755.1-RS_2024_08/GCF_009914755.1-RS_2024_08_annotation_report.xml"
            ]
        ]
    }
    return annotations_urls


def _report_info(xml_content):
    tree = ET.ElementTree(ET.fromstring(xml_content))
    root = tree.getroot()
    return {
        "freeze_date_id": root.find("./BuildInfo/FreezeDateId").text,
        "full_assembly_name": root.find("./AssembliesReport/FullAssembly/Name").text,
        "full_assembly_accession": root.find(
            "./AssembliesReport/FullAssembly/Accession"
        ).text,
    }


def download_annotation_releases(urls, directory="./ncbi_annotation_releases", assembly_id_start=None):
    """
    Download the annotations (GFF3) report (XML) files in the provided directory.
    """
    print("Downloading assembly releases:")
    path_dir = Path(directory)
    path_dir.mkdir(parents=True, exist_ok=True)

    metadata = {}

    for assembly_id in urls:
        if assembly_id_start is not None and not assembly_id.startswith(assembly_id_start):
            continue
        print(f" - assembly: {assembly_id}")
        metadata[assembly_id] = {}
        for url_gff3, url_report in urls[assembly_id]:
            response_xml = requests.get(url_report)
            if response_xml.status_code == 200:
                freeze_date_id = _report_info(response_xml.content)["freeze_date_id"]
                response_gff3 = requests.get(url_gff3)
                if response_gff3.status_code == 200:
                    path_dir = Path(directory) / assembly_id / freeze_date_id
                    print(f"   - dir: {path_dir}")

                    file_name_gff3 = url_gff3.split("/")[-1]
                    file_name_report = url_report.split("/")[-1]
                    path_dir.mkdir(parents=True, exist_ok=True)

                    path_file_xml = path_dir / file_name_report
                    open(path_file_xml, "wb").write(response_xml.content)

                    path_file_gff3 = path_dir / file_name_gff3
                    open(path_file_gff3, "wb").write(response_gff3.content)

                    metadata[assembly_id][freeze_date_id] = {
                        "xml": str(path_file_xml),
                        "gff3": str(path_file_gff3),
                    }
    return metadata


def _get_gene(g_id, model):
    if model.get("features"):
        for gene in model["features"]:
            if gene["id"] == g_id:
                return gene


def _get_gene_i(g_id, model):
    if model.get("features"):
        for i, gene in enumerate(model["features"]):
            if gene["id"] == g_id:
                return i


def _get_gene_transcript_ids(gene):
    transcripts = []
    if gene.get("features"):
        for feature in gene["features"]:
            transcripts.append(feature["id"])
    return transcripts


def _get_transcripts_mappings(model):
    transcripts = {}
    if model.get("features"):
        for i_g, gene in enumerate(model["features"]):
            if gene.get("features"):
                for i_t, transcript in enumerate(gene["features"]):
                    if transcript.get("id") is None:
                        continue
                    elif transcript["id"] in transcripts:
                        raise Exception(
                            f"Multiple transcripts with same id ({transcript['id']}) in model."
                        )
                    else:
                        transcripts[transcript["id"]] = {
                            "i_g": i_g,
                            "gene_id": gene["id"],
                            "i_t": i_t,
                        }

    return transcripts


def _added_from(feature, model):
    if feature.get("qualifiers") is None:
        feature["qualifiers"] = {"annotation_added_from": {}}
    if feature.get("qualifiers").get("annotation_added_from") is None:
        feature["qualifiers"]["annotation_added_from"] = {}
    if feature["qualifiers"]["annotation_added_from"].get("freeze_date_id") is None:
        feature["qualifiers"]["annotation_added_from"]["freeze_date_id"] = model[
            "qualifiers"
        ]["annotations"]["freeze_date_id"]
    if feature["qualifiers"]["annotation_added_from"].get("full_assembly_name") is None:
        feature["qualifiers"]["annotation_added_from"]["full_assembly_name"] = model[
            "qualifiers"
        ]["annotations"]["full_assembly_name"]
    if (
        feature["qualifiers"]["annotation_added_from"].get("full_assembly_accession")
        is None
    ):
        feature["qualifiers"]["annotation_added_from"]["full_assembly_accession"] = (
            model["qualifiers"]["annotations"]["full_assembly_accession"]
        )


def _gene_added_from(gene, model):
    _added_from(gene, model)
    if gene.get("features"):
        for transcript in gene["features"]:
            _added_from(transcript, model)


def _merge(new, old):
    ts_new = _get_transcripts_mappings(new)
    ts_old = _get_transcripts_mappings(old)

    ts_not_in = set(ts_old.keys()) - set(ts_new.keys())

    for t_not_in_id in ts_not_in:
        if t_not_in_id in ts_new:
            continue
        gene_new = _get_gene(ts_old[t_not_in_id]["gene_id"], new)
        if not gene_new:
            gene_old = deepcopy(_get_gene(ts_old[t_not_in_id]["gene_id"], old))
            gene_ts = _get_gene_transcript_ids(gene_old)
            gene_ts_already_in = []
            for i, t in enumerate(gene_ts):
                if t in ts_new:
                    gene_ts_already_in.append(i)
            for i in gene_ts_already_in[::-1]:
                gene_old["features"].pop(i)
            _gene_added_from(gene_old, old)
            if new.get("features") is None:
                new["features"] = []
            new["features"].append(gene_old)
            for t in set(gene_ts) - set(gene_ts_already_in):
                ts_new[t] = {"i_g": len(new["features"]), "gene_id": gene_old["id"]}
        else:
            transcript = old["features"][ts_old[t_not_in_id]["i_g"]]["features"][
                ts_old[t_not_in_id]["i_t"]
            ]
            _added_from(transcript, old)
            if gene_new.get("features") is None:
                gene_new["features"] = []
            gene_new["features"].append(deepcopy(transcript))
            ts_new[t_not_in_id] = {
                "i_g": _get_gene_i(ts_old[t_not_in_id]["gene_id"], new),
                "gene_id": gene_new["id"],
            }

    ts_not_in = set(ts_old.keys()) - set(ts_new.keys())
    if len(ts_not_in) != 0:
        raise Exception("Not all the transcripts were added.")


def _directory_metadata(directory="./ncbi_annotation_releases"):
    m = {}
    for path_assembly in Path(directory).iterdir():
        assembly = str(path_assembly).split("/")[-1]
        m[assembly] = {}
        for path_date in Path(path_assembly).iterdir():
            date = str(path_date).split("/")[-1]
            m[assembly][date] = {}
            for path_file in Path(path_date).iterdir():
                f = str(path_file).split("/")[-1]
                if f.endswith(".gff.gz"):
                    m[assembly][date]["gff3"] = str(path_file)
                elif f.endswith(".xml"):
                    m[assembly][date]["xml"] = str(path_file)
    return m


def get_annotation_models(
    directory_input="./ncbi_annotation_releases",
    assembly_id_start=None,
    ref_id_start=None,
):
    print("Get annotation models:")
    metadata = _directory_metadata(directory_input)
    out = {}
    for assembly in metadata:
        models = {}
        assemblies = []
        ref_ids = []
        if assembly_id_start is not None and not assembly.startswith(assembly_id_start):
            continue
        for freeze_date_id in sorted(metadata[assembly]):
            print(f"- get from: {assembly}, date: {freeze_date_id}")
            assembly_details = _report_info(open(metadata[assembly][freeze_date_id]["xml"]).read().strip())
            assemblies.append(assembly_details)

            with gzip.open(metadata[assembly][freeze_date_id]["gff3"], "rb") as f:
                current_id = ""
                current_content = ""
                extras = ""
                for line in f:
                    s_line = line.decode()
                    if s_line.startswith("#!"):
                        extras += s_line
                    elif s_line.startswith("##sequence-region"):
                        if current_id and (ref_id_start is None or current_id.startswith(ref_id_start)):
                            current_model = parse(current_content, "gff3")
                            current_model["qualifiers"]["annotations"] = assembly_details
                            print(f"  - {current_id}")
                            if current_id not in models:
                                models[current_id] = current_model
                            else:
                                _merge(current_model, models[current_id])
                                models[current_id] = current_model
                            ref_ids.append(current_id)

                        current_id = s_line.split(" ")[1]
                        current_content = f"##gff-version 3\n{extras}{s_line}"
                    elif s_line.startswith("##species") or s_line.startswith(
                        current_id
                    ):
                        current_content += s_line
        for ref_id in ref_ids:
            models[ref_id]["qualifiers"]["annotations"] = assemblies

        out[assembly] = models

    return out


def annotations_summary(models_directory, ref_id_start=None):
    """
    Print information about how many genes and transcripts are present
    in the models, including how many transcripts were added
    from older releases.

    :param models_directory: Directory with the reference model files.
    :param ref_id_start: Limit to specific reference(s) ID.
    """

    def _per_model():
        output = {}
        for file in Path(models_directory).glob(glob):
            model = json.load(open(file))
            summary = {"genes": 0, "transcripts": 0, "added": 0}
            if model.get("features"):
                summary["genes"] += len(model["features"])
                for gene in model["features"]:
                    if gene.get("features"):
                        summary["transcripts"] += len(gene)
                        for transcript in gene["features"]:
                            if transcript.get("qualifiers") and transcript["qualifiers"].get("annotation_added_from"):
                                summary["added"] += 1
            output[model["id"]] = summary
        total_genes = sum([output[ref_id]["genes"] for ref_id in output])
        total_transcripts = sum([output[ref_id]["transcripts"] for ref_id in output])
        total_added = sum([output[ref_id]["added"] for ref_id in output])

        header = f"{'Reference ID':15} {'genes':>10}{'transcripts':>15}{'added':>10}"
        print(f"\n{header}\n{'-' * len(header)}")
        for ref_id in sorted(output):
            genes = f"{output[ref_id]['genes']:>10}"
            transcripts = f"{output[ref_id]['transcripts']:>15}"
            added = f"{output[ref_id]['added']:>10}"
            print(f"{ref_id:15} {genes}{transcripts}{added}")
        total = (
            f"{'Total':15} {total_genes:>10}{total_transcripts:>15}{total_added:>10}"
        )
        print(f"{'-' * len(header)}\n{total}\n")

    glob = "*"
    if ref_id_start is not None:
        glob = f"{ref_id_start}{glob}"

    _per_model()


def retrieve_assemblies(
    directory_input="./ncbi_annotation_releases",
    directory_output="./ncbi_annotation_models",
    assembly_id_start=None,
    ref_id_start=None,
    downloaded=False,
    include_sequence=False,
):
    if not downloaded:
        download_annotation_releases(_annotations_urls(), directory_input, assembly_id_start)
    else:
        print(f"Using downloaded releases from:\n {directory_input}")
    models = get_annotation_models(
        directory_input=directory_input,
        assembly_id_start=assembly_id_start,
        ref_id_start=ref_id_start,
    )

    Path(directory_output).mkdir(parents=True, exist_ok=True)
    for assembly_id in models:
        for r_id in models[assembly_id]:
            file_path = f"{directory_output}/{r_id}"
            print(f"- writing {file_path}.annotations")
            open(f"{file_path}.annotations", "w").write(json.dumps(models[assembly_id][r_id]))

    if include_sequence:
        print("Downloading the sequences:")
        for assembly_id in models:
            for r_id in models[assembly_id]:
                file_path = f"{directory_output}/{r_id}"
                print(f"- get the sequence for {r_id}")
                seq = parse(fetch(r_id, "fasta")[0], "fasta")["seq"]
                print(f"- writing {file_path}.sequence")
                open(f"{file_path}.sequence", "w").write(seq)

    print("\n")
