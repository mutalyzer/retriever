import gzip
import json
import xml.etree.ElementTree as ET
from copy import deepcopy
from io import BytesIO
from pathlib import Path

import requests

from ..parser import parse
from .ncbi import fetch


def _get_gene(g_id, model):
    features = model.get("features", [])
    return next((gene for gene in features if gene.get("id") == g_id), None)


def _get_gene_i(g_id, model):
    features = model.get("features", [])
    return next((i for i, gene in enumerate(features) if gene.get("id") == g_id), None)


def _get_gene_transcript_ids(gene):
    return [feature.get("id") for feature in gene.get("features", []) if "id" in feature]


def _get_transcripts_mappings(model):
    transcripts = {}

    for i_g, gene in enumerate(model.get("features", [])):
        gene_id = gene.get("id")
        if not gene_id:
            continue  # Skip genes without an ID

        for i_t, transcript in enumerate(gene.get("features", [])):
            transcript_id = transcript.get("id")
            if transcript_id is None:
                continue  # Skip transcripts without an ID

            if transcript_id in transcripts:
                raise ValueError(f"Duplicate transcript ID detected: {transcript_id}")

            transcripts[transcript_id] = {
                "i_g": i_g,
                "gene_id": gene_id,
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
    """
    Include into the new reference model the missing genes and transcripts found only in the old
    model.

    :param new: Reference model with a newer freeze date.
    :param old: Reference model with an older freeze date.
    """
    ts_new = _get_transcripts_mappings(new)
    ts_old = _get_transcripts_mappings(old)

    ts_not_in = sorted(set(ts_old) - set(ts_new))
    gene_cache = {}

    for t_not_in_id in ts_not_in:
        if t_not_in_id in ts_new:
            continue
        old_gene_id = ts_old[t_not_in_id]["gene_id"]
        if old_gene_id not in gene_cache:
            gene_cache[old_gene_id] = _get_gene(old_gene_id, new)
        gene_new = gene_cache[old_gene_id]
        if not gene_new:
            gene_old = deepcopy(_get_gene(old_gene_id, old))
            gene_ts = set(_get_gene_transcript_ids(gene_old))
            gene_old["features"] = [t for t in gene_old["features"] if t["id"] not in ts_new]
            _gene_added_from(gene_old, old)
            new.setdefault("features", []).append(gene_old)
            gene_index = len(new["features"]) - 1  # Store the new gene index
            for t in gene_ts - ts_new.keys():
                ts_new[t] = {"i_g": gene_index, "gene_id": gene_old["id"]}
        else:
            transcript = deepcopy(old["features"][ts_old[t_not_in_id]["i_g"]]["features"][ts_old[t_not_in_id]["i_t"]])
            _added_from(transcript, old)
            gene_new.setdefault("features", []).append(transcript)
            if old_gene_id not in gene_cache:
                gene_cache[old_gene_id] = _get_gene_i(old_gene_id, new)
            ts_new[t_not_in_id] = {
                "i_g": gene_cache[old_gene_id],
                "gene_id": gene_new["id"],
            }

    # Final check for any remaining missing transcripts
    if set(ts_old) - set(ts_new):
        raise ValueError(f"Not all transcripts were added: {set(ts_old) - set(ts_new)}")


def _common_url():
    return "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/"


def _annotations_urls():
    common = _common_url()
    annotations_urls = {
        "GRCh37": [
            [
                common
                + "105.20190906/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gff.gz",
                common + "105.20190906/Homo_sapiens_AR105_annotation_report.xml",
            ],
            [
                common
                + "105.20220307/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.gff.gz",
                common
                + "105.20220307/Homo_sapiens_AR105.20220307_annotation_report.xml",
            ],
            [
                common
                + "GCF_000001405.25-RS_2024_09/GCF_000001405.25_GRCh37.p13_genomic.gff.gz",
                common
                + "GCF_000001405.25-RS_2024_09/GCF_000001405.25-RS_2024_09_annotation_report.xml",
            ],
        ],
        "GRCh38": [
            [
                common
                + "109/GCF_000001405.38_GRCh38.p12/GCF_000001405.38_GRCh38.p12_genomic.gff.gz",
                common
                + "109/GCF_000001405.38_GRCh38.p12/Homo_sapiens_AR109_annotation_report.xml",
            ],
            [
                common
                + "109.20190607/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz",
                common
                + "109.20190607/GCF_000001405.39_GRCh38.p13/Homo_sapiens_AR109.20190607_annotation_report.xml",
            ],
            [
                common
                + "109.20190905/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz",
                common
                + "109.20190905/GCF_000001405.39_GRCh38.p13/Homo_sapiens_AR109.20190905_annotation_report.xml",
            ],
            [
                common
                + "109.20191205/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz",
                common
                + "109.20191205/GCF_000001405.39_GRCh38.p13/Homo_sapiens_AR109.20191205_annotation_report.xml",
            ],
            [
                common
                + "109.20200228/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz",
                common
                + "109.20200228/GCF_000001405.39_GRCh38.p13/Homo_sapiens_AR109.20200228_annotation_report.xml",
            ],
            [
                common
                + "109.20200522/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz",
                common
                + "109.20200522/GCF_000001405.39_GRCh38.p13/Homo_sapiens_AR109.20200522_annotation_report.xml",
            ],
            [
                common
                + "109.20200815/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz",
                common
                + "109.20200815/GCF_000001405.39_GRCh38.p13/Homo_sapiens_AR109.20200815_annotation_report.xml",
            ],
            [
                common
                + "109.20201120/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz",
                common
                + "109.20201120/GCF_000001405.39_GRCh38.p13/Homo_sapiens_AR109.20201120_annotation_report.xml",
            ],
            [
                common
                + "109.20210226/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz",
                common
                + "109.20210226/GCF_000001405.39_GRCh38.p13/Homo_sapiens_AR109.20210226_annotation_report.xml",
            ],
            [
                common
                + "109.20210514/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz",
                common
                + "109.20210514/GCF_000001405.39_GRCh38.p13/Homo_sapiens_AR109.20210514_annotation_report.xml",
            ],
            [
                common
                + "109.20211119/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz",
                common
                + "109.20211119/GCF_000001405.39_GRCh38.p13/Homo_sapiens_AR109.20211119_annotation_report.xml",
            ],
            [
                common
                + "110/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz",
                common + "110/Homo_sapiens_AR110_annotation_report.xml",
            ],
            [
                common
                + "GCF_000001405.40-RS_2023_03/GCF_000001405.40_GRCh38.p14_genomic.gff.gz",
                common
                + "GCF_000001405.40-RS_2023_03/GCF_000001405.40-RS_2023_03_annotation_report.xml",
            ],
            [
                common
                + "GCF_000001405.40-RS_2023_10/GCF_000001405.40_GRCh38.p14_genomic.gff.gz",
                common
                + "GCF_000001405.40-RS_2023_10/GCF_000001405.40-RS_2023_10_annotation_report.xml",
            ],
            [
                common
                + "GCF_000001405.40-RS_2024_08/GCF_000001405.40_GRCh38.p14_genomic.gff.gz",
                common
                + "GCF_000001405.40-RS_2024_08/GCF_000001405.40-RS_2024_08_annotation_report.xml",
            ],
            [
                common
                + "GCF_000001405.40-RS_2025_08/GCF_000001405.40_GRCh38.p14_genomic.gff.gz",
                common
                + "GCF_000001405.40-RS_2025_08/GCF_000001405.40-RS_2025_08_annotation_report.xml",
            ],
        ],
        "T2T-CHM13v2": [
            [
                common
                + "GCF_009914755.1-RS_2023_03/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz",
                common
                + "GCF_009914755.1-RS_2023_03/GCF_009914755.1-RS_2023_03_annotation_report.xml",
            ],
            [
                common
                + "GCF_009914755.1-RS_2023_10/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz",
                common
                + "GCF_009914755.1-RS_2023_10/GCF_009914755.1-RS_2023_10_annotation_report.xml",
            ],
            [
                common
                + "GCF_009914755.1-RS_2024_08/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz",
                common
                + "GCF_009914755.1-RS_2024_08/GCF_009914755.1-RS_2024_08_annotation_report.xml",
            ],
            [
                common
                + "GCF_009914755.1-RS_2025_08/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz",
                common
                + "GCF_009914755.1-RS_2025_08/GCF_009914755.1-RS_2025_08_annotation_report.xml",
            ],
        ],
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


def download_ncbi_releases_file(url):
    response = requests.get(url, timeout=10)
    if response.status_code == 200:
        return response.content
    raise response.raise_for_status()


def _file_name(url):
    # return url.replace(_common_url(), "").replace("/", "__")
    return url.split("/")[-1]


def download_annotation_releases(urls, assembly_id_start=None):
    """
    Download the annotations (GFF3) report (XML) files in the provided directory.
    """
    print("\nDownloading assembly releases:")
    metadata = {}
    for assembly_id in urls:
        if assembly_id_start is not None and not assembly_id.startswith(
            assembly_id_start
        ):
            continue
        print(f" - assembly: {assembly_id}")
        metadata[assembly_id] = {}
        for i, (url_gff3, url_xml) in enumerate(urls[assembly_id]):
            print(f"  - {i+1}/{len(urls[assembly_id])}:")
            print(f"   - downloading xml: {url_xml}")
            response_xml = download_ncbi_releases_file(url_xml)
            freeze_date_id = _report_info(response_xml)["freeze_date_id"]
            print(f"   - freeze date: {freeze_date_id}")
            print(f"   - downloading gff3: {url_gff3}")
            response_gff3 = download_ncbi_releases_file(url_gff3)
            metadata[assembly_id][freeze_date_id] = {
                "file_name_gff3": _file_name(url_gff3),
                "file_name_xml": _file_name(url_xml),
                "xml": response_xml,
                "gff3": response_gff3,
            }
    return metadata


def write_annotations_releases(metadata, directory="./ncbi_annotation_releases"):
    print("\nWriting the files:")
    path_dir = Path(directory)
    path_dir.mkdir(parents=True, exist_ok=True)
    for assembly_id, freeze_dates in metadata.items():
        for freeze_date_id, file_data in freeze_dates.items():
            file_name_xml = f"{assembly_id}__{freeze_date_id}__{metadata[assembly_id][freeze_date_id]['file_name_xml']}"
            file_name_gff3 = f"{assembly_id}__{freeze_date_id}__{metadata[assembly_id][freeze_date_id]['file_name_gff3']}"
            path_xml = Path(directory) / file_name_xml
            path_gff3 = Path(directory) / file_name_gff3
            print(f" - writing xml: {path_xml}")
            with open(path_xml, "wb") as xml_file:
                xml_file.write(file_data["xml"])
            print(f" - writing gff3: {path_gff3}")
            with open(path_gff3, "wb") as gff3_file:
                gff3_file.write(file_data["gff3"])


def get_annotation_models(metadata, ref_id_start=None):
    print("\nParse and extract the annotation models:")
    out = {}

    for assembly in metadata:
        print(f"- assembly: {assembly}")
        models = {}
        assemblies = []
        processed_ref_ids = set()

        for freeze_date_id in sorted(metadata[assembly]):
            print(f"  - freeze date: {freeze_date_id}")
            assembly_details = _report_info(metadata[assembly][freeze_date_id]["xml"])
            assemblies.append(assembly_details)

            gff3_data = metadata[assembly][freeze_date_id]["gff3"]
            _process_gff3_file(gff3_data, models, assembly_details, ref_id_start, processed_ref_ids)

        # Update all models with complete assembly list
        for ref_id in processed_ref_ids:
            models[ref_id]["qualifiers"]["annotations"] = assemblies

        out[assembly] = models

    return out


def _process_gff3_file(gff3_data, models, assembly_details, ref_id_start, processed_ref_ids):
    """Process a GFF3 file and update models dictionary."""
    with gzip.open(BytesIO(gff3_data), 'rt') as f:
        current_id = ""
        current_content = ""
        header_lines = ""

        for line in f:
            if line.startswith("#!"):
                header_lines += line
            elif line.startswith("##sequence-region"):
                _finalize_current_model(
                    current_id, current_content, models, assembly_details, ref_id_start, processed_ref_ids
                )

                current_id = line.split()[1]
                current_content = f"##gff-version 3\n{header_lines}{line}"
            elif line.startswith("##species") or line.startswith(current_id):
                current_content += line

        # Process final entry
        _finalize_current_model(
            current_id, current_content, models, assembly_details, ref_id_start, processed_ref_ids
        )


def _finalize_current_model(ref_id, content, models, assembly_details, ref_id_start, processed_ref_ids):
    """Finalize and store the current model if it should be processed."""
    if _should_process(ref_id, ref_id_start):
        return

    print(f"   - parsing reference id: {ref_id}")
    current_model = parse(content, "gff3")
    current_model["qualifiers"]["annotations"] = assembly_details

    if ref_id in models:
        print(f"   - merging: {ref_id}")
        _merge(current_model, models[ref_id])

    models[ref_id] = current_model
    processed_ref_ids.add(ref_id)


def _should_process(ref_id, ref_id_start):
    """Check if a reference ID should be processed based on the filter."""
    return not ref_id or not (ref_id_start is None or ref_id.startswith(ref_id_start))


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
            with open(file, encoding="utf-8") as json_file:
                model = json.load(json_file)
            summary = {"genes": 0, "transcripts": 0, "added": 0}
            if model.get("features"):
                summary["genes"] += len(model["features"])
                for gene in model["features"]:
                    if gene.get("features"):
                        summary["transcripts"] += len(gene)
                        for transcript in gene["features"]:
                            if transcript.get("qualifiers") and transcript[
                                "qualifiers"
                            ].get("annotation_added_from"):
                                summary["added"] += 1
            output[model["id"]] = summary
        total_genes = sum(model["genes"] for ref_id, model in output.items())
        total_transcripts = sum(model["transcripts"] for ref_id, model in output.items())
        total_added = sum(model["added"] for ref_id, model in output.items())

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


def load_metadata(directory="./ncbi_annotation_releases", assembly_id_start=None):
    print("\nLoading annotation releases:")
    metadata = {}
    for path_file in Path(directory).iterdir():
        file_name = str(path_file).rsplit("/", maxsplit=1)[-1]
        assembly, freeze_date, file_path = str(file_name).split("__")
        if assembly_id_start is None or (
            assembly_id_start and assembly.startswith(assembly_id_start)
        ):
            print(f" - assembly: {assembly}")
            if assembly not in metadata:
                metadata[assembly] = {}
            if freeze_date not in metadata[assembly]:
                metadata[assembly][freeze_date] = {}
            print(f"  - date: {freeze_date}")
            if file_path.endswith("gz"):
                print(f"   - load gff: {file_name}")
                metadata[assembly][freeze_date]["file_name_gff3"] = file_name.split("__")[-1]
                with open(path_file, "rb") as gff_file:
                    metadata[assembly][freeze_date]["gff3"] = gff_file.read()
            if file_path.endswith("xml"):
                print(f"   - load xml: {file_name}")
                metadata[assembly][freeze_date]["file_name_xml"] = file_name.split("__")[-1]
                with open(path_file, "rb") as xml_file:
                    metadata[assembly][freeze_date]["xml"] = xml_file.read()
    return metadata


def retrieve_assemblies(
    directory_input="./ncbi_annotation_releases",
    directory_output="./ncbi_annotation_models",
    assembly_id_start=None,
    ref_id_start=None,
    downloaded=False,
    write_downloaded=False,
    include_sequence=False,
):
    """
    Retrieve the models (including historical transcripts) from the NCBI FTP locations.

    :param directory_input: Where to write the release files.
    :param directory_output: Where to write the annotation (and sequence) files.
    :param assembly_id_start: Restrict only to certain assemblies.
    :param ref_id_start: Restrict only to certain reference ids.
    :param downloaded: The release files are already downloaded in the directory_input.
    :param write_downloaded: Write the release files locally in the in directory_input.
    :param include_sequence: Download also the sequence and write it in the directory_output.
    """
    if downloaded:
        metadata = load_metadata(directory_input, assembly_id_start)
    else:
        metadata = download_annotation_releases(_annotations_urls(), assembly_id_start)
        if write_downloaded:
            write_annotations_releases(metadata, directory_input)

    models = get_annotation_models(metadata, ref_id_start=ref_id_start)

    print("\nWriting the models:")
    Path(directory_output).mkdir(parents=True, exist_ok=True)
    for assembly_id, references in models.items():
        print(f" - assembly: {assembly_id}")
        for r_id, r_model in references.items():
            file_path = f"{directory_output}/{r_id}"
            print(f"  - writing: {file_path}.annotations")
            with open(f"{file_path}.annotations", "w", encoding="utf-8") as annotations_file:
                annotations_file.write(json.dumps(r_model))

    if include_sequence:
        print("\nDownloading and writing the sequences:")
        for assembly_id, references in models.items():
            for r_id in references:
                file_path = f"{directory_output}/{r_id}"
                print(f" - download the sequence for: {r_id}")
                seq = parse(fetch(r_id, "fasta")[0], "fasta")["seq"]
                print(f"  - writing: {file_path}.sequence")
                with open(f"{file_path}.sequence", "w", encoding="utf-8") as sequence_file:
                    sequence_file.write(seq)
