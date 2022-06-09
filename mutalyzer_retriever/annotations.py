import gzip
import io
import json
import xml.etree.ElementTree as ET
from copy import deepcopy
from ftplib import FTP, error_perm
from pathlib import Path

import requests

from mutalyzer_retriever.parser import parse
from mutalyzer_retriever.retriever import retrieve_raw


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
                    if transcript["id"] in transcripts:
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


def _merge(new, old):
    ts_new = _get_transcripts_mappings(new)
    ts_old = _get_transcripts_mappings(old)

    ts_not_in = set(ts_old.keys()) - set(ts_new.keys())
    # print("- before")
    # print(len(ts_new), len(ts_old))
    # print(len(ts_not_in))

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
            new["features"].append(gene_old)
            for t in set(gene_ts) - set(gene_ts_already_in):
                ts_new[t] = {"i_g": len(new["features"]), "gene_id": gene_old["id"]}
        else:
            # print(ts_old[t_not_in_id])
            transcript = old["features"][ts_old[t_not_in_id]["i_g"]]["features"][
                ts_old[t_not_in_id]["i_t"]
            ]
            # print(" - add:", transcript)
            if gene_new.get("features") is None:
                gene_new["features"] = []
            gene_new["features"].append(deepcopy(transcript))
            ts_new[t_not_in_id] = {
                "i_g": _get_gene_i(ts_old[t_not_in_id]["gene_id"], new),
                "gene_id": gene_new["id"],
            }

    ts_not_in = set(ts_old.keys()) - set(ts_new.keys())
    ts_new = _get_transcripts_mappings(new)
    ts_old = _get_transcripts_mappings(old)
    # print("- after")
    # print(len(ts_new), len(ts_old))
    # print(len(ts_not_in))
    if len(ts_not_in) != 0:
        raise Exception("Not all the transcripts were added.")


def group_by_accession(annotations):
    groups = {}
    sorted_annotations = sorted(annotations, key=lambda d: d["freeze_date_id"])
    for annotation in sorted_annotations:
        assembly = annotation["assembly_name"].split(".")[0]
        if assembly not in groups:
            groups[assembly] = []
        groups[assembly].append(annotation)
    return groups


class Annotations:
    def __init__(
        self,
        path_input="./downloads",
        path_output="./models",
        downloaded=False,
        ref_id=None,
    ):
        self.ftp_url = "ftp.ncbi.nlm.nih.gov"
        self.ftp_dir = (
            "genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases"
        )

        self.local_input_dir = path_input
        self.local_output_dir = path_output

        self.ref_id_start = ref_id

        if downloaded:
            self.annotations = json.loads(open(self._metadata_path(), "r").read())
        else:
            self._raw_start()

        self._get_models()

    def _raw_start(self):
        self.annotations = self._get_ftp_locations()
        self._input_directory_setup()
        self._retrieve_files()
        self._update_dates()
        open(self._metadata_path(), "w").write(json.dumps(self.annotations, indent=2))

    def _get_ftp_locations(self):
        print("- get ftp locations")
        locations = []
        with FTP(self.ftp_url) as ftp:
            ftp.login()
            ftp.cwd(self.ftp_dir)
            for d_a in ftp.nlst():
                try:
                    ftp.cwd(d_a)
                except error_perm:
                    continue
                annotation = {"id": d_a}
                for d_d in ftp.nlst():
                    if d_d.endswith("annotation_report.xml"):
                        annotation["annotation_report"] = d_d
                    if d_d.startswith("GCF_") and "GRCh" in d_d:
                        annotation["dir"] = d_d
                        try:
                            ftp.cwd(d_d)
                        except error_perm:
                            continue
                        for d_f in ftp.nlst():
                            if d_f.endswith("_genomic.gff.gz"):
                                annotation["file"] = d_f
                        ftp.cwd("..")
                ftp.cwd("..")
                locations.append(annotation)
        print("  done")
        return locations

    def _input_directory_setup(self):
        print(f"- local input directory set up to {self.local_input_dir}")
        local_dir_path = Path(self.local_input_dir)

        if not local_dir_path.is_dir():
            print("  created")
            local_dir_path.mkdir()
        print("  done")

    def _output_directory_setup(self):
        print(f"- local output directory set up to {self.local_output_dir}")
        local_dir_path = Path(self.local_output_dir)

        if not local_dir_path.is_dir():
            print("  created")
            local_dir_path.mkdir()
        print("  done")

    def _retrieve_files(self):
        print("- retrieve files")

        common_url = "https://" + self.ftp_url + "/" + self.ftp_dir
        for annotation in self.annotations:
            url = f"{common_url}/{annotation['id']}/{annotation['dir']}/{annotation['file']}"
            # print(url)
            r = requests.get(url)
            open(self._gff_file_name(annotation), "wb").write(r.content)

            url = f"{common_url}/{annotation['id']}/{annotation['annotation_report']}"
            # print(url)
            r = requests.get(url)
            open(self._report_file_name(annotation), "wb").write(r.content)
        print("  done")

    def _gff_file_name(self, location):
        return self.local_input_dir + "/" + location["id"] + "_" + location["file"]

    def _report_file_name(self, annotation):
        return (
            self.local_input_dir
            + "/"
            + annotation["id"]
            + "_"
            + annotation["annotation_report"]
        )

    def _metadata_path(self):
        return self.local_input_dir + "/" + "metadata.json"

    def _update_dates(self):
        print("- update dates")
        for annotation in self.annotations:
            annotation.update(self._report_info(self._report_file_name(annotation)))
        print("  done")

    @staticmethod
    def _report_info(report_file):
        tree = ET.parse(report_file)
        root = tree.getroot()
        return {
            "freeze_date_id": root.find("./BuildInfo/FreezeDateId").text,
            "assembly_name": root.find("./AssembliesReport/FullAssembly/Name").text,
            "assembly_accession": root.find(
                "./AssembliesReport/FullAssembly/Accession"
            ).text,
        }

    def _get_models(self):
        self._output_directory_setup()

        assemblies = group_by_accession(self.annotations)
        for assembly in assemblies:
            self.get_assembly_model(assemblies[assembly])

    def get_assembly_model(self, annotations):
        out = {}
        for annotation in annotations:
            print(
                f"- processing {annotation['id']} from {annotation['freeze_date_id']}, ({annotation['assembly_name']}, {annotation['assembly_accession']})"
            )

            with gzip.open(self._gff_file_name(annotation), "rb") as f:
                current_id = ""
                current_content = ""
                extras = ""
                for line in f:
                    s_line = line.decode()
                    if s_line.startswith("#!"):
                        extras += s_line
                    elif s_line.startswith("##sequence-region"):
                        if current_id and (
                            self.ref_id_start is None
                            or current_id.startswith(self.ref_id_start)
                        ):
                            current_model = parse(current_content, "gff3")
                            print(f"  - {current_id}")
                            # print(current_model["qualifiers"])
                            if current_id not in out:
                                out[current_id] = current_model
                            else:
                                _merge(current_model, out[current_id])
                                out[current_id] = current_model

                        current_id = s_line.split(" ")[1]
                        current_content = f"##gff-version 3\n{extras}{s_line}"
                    elif s_line.startswith("##species") or s_line.startswith(
                        current_id
                    ):
                        current_content += s_line

        for r_id in out:
            print(f"- writing {r_id}")
            fasta = retrieve_raw(r_id, "ncbi", "fasta", timeout=10)
            model = {"annotations": out[r_id], "sequence": parse(fasta[0], "fasta")}
            open(self.local_output_dir + r_id + ".json", "w").write(json.dumps(model))
        print("\n")


def retrieve_annotations(path_input, path_output, downloaded, ref_id_start=None):
    print(path_input, path_output, downloaded, ref_id_start)

    Annotations(path_input, path_output, downloaded, ref_id_start)
