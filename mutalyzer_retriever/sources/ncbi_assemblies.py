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


def _added_from(feature, model):
    if feature.get("qualifiers") is None:
        feature["qualifiers"] = {}
    if feature["qualifiers"].get("added_freeze_date_id") is None:
        feature["qualifiers"]["added_freeze_date_id"] = model["qualifiers"][
            "freeze_date_id"
        ]
    if feature["qualifiers"].get("added_annotation_id") is None:
        feature["qualifiers"]["added_annotation_id"] = model["qualifiers"][
            "annotation_id"
        ]


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


def group_by_accession(annotations):
    groups = {}
    sorted_annotations = sorted(annotations, key=lambda d: d["freeze_date_id"])
    for annotation in sorted_annotations:
        assembly = annotation["assembly_name"].split(".")[0]
        if assembly not in groups:
            groups[assembly] = []
        groups[assembly].append(annotation)
    return groups


class Assemblies:
    """
    Retrieve reference models for human chromosomes.
    """
    def __init__(
        self,
        path_input=None,
        path_output=None,
        downloaded=False,
        ref_id=None,
        split=False,
        only_annotations=False,
    ):
        self.ftp_url = "ftp.ncbi.nlm.nih.gov"
        self.ftp_dir = "genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases"

        self.local_input_dir = path_input if path_input else "./downloads"
        self.local_output_dir = path_output if path_output else "./models"

        self.ref_id_start = ref_id
        self.split = split
        self.only_annotations = only_annotations

        if downloaded:
            self.annotations = json.loads(open(self._metadata_path(), "r").read())
        else:
            self._raw_start()

        self._get_models()

    def _raw_start(self):
        self.annotations = self._get_ftp_locations()
        self._input_directory_setup()
        self._retrieve_gff_files()
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
                                annotation["file_gff"] = d_f
                            elif d_f.endswith("_genomic.fna.gz"):
                                annotation["file_fasta"] = d_f
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

    def _retrieve_gff_files(self):
        print("- retrieve gff files")

        common_url = "https://" + self.ftp_url + "/" + self.ftp_dir
        for annotation in self.annotations:
            url = f"{common_url}/{annotation['id']}/{annotation['dir']}/{annotation['file_gff']}"
            # print(url)
            r = requests.get(url)
            open(self._file_name_gff(annotation), "wb").write(r.content)

            url = f"{common_url}/{annotation['id']}/{annotation['annotation_report']}"
            # print(url)
            r = requests.get(url)
            open(self._file_name_report(annotation), "wb").write(r.content)
        print("  done")

    def _file_name_gff(self, location):
        return self.local_input_dir + "/" + location["id"] + "_" + location["file_gff"]

    def _file_name_fasta(self, location):
        return (
            self.local_input_dir + "/" + location["id"] + "_" + location["file_fasta"]
        )

    def _file_name_report(self, location):
        return (
            self.local_input_dir
            + "/"
            + location["id"]
            + "_"
            + location["annotation_report"]
        )

    def _metadata_path(self):
        return self.local_input_dir + "/" + "metadata.json"

    def _update_dates(self):
        print("- update dates")
        for annotation in self.annotations:
            annotation.update(self._report_info(self._file_name_report(annotation)))
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

    @staticmethod
    def _add_annotations_details(model, annotation):
        if model.get("qualifiers") is None:
            model["qualifiers"] = {}
        model["qualifiers"]["freeze_date_id"] = annotation["freeze_date_id"]
        model["qualifiers"]["annotation_id"] = annotation["id"]
        model["qualifiers"]["assembly_name"] = annotation["assembly_name"]
        model["qualifiers"]["assembly_accession"] = annotation["assembly_accession"]

    def get_assembly_model(self, annotations):
        out = {}
        for annotation in annotations:
            print(
                f"- processing {annotation['id']} from {annotation['freeze_date_id']}, ({annotation['assembly_name']}, {annotation['assembly_accession']})"
            )

            with gzip.open(self._file_name_gff(annotation), "rb") as f:
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
                            open(f"{self.local_output_dir}/{current_id}.gff3", "w").write(current_content)
                            self._add_annotations_details(current_model, annotation)
                            print(f"  - {current_id}")
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
            file_path = f"{self.local_output_dir}/{r_id}"
            if self.only_annotations:
                print(f"- writing {file_path}.annotations")
                open(f"{file_path}.annotations", "w").write(json.dumps(out[r_id]))
            else:
                fasta = retrieve_raw(r_id, "ncbi", "fasta", timeout=10)
                seq = parse(fasta[0], "fasta")
                if self.split:
                    print(f"- writing {file_path}.annotations")
                    open(f"{file_path}.annotations", "w").write(json.dumps(out[r_id]))
                    print(f"- writing {file_path}.sequence")
                    open(f"{file_path}.sequence", "w").write(seq["seq"])
                else:
                    print(f"- writing {self.local_output_dir}/{r_id}")
                    model = {"annotations": out[r_id], "sequence": seq}
                    open(file_path, "w").write(json.dumps(model))
        print("\n")


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
            model = json.load(open(file))["annotations"]
            summary = {"genes": 0, "transcripts": 0, "added": 0}
            if model.get("features"):
                summary["genes"] += len(model["features"])
                for gene in model["features"]:
                    if gene.get("features"):
                        summary["transcripts"] += len(gene)
                        for transcript in gene["features"]:
                            if transcript.get("qualifiers") and transcript[
                                "qualifiers"
                            ].get("added_freeze_date_id"):
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
