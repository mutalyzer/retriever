GRCH38 = {
    "1": "NC_000001.11",
    "2": "NC_000002.12",
    "3": "NC_000003.12",
    "4": "NC_000004.12",
    "5": "NC_000005.10",
    "6": "NC_000006.12",
    "7": "NC_000007.14",
    "8": "NC_000008.11",
    "9": "NC_000009.12",
    "10": "NC_000010.11",
    "11": "NC_000011.10",
    "12": "NC_000012.12",
    "13": "NC_000013.11",
    "14": "NC_000014.9",
    "15": "NC_000015.10",
    "16": "NC_000016.10",
    "17": "NC_000017.11",
    "18": "NC_000018.10",
    "19": "NC_000019.10",
    "20": "NC_000020.11",
    "21": "NC_000021.9",
    "22": "NC_000022.11",
    "23": "NC_000023.11",
    "X": "NC_000023.11",
    "24": "NC_000024.10",
    "Y": "NC_000024.10",
}

GRCH37 = {
    "1": "NC_000001.10",
    "2": "NC_000002.11",
    "3": "NC_000003.11",
    "4": "NC_000004.11",
    "5": "NC_000005.9",
    "6": "NC_000006.11",
    "7": "NC_000007.13",
    "8": "NC_000008.10",
    "9": "NC_000009.11",
    "10": "NC_000010.10",
    "11": "NC_000011.9",
    "12": "NC_000012.11",
    "13": "NC_000013.10",
    "14": "NC_000014.8",
    "15": "NC_000015.9",
    "16": "NC_000016.9",
    "17": "NC_000017.10",
    "18": "NC_000018.9",
    "19": "NC_000019.9",
    "20": "NC_000020.10",
    "21": "NC_000021.8",
    "22": "NC_000022.10",
    "23": "NC_000023.10",
    "X": "NC_000023.10",
    "24": "NC_000024.9",
    "Y": "NC_000024.9",
}

ASSEMBLY_ALIASES = {
    "HG38": "GRCH38",
    "HG19": "GRCH37",
}

ASSEMBLIES = {
    "GRCH38": GRCH38,
    "GRCH37": GRCH37,
}


def get_model_qualifier(model, qualifier):
    if (
        model.get("annotations")
        and model["annotations"].get("qualifiers")
        and model["annotations"]["qualifiers"]
    ):
        return model["annotations"]["qualifiers"].get(qualifier)


def get_reference_mol_type(model):
    return get_model_qualifier(model, "mol_type")


def get_chromosome_accession_from_mrna_model(ref_id, model):
    if get_reference_mol_type(model) == "mRNA" and ref_id.startswith("NM"):
        chromosome_number = get_model_qualifier(model, "chromosome")
        if chromosome_number is not None:
            chromosome_accessions = []
            for assembly in ASSEMBLIES:
                chromosome_accession = ASSEMBLIES[assembly].get(chromosome_number)
                if chromosome_accession is not None:
                    chromosome_accessions.append((assembly, chromosome_accession))
            if chromosome_accessions:
                return chromosome_accessions


def get_assembly_id(assembly_id):
    if assembly_id is not None:
        a_id = assembly_id.upper()
    else:
        return
    if a_id in ASSEMBLY_ALIASES:
        a_id = ASSEMBLY_ALIASES[a_id]
    if a_id in ASSEMBLIES:
        return a_id
    return None


def get_assembly_chromosome_accession(r_id, s_id):
    if r_id is not None:
        assembly = r_id.upper()
    else:
        return
    if assembly in ASSEMBLY_ALIASES:
        assembly = ASSEMBLY_ALIASES[assembly]
    elif assembly not in ASSEMBLIES:
        return None

    if s_id is not None:
        if s_id.upper().startswith("CHR"):
            chr_number = s_id.upper().split("CHR")[1]
        else:
            chr_number = s_id.upper()
    else:
        return None

    if chr_number in ASSEMBLIES[assembly]:
        return ASSEMBLIES[assembly][chr_number]
