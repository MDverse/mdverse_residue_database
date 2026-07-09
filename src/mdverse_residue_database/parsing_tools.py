import json
import re
from pathlib import Path

import click
from loguru import logger
from pydantic import BaseModel


class ResidueEntry(BaseModel):
    name: str
    atoms: list[str]
    type: str
    force_fields: str
    source_files: list[str]


RTP_FILENAME_TO_TYPE = {
    "aminoacids": "PROTEIN",
    "dna": "NUCLEIC_ACID",
    "rna": "NUCLEIC_ACID",
    "ions": "ION",
    "tip3p": "WATER",
    "tip4pew": "WATER",
    "tip4p": "WATER",
    "tip5p": "WATER",
    "spc": "WATER",
    "spce": "WATER",
}


def parse_rtp(
    rtp_filepath: Path, force_field: str, molecule_type: str
) -> list[ResidueEntry]:
    """
    Parse a single GROMACS .rtp file and extract all residues as ResidueEntry objects.

    Returns
    -------
    list[ResidueEntry]
        list of ResidueEntry objects extracted from the file.
    """
    residues: list[ResidueEntry] = []
    current_name: str | None = None
    current_atoms: list[str] = []
    in_atoms = False

    with open(rtp_filepath, "r") as f:
        for line in f:
            line = line.strip()

            if line and not line.startswith(";"):
                section_match = re.match(r"^\[\s*(\w+)\s*\]$", line)

                if section_match:
                    logger.success(f"Match found in line {line}")
                    section_name = section_match.group(1).strip()
                    logger.info(f"Name found : {section_name}")

                    if section_name == "atoms":
                        in_atoms = True

                    elif section_name.isupper():
                        logger.info(f"Residue Name : {section_name}")
                        if current_name and current_atoms:
                            logger.info(
                                f"Creating complete residue entry for {current_name}"
                            )
                            entry = ResidueEntry(
                                name=current_name,
                                atoms=current_atoms.copy(),
                                type=molecule_type,
                                force_fields=force_field,
                                source_files=[rtp_filepath.name],
                            )
                            residues.append(entry)

                        current_name = section_name
                        current_atoms = []
                        in_atoms = False

                    else:
                        in_atoms = False

                elif in_atoms and current_name:
                    logger.info(f"Adding atoms composing {current_atoms}")
                    atom_parts = line.split()
                    if atom_parts:
                        current_atoms.append(atom_parts[0])
        logger.info(f"Added {len(current_atoms)} atoms for {current_name}")

    if current_name and current_atoms:
        logger.info(f"Creating an entry for the last residue: {current_name}")
        entry = ResidueEntry(
            name=current_name,
            atoms=current_atoms,
            type=molecule_type,
            force_fields=force_field,
            source_files=[rtp_filepath.name],
        )
        residues.append(entry)

    if not residues:
        logger.warning(f"No residues extracted from {rtp_filepath.name}")
    else:
        logger.debug(f"{len(residues)} residues extracted from {rtp_filepath.name}")

    return residues


def parse_itp(
    itp_filepath: Path, force_field: str, molecule_type: str
) -> list[ResidueEntry]:
    """
    Parse a single GROMACS .itp file and extract all residues as ResidueEntry objects.

    Returns
    -------
    list[ResidueEntry]
        list of ResidueEntry objects extracted from the file.
    """
    residues: list[ResidueEntry] = []
    current_name: str | None = None
    current_atoms: list[str] = []
    in_moleculetype = False
    in_atoms = False

    with open(itp_filepath, "r") as f:
        for line in f:
            line = line.strip()

            if line and not line.startswith(";"):
                section_match = re.match(r"^\[\s*(\w+)\s*\]$", line)

                if section_match:
                    logger.success(f"Match found in line {line}")
                    section_name = section_match.group(1).strip()
                    logger.info(f"Name found : {section_name}")

                    if section_name == "moleculetype":
                        in_moleculetype = True
                        in_atoms = False

                    elif section_name == "atoms":
                        in_atoms = True
                        in_moleculetype = False

                    else:
                        in_moleculetype = False
                        in_atoms = False

                elif in_moleculetype:
                    molecule_parts = line.split()
                    logger.info(f"Information of the molecules : {molecule_parts}")
                    if molecule_parts:
                        if current_name and current_atoms:
                            logger.info(
                                f"Creating complete residue entry for {current_name}"
                            )
                            entry = ResidueEntry(
                                name=current_name,
                                atoms=current_atoms.copy(),
                                type=molecule_type,
                                force_fields=force_field,
                                source_files=[itp_filepath.name],
                            )
                            residues.append(entry)

                        current_name = molecule_parts[0]
                        current_atoms = []
                    in_moleculetype = False

                elif in_atoms and current_name:
                    logger.info(f"Adding atoms composing {current_atoms}")
                    atom_parts = line.split()
                    if len(atom_parts) >= 5:
                        current_atoms.append(atom_parts[4])

        logger.info(f"Added {len(current_atoms)} atoms for {current_name}")

    if current_name and current_atoms:
        logger.info(f"Creating an entry for the last residue: {current_name}")
        entry = ResidueEntry(
            name=current_name,
            atoms=current_atoms,
            type=molecule_type,
            force_fields=force_field,
            source_files=[itp_filepath.name],
        )
        residues.append(entry)

    if not residues:
        logger.warning(f"No residues extracted from {itp_filepath.name}")
    else:
        logger.debug(f"{len(residues)} residues extracted from {itp_filepath.name}")

    return residues


def apply_tdb_patches_to_residue(
    residue: ResidueEntry,
    tdb_filepath: Path,
) -> list[ResidueEntry]:
    """
    Read a .tdb file line by line an applying modifications.
    Returns
    -------
    list[ResidueEntry]
        One ResidueEntry variant per completed patch.
    """
    modified_residues: list[ResidueEntry] = []

    atoms = residue.atoms.copy()
    mode: str | None = None

    with open(tdb_filepath, "r") as f:
        lines = [line.strip() for line in f.readlines()]

    i = 0
    n = len(lines)
    while i < n:
        line = lines[i]
        i += 1

        if line and not line.startswith(";"):
            if line == "[ replace ]":
                mode = "replace"

            elif line == "[ add ]":
                mode = "add"

            elif line == "[ delete ]":
                mode = "delete"

            elif line.startswith("[") and line.endswith("]"):
                if mode == "delete":
                    logger.info(f"Patch complete, saving variant for {residue.name}")
                    modified_residues.append(
                        residue.model_copy(
                            update={
                                "atoms": atoms,
                                "source_files": residue.source_files
                                + [tdb_filepath.name],
                            }
                        )
                    )
                    atoms = residue.atoms.copy()
                mode = None

            elif mode == "replace":
                old_name, new_name = line.split()[:2]
                if old_name in atoms:
                    atoms[atoms.index(old_name)] = new_name
                    logger.info(
                        f"Replaced {old_name} with {new_name} in {residue.name}"
                    )
                else:
                    logger.warning(f"{old_name} to replace not found in {residue.name}")

            elif mode == "add":
                nadd = int(line.split()[0])
                for _ in range(nadd):
                    atom_line = lines[i]
                    i += 1
                    new_atom_name = atom_line.split()[0]
                    atoms.append(new_atom_name)
                    logger.info(f"Added {new_atom_name} to {residue.name}")

            elif mode == "delete":
                for atom_name in line.split():
                    if atom_name in atoms:
                        atoms.remove(atom_name)
                        logger.info(f"Deleted {atom_name} from {residue.name}")
                    else:
                        logger.warning(
                            f"{atom_name} to delete not found in {residue.name}"
                        )

    if mode == "delete":
        logger.info(f"Patch complete (end of file), saving variant for {residue.name}")
        modified_residues.append(
            residue.model_copy(
                update={
                    "atoms": atoms,
                    "source_files": residue.source_files + [tdb_filepath.name],
                }
            )
        )

    return modified_residues


def apply_terminaison_modification_all_res(
    list_entry_residue: list[ResidueEntry], tdb_filepath: Path
) -> dict[str, list[ResidueEntry]]:
    """
    Apply terminus patches from a .tdb file to a list of residues and
    group the resulting variants by residue name.

    Returns
    -------
    dict[str, list[ResidueEntry]]
        Mapping of residue name to the list of its ResidueEntry variants
        (one per patch applied).
    """
    residue_database: dict[str, list[ResidueEntry]] = {}

    for residue in list_entry_residue:
        variant_list = apply_tdb_patches_to_residue(residue, tdb_filepath)

        # [WIP: Changes are needed so that even if a residue is already present withe teh same name and atoms but different force field, it should be added to the database]
        if residue.name not in residue_database:
            residue_database[residue.name] = []
        residue_database[residue.name].extend(variant_list)

        logger.info(f"{len(variant_list)} variants added for {residue.name}")

    return residue_database


@click.command()
@click.option(
    "--input-filepath",
    type=click.Path(exists=True, path_type=Path),
    help="rtp or itp file conatinaing the residue or wter molecules",
)
@click.option(
    "--force-field",
    required=True,
    help="Name of the force field from where the file was extracted.",
)
@click.option(
    "--tdb-filepath",
    type=click.Path(exists=True, path_type=Path),
    help="Optional.tdb file to apply terminaison patches",
)
@click.option(
    "--residue-db-json",
    type=click.Path(path_type=Path),
    help="Ouput JSON filepath that will contain the different residue entry and their informations",
)
def main(
    input_filepath: Path,
    force_field: str,
    tdb_filepath: Path | None,
    residue_db_json: Path | None,
):
    file_tye = input_filepath.suffix.lower()

    molecule_type = RTP_FILENAME_TO_TYPE.get(input_filepath.stem.lower(), "RESIDUE")
    logger.info(f"Molecule type : {molecule_type}")

    if file_tye == ".rtp":
        residue_entry_list = parse_rtp(input_filepath, force_field, molecule_type)
    elif file_tye == ".itp":
        residue_entry_list = parse_itp(input_filepath, force_field, molecule_type)

    logger.debug(f"{len(residue_entry_list)} residues parsed")
    ala_count = sum(1 for r in residue_entry_list if r.name == "ALA")
    logger.debug(f"{ala_count} ALA entries found before patching")

    if tdb_filepath:
        residue_database = apply_terminaison_modification_all_res(
            residue_entry_list, tdb_filepath
        )
    else:
        residue_database = {}
        for residue in residue_entry_list:
            residue_database.setdefault(residue.name, []).append(residue)

    residue_db = {
        name: [entry.model_dump() for entry in entries]
        for name, entries in residue_database.items()
    }

    if residue_db_json:
        with open(residue_db_json, "w") as f:
            json.dump(residue_db, f, indent=2)
        logger.success(f"Residue database written to {residue_db_json}")
    else:
        click.echo(json.dumps(residue_db, indent=2))


if __name__ == "__main__":
    main()
