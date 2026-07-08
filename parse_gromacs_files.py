import json
from datetime import UTC, datetime
from pathlib import Path

import click
from loguru import logger
from pydantic import BaseModel

from mdverse_residue_database.parsing_tools import (
    RTP_FILENAME_TO_TYPE,
    ResidueEntry,
    apply_terminaison_modification_all_res,
    parse_itp,
    parse_rtp,
)


class ResidueDatabase(BaseModel):
    version: str
    residues: dict[str, list[ResidueEntry]]


def find_force_field_directories(top_dir: Path) -> list[Path]:
    """
    Find all force field directories under the GROMACS
    top directory.

    Returns
    -------
    list[Path]
        list of paths to force field directories.
    """
    ff_directory_list: list[Path] = []

    for folders in sorted(top_dir.iterdir()):
        if folders.is_dir() and folders.suffix == ".ff":
            ff_directory_list.append(folders)
            logger.info(f"Force field directory found : {folders.name}")

    if not ff_directory_list:
        logger.warning(f"No force field directory found in {top_dir}")
    else:
        logger.debug(f"{len(ff_directory_list)} force field directories found")

    return ff_directory_list


def apply_all_terminus_patches(
    residue_entry_list: list[ResidueEntry], n_tdb_filepath: Path, c_tdb_filepath: Path
) -> list[ResidueEntry]:
    """
    Apply N-terminus and C-terminus patches to a list of residues.

    Returns
    -------
    list[ResidueEntry]
        Original residues plus all patched variants found.
    """
    patched_residue_list: list[ResidueEntry] = list(residue_entry_list)

    if n_tdb_filepath.exists():
        logger.info(f"Applying N-terminus patches from {n_tdb_filepath.name}")
        n_patch_database = apply_terminaison_modification_all_res(
            residue_entry_list, n_tdb_filepath
        )
        for residue_name in n_patch_database:
            patched_residue_list.extend(n_patch_database[residue_name])

    if c_tdb_filepath.exists():
        logger.info(f"Applying C-terminus patches from {c_tdb_filepath.name}")
        c_patch_database = apply_terminaison_modification_all_res(
            residue_entry_list, c_tdb_filepath
        )
        for residue_name in c_patch_database:
            patched_residue_list.extend(c_patch_database[residue_name])

    return patched_residue_list


def process_one_force_field_directory(ff_directory: Path) -> list[ResidueEntry]:
    """
    Process all .rtp and .itp files in one force field directory and check terminus patches

    Returns
    -------
    list[ResidueEntry]
        list of ResidueEntry objects, including terminus patch variants.
    """
    force_field = ff_directory.name.replace(".ff", "")
    logger.info(f"Processing force field : {force_field}")

    all_residue_list: list[ResidueEntry] = []

    rtp_filepath_list: list[Path] = []
    itp_filepath_list: list[Path] = []

    for files in sorted(ff_directory.iterdir()):
        if files.is_file() and files.suffix == ".rtp":
            rtp_filepath_list.append(files)
        elif files.is_file() and files.suffix == ".itp":
            itp_filepath_list.append(files)

    for rtp_filepath in rtp_filepath_list:
        molecule_type = RTP_FILENAME_TO_TYPE.get(rtp_filepath.stem.lower(), "Unlisted")
        logger.info(f"Parsing {rtp_filepath.name} as {molecule_type}")

        residue_entry_list = parse_rtp(rtp_filepath, force_field, molecule_type)

        n_tdb_filepath = ff_directory / f"{rtp_filepath.stem}.n.tdb"
        c_tdb_filepath = ff_directory / f"{rtp_filepath.stem}.c.tdb"

        if n_tdb_filepath.exists() or c_tdb_filepath.exists():
            residue_entry_list = apply_all_terminus_patches(
                residue_entry_list, n_tdb_filepath, c_tdb_filepath
            )

        all_residue_list.extend(residue_entry_list)

    for itp_filepath in itp_filepath_list:
        molecule_type = RTP_FILENAME_TO_TYPE.get(itp_filepath.stem.lower(), "Unlisted")
        logger.info(f"Parsing {itp_filepath.name} as {molecule_type}")

        residue_entry_list = parse_itp(itp_filepath, force_field, molecule_type)
        all_residue_list.extend(residue_entry_list)

    logger.debug(f"{len(all_residue_list)} residue entries collected for {force_field}")

    return all_residue_list


def build_residue_database(top_dir: Path) -> ResidueDatabase:
    """
    Build the complete residue database by parsing every force field directory.

    Returns
    -------
    ResidueDatabase
        The complete residue database, grouped by residue name.
    """
    residues: dict[str, list[ResidueEntry]] = {}
    force_fields_processed: list[str] = []

    ff_directory_list = find_force_field_directories(top_dir)

    for ff_directory in ff_directory_list:
        force_field = ff_directory.name.replace(".ff", "")
        residue_entry_list = process_one_force_field_directory(ff_directory)

        for residue_entry in residue_entry_list:
            if residue_entry.name not in residues:
                residues[residue_entry.name] = []
            residues[residue_entry.name].append(residue_entry)

        force_fields_processed.append(force_field)

    total_residues = 0
    for residue_name in residues:
        total_residues = total_residues + len(residues[residue_name])

    logger.success(
        f"Residue database built : {total_residues} entries across {len(force_fields_processed)} force fields"
    )

    residue_database = ResidueDatabase(
        version=datetime.now(UTC).isoformat(),
        residues=residues,
    )

    return residue_database


@click.command()
@click.option(
    "--top-dir",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    help="Path to the GROMACS top directory containing all *.ff force field folders.",
)
@click.option(
    "--database-res-json",
    type=click.Path(path_type=Path),
    help="Output JSON filepath for the complete residue database.",
)
def main(top_dir: Path, database_res_json: Path):
    residue_database = build_residue_database(top_dir)

    with open(database_res_json, "w") as f:
        json.dump(residue_database.model_dump(), f, indent=2)

    logger.success(f"Residue database written to {database_res_json}")


if __name__ == "__main__":
    main()
