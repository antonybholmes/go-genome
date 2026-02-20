# -*- coding: utf-8 -*-
"""
Encode read counts per base in 2 bytes

@author: Antony Holmes
"""
import argparse
import os
import sqlite3

import uuid_utils as uuid

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dir", help="sample name")
args = parser.parse_args()

dir = args.dir  # sys.argv[1]

data = []


db = os.path.join(dir, "genome.db")


if os.path.exists(db):
    os.remove(db)

conn = sqlite3.connect(db)
cursor = conn.cursor()

cursor.execute("PRAGMA journal_mode = WAL;")
cursor.execute("PRAGMA foreign_keys = ON;")

cursor.execute(
    f"""
    CREATE TABLE genomes (
        id INTEGER PRIMARY KEY,
        public_id TEXT NOT NULL UNIQUE,
        name TEXT NOT NULL,
        scientific_name TEXT NOT NULL,
        UNIQUE(name, scientific_name));
    """,
)

cursor.execute("CREATE INDEX idx_genomes_name_id ON genomes(LOWER(name));")

cursor.execute(
    f"INSERT INTO genomes (id, public_id, name, scientific_name) VALUES (1, '{uuid.uuid7()}', 'Human', 'Homo sapiens');"
)
cursor.execute(
    f"INSERT INTO genomes (id, public_id, name, scientific_name) VALUES (2, '{uuid.uuid7()}', 'Mouse', 'Mus musculus');"
)

cursor.execute(
    f"""
    CREATE TABLE assemblies (
        id INTEGER PRIMARY KEY,
        public_id TEXT NOT NULL UNIQUE,
        genome_id INTEGER NOT NULL,
        name TEXT NOT NULL UNIQUE,
        FOREIGN KEY (genome_id) REFERENCES genomes(id) ON DELETE CASCADE);
    """,
)

cursor.execute("CREATE INDEX idx_assemblies_name_id ON assemblies(LOWER(name));")
cursor.execute("CREATE INDEX idx_assemblies_genome_id ON assemblies(genome_id);")

cursor.execute(
    f"INSERT INTO assemblies (id, public_id, genome_id, name) VALUES (1, '{uuid.uuid7()}', 1, 'GRCh37');"
)
cursor.execute(
    f"INSERT INTO assemblies (id, public_id, genome_id, name) VALUES (2, '{uuid.uuid7()}', 1, 'GRCh38');"
)
cursor.execute(
    f"INSERT INTO assemblies (id, public_id, genome_id, name) VALUES (3, '{uuid.uuid7()}', 2, 'GRCm38');"
)

cursor.execute(
    f"INSERT INTO assemblies (id, public_id, genome_id, name) VALUES (4, '{uuid.uuid7()}', 2, 'GRCm39');"
)

cursor.execute(
    f"""
    CREATE TABLE assembly_aliases (
        id INTEGER PRIMARY KEY,
        assembly_id INTEGER NOT NULL,
        alias TEXT NOT NULL UNIQUE,
        FOREIGN KEY (assembly_id) REFERENCES assemblies(id) ON DELETE CASCADE);
    """,
)

cursor.execute(f"INSERT INTO assembly_aliases (assembly_id, alias) VALUES (1, 'hg19');")
cursor.execute(
    f"INSERT INTO assembly_aliases (assembly_id, alias) VALUES (1, 'GRCh37');"
)

cursor.execute(f"INSERT INTO assembly_aliases (assembly_id, alias) VALUES (2, 'hg38');")
cursor.execute(
    f"INSERT INTO assembly_aliases (assembly_id, alias) VALUES (2, 'GRCh38');"
)

cursor.execute(f"INSERT INTO assembly_aliases (assembly_id, alias) VALUES (3, 'mm10');")
cursor.execute(
    f"INSERT INTO assembly_aliases (assembly_id, alias) VALUES (3, 'GRCm38');"
)

cursor.execute(f"INSERT INTO assembly_aliases (assembly_id, alias) VALUES (4, 'mm39');")
cursor.execute(
    f"INSERT INTO assembly_aliases (assembly_id, alias) VALUES (4, 'GRCm39');"
)


assembly_map = {
    "hg19": 1,
    "grch37": 1,
    "grch38": 2,
    "mm10": 3,
    "grcm38": 3,
    "mm39": 4,
    "grcm39": 4,
}

cursor.execute(
    f"""
    CREATE TABLE annotation_types (
        id INTEGER PRIMARY KEY,
        public_id TEXT NOT NULL UNIQUE,
        name TEXT NOT NULL UNIQUE);
    """,
)

cursor.execute(
    "CREATE INDEX idx_annotation_types_name_id ON annotation_types(LOWER(name));"
)

cursor.execute(
    f"INSERT INTO annotation_types (id, public_id, name) VALUES (1, '{uuid.uuid7()}', 'GTF');"
)

cursor.execute(
    f""" CREATE TABLE annotations (
	id INTEGER PRIMARY KEY,
    public_id TEXT NOT NULL UNIQUE,
    assembly_id INTEGER NOT NULL,
    annotation_type_id INTEGER NOT NULL,
	name TEXT NOT NULL UNIQUE,
    url TEXT NOT NULL,
	description TEXT NOT NULL DEFAULT '',
	tags TEXT NOT NULL DEFAULT '',
    FOREIGN KEY (assembly_id) REFERENCES assemblies(id),
    FOREIGN KEY (annotation_type_id) REFERENCES annotation_types(id)
);
"""
)

type_map = {
    "gtf": 1,
}

for root, dirs, files in os.walk(dir):
    if "trash" in root:
        continue

    for filename in files:
        if "gtf" in filename and filename.endswith(".db"):
            # relative_dir = root.replace(dir, "")[1:]

            # genome, platform, dataset = relative_dir.split("/")

            # filepath = os.path.join(root, filename)
            # print(root, filename, relative_dir, platform, genome, dataset,)

            # path = os.path.join(root, filename)

            # gex_path = os.path.join(relative_dir, "gex")

            conn2 = sqlite3.connect(os.path.join(root, filename))
            conn2.row_factory = sqlite3.Row

            print(filename)

            # Create a cursor object
            cursor2 = conn2.cursor()

            cursor2.execute("SELECT public_id, name, assembly FROM info")

            results = cursor2.fetchone()

            print(results["assembly"])
            # genome = genome_map[results["genome"]]
            assembly = assembly_map[results["assembly"]]

            cursor.execute(
                f"""INSERT INTO annotations (public_id, assembly_id, annotation_type_id, name, url) VALUES (
                '{results["public_id"]}',
                '{assembly}',
                '{type_map["gtf"]}',
                '{results["name"]}',
                '{filename}');"""
            )

            conn2.close()


cursor.execute(
    f""" CREATE INDEX idx_annotations_name ON annotations (LOWER(name));
"""
)

cursor.execute(
    f""" CREATE INDEX idx_annotations_assembly_id ON annotations (assembly_id);
"""
)

cursor.execute(
    f""" CREATE INDEX idx_annotations_annotation_type_id ON annotations (annotation_type_id);
"""
)

cursor.execute(
    f""" CREATE INDEX idx_assembly_aliases_alias ON assembly_aliases (LOWER(alias));
"""
)

conn.commit()
conn.close()
