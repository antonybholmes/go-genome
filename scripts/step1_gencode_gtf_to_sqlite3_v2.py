# -*- coding: utf-8 -*-
"""
Generate a tss distribution for a region file

Created on Thu Jun 26 10:35:40 2014

@author: Antony Holmes
"""

import collections
import gzip
import json
import os
import re
import sqlite3
import sys

import numpy as np
import pandas as pd
import uuid_utils as uuid

tables = ["genes"]

level_map = {"gene": 1, "transcript": 2, "exon": 3}

with open("files.json", "r") as f:
    files = json.load(f)


INFO_SQL = """CREATE TABLE info (id 
    INTEGER PRIMARY KEY,
    public_id TEXT NOT NULL UNIQUE,
    genome TEXT NOT NULL DEFAULT '',
    assembly TEXT NOT NULL DEFAULT '',
    name TEXT NOT NULL DEFAULT '',
    file TEXT NOT NULL DEFAULT ''
);"""

GENE_TYPES_SQL = """CREATE TABLE gene_types (
    id INTEGER PRIMARY KEY,
    public_id TEXT NOT NULL UNIQUE,
    name TEXT NOT NULL UNIQUE);"""

GENES_SQL = """CREATE TABLE genes (
    id INTEGER PRIMARY KEY,
    gene_id TEXT NOT NULL DEFAULT '',
    gene_symbol TEXT NOT NULL DEFAULT '',
    chr TEXT NOT NULL DEFAULT 'chr1',
    start INT NOT NULL DEFAULT 1,
    end INT NOT NULL DEFAULT 1,
    strand TEXT NOT NULL DEFAULT '+',
    gene_type_id INT NOT NULL,
    FOREIGN KEY (gene_type_id) REFERENCES gene_types(id)
);"""

TRANSCRIPT_TYPES_SQL = """CREATE TABLE transcript_types (
    id INTEGER PRIMARY KEY,
    public_id TEXT NOT NULL UNIQUE,
    name TEXT NOT NULL UNIQUE);"""

TRANSCRIPTS_SQL = """CREATE TABLE transcripts (
    id INTEGER PRIMARY KEY,
    gene_id INT NOT NULL,
    transcript_id TEXT NOT NULL DEFAULT '',
    start INT NOT NULL DEFAULT 1,
    end INT NOT NULL DEFAULT 1,
    is_canonical INT NOT NULL DEFAULT 0,
    is_longest INT NOT NULL DEFAULT 0,
    transcript_type_id INT NOT NULL,
    FOREIGN KEY (gene_id) REFERENCES genes(id),
    FOREIGN KEY (transcript_type_id) REFERENCES transcript_types(id)
);"""

# exon ids are stored in a separate table to save space,
# as they are often repeated across exons, cds and utrs
# EXONS_IDS_SQL = """CREATE TABLE exon_ids (
#     id INTEGER PRIMARY KEY,
#     name TEXT NOT NULL UNIQUE
# );"""

FEATURE_TYPES_SQL = """CREATE TABLE feature_types (
    id INTEGER PRIMARY KEY,
    public_id TEXT NOT NULL UNIQUE,
    name TEXT NOT NULL UNIQUE);"""


EXONS_SQL = """CREATE TABLE exons (
    id INTEGER PRIMARY KEY,
    transcript_id INT NOT NULL,
    exon_id TEXT NOT NULL DEFAULT '',
    exon_number INT NOT NULL DEFAULT 1,
    UNIQUE(transcript_id, exon_id, exon_number),
    FOREIGN KEY (transcript_id) REFERENCES transcripts(id)
);"""

FEATURES_SQL = """CREATE TABLE features (
    id INTEGER PRIMARY KEY,
    transcript_id INT NOT NULL,
    exon_id INTEGER NOT NULL,
    feature_type_id INT NOT NULL,
    start INT NOT NULL DEFAULT 1,
    end INT NOT NULL DEFAULT 1,
    UNIQUE(transcript_id, exon_id, feature_type_id, start, end),
    FOREIGN KEY (transcript_id) REFERENCES transcripts(id),
    FOREIGN KEY (exon_id) REFERENCES exons(id),
    FOREIGN KEY (feature_type_id) REFERENCES feature_types(id)
);"""

feature_type_map = {
    "exon": 1,
    "CDS": 2,
    "UTR": 3,
    "start_codon": 4,
    "stop_codon": 5,
    "five_prime_utr": 6,
    "three_prime_utr": 7,
}

for file_desc in files:
    print(file_desc)

    db = f"../data/modules/genome/gtf_{file_desc['version']}.db"
    if os.path.exists(db):
        os.remove(db)

    # Connect to the SQLite database
    conn = sqlite3.connect(db)
    cursor = conn.cursor()

    cursor.execute("PRAGMA journal_mode = WAL;")
    cursor.execute("PRAGMA foreign_keys = ON;")

    cursor.execute(INFO_SQL)
    cursor.execute(GENE_TYPES_SQL)
    cursor.execute(TRANSCRIPT_TYPES_SQL)
    cursor.execute(GENES_SQL)
    cursor.execute(TRANSCRIPTS_SQL)
    # cursor.execute(EXONS_IDS_SQL)
    cursor.execute(FEATURE_TYPES_SQL)
    cursor.execute(EXONS_SQL)
    cursor.execute(FEATURES_SQL)

    cursor.execute(
        f"INSERT INTO info (public_id, genome, assembly, name, file) VALUES('{uuid.uuid7()}', '{file_desc['genome']}', '{file_desc['assembly']}', '{file_desc['version']}', '{file_desc['file']}');"
    )

    cursor.execute(
        f"INSERT INTO feature_types (id, public_id, name) VALUES (1, '{uuid.uuid7()}', 'exon');"
    )
    cursor.execute(
        f"INSERT INTO feature_types (id, public_id, name) VALUES (2, '{uuid.uuid7()}', 'cds');"
    )
    cursor.execute(
        f"INSERT INTO feature_types (id, public_id, name) VALUES (3, '{uuid.uuid7()}', 'utr');"
    )

    cursor.execute(
        f"INSERT INTO feature_types (id, public_id, name) VALUES (4, '{uuid.uuid7()}', 'start_codon');"
    )

    cursor.execute(
        f"INSERT INTO feature_types (id, public_id, name) VALUES (5, '{uuid.uuid7()}', 'stop_codon');"
    )

    cursor.execute(
        f"INSERT INTO feature_types (id, public_id, name) VALUES (6, '{uuid.uuid7()}', 'five_prime_utr');"
    )
    cursor.execute(
        f"INSERT INTO feature_types (id, public_id, name) VALUES (7, '{uuid.uuid7()}', 'three_prime_utr');"
    )

    record = 1

    gene_id = ""
    gene_name = ""
    gene_type = ""
    transcript_id = ""
    transcript_type = ""
    exon_id = ""
    chr = ""
    start = 0
    end = 0
    stranded_start = 0
    stranded_end = 0

    gene_map = {}
    gene_types = {}
    exon_map = {}
    transcript_map = {}
    transcript_types = {}

    with gzip.open(
        file_desc["file"],
        "rt",
    ) as f:
        for line in f:
            if line.startswith("#"):
                continue

            tokens = line.strip().split("\t")

            level = tokens[2]

            # if 'tag "basic"' not in line:
            #    continue

            # print(line)

            tags = set()

            if (
                "Ensembl_canonical" in line
                or "appris_principal" in line
                or "MANE_select" in line
            ):
                tags.add("canonical:true")
                is_canonical = 1
            else:
                is_canonical = 0

            # gene
            matcher = re.search(r'gene_id "(.+?)";', tokens[8])

            if matcher:
                # remove version
                gene_id = re.sub(r"\..+", "", matcher.group(1))

                if gene_id not in gene_map:
                    gene_map[gene_id] = len(gene_map) + 1

            matcher = re.search(r'gene_name "(.+?)";', tokens[8])

            if matcher:
                gene_name = matcher.group(1)

            # transcript_type
            matcher = re.search(r'gene_type "(.+?)";', tokens[8])

            if matcher:
                gene_type = re.sub(r"\..+", "", matcher.group(1))
                tags.add(f"gene_type:{gene_type}")

            # transcript
            matcher = re.search(r'transcript_id "(.+?)";', tokens[8])

            if matcher:
                transcript_id = re.sub(r"\..+", "", matcher.group(1))

                if transcript_id not in transcript_map:
                    transcript_map[transcript_id] = len(transcript_map) + 1

            matcher = re.search(r'transcript_type "(.+?)";', tokens[8])

            if matcher:
                transcript_type = re.sub(r"\..+", "", matcher.group(1))
                tags.add(f"transcript_type:{transcript_type}")

            # exon

            matcher = re.search(r"exon_number (\d+);", tokens[8])

            if matcher:
                exon_number = int(matcher.group(1))

            matcher = re.search(r'exon_id "(.+?)";', tokens[8])

            if matcher:
                exon_id = re.sub(r"\..+", "", matcher.group(1))

                if exon_id not in exon_map:
                    exon_map[exon_id] = len(exon_map) + 1

            if gene_type not in gene_types:
                gene_types[gene_type] = len(gene_types) + 1

                cursor.execute(
                    "INSERT INTO gene_types (id, public_id, name) VALUES (?, ?, ?)",
                    (
                        gene_types[gene_type],
                        str(uuid.uuid7()),
                        gene_type,
                    ),
                )

            if transcript_type not in transcript_types:
                transcript_types[transcript_type] = len(transcript_types) + 1

                cursor.execute(
                    "INSERT INTO transcript_types (id, public_id, name) VALUES (?, ?, ?)",
                    (
                        transcript_types[transcript_type],
                        str(uuid.uuid7()),
                        transcript_type,
                    ),
                )

            chr = tokens[0]
            start = int(tokens[3])
            end = int(tokens[4])
            # mid = int((end + start) / 2)
            strand = tokens[6]

            # invert coordinates to make searching easier
            if strand == "-":
                stranded_start = end
                # stranded_end = start
            else:
                stranded_start = start
                # stranded_end = end

            tag_str = ",".join(sorted(tags))

            if level == "gene":

                cursor.execute(
                    "INSERT INTO genes (id, gene_id, gene_symbol, chr, start, end, strand, gene_type_id) VALUES (?, ?, ?, ?, ?, ?, ?, ?) ON CONFLICT DO NOTHING",
                    (
                        gene_map[gene_id],
                        gene_id,
                        gene_name,
                        chr,
                        start,
                        end,
                        strand,
                        gene_types[gene_type],
                    ),
                )
                # record = cursor.lastrowid
            elif level == "transcript":

                cursor.execute(
                    "INSERT INTO transcripts (id, gene_id, transcript_id, start, end, is_canonical, transcript_type_id) VALUES (?, ?, ?, ?, ?, ?, ?) ON CONFLICT DO NOTHING",
                    (
                        transcript_map[transcript_id],
                        gene_map[gene_id],
                        transcript_id,
                        start,
                        end,
                        is_canonical,
                        transcript_types[transcript_type],
                    ),
                )
                # record = cursor.lastrowid
            else:

                cursor.execute(
                    "INSERT INTO exons (id, transcript_id, exon_id, exon_number) VALUES (?, ?, ?, ?) ON CONFLICT DO NOTHING",
                    (
                        exon_map[exon_id],
                        transcript_map[transcript_id],
                        exon_id,
                        exon_number,
                    ),
                )

                if level in feature_type_map:

                    cursor.execute(
                        "INSERT INTO features (transcript_id, exon_id, feature_type_id, start, end) VALUES (?, ?, ?, ?, ?) ON CONFLICT DO NOTHING",
                        (
                            transcript_map[transcript_id],
                            exon_map[exon_id],
                            feature_type_map[level],
                            start,
                            end,
                        ),
                    )

            record += 1

    cursor.execute("CREATE INDEX idx_genes_gene_id ON genes(gene_id);")
    cursor.execute("CREATE INDEX idx_genes_gene_symbol ON genes(gene_symbol);")
    cursor.execute("CREATE INDEX idx_genes_gene_type_id ON genes(gene_type_id);")
    cursor.execute(
        "CREATE INDEX idx_genes_chr_start_end_strand ON genes(chr, start, end, strand);"
    )

    cursor.execute("CREATE INDEX idx_transcripts_gene_id ON transcripts(gene_id);")
    cursor.execute(
        "CREATE INDEX idx_transcripts_transcript_id ON transcripts(transcript_id);"
    )
    cursor.execute(
        "CREATE INDEX idx_transcripts_transcript_type_id ON transcripts(transcript_type_id);"
    )
    cursor.execute("CREATE INDEX idx_transcripts_start_end ON transcripts(start, end);")
    cursor.execute(
        "CREATE INDEX idx_transcripts_is_canonical ON transcripts(is_canonical);"
    )
    cursor.execute(
        "CREATE INDEX idx_transcripts_is_longest ON transcripts(is_longest);"
    )

    cursor.execute("CREATE INDEX idx_exons_exon_id ON exons(exon_id);")

    cursor.execute("CREATE INDEX idx_exons_transcript_id ON exons(transcript_id);")

    cursor.execute(
        "CREATE INDEX idx_features_feature_type_id_start_end ON features(feature_type_id, start, end);"
    )

    # work out who is longest transcript per gene
    print("Finding longest transcripts...")

    cursor.execute(
        f"""
        WITH max_lengths AS (
            SELECT 
                t.gene_id,
                t.transcript_id,
                ROW_NUMBER() OVER (PARTITION BY t.gene_id ORDER BY ABS(t.end - t.start) DESC) AS rank
            FROM transcripts t
        )
        SELECT gene_id, transcript_id FROM max_lengths
        WHERE rank = 1
        ORDER BY gene_id, transcript_id;
        """,
    )

    rows = cursor.fetchall()

    queries = []

    for row in rows:
        # print(row)
        queries.append({"gene_id": row[0], "transcript_id": row[1]})

    cursor.executemany(
        "UPDATE transcripts SET is_longest = 1 WHERE gene_id = :gene_id AND transcript_id = :transcript_id",
        queries,
    )

    print(len(queries), "transcripts to set as longest")

    # Commit and close
    conn.commit()
    conn.close()
