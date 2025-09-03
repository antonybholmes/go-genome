# -*- coding: utf-8 -*-
"""
Generate a tss distribution for a region file

Created on Thu Jun 26 10:35:40 2014

@author: Antony Holmes
"""

import json
import sys
import collections
import re
import os
from attrs import fields
import pandas as pd
import numpy as np
import gzip
import sqlite3

tables = ["genes"]

feature_map = {"gene": 1, "transcript": 2, "exon": 3}

files = [
    [
        ["hg19", "Gencode GRCh37 v48"],
        "/home/antony/development/gencode.v48lift37.basic.annotation.gtf.gz",
    ],
    [
        ["grch38", "Gencode GRCh38 v48"],
        "/home/antony/development/gencode.v48.basic.annotation.gtf.gz",
    ],
    [
        ["mm10", "Gencode mm10 vM25"],
        "/home/antony/development/gencode.vM25.basic.annotation.gtf.gz",
    ],
]

FEATURE_SQL = """CREATE TABLE features (
    id INTEGER PRIMARY KEY,
    name TEXT NOT NULL
);
"""


GTF_SQL = """CREATE TABLE gtf (
    id INTEGER PRIMARY KEY,
    seqname TEXT NOT NULL DEFAULT '',
    source TEXT NOT NULL DEFAULT '',
    feature TEXT NOT NULL DEFAULT 'gene',
    start INTEGER NOT NULL DEFAULT 1,
    end INTEGER NOT NULL DEFAULT 1,
    tss INTEGER NOT NULL DEFAULT 1,
    score TEXT NOT NULL DEFAULT '.',
    strand TEXT NOT NULL DEFAULT '+',
    frame TEXT NOT NULL DEFAULT '.',
    gene_id TEXT NOT NULL DEFAULT '',
    gene_name TEXT NOT NULL DEFAULT '',
    transcript_id TEXT NOT NULL DEFAULT '',
    transcript_name TEXT NOT NULL DEFAULT '',
    exon_number INTEGER NOT NULL DEFAULT 1,
    type TEXT NOT NULL DEFAULT '',
    is_canonical BOOLEAN NOT NULL DEFAULT 0,
    attributes TEXT NOT NULL DEFAULT ''
);
"""

GTF_ATTRIBUTES_SQL = """CREATE TABLE gtf_attributes (
    id INTEGER PRIMARY KEY,
    gtf_id INTEGER,
    key TEXT NOT NULL,
    value TEXT NOT NULL DEFAULT '',
    FOREIGN KEY (gtf_id) REFERENCES gtf(id)
);
"""

VALID_GENE_TYPES = {"protein_coding"}  # , "lincRNA", "antisense", "pseudogene"}


def parse_gtf_line(line):
    if line.startswith("#"):
        return None

    fields = line.strip().split("\t")
    if len(fields) != 9:
        return None

    source = fields[1]
    feature = fields[2]

    if not re.match(r"^(gene|transcript|exon)$", feature):
        return None

    is_canonical = "Ensembl_canonical" in line or "appris_principal" in line

    gtf = {
        "seqname": fields[0],
        "source": source,
        "feature": feature,
        # "feature_id": feature_map.get(feature, 1),
        "start": int(fields[3]),
        "end": int(fields[4]),
        "score": fields[5],
        "strand": fields[6],
        "frame": fields[7],
        "is_canonical": is_canonical,
    }

    gtf["tss"] = int(fields[3]) if gtf["strand"] == "+" else int(fields[4])

    attr_dict = {}
    attributes = fields[8].strip().split(";")
    for attr in attributes:
        if attr.strip():
            key, value = attr.strip().split(" ", 1)
            attr_dict[key] = value.strip('"').strip()

    atts = json.dumps(attr_dict, sort_keys=True, separators=(",", ":"))

    # if attr_dict.get("gene_type", "") not in VALID_GENE_TYPES:
    #    return None

    # Add commonly used attributes (add more if needed)
    gtf["gene_id"] = attr_dict.get("gene_id", "")
    gtf["type"] = attr_dict.get("gene_type", "")
    gtf["transcript_id"] = attr_dict.get("transcript_id", "")
    gtf["gene_name"] = attr_dict.get("gene_name", "")
    gtf["transcript_name"] = attr_dict.get("transcript_name", "")
    gtf["exon_number"] = int(attr_dict.get("exon_number", "1"))
    gtf["attributes"] = ""  # atts

    return gtf


def insert_gtf_line(cur, insert_row, used_attributes):
    cur.execute(
        f"""
        INSERT INTO gtf ({', '.join(insert_row.keys())})
        VALUES ({', '.join(':' + k for k in insert_row.keys())})
    """,
        insert_row,
    )


for file_desc in files:
    print(file_desc)

    db = f"data/modules/genome/gtf_{file_desc[0][0]}.db"
    if os.path.exists(db):
        os.remove(db)

    # Connect to the SQLite database
    conn = sqlite3.connect(db)
    cursor = conn.cursor()

    cursor.execute("PRAGMA journal_mode = WAL;")
    cursor.execute("PRAGMA foreign_keys = ON;")
    cursor.execute("BEGIN TRANSACTION;")

    cursor.execute(GTF_SQL)

    cursor.execute("END TRANSACTION;")

    cursor.execute("BEGIN TRANSACTION;")

    used_attributes = {}

    with gzip.open(
        file_desc[1],
        "rt",
    ) as f:
        for line in f:
            parsed = parse_gtf_line(line)
            if parsed:
                insert_gtf_line(cursor, parsed, used_attributes)

    cursor.execute("END TRANSACTION;")

    cursor.execute("BEGIN TRANSACTION;")

    cursor.execute(
        "CREATE INDEX idx_gtf_feature ON gtf(feature, seqname, tss, strand);"
    )
    cursor.execute("CREATE INDEX idx_gtf_coords ON gtf(seqname, start, end);")

    cursor.execute("CREATE INDEX idx_gtf_type ON gtf(type);")

    cursor.execute("CREATE INDEX idx_gtf_is_canonical ON gtf(is_canonical);")

    cursor.execute(
        f"CREATE TABLE info (id INTEGER PRIMARY KEY ASC, genome TEXT NOT NULL, version TEXT NOT NULL);",
    )

    cursor.execute(
        f"INSERT INTO info (genome, version) VALUES('{file_desc[0][0]}', '{file_desc[0][1]}');",
    )

    cursor.execute("END TRANSACTION;")

    # Commit and close
    conn.commit()
    conn.close()
