# -*- coding: utf-8 -*-
"""
Generate a tss distribution for a region file

Created on Thu Jun 26 10:35:40 2014

@author: Antony Holmes
"""

import sys
import collections
import re
import os
import pandas as pd
import numpy as np
import gzip
import sqlite3


files = [
    [
        ["grch38", "GCB TAD v2"],
        "/ifs/archive/cancer/Lab_RDF/scratch_Lab_RDF/ngs/scripts/python/gal/assets/human/grch38/tads_20250902.tsv",
    ],
]

TAD_SQL = """CREATE TABLE tads (
    id INTEGER PRIMARY KEY ASC,
    chr TEXT NOT NULL DEFAULT 'chr1',
    start INT NOT NULL DEFAULT 1,
    end INT NOT NULL DEFAULT 1,
    gene_ids TEXT NOT NULL DEFAULT '',
    gene_names TEXT NOT NULL DEFAULT ''
);"""


for file_desc in files:
    db = f"data/modules/genome/tads_{file_desc[0][0]}.db"
    if os.path.exists(db):
        os.remove(db)

    # Connect to the SQLite database
    conn = sqlite3.connect(db)
    cursor = conn.cursor()

    cursor.execute("PRAGMA journal_mode = WAL;")
    cursor.execute("PRAGMA foreign_keys = ON;")
    cursor.execute("BEGIN TRANSACTION;")

    cursor.execute(TAD_SQL)

    cursor.execute("END TRANSACTION;")

    cursor.execute("BEGIN TRANSACTION;")

    df = pd.read_csv(
        file_desc[1],
        sep="\t",
        header=0,
        keep_default_na=False,
    )

    for _, row in df.iterrows():

        cursor.execute(
            "INSERT INTO tads (chr, start, end, gene_ids, gene_names) VALUES (?, ?, ?, ?, ?)",
            (row["chr"], row["start"], row["end"], row["gene_ids"], row["gene_names"]),
        )

    cursor.execute("END TRANSACTION;")

    cursor.execute("BEGIN TRANSACTION;")

    cursor.execute("CREATE INDEX idx_tads_chr_start_end ON tads(chr, start, end);")

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
