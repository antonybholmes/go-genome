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

with open("blacklist.json", "r") as f:
    files = json.load(f)


REGIONS_SQL = """CREATE TABLE regions 
    (id INTEGER PRIMARY KEY ASC,
    chr TEXT NOT NULL DEFAULT 'chr1',
    start INT NOT NULL DEFAULT 1,
    end INT NOT NULL DEFAULT 1,
    notes TEXT NOT NULL DEFAULT '');"""


for file_desc in files:
    print(file_desc)

    db = f"../data/modules/genome/blacklist_{file_desc['assembly']}.db"
    if os.path.exists(db):
        os.remove(db)

    # Connect to the SQLite database
    conn = sqlite3.connect(db)
    cursor = conn.cursor()

    cursor.execute("PRAGMA journal_mode = WAL;")
    cursor.execute("PRAGMA foreign_keys = ON;")
    cursor.execute("BEGIN TRANSACTION;")

    cursor.execute(REGIONS_SQL)

    cursor.execute("END TRANSACTION;")

    cursor.execute("BEGIN TRANSACTION;")

    with open(
        file_desc["file"],
        "r",
    ) as f:
        for line in f:
            if line.startswith("#"):
                continue

            tokens = line.strip().split("\t")

            chromosome = tokens[0]
            start = int(tokens[1]) + 1
            end = int(tokens[2])
            notes = tokens[3] if len(tokens) > 3 else ""

            cursor.execute(
                "INSERT INTO regions (chr, start, end, notes) VALUES (?, ?, ?, ?)",
                (chromosome, start, end, notes),
            )

    cursor.execute("END TRANSACTION;")

    cursor.execute("BEGIN TRANSACTION;")

    cursor.execute(
        "CREATE INDEX idx_regions_chr_start_end_notes ON regions(chr, start, end, notes);"
    )

    cursor.execute(
        f"CREATE TABLE info (id INTEGER PRIMARY KEY, public_id TEXT NOT NULL UNIQUE, genome TEXT NOT NULL UNIQUE, version TEXT NOT NULL);",
    )

    cursor.execute(
        f"INSERT INTO info (public_id, genome, version) VALUES('{uuid.uuid7()}', '{file_desc['assembly']}', '{file_desc['name']}');",
    )

    cursor.execute("END TRANSACTION;")

    # Commit and close
    conn.commit()
    conn.close()
