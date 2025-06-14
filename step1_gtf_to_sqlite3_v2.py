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

tables = ["genes"]

level_map = {"gene": 1, "transcript": 2, "exon": 3}

files = [
    [
        ["hg19", "Gencode v38"],
        "/home/antony/development/gencode.v48lift37.basic.annotation.gtf.gz",
    ],
    [
        ["grch38", "Gencode v44"],
        "/home/antony/development/gencode.v48.basic.annotation.gtf.gz",
    ],
    [
        ["mm10", "Gencode vM25"],
        "/home/antony/development/gencode.vM25.basic.annotation.gtf.gz",
    ],
]

GENES_SQL = """CREATE TABLE genes 
    (id INTEGER PRIMARY KEY ASC,
    gene_id TEXT NOT NULL DEFAULT '',
    gene_symbol TEXT NOT NULL DEFAULT '',
    chr TEXT NOT NULL DEFAULT 'chr1',
    start INT NOT NULL DEFAULT 1,
    end INT NOT NULL DEFAULT 1,
    strand CHAR(1) NOT NULL DEFAULT '+',
    is_canonical INT NOT NULL DEFAULT 0,
    gene_type TEXT NOT NULL DEFAULT '');"""

TRANSCRIPTS_SQL = """CREATE TABLE transcripts 
    (id INTEGER PRIMARY KEY ASC,
    gene_id INT NOT NULL REFERENCES genes(id),
    transcript_id TEXT NOT NULL DEFAULT '',
    start INT NOT NULL DEFAULT 1,
    end INT NOT NULL DEFAULT 1);"""

EXONS_SQL = """CREATE TABLE exons 
    (id INTEGER PRIMARY KEY ASC,
    transcript_id INT NOT NULL REFERENCES transcripts(id),
    exon_id TEXT NOT NULL DEFAULT '',
    exon_number INT NOT NULL DEFAULT 1,
    start INT NOT NULL DEFAULT 1,
    end INT NOT NULL DEFAULT 1);"""

for file_desc in files:
    print(file_desc)

    db = f"data/modules/genome/{file_desc[0][0]}.db"
    if os.path.exists(db):
        os.remove(db)

    # Connect to the SQLite database
    conn = sqlite3.connect(f"data/modules/genome/{file_desc[0][0]}.db")
    cursor = conn.cursor()

    cursor.execute("PRAGMA journal_mode = WAL;")
    cursor.execute("PRAGMA foreign_keys = ON;")
    cursor.execute(GENES_SQL)
    cursor.execute(TRANSCRIPTS_SQL)
    cursor.execute(EXONS_SQL)
    cursor.execute("CREATE INDEX idx_genes_gene_id ON genes(gene_id);")
    cursor.execute(
        "CREATE INDEX idx_transcripts_transcript_id ON transcripts(transcript_id);"
    )
    cursor.execute("CREATE INDEX idx_exons_exon_id ON exons(exon_id);")
    cursor.execute(
        "CREATE INDEX idx_genes_chr_start_end_strand ON genes(chr, start, end, strand);"
    )

    cursor.execute(
        f"CREATE TABLE info (id INTEGER PRIMARY KEY ASC, genome TEXT NOT NULL, version TEXT NOT NULL);",
    )

    cursor.execute(
        f"INSERT INTO info (genome, version) VALUES('{file_desc[0][0]}', '{file_desc[0][1]}');",
    )

    cursor.execute("END TRANSACTION;")

    cursor.execute("BEGIN TRANSACTION;")

    record = 1
    gene_record_id = 0
    transcript_record_id = 0

    with gzip.open(
        file_desc[1],
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

            is_canonical = 0
            tags = set()
            parent_record_id = -1
            gene_id = "n/a"
            gene_name = "n/a"
            gene_type = "n/a"
            transcript_id = "n/a"
            exon_id = "n/a"
            chr = ""
            start = 0
            end = 0
            stranded_start = 0
            stranded_end = 0

            if level == "gene":
                parent_record_id = -1
                gene_record_id += 1

            if level == "transcript":
                parent_record_id = gene_record_id
                transcript_record_id += 1
                exon_number = 1

            if level == "exon":
                parent_record_id = transcript_record_id
                exon_record_id = record

            # gene
            matcher = re.search(r'gene_id "(.+?)";', tokens[8])

            if matcher:
                # remove version
                gene_id = re.sub(r"\..+", "", matcher.group(1))

            matcher = re.search(r'gene_name "(.+?)";', tokens[8])

            if matcher:
                gene_name = matcher.group(1)

            # transcript_type
            matcher = re.search(r'gene_type "(.+?)";', tokens[8])

            if matcher:
                gene_type = re.sub(r"\..+", "", matcher.group(1))
                tags.add(f"type:{gene_type}")

            # transcript
            matcher = re.search(r'transcript_id "(.+?)";', tokens[8])

            if matcher:
                transcript_id = re.sub(r"\..+", "", matcher.group(1))

            if "Ensembl_canonical" in line:
                tags.add("canonical:true")
                is_canonical = 1

            # exon
            matcher = re.search(r'exon_id "(.+?)";', tokens[8])

            if matcher:
                exon_id = re.sub(r"\..+", "", matcher.group(1))

            chr = tokens[0]
            start = int(tokens[3])
            end = int(tokens[4])
            # mid = int((end + start) / 2)
            strand = tokens[6]

            # invert coordinates to make searching easier
            # if strand == "-":
            #     stranded_start = end
            #     stranded_end = start
            # else:
            #     stranded_start = start
            #     stranded_end = end

            tag_str = ",".join(sorted(tags))

            # gene_id = id_map.get(gene_id, 1)
            # gene_name = id_map.get(gene_name, 1)
            # transcript_id = id_map.get(transcript_id, 1)
            # exon_id = id_map.get(exon_id, 1)
            # gene_type = id_map.get(gene_type, 1)

            if level == "gene":
                cursor.execute(
                    "INSERT INTO genes (gene_id, gene_symbol, chr, start, end, strand, is_canonical, gene_type) VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                    (
                        gene_id,
                        gene_name,
                        chr,
                        start,
                        end,
                        strand,
                        is_canonical,
                        gene_type,
                    ),
                )
                # record = cursor.lastrowid
            elif level == "transcript":
                cursor.execute(
                    "INSERT INTO transcripts (transcript_id, gene_id, start, end) VALUES (?, ?, ?, ?)",
                    (
                        transcript_id,
                        gene_record_id,
                        start,
                        end,
                    ),
                )
                # record = cursor.lastrowid
            elif level == "exon":

                cursor.execute(
                    "INSERT INTO exons (exon_id, transcript_id, exon_number, start, end) VALUES (?, ?, ?, ?, ?)",
                    (
                        exon_id,
                        transcript_record_id,
                        exon_number,
                        start,
                        end,
                    ),
                )
                exon_number += 1
                # record = cursor.lastrowid
            else:
                # print(f"Unknown level: {level}")
                pass
            # break

            record += 1

    # cursor.execute("COMMIT;")

    # Commit and close
    conn.commit()
    conn.close()
