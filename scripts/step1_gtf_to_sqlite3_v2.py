# -*- coding: utf-8 -*-
"""
Generate a tss distribution for a region file

Created on Thu Jun 26 10:35:40 2014

@author: Antony Holmes
"""

import collections
import gzip
import os
import re
import sqlite3
import sys

import numpy as np
import pandas as pd

tables = ["genes"]

level_map = {"gene": 1, "transcript": 2, "exon": 3}

files = [
    [
        ["hg19", "Gencode GRCh37 v48"],
        "/ifs/archive/cancer/Lab_RDF/scratch_Lab_RDF/ngs/references/gencode/grch37/gencode.v48lift37.basic.annotation.gtf.gz",
    ],
    [
        ["grch38", "Gencode GRCh38 v48"],
        "/ifs/archive/cancer/Lab_RDF/scratch_Lab_RDF/ngs/references/gencode/grch38/gencode.v48.basic.annotation.gtf.gz",
    ],
    [
        ["mm10", "Gencode mm10 vM25"],
        "/ifs/archive/cancer/Lab_RDF/scratch_Lab_RDF/ngs/references/gencode/mm10/gencode.vM25.basic.annotation.gtf.gz",
    ],
]

GENE_TYPES_SQL = """CREATE TABLE gene_types 
    (id INTEGER PRIMARY KEY ASC,
    name TEXT NOT NULL UNIQUE DEFAULT '');"""

GENES_SQL = """CREATE TABLE genes 
    (id INTEGER PRIMARY KEY ASC,
    gene_id TEXT NOT NULL DEFAULT '',
    gene_symbol TEXT NOT NULL DEFAULT '',
    chr TEXT NOT NULL DEFAULT 'chr1',
    start INT NOT NULL DEFAULT 1,
    end INT NOT NULL DEFAULT 1,
    strand TEXT NOT NULL DEFAULT '+',
    gene_type_id INT NOT NULL,
    FOREIGN KEY (gene_type_id) REFERENCES gene_types(id)
);"""

TRANSCRIPT_TYPES_SQL = """CREATE TABLE transcript_types 
    (id INTEGER PRIMARY KEY ASC,
    name TEXT NOT NULL UNIQUE DEFAULT '');"""

TRANSCRIPTS_SQL = """CREATE TABLE transcripts 
    (id INTEGER PRIMARY KEY ASC,
    gene_id INT NOT NULL,
    transcript_id TEXT NOT NULL DEFAULT '',
    chr TEXT NOT NULL DEFAULT 'chr1',
    start INT NOT NULL DEFAULT 1,
    end INT NOT NULL DEFAULT 1,
    strand TEXT NOT NULL DEFAULT '+',
    is_canonical INT NOT NULL DEFAULT 0,
    is_longest INT NOT NULL DEFAULT 0,
    transcript_type_id INT NOT NULL,
    FOREIGN KEY (gene_id) REFERENCES genes(id),
    FOREIGN KEY (transcript_type_id) REFERENCES transcript_types(id)
);"""


EXONS_SQL = """CREATE TABLE exons 
    (id INTEGER PRIMARY KEY ASC,
    transcript_id INT NOT NULL,
    exon_id TEXT NOT NULL DEFAULT '',
    chr TEXT NOT NULL DEFAULT 'chr1',
    start INT NOT NULL DEFAULT 1,
    end INT NOT NULL DEFAULT 1,
    strand TEXT NOT NULL DEFAULT '+',
    exon_number INT NOT NULL DEFAULT 1,
    FOREIGN KEY (transcript_id) REFERENCES transcripts(id)
);"""


for file_desc in files:
    print(file_desc)

    db = f"../data/modules/genome/{file_desc[0][0]}.db"
    if os.path.exists(db):
        os.remove(db)

    # Connect to the SQLite database
    conn = sqlite3.connect(db)
    cursor = conn.cursor()

    cursor.execute("PRAGMA journal_mode = WAL;")
    cursor.execute("PRAGMA foreign_keys = ON;")
    cursor.execute("BEGIN TRANSACTION;")

    cursor.execute(GENE_TYPES_SQL)
    cursor.execute(TRANSCRIPT_TYPES_SQL)
    cursor.execute(GENES_SQL)
    cursor.execute(TRANSCRIPTS_SQL)
    cursor.execute(EXONS_SQL)

    cursor.execute("END TRANSACTION;")

    cursor.execute("BEGIN TRANSACTION;")

    record = 1
    gene_record_id = 0
    transcript_record_id = 0

    gene_types = {}
    transcript_types = {}

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

            tags = set()
            parent_record_id = -1
            gene_id = "n/a"
            gene_name = "n/a"
            gene_type = "n/a"
            transcript_id = "n/a"
            transcript_type = "n/a"
            exon_id = "n/a"
            chr = ""
            start = 0
            end = 0
            stranded_start = 0
            stranded_end = 0

            if "Ensembl_canonical" in line or "appris_principal" in line:
                tags.add("canonical:true")
                is_canonical = 1
            else:
                is_canonical = 0

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
                tags.add(f"gene_type:{gene_type}")

            # transcript
            matcher = re.search(r'transcript_id "(.+?)";', tokens[8])

            if matcher:
                transcript_id = re.sub(r"\..+", "", matcher.group(1))

            matcher = re.search(r'transcript_type "(.+?)";', tokens[8])

            if matcher:
                transcript_type = re.sub(r"\..+", "", matcher.group(1))
                tags.add(f"transcript_type:{transcript_type}")

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
            if strand == "-":
                stranded_start = end
                # stranded_end = start
            else:
                stranded_start = start
                # stranded_end = end

            tag_str = ",".join(sorted(tags))

            # gene_id = id_map.get(gene_id, 1)
            # gene_name = id_map.get(gene_name, 1)
            # transcript_id = id_map.get(transcript_id, 1)
            # exon_id = id_map.get(exon_id, 1)
            # gene_type = id_map.get(gene_type, 1)

            if gene_type not in gene_types:
                cursor.execute("INSERT INTO gene_types (name) VALUES (?)", (gene_type,))
                gene_types[gene_type] = len(gene_types) + 1
            gene_type_id = gene_types[gene_type]

            if transcript_type not in transcript_types:
                cursor.execute(
                    "INSERT INTO transcript_types (name) VALUES (?)", (transcript_type,)
                )
                transcript_types[transcript_type] = len(transcript_types) + 1

            transcript_type_id = transcript_types[transcript_type]

            if level == "gene":
                cursor.execute(
                    "INSERT INTO genes (gene_id, gene_symbol, chr, start, end, strand, gene_type_id) VALUES (?, ?, ?, ?, ?, ?, ?)",
                    (
                        gene_id,
                        gene_name,
                        chr,
                        start,
                        end,
                        strand,
                        gene_type_id,
                    ),
                )
                # record = cursor.lastrowid
            elif level == "transcript":
                cursor.execute(
                    "INSERT INTO transcripts (gene_id, transcript_id, chr, start, end, strand, is_canonical, transcript_type_id) VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                    (
                        gene_record_id,
                        transcript_id,
                        chr,
                        start,
                        end,
                        strand,
                        is_canonical,
                        transcript_type_id,
                    ),
                )
                # record = cursor.lastrowid
            elif level == "exon":
                cursor.execute(
                    "INSERT INTO exons (transcript_id, exon_id, chr, start, end, strand, exon_number) VALUES (?, ?, ?, ?, ?, ?, ?)",
                    (
                        transcript_record_id,
                        exon_id,
                        chr,
                        start,
                        end,
                        strand,
                        exon_number,
                    ),
                )
                exon_number += 1
                # record = cursor.lastrowid
            else:
                # print(f"Unknown level: {level}")
                pass
            # break

            record += 1

    cursor.execute("END TRANSACTION;")

    cursor.execute("BEGIN TRANSACTION;")

    cursor.execute("CREATE INDEX idx_genes_gene_id ON genes(gene_id);")
    cursor.execute("CREATE INDEX idx_genes_gene_symbol ON genes(gene_symbol);")
    cursor.execute("CREATE INDEX idx_genes_gene_type_id ON genes(gene_type_id);")
    # cursor.execute("CREATE INDEX idx_genes_tss ON genes(tss);")
    cursor.execute(
        "CREATE INDEX idx_genes_chr_start_end_strand ON genes(chr, start, end, strand);"
    )

    cursor.execute(
        "CREATE INDEX idx_transcripts_transcript_id ON transcripts(transcript_id);"
    )
    cursor.execute("CREATE INDEX idx_transcripts_gene_id ON transcripts(gene_id);")
    cursor.execute(
        "CREATE INDEX idx_transcripts_chr_start_end_strand ON transcripts(chr, start, end, strand);"
    )
    cursor.execute(
        "CREATE INDEX idx_transcripts_is_canonical ON transcripts(is_canonical);"
    )
    cursor.execute(
        "CREATE INDEX idx_transcripts_is_longest ON transcripts(is_longest);"
    )

    cursor.execute(
        "CREATE INDEX idx_transcripts_transcript_type_id ON transcripts(transcript_type_id);"
    )

    cursor.execute("CREATE INDEX idx_exons_exon_id ON exons(exon_id);")
    cursor.execute("CREATE INDEX idx_exons_transcript_id ON exons(transcript_id);")
    cursor.execute(
        "CREATE INDEX idx_exons_chr_start_end_strand ON exons(chr, start, end, strand);"
    )

    cursor.execute(
        f"CREATE TABLE info (id INTEGER PRIMARY KEY ASC, genome TEXT NOT NULL, version TEXT NOT NULL);",
    )

    cursor.execute(
        f"INSERT INTO info (genome, version) VALUES('{file_desc[0][0]}', '{file_desc[0][1]}');",
    )

    cursor.execute("END TRANSACTION;")

    # work out who is longest transcript per gene
    print("Finding longest transcripts...")

    cursor.execute(
        f"""
        WITH max_lengths AS (
            SELECT gene_id,
            transcript_id,
            ROW_NUMBER() OVER (PARTITION BY gene_id ORDER BY ABS(end - start) DESC) AS rank
            FROM transcripts
        )
        SELECT gene_id, transcript_id FROM max_lengths
        WHERE rank = 1
        ORDER BY gene_id, transcript_id;
        """,
    )

    rows = cursor.fetchall()

    cursor.execute("BEGIN TRANSACTION;")

    queries = []

    for row in rows:
        # print(row)
        queries.append({"gene_id": row[0], "transcript_id": row[1]})

    cursor.executemany(
        "UPDATE transcripts SET is_longest = 1 WHERE gene_id = :gene_id AND transcript_id = :transcript_id",
        queries,
    )

    cursor.execute("END TRANSACTION;")

    print(len(queries), "transcripts to set as longest")

    # Commit and close
    conn.commit()
    conn.close()
