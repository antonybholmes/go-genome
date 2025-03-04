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


tables = ["genes"]

level_map = {"gene": 1, "transcript": 2, "exon": 3}

files = [
    [
        "hg19.sql",
        ["hg19", "Gencode v38"],
        "/ifs/scratch/cancer/Lab_RDF/ngs/references/gencode/grch37/gencode.v38lift37.annotation.gtf.gz",
    ],
    [
        "grch38.sql",
        ["grch38", "Gencode v44"],
        "/ifs/scratch/cancer/Lab_RDF/ngs/references/gencode/grch38/gencode.v44.annotation.gtf.gz",
    ],
    [
        "mm10.sql",
        ["mm10", "Gencode vM25"],
        "/ifs/scratch/cancer/Lab_RDF/ngs/references/gencode/mm10/gencode.vM25.annotation.gtf.gz",
    ],
]

table = "genes"

for f in files:
    print(f)
    with open(
        f"data/modules/genome/{f[0]}",
        "w",
    ) as out:
        print(
            f"CREATE TABLE info (id INTEGER PRIMARY KEY ASC, genome TEXT NOT NULL, version TEXT NOT NULL);",
            file=out,
        )

        print(
            f"CREATE TABLE {table} (id INTEGER PRIMARY KEY ASC, parent_id INTEGER NOT NULL, level INTEGER NOT NULL, chr TEXT NOT NULL, start INTEGER NOT NULL, end INTEGER NOT NULL, tss INTEGER NOT NULL, strand TEXT NOT NULL, gene_id TEXT NOT NULL DEFAULT '', gene_symbol TEXT NOT NULL DEFAULT '', transcript_id TEXT NOT NULL DEFAULT '', exon_id TEXT NOT NULL DEFAULT '', is_canonical BOOLEAN NOT NULL DEFAULT 0);",
            file=out,
        )
        # print(f"CREATE INDEX {table}_level ON {table} (level);", file=out)
        # print(f"CREATE INDEX {table}_chr ON {table} (chr);", file=out)
        # print(f"CREATE INDEX {table}_start ON {table} (start);", file=out)
        # print(f"CREATE INDEX {table}_end ON {table} (end);", file=out)

        print(
            f"CREATE INDEX {table}_level_chr_start_end_strand_idx ON {table} (level, chr, start, end, strand);",
            file=out,
        )

        print(
            f"CREATE INDEX {table}_transcripts_idx ON {table} (level, gene_id, is_canonical);",
            file=out,
        )
        print(
            f"CREATE INDEX {table}_exons_idx ON {table} (level, transcript_id, is_canonical);",
            file=out,
        )

        print(f"CREATE INDEX {table}_gene_id_idx ON {table} (gene_id);", file=out)

        print(
            f"CREATE INDEX {table}_gene_symbol_idx ON {table} (gene_symbol);",
            file=out,
        )

        # print(
        #     f"CREATE INDEX {table}_transcript_id_idx ON {table} (transcript_id);",
        #     file=out,
        # )

        # print(
        #     f"CREATE INDEX {table}_exon_id_idx ON {table} (exon_id);",
        #     file=out,
        # )

        # print(
        #    f"CREATE INDEX {table}_level_chr_stranded_start_stranded_end ON {table} (level, chr, stranded_start, stranded_end);",
        #    file=out,
        # )

        print(
            f"INSERT INTO info (genome, version) VALUES('{f[1][0]}', '{f[1][1]}');",
            file=out,
        )

        print()

        levels = {"gene", "transcript", "exon"}

        print("BEGIN TRANSACTION;", file=out)

        record = 1
        gene_record_id = -1
        transcript_record_id = -1
        exon_record_id = -1
        tags = set()
        is_canonical = 0

        with gzip.open(
            f[2],
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

                if level not in levels:
                    continue

                if level == "gene":
                    parent_record_id = -1
                    gene_record_id = record

                if level == "transcript":
                    parent_record_id = gene_record_id
                    transcript_record_id = record
                    tags = set()
                    is_canonical = 0

                if level == "exon":
                    parent_record_id = transcript_record_id
                    exon_record_id = record

                # gene
                matcher = re.search(r'gene_id "(.+?)";', tokens[8])

                if matcher:
                    # remove version
                    gene_id = re.sub(r"\..+", "", matcher.group(1))
                else:
                    gene_id = ""

                matcher = re.search(r'gene_name "(.+?)";', tokens[8])

                if matcher:
                    gene_name = matcher.group(1)
                else:
                    gene_name = ""

                # transcript
                matcher = re.search(r'transcript_id "(.+?)";', tokens[8])

                if matcher:
                    transcript_id = re.sub(r"\..+", "", matcher.group(1))
                else:
                    transcript_id = ""

                if "Ensembl_canonical" in line:
                    tags.add("canonical")
                    is_canonical = 1

                # exon
                matcher = re.search(r'exon_id "(.+?)";', tokens[8])

                if matcher:
                    exon_id = re.sub(r"\..+", "", matcher.group(1))
                else:
                    exon_id = ""

                chr = tokens[0]
                start = int(tokens[3])
                end = int(tokens[4])
                # mid = int((end + start) / 2)
                strand = tokens[6]

                # invert coordinates to make searching easier
                if strand == "-":
                    stranded_start = end
                    stranded_end = start
                else:
                    stranded_start = start
                    stranded_end = end

                tag_str = ",".join(sorted(tags))

                print(
                    f"INSERT INTO genes (parent_id, level, chr, start, end, tss, strand, gene_id, gene_symbol, transcript_id, exon_id, is_canonical) VALUES ({parent_record_id}, {level_map[level]}, '{chr}', {start}, {end}, {stranded_start}, '{strand}', '{gene_id}', '{gene_name}', '{transcript_id}', '{exon_id}', {is_canonical});",
                    file=out,
                )

                # break

                record += 1

        print("COMMIT;", file=out)
