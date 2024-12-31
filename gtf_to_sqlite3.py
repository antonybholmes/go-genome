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
        "/ifs/scratch/cancer/Lab_RDF/ngs/references/gencode/grch37/gencode.v38lift37.annotation.gtf.gz",
    ],
    [
        "grch38.sql",
        "/ifs/scratch/cancer/Lab_RDF/ngs/references/gencode/grch38/gencode.v44.annotation.gtf.gz",
    ],
    [
        "mm10.sql",
        "/ifs/scratch/cancer/Lab_RDF/ngs/references/gencode/mm10/gencode.vM25.annotation.gtf.gz",
    ],
]


for f in files:
    print(f)
    with open(
        f"data/modules/genes/{f[0]}",
        "w",
    ) as out:
        for table in tables:
            print(
                f"CREATE TABLE {table} (id INTEGER PRIMARY KEY ASC, parent_id INTEGER NOT NULL, level INTEGER NOT NULL, chr TEXT NOT NULL, start INTEGER NOT NULL, end INTEGER NOT NULL, tss INTEGER NOT NULL, strand TEXT NOT NULL, gene_id TEXT NOT NULL DEFAULT '', gene_symbol TEXT NOT NULL DEFAULT '', transcript_id TEXT NOT NULL DEFAULT '', exon_id TEXT NOT NULL DEFAULT '', tags NOT NULL DEFAULT '');",
                file=out,
            )
            # print(f"CREATE INDEX {table}_level ON {table} (level);", file=out)
            # print(f"CREATE INDEX {table}_chr ON {table} (chr);", file=out)
            # print(f"CREATE INDEX {table}_start ON {table} (start);", file=out)
            # print(f"CREATE INDEX {table}_end ON {table} (end);", file=out)
            print(f"CREATE INDEX {table}_tags_idx ON {table} (tags);", file=out)
            print(f"CREATE INDEX {table}_gene_id_idx ON {table} (gene_id);", file=out)
            print(
                f"CREATE INDEX {table}_gene_symbol_idx ON {table} (gene_symbol);",
                file=out,
            )
            print(
                f"CREATE INDEX {table}_level_chr_start_end_strand_idx ON {table} (level, chr, start, end, strand);",
                file=out,
            )
            print(
                f"CREATE INDEX {table}_transcript_id_idx ON {table} (transcript_id);",
                file=out,
            )
            print(
                f"CREATE INDEX {table}_exon_id_idx ON {table} (exon_id);",
                file=out,
            )

            # print(
            #    f"CREATE INDEX {table}_level_chr_stranded_start_stranded_end ON {table} (level, chr, stranded_start, stranded_end);",
            #    file=out,
            # )

            print()

        levels = {"gene", "transcript", "exon"}

        print("BEGIN TRANSACTION;", file=out)

        record = 1
        gene_record_id = -1
        transcript_record_id = -1
        exon_record_id = -1
        tags = set()
        
        with gzip.open(
            f[1],
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

                if level == "transcript":
                    parent_record_id = gene_record_id
                    transcript_record_id = record

                    matcher = re.search(r'transcript_id "(.+?)";', tokens[8])

                    if matcher:
                        transcript_id = re.sub(r"\..+", "", matcher.group(1))
                    else:
                        transcript_id = ""

                    tags = set()

                    if "Ensembl_canonical" in line:
                        tags.add("canonical")

                if level == "exon":
                    parent_record_id = transcript_record_id
                    exon_record_id = record

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

                tags = ",".join(sorted(tags))

                print(
                    f"INSERT INTO genes (parent_id, level, chr, start, end, tss, strand, gene_id, gene_symbol, transcript_id, exon_id) VALUES ({parent_record_id}, {level_map[level]}, '{chr}', {start}, {end}, {stranded_start}, '{strand}', '{gene_id}', '{gene_name}', '{transcript_id}', '{exon_id}', '{tags}');",
                    file=out,
                )

                # break

                record += 1

        print("COMMIT;", file=out)
