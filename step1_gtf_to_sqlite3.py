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

GENES_SQL = """CREATE TABLE genes 
    (id INTEGER PRIMARY KEY ASC, 
    parent_id INTEGER NOT NULL, 
    level INTEGER NOT NULL, 
    chr TEXT NOT NULL, 
    start INTEGER NOT NULL, 
    end INTEGER NOT NULL, 
    tss INTEGER NOT NULL, 
    strand TEXT NOT NULL, 
    gene_id TEXT NOT NULL DEFAULT "", 
    gene_symbol TEXT NOT NULL DEFAULT "", 
    transcript_id TEXT NOT NULL DEFAULT "", 
    exon_id TEXT NOT NULL DEFAULT "", 
    is_canonical INTEGER NOT NULL DEFAULT 0, 
    gene_type TEXT NOT NULL DEFAULT "");"""

for file_desc in files:
    print(file_desc)
    with open(
        f"data/modules/genome/{file_desc[0]}",
        "w",
    ) as out:
        print(
            f"CREATE TABLE info (id INTEGER PRIMARY KEY ASC, genome TEXT NOT NULL, version TEXT NOT NULL);",
            file=out,
        )

        # creat the gene id table
        # print(
        #     f"CREATE TABLE ids (id INTEGER PRIMARY KEY ASC, name TEXT NOT NULL);",
        #     file=out,
        # )
        # print(
        #     f"CREATE INDEX ids_name_idx ON ids (name);",
        #     file=out,
        # )

        # create the transcript type table
        # print(
        #     f"CREATE TABLE transcript_types (id INTEGER PRIMARY KEY ASC, name TEXT NOT NULL);",
        #     file=out,
        # )

        # create the genes table
        print(
            GENES_SQL,
            file=out,
        )
        # print(f"CREATE INDEX genes_level ON genes (level);", file=out)
        # print(f"CREATE INDEX genes_chr ON genes (chr);", file=out)
        # print(f"CREATE INDEX genes_start ON genes (start);", file=out)
        # print(f"CREATE INDEX genes_end ON genes (end);", file=out)

        print(
            f"CREATE INDEX genes_level_chr_start_end_strand_idx ON genes (level, chr, start, end, strand);",
            file=out,
        )

        print(
            f"CREATE INDEX genes_transcripts_idx ON genes (level, gene_id, is_canonical);",
            file=out,
        )
        print(
            f"CREATE INDEX genes_exons_idx ON genes (level, transcript_id, is_canonical);",
            file=out,
        )

        print(
            f"CREATE INDEX genes_gene_types_idx ON genes (level, gene_id, gene_type);",
            file=out,
        )
        print(
            f"CREATE INDEX genes_transcript_types_idx ON genes (level, transcript_id, gene_type);",
            file=out,
        )
        print(
            f"CREATE INDEX genes_exons_types_idx ON genes (level, exon_id, gene_type);",
            file=out,
        )

        print(f"CREATE INDEX genes_gene_id_idx ON genes (gene_id);", file=out)

        print(
            f"CREATE INDEX genes_gene_symbol_idx ON genes (gene_symbol);",
            file=out,
        )

        # print(
        #     f"CREATE INDEX genes_transcript_id_idx ON genes (transcript_id);",
        #     file=out,
        # )

        # print(
        #     f"CREATE INDEX genes_exon_id_idx ON genes (exon_id);",
        #     file=out,
        # )

        # print(
        #    f"CREATE INDEX genes_level_chr_stranded_start_stranded_end ON genes (level, chr, stranded_start, stranded_end);",
        #    file=out,
        # )

        print(
            f"INSERT INTO info (genome, version) VALUES('{file_desc[1][0]}', '{file_desc[1][1]}');",
            file=out,
        )

        print()

        levels = {"gene", "transcript", "exon"}

        record = 1
        gene_record_id = -1
        transcript_record_id = -1
        exon_record_id = -1
        tags = set()
        is_canonical = 0
        id_map = {"n/a": 1}
        # transcript_types = {}

        # read the GTF file for building id maps

        with gzip.open(
            file_desc[2],
            "rt",
        ) as f:
            for line in f:
                if line.startswith("#"):
                    continue

                tokens = line.strip().split("\t")

                # gene
                matcher = re.search(r'gene_id "(.+?)";', tokens[8])

                if matcher:
                    # remove version
                    gene_id = re.sub(r"\..+", "", matcher.group(1))

                    if gene_id not in id_map:
                        id_map[gene_id] = len(id_map) + 1

                matcher = re.search(r'gene_name "(.+?)";', tokens[8])

                if matcher:
                    gene_name = matcher.group(1)

                    if gene_name not in id_map:
                        id_map[gene_name] = len(id_map) + 1

                matcher = re.search(r'transcript_id "(.+?)";', tokens[8])

                if matcher:
                    transcript_id = matcher.group(1)

                    if transcript_id not in id_map:
                        id_map[transcript_id] = len(id_map) + 1

                matcher = re.search(r'exon_id "(.+?)";', tokens[8])

                if matcher:
                    exon_id = matcher.group(1)

                    if exon_id not in id_map:
                        id_map[exon_id] = len(id_map) + 1

                # transcript_type
                matcher = re.search(r'gene_type "(.+?)";', tokens[8])

                if matcher:
                    gene_type = re.sub(r"\..+", "", matcher.group(1))

                    if gene_type not in id_map:
                        id_map[gene_type] = len(id_map) + 1

        # sort the id_map by id
        ids = list(sorted(id_map.items(), key=lambda item: item[1]))

        # print("BEGIN TRANSACTION;", file=out)

        # for name, id in ids:
        #     print(
        #         f"INSERT INTO ids (name) VALUES ('{name}');",
        #         file=out,
        #     )

        # print("COMMIT;", file=out)
        # print()

        # parse file again to generate SQL insert statements

        print("BEGIN TRANSACTION;", file=out)

        with gzip.open(
            file_desc[2],
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
                    gene_record_id = record

                if level == "transcript":
                    parent_record_id = gene_record_id
                    transcript_record_id = record

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
                if strand == "-":
                    stranded_start = end
                    stranded_end = start
                else:
                    stranded_start = start
                    stranded_end = end

                tag_str = ",".join(sorted(tags))

                # gene_id = id_map.get(gene_id, 1)
                # gene_name = id_map.get(gene_name, 1)
                # transcript_id = id_map.get(transcript_id, 1)
                # exon_id = id_map.get(exon_id, 1)
                # gene_type = id_map.get(gene_type, 1)

                print(
                    f"INSERT INTO genes (parent_id, level, chr, start, end, tss, strand, gene_id, gene_symbol, transcript_id, exon_id, is_canonical, gene_type) VALUES ({parent_record_id}, {level_map[level]}, '{chr}', {start}, {end}, {stranded_start}, '{strand}', {gene_id}, {gene_name}, {transcript_id}, {exon_id}, {is_canonical}, {gene_type});",
                    file=out,
                )

                # break

                record += 1

        print("COMMIT;", file=out)
