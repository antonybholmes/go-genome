CREATE TABLE genes (id INTEGER PRIMARY KEY ASC, parent_id INTEGER NOT NULL, level INTEGER NOT NULL, chr TEXT NOT NULL, start INTEGER NOT NULL, end INTEGER NOT NULL, tss INTEGER NOT NULL, strand TEXT NOT NULL, gene_id TEXT, gene_symbol TEXT, transcript_id TEXT, exon_id TEXT);
CREATE INDEX genes_gene_id_idx ON genes (gene_id);
CREATE INDEX genes_gene_symbol_idx ON genes (gene_symbol);
CREATE INDEX genes_level_chr_start_end_strand_idx ON genes (level, chr, start, end, strand);
CREATE INDEX genes_transcript_id_idx ON genes (transcript_id);
CREATE INDEX genes_exon_id_idx ON genes (exon_id);
