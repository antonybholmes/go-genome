use std::{
    cmp,
    collections::{BTreeMap, BTreeSet, HashMap},
    error::Error,
};

use dna::Location;
use loctogene::{GenomicFeature, TSSRegion};
use serde::Serialize;

mod tests;

pub const NA: &str = "n/a";
pub const PROMOTER: &str = "promoter";
pub const EXONIC: &str = "exonic";
pub const INTRONIC: &str = "intronic";
pub const INTERGENIC: &str = "intergenic";

//const ERROR_FEATURES:Features= Features{location: dna::EMPTY_STRING, level: dna::EMPTY_STRING, features: [].to_vec()};

#[derive(Serialize, Clone)]
pub struct ClosestGene {
    pub gene_id: String,
    pub gene_symbol: String,
    pub prom_label: String,
    pub tss_dist: i32,
}

#[derive(Serialize)]
pub struct GeneAnnotation {
    pub gene_ids: String,
    pub gene_symbols: String,
    pub prom_labels: String,
    pub tss_dists: String,
    pub closest_genes: Vec<ClosestGene>,
}

struct GeneProm {
    is_promoter: bool,
    is_intronic: bool,
    is_exon: bool,
    abs_d: i32,
    d: i32,
}

pub struct Annotate {
    genesdb: loctogene::Loctogene,
    tss_region: TSSRegion,
    n: u16,
}

 

impl Annotate {
    pub fn new(genesdb: loctogene::Loctogene, tss_region: TSSRegion, n: u16) -> Self {
        return Annotate {
            genesdb,
            tss_region,
            n,
        };
    }

    pub fn annotate(&self, location: &Location) -> Result<GeneAnnotation, Box<dyn Error>> {
        let mid: i32 = location.mid();

        // extend search area to account  for promoter
        // let search_loc: Location = Location::new(
        //     &location.chr,
        //     location.start - self.tss_region.offset_5p.abs(),
        //     location.end + self.tss_region.offset_5p.abs(),
        // )?;

        let genes_within: Vec<loctogene::GenomicFeature> = self.genesdb.get_genes_within_promoter(
            &location,
            loctogene::Level::Transcript,
            cmp::max(self.tss_region.offset_5p(), self.tss_region.offset_3p()),
        )?;

        // we need the unique ids to symbols
        let mut id_map: HashMap<&str, &str> = HashMap::new();
        let mut promoter_map: HashMap<&str, GeneProm> = HashMap::new();
        //let mut dist_map: HashMap<&str, bool> = HashMap::new();

        for gene in genes_within.iter() {
            let id: &str = &gene.gene_id;

            // println!(
            //     "{} {} {} {} {}",
            //     gene.gene_id, gene.gene_symbol, gene.start, gene.end, gene.strand
            // );

            id_map.insert(id, &gene.gene_symbol);

            //let labels = self.classify_location(location, gene);

            let exons: Vec<loctogene::GenomicFeature> = self.genesdb.in_exon(&location, id)?;

            let is_exon: bool = exons.len() > 0;

            let is_promoter: bool = (gene.strand == "+"
                && mid >= gene.start - self.tss_region.offset_5p()
                && mid <= gene.start + self.tss_region.offset_3p())
                || (gene.strand == "-"
                    && mid >= gene.end - self.tss_region.offset_3p()
                    && mid <= gene.end + self.tss_region.offset_5p());

            let is_intronic = mid >= gene.start && mid <= gene.end;

            let d: i32 = if gene.strand == "+" {
                gene.start - mid
            } else {
                gene.end - mid
            };

            //println!("{} {} {}", gene.end - mid, gene.end, mid);

            // update by inserting default case and then updating
            promoter_map
                .entry(id)
                .and_modify(|v| {
                    v.is_intronic = v.is_intronic || is_intronic;
                    v.is_promoter = v.is_promoter || is_promoter;
                    v.is_exon = v.is_exon || exons.len() > 0;

                    let abs_d: i32 = d.abs();

                    if abs_d < v.abs_d {
                        v.d = d;
                        v.abs_d = abs_d;
                    }
                })
                .or_insert(GeneProm {
                    is_promoter,
                    is_intronic,
                    is_exon,
                    d,
                    abs_d: d.abs(),
                });
        }

        // sort the ids by distance
        let mut dist_map: BTreeMap<i32, BTreeSet<&str>> = BTreeMap::new();

        for id in id_map.keys().map(|k| *k) {
            dist_map
                .entry(promoter_map.get(id).unwrap().abs_d)
                .and_modify(|v| {
                    v.insert(id);
                })
                .or_insert({
                    let mut dids: BTreeSet<&str> = BTreeSet::new();
                    dids.insert(id);
                    dids
                });
        }

        // now put the ids in distance order
        let mut ids: Vec<String> = Vec::with_capacity(id_map.len());

        for d in dist_map.keys() {
            for id in dist_map.get(d).unwrap() {
                ids.push(id.to_string());
            }
        }

        println!("within d {}", ids.len());


        // make a list of the symbols in distance order
        let gene_symbols : Vec<String> = ids
            .iter()
            .map(|id:&str| id_map.get(id).unwrap().to_string()) 
            .collect::<Vec<String>>();
       

        let prom_labels: Vec<String> = ids
            .iter()
            .map(|id: &str| {
                let p = promoter_map.get(id).unwrap();
                make_label(p.is_promoter, p.is_exon, p.is_intronic)
            })
            .collect::<Vec<String>>();
        

        let tss_dists: Vec<String> = ids
            .iter()
            .map(|id: &str| promoter_map.get(id).unwrap().d.to_string())
            .collect::<Vec<String>>();

        if ids.len() == 0 {
            ids.push(NA.to_owned());
            gene_symbols.push(NA.to_owned());
            tss_dists.push(NA.to_owned());
        }
        
        println!("{} geneids", ids.join(";"));
        println!("{}", gene_symbols);
        println!("{}", prom_labels);
        println!("{}", tss_dists);

        let closest_genes: Vec<loctogene::GenomicFeature> =
            self.genesdb
                .get_closest_genes(&location, self.n, loctogene::Level::Gene)?;

        let annotation: GeneAnnotation = GeneAnnotation {
            gene_ids: ids.join(";"),
            gene_symbols,
            prom_labels,
            tss_dists,
            closest_genes: closest_genes
                .iter()
                .map(|cg| ClosestGene {
                    gene_id: cg.gene_id.to_owned(),
                    gene_symbol: cg.gene_symbol.to_owned(),
                    tss_dist: cg.dist,
                    prom_label: self.classify_location(location, cg),
                })
                .collect(),
        };

        Ok(annotation)
    }

    fn classify_location(&self, location: &Location, feature: &GenomicFeature) -> String {
        let mid = location.mid();

        let s = if feature.strand == "+" {
            feature.start - self.tss_region.offset_5p()
        } else {
            feature.start
        };

        let e = if feature.strand == "-" {
            feature.end + self.tss_region.offset_5p()
        } else {
            feature.end
        };

        if location.start > e || location.end < s {
            return INTERGENIC.to_string();
        }

        let is_promoter: bool = (feature.strand == "+"
            && mid >= s
            && mid <= feature.start + self.tss_region.offset_3p())
            || (feature.strand == "-"
                && mid >= feature.end - self.tss_region.offset_3p()
                && mid <= e);

        let exons: Vec<loctogene::GenomicFeature> =
            match self.genesdb.in_exon(&location, &feature.gene_id) {
                Ok(exons) => exons,
                Err(_) => vec![],
            };

        let is_exon = exons.len() > 0;

        let is_intronic = mid >= feature.start && mid <= feature.end;

        return make_label(is_promoter, is_exon, is_intronic);
    }
}

fn make_label(is_promoter: bool, is_exon: bool, is_intronic: bool) -> String {
    let mut labels: Vec<&str> = Vec::with_capacity(2);

    if is_promoter {
        labels.push(PROMOTER);
    }

    if is_exon {
        labels.push(EXONIC);
    } else {
        if is_intronic {
            labels.push(INTRONIC);
        }
    }

    return labels.join(",");
}
