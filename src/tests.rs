#[cfg(test)]
use dna::Location;
#[cfg(test)]
use loctogene::{Loctogene,TSSRegion};
#[cfg(test)]
use std::error::Error;

#[cfg(test)]
use serde_json::json;

#[cfg(test)]
use crate::{Annotate, GeneAnnotation};


#[cfg(test)]

#[test]
fn test_annotation() ->Result<(), Box<dyn Error>>{
    

    //let loc: Location = Location::parse("chr3:187721370-187733550")?;


    let loc: Location = Location::parse("chr3:187745448-187745468")?;

    let genesdb: Loctogene = Loctogene::new("../docker-rust-edb-api/data/loctogene/grch38.db")?;


    let annotatedb: Annotate = Annotate::new(genesdb, TSSRegion::default(), 10);

    let annotation: GeneAnnotation = annotatedb.annotate(&loc)?;

    //let js: serde_json::Value = json!(records);

    //println!("{}", js);

    Ok(())
}

