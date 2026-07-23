//! Real-data fastSA lookup driver used by benchmark_duckvep_fastvep.Rmd.
//!
//! Copy this file into the pinned FastVEP checkout as
//! crates/fastvep-cli/examples/real_lookup.rs, then run:
//!
//! cargo run --release -p fastvep-cli --example real_lookup -- \
//!   CLINVAR.osa GIAB_ALLELES.tsv THREADS REPEATS

use anyhow::{Context, Result};
use fastvep_cache::annotation::{AnnotationProvider, AnnotationValue};
use fastvep_sa::reader::SaReader;
use rayon::prelude::*;
use rayon::ThreadPoolBuilder;
use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::time::Instant;

const BATCH_SIZE: usize = 1024;

#[derive(Debug)]
struct Query {
    chrom: String,
    pos: u64,
    reference: String,
    alternate: String,
}

fn read_queries(path: &Path) -> Result<Vec<Query>> {
    let file = File::open(path)
        .with_context(|| format!("opening query relation {}", path.display()))?;
    let mut queries = Vec::new();
    for (line_number, line) in BufReader::new(file).lines().enumerate() {
        let line = line.with_context(|| {
            format!("reading line {} from {}", line_number + 1, path.display())
        })?;
        let mut fields = line.split('\t');
        let chrom = fields.next();
        let pos = fields.next();
        let reference = fields.next();
        let alternate = fields.next();
        if chrom.is_none()
            || pos.is_none()
            || reference.is_none()
            || alternate.is_none()
            || fields.next().is_some()
        {
            anyhow::bail!(
                "{}:{}: expected exactly four tab-separated fields",
                path.display(),
                line_number + 1
            );
        }
        queries.push(Query {
            chrom: chrom.unwrap().to_owned(),
            pos: pos
                .unwrap()
                .parse()
                .with_context(|| format!("parsing position on line {}", line_number + 1))?,
            reference: reference.unwrap().to_owned(),
            alternate: alternate.unwrap().to_owned(),
        });
    }
    Ok(queries)
}

fn annotation_bytes(value: AnnotationValue) -> usize {
    match value {
        AnnotationValue::Json(value) | AnnotationValue::Positional(value) => value.len(),
        AnnotationValue::Interval(values) => values.iter().map(String::len).sum(),
    }
}

fn run_once(reader: &SaReader, queries: &[Query]) -> Result<(u64, u64)> {
    let mut hits = 0u64;
    let mut payload_bytes = 0u64;
    for batch in queries.chunks(BATCH_SIZE) {
        let mut positions: HashMap<&str, Vec<u64>> = HashMap::new();
        for query in batch {
            positions
                .entry(query.chrom.as_str())
                .or_default()
                .push(query.pos);
        }
        for (chrom, positions) in positions {
            reader.preload(chrom, &positions)?;
        }

        let batch_result = batch
            .par_iter()
            .try_fold(
                || (0u64, 0u64),
                |(hits, bytes), query| {
                    let value = reader.annotate_position(
                        &query.chrom,
                        query.pos,
                        &query.reference,
                        &query.alternate,
                    )?;
                    Ok::<_, anyhow::Error>(match value {
                        Some(value) => (
                            hits + 1,
                            bytes + u64::try_from(annotation_bytes(value)).unwrap_or(u64::MAX),
                        ),
                        None => (hits, bytes),
                    })
                },
            )
            .try_reduce(
                || (0u64, 0u64),
                |left, right| {
                    Ok::<_, anyhow::Error>((
                        left.0.saturating_add(right.0),
                        left.1.saturating_add(right.1),
                    ))
                },
            )?;
        hits = hits.saturating_add(batch_result.0);
        payload_bytes = payload_bytes.saturating_add(batch_result.1);
    }
    Ok((hits, payload_bytes))
}

fn main() -> Result<()> {
    let args: Vec<String> = env::args().collect();
    if args.len() != 5 {
        anyhow::bail!(
            "usage: {} CLINVAR.osa GIAB_ALLELES.tsv THREADS REPEATS",
            args.first().map(String::as_str).unwrap_or("real_lookup")
        );
    }
    let threads: usize = args[3].parse().context("parsing THREADS")?;
    let repeats: usize = args[4].parse().context("parsing REPEATS")?;
    if threads == 0 || repeats == 0 {
        anyhow::bail!("THREADS and REPEATS must be positive");
    }

    let queries = read_queries(Path::new(&args[2]))?;
    let reader = SaReader::open(Path::new(&args[1]))?;
    let pool = ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .context("building Rayon thread pool")?;

    println!("run,threads,queries,hits,payload_bytes,seconds,queries_per_second");
    for run in 1..=repeats {
        let started = Instant::now();
        let (hits, payload_bytes) =
            pool.install(|| run_once(&reader, &queries))?;
        let seconds = started.elapsed().as_secs_f64();
        println!(
            "{run},{threads},{},{hits},{payload_bytes},{seconds:.9},{:.3}",
            queries.len(),
            queries.len() as f64 / seconds
        );
    }
    Ok(())
}
