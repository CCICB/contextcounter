use anyhow::Context;
use clap::Parser;
use contextcounter::counts::{CountsDi, CountsPenta, CountsTri};
use fern::colors::ColoredLevelConfig;
use log::info;
use noodles::fasta;
use std::{
    collections::HashSet,
    fs::{self, File},
    io::BufReader,
    path::{Path, PathBuf},
    time::SystemTime,
};

#[derive(Parser, Debug)]
#[command(
    author,
    version = "0.0.1",
    about = "Count frequency of di/tri/penta nucleotide sequences in a fasta file"
)]
struct Cli {
    /// Path to the input FASTA file
    fasta: PathBuf,

    /// Folder to write count files
    #[arg(short, long, value_name = "OUTFOLDER", default_value = "contexts")]
    outdir: PathBuf,

    /// Folder to write count files
    #[arg(short, long, default_value_t = false)]
    print_counts: bool,
}

fn setup_logger() -> Result<(), fern::InitError> {
    let colors = ColoredLevelConfig::new().info(fern::colors::Color::Green);

    fern::Dispatch::new()
        .format(move |out, message, record| {
            out.finish(format_args!(
                "[{} {} {}] {}",
                humantime::format_rfc3339_seconds(SystemTime::now()),
                colors.color(record.level()),
                record.target(),
                message
            ))
        })
        .level(log::LevelFilter::Info)
        .chain(std::io::stderr())
        // .chain(fern::log_file(logfile)?)
        .apply()?;
    Ok(())
}

pub fn log_error_and_panic(err: anyhow::Error) {
    let causes = err
        .chain()
        .skip(1)
        .map(|cause| format!("Caused by: {}", cause))
        .collect::<Vec<_>>()
        .join(" <= ");

    log::error!("Fatal Error: {}. {}", err, causes);
    std::process::exit(1)
}

fn main() {
    // Set up logging
    setup_logger().expect("Error setting up logging");

    // Log any bubbled up errors and panic
    if let Err(err) = run() {
        log_error_and_panic(err);
        std::process::exit(1)
    }
}

fn run() -> Result<(), anyhow::Error> {
    // let fasta = std::path::PathBuf::from("testfiles/test.fasta");

    // Parce CLI arguments
    let cli = Cli::parse();
    let fasta = cli.fasta;
    let outdir = cli.outdir;
    let print_counts = cli.print_counts;
    let skip: HashSet<String> = HashSet::new();

    // Create output directory if it doesn't exist
    fs::create_dir_all(&outdir)
        .with_context(|| format!("Failed to create output directory: {}", outdir.display()))?;

    let stem = fasta
        .file_stem()
        .and_then(|s| s.to_str())
        .ok_or_else(|| anyhow::anyhow!("Invalid fasta file stem"))?;

    let prefix = outdir.join(stem);

    let trinucleotides = count_trinucleotides(&fasta, &skip, print_counts)?;
    let pentanucleotides = count_pentanucleotides(&fasta, &skip, print_counts)?;
    let dinucleotides = count_dinucleotides(&fasta, &skip, print_counts)?;

    // Write output
    info!("Writing files to: {}", outdir.canonicalize()?.display());
    let _ = write_context_file("trinucleotide", &prefix, trinucleotides.to_string());
    let _ = write_context_file("dinucleotide", &prefix, dinucleotides.to_string());
    let _ = write_context_file("pentanucleotide", &prefix, pentanucleotides.to_string());
    Ok(())
}

fn write_context_file(
    context_type: &str,
    prefix: &Path,
    content: String,
) -> Result<(), anyhow::Error> {
    let filename = prefix.with_file_name(format!(
        "{}_{}.tsv",
        prefix.file_name().unwrap().to_string_lossy(),
        context_type
    ));
    Ok(fs::write(filename, content)?)
}

/// Count trinucleotide contexts (pyrimidine centered)
fn count_trinucleotides(
    fasta: &PathBuf,
    skip: &HashSet<String>,
    print_counts: bool,
) -> Result<CountsTri, anyhow::Error> {
    // Configure windowsize
    let windowsize: usize = 3;

    // Initialise counts of each trinucleotide to zero
    let mut tnc_counts = contextcounter::counts::CountsTri::default();

    info!("Contig: [{}]", fasta.display());

    // Create FASTA file reader
    let mut reader = File::open(fasta)
        .map(BufReader::new)
        .map(fasta::io::Reader::new)
        .context("Failed to read TNC counts")?;

    // Read each record (contains references to sequence names & info)
    for result in reader.records() {
        let record = result?;
        let contig_name = record.definition().name();

        // Check if contig should be skipped (often used to exclude sex chromosomes from counts
        if !skip.is_empty() && skip.contains(contig_name) {
            info!("Contig: {} (skipped: in blacklist)", contig_name);
            continue;
        }

        // Otherwise, proceed with TNC counting
        info!("Contig: {}", contig_name);

        let seq_bytes = record.sequence().as_ref();

        // Skip sequences shorter than the window size
        if seq_bytes.len() < windowsize {
            continue;
        }

        // Sliding window of size 3
        for window in seq_bytes.windows(windowsize) {
            // Fetch trinucleotide sequence
            let tri = String::from_utf8(Vec::from(window)).unwrap();
            tnc_counts.increment(&tri);
        }
    }

    // Display count matrix
    if print_counts {
        println!("{}", tnc_counts)
    };

    Ok(tnc_counts)
}

/// Count pentanucleotide contexts (pyrimidine centered)
fn count_pentanucleotides(
    fasta: &PathBuf,
    skip: &HashSet<String>,
    print_counts: bool,
) -> Result<CountsPenta, anyhow::Error> {
    // Configure windowsize
    let windowsize: usize = 5;

    // Initialise counts of each trinucleotide to zero
    let mut penta_counts = contextcounter::counts::CountsPenta::default();

    info!("Counting pentanucleotide contexts in [{}]", fasta.display());

    // Create FASTA file reader
    let mut reader = File::open(fasta)
        .map(BufReader::new)
        .map(fasta::io::Reader::new)
        .context("Failed to read pentanucleotide counts")?;

    // Read each record (contains references to sequence names & info)
    for result in reader.records() {
        let record = result?;
        let contig_name = record.definition().name();

        // Check if contig should be skipped (often used to exclude sex chromosomes from counts
        if !skip.is_empty() && skip.contains(contig_name) {
            info!("Contig: {} (skipped: in blacklist)", contig_name);
            continue;
        }

        // Otherwise, proceed with TNC counting
        info!("Contig: {}", contig_name);

        let seq_bytes = record.sequence().as_ref();

        // Skip sequences shorter than the window size
        if seq_bytes.len() < windowsize {
            continue;
        }

        // Sliding window along fasta entry
        for window in seq_bytes.windows(windowsize) {
            // Fetch trinucleotide sequence
            let penta = String::from_utf8(Vec::from(window)).unwrap();
            penta_counts.increment(&penta);
        }
    }

    // Display count matrix
    if print_counts {
        println!("{}", penta_counts)
    };

    Ok(penta_counts)
}

fn count_dinucleotides(
    fasta: &PathBuf,
    skip: &HashSet<String>,
    print_counts: bool,
) -> Result<CountsDi, anyhow::Error> {
    // Configure windowsize
    let windowsize: usize = 2;

    // Initialise counts of each trinucleotide to zero
    let mut dinucleotide_counts = contextcounter::counts::CountsDi::default();

    info!("Counting dinucleotide contexts in [{}]", fasta.display());

    // Create FASTA file reader
    let mut reader = File::open(fasta)
        .map(BufReader::new)
        .map(fasta::io::Reader::new)
        .context("Failed to read dinucleotide counts")?;

    // Read each record (contains references to sequence names & info)
    for result in reader.records() {
        let record = result?;
        let contig_name = record.definition().name();

        // Check if contig should be skipped (often used to exclude sex chromosomes from counts
        if !skip.is_empty() && skip.contains(contig_name) {
            info!("Contig: {} (skipped: in blacklist)", contig_name);
            continue;
        }

        // Otherwise, proceed with TNC counting
        info!("Contig: {}", contig_name);

        let seq_bytes = record.sequence().as_ref();

        // Skip sequences shorter than the window size
        if seq_bytes.len() < windowsize {
            continue;
        }

        // Sliding window along fasta entry
        for window in seq_bytes.windows(windowsize) {
            // Fetch trinucleotide sequence
            let dinucleotide = String::from_utf8(Vec::from(window)).unwrap();
            dinucleotide_counts.increment(&dinucleotide);
        }
    }

    // Display count matrix
    if print_counts {
        println!("{}", dinucleotide_counts);
    }

    Ok(dinucleotide_counts)
}
