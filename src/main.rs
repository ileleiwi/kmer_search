use clap::{Arg, Command};
use std::fs;
use std::path::PathBuf;
pub mod file_utils;
pub mod sequence_utils;

const KMER_SIZES: [usize; 5] = [1, 2, 3, 4, 5];
const RMER_SIZES: [usize; 5] = [6, 7, 8, 9, 10];
const EXTENSIONS: [&str; 3] = ["fa", "fasta", "fna"];

fn main() {
    // Get command-line arguments
    let matches: clap::ArgMatches = Command::new("kmer_search")
        .arg(
            Arg::new("input_dir")
                .short('i')
                .long("input_dir")
                .value_name("DIRECTORY")
                .help("Path to directory of input FASTA files")
                .required(true),
        )
        .arg(
            Arg::new("output_dir")
                .short('o')
                .long("output_dir")
                .value_name("DIRECTORY")
                .help("Output directory")
                .required(false),
        )
        .get_matches();

    // Get input directory
    let input_dir = matches.get_one::<String>("input_dir").unwrap();

    // Use default output file names if not provided
    let default_output_dir = String::from("./kmer_search_output");
    let output_dir = matches
        .get_one::<String>("output_dir")
        .map_or_else(|| PathBuf::from(&default_output_dir), PathBuf::from);

    // Create output directory if it doesn't exist
    match fs::create_dir_all(&output_dir) {
        Ok(_) => println!("Output directory created"),
        Err(e) => {
            println!("Error creating output directory: {}", e);
        }
    }

    let full_out_path = output_dir.join("mer_counts_output.tsv");
    let full_out_path_stats = output_dir.join("sequence_stats.tsv");

    // Process files and write output
    let fasta_objects = file_utils::instantiate_fastaseq_objects(input_dir);
    let mer_hash_vec = file_utils::process_files_with_extensions(&fasta_objects);
    let column_names = file_utils::collect_headers(&fasta_objects);

    match file_utils::write_stats_tsv(&fasta_objects, &full_out_path_stats) {
        Ok(()) => println!("Stats file successfully written"),
        Err(err) => eprintln!("Error writing Counts file: {}", err),
    }
    match file_utils::write_mer_tsv(&mer_hash_vec, &full_out_path, &column_names) {
        Ok(()) => println!("Counts file successfully written"),
        Err(err) => eprintln!("Error writing Counts file: {}", err),
    }
}
