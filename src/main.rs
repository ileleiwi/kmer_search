use clap::{Arg, Command};
use std::collections::HashMap;
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

    // Process files and write output
    let mer_hash_vec = file_utils::process_files_with_extensions(input_dir);
    let column_names = file_utils::collect_headers(input_dir);
    match file_utils::write_tsv(&mer_hash_vec, &full_out_path, &column_names) {
        Ok(()) => println!("TSV file successfully written"),
        Err(err) => eprintln!("Error writing TSV file: {}", err),
    }
}

fn collapse_hashes(mer_hash_vec: &Vec<HashMap<String, usize>>) -> HashMap<String, usize> {
    let mut collapsed_hash: HashMap<String, usize> = HashMap::new();
    for hash in mer_hash_vec {
        for (id, count) in hash {
            let mer = collapsed_hash.entry(id.to_string()).or_insert(0);
            *mer += count;
        }
    }
    collapsed_hash
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_collect_headers() {
        let test_dir = "/Users/leleiwi1/Desktop/LLNL_postdoc/kmer_tool_16S/test_input";
        let file_names = file_utils::collect_headers(&test_dir.to_string());
        assert_eq!(file_names.len(), 4);
        assert_eq!(file_names[0], "t.fasta");
        assert_eq!(file_names[1], "t.fna");
        assert_eq!(file_names[2], "test1.fa");
        assert_eq!(file_names[3], "test2.fa");
    }
}
