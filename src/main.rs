use clap::{Arg, Command};
use std::collections::HashMap;
use std::collections::HashSet;
use std::fs;
use std::fs::File;
use std::hash;
use std::io::Write;
use std::path::PathBuf;
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
    let mer_hash_vec = process_files_with_extensions(input_dir);
    let column_names = collect_headers(input_dir);
    match write_tsv(&mer_hash_vec, &full_out_path, &column_names) {
        Ok(()) => println!("TSV file successfully written"),
        Err(err) => eprintln!("Error writing TSV file: {}", err),
    }
}

fn process_files_with_extensions(folder_path: &str) -> Vec<HashMap<String, usize>> {
    // Read the contents of the directory
    let dir = match fs::read_dir(folder_path) {
        Ok(dir) => dir,
        Err(e) => {
            eprintln!("Error reading directory: {}", e);
            std::process::exit(1);
        }
    };

    // Iterate over the files in the directory
    let mut all_file_hashes: Vec<HashMap<String, usize>> = Vec::new();
    for f in dir {
        let f = match f {
            Ok(f) => f,
            Err(e) => {
                eprintln!("Error reading file: {}", e);
                continue;
            }
        };
        let file_path = f.path();
        if check_extension(&file_path, &EXTENSIONS) {
            let fasta: sequence_utils::FastaSeq =
                sequence_utils::fasta_seq_construct(file_path.to_str().unwrap());

            // Count k-mers
            let mut kmer_counts: HashMap<String, usize> = HashMap::new();
            let mut rev_comp_kmer_counts: HashMap<String, usize> = HashMap::new();
            let mut rmer_counts: HashMap<String, usize> = HashMap::new();
            let mut rev_comp_rmer_counts: HashMap<String, usize> = HashMap::new();

            for k_size in KMER_SIZES.iter() {
                let add_hash = sequence_utils::count_kmers(&fasta.sequence, k_size);
                kmer_counts.extend(add_hash);
            }

            for k_size in KMER_SIZES.iter() {
                let add_hash = sequence_utils::count_kmers(&fasta.revcomp_seq(), k_size);
                rev_comp_kmer_counts.extend(add_hash);
            }

            for r_size in RMER_SIZES.iter() {
                let add_hash = sequence_utils::count_kmers(&fasta.yr_seq(), r_size);
                rmer_counts.extend(add_hash);
            }

            for r_size in RMER_SIZES.iter() {
                let add_hash = sequence_utils::count_kmers(&fasta.rev_comp_yr(), r_size);
                rev_comp_rmer_counts.extend(add_hash);
            }

            let combined_kr_hashes: Vec<HashMap<String, usize>> = vec![
                kmer_counts,
                rev_comp_kmer_counts,
                rmer_counts,
                rev_comp_rmer_counts,
            ];

            let collapsed_combined_kr_hashes = collapse_hashes(&combined_kr_hashes);
            all_file_hashes.push(collapsed_combined_kr_hashes);
        };
    }
    all_file_hashes
}

fn check_extension(file_path: &PathBuf, extensions: &[&str]) -> bool {
    let ext = file_path.extension().unwrap(); // Get extension
    println!("Extension: {:?}", ext);
    let ext_str = ext.to_string_lossy().to_lowercase(); // Convert to lowercase for case insensitivity
    println!("Extension string: {}", ext_str);
    if extensions.iter().any(|&e| e == ext_str) {
        println!("Processing file: {:?}", file_path);
        true
    } else {
        false
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

fn write_tsv(
    vec_of_hashmaps: &Vec<HashMap<String, usize>>,
    out_file: &PathBuf,
    file_names: &[String],
) -> Result<(), std::io::Error> {
    // Collect all unique keys
    let mut unique_keys: HashSet<&String> = HashSet::new();
    for hashmap in vec_of_hashmaps {
        for key in hashmap.keys() {
            unique_keys.insert(key);
        }
    }

    // Open the file for writing
    let mut file = File::create(out_file)?;

    // Write column names
    let mut col_names = Vec::new();
    col_names.push(String::from("mer"));
    for name in file_names {
        col_names.push(name.to_string());
    }
    writeln!(file, "{}", col_names.join("\t"))?;

    // Write data
    let mut row_values = Vec::new();
    for key in &unique_keys {
        row_values.clear();
        row_values.push(key.to_string());
        for hashmap in vec_of_hashmaps {
            let value = hashmap.get(*key).unwrap_or(&0).to_string(); // If key not found, default to 0
            row_values.push(value);
        }
        writeln!(file, "{}", row_values.join("\t"))?;
    }

    Ok(())
}

fn collect_headers(folder_path: &String) -> Vec<String> {
    // Read the contents of the directory
    let dir = match fs::read_dir(folder_path) {
        Ok(dir) => dir,
        Err(e) => {
            eprintln!("Error reading directory: {}", e);
            std::process::exit(1);
        }
    };

    // Iterate over the files in the directory
    let mut file_names: Vec<String> = Vec::new();
    for f in dir {
        let f = match f {
            Ok(f) => f,
            Err(e) => {
                eprintln!("Error reading file: {}", e);
                continue;
            }
        };
        let file_path = f.path();
        println!("Processing file: {:?}", file_path);

        if check_extension(&file_path, &EXTENSIONS) {
            let file_name = file_path.file_name().unwrap().to_str().unwrap();
            file_names.push(file_name.to_string());
        }
    }
    file_names
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_collect_headers() {
        let test_dir = "/Users/leleiwi1/Desktop/LLNL_postdoc/kmer_tool_16S/test_input";
        let file_names = collect_headers(&test_dir.to_string());
        assert_eq!(file_names.len(), 4);
        assert_eq!(file_names[0], "t.fasta");
        assert_eq!(file_names[1], "t.fna");
        assert_eq!(file_names[2], "test1.fa");
        assert_eq!(file_names[3], "test2.fa");
    }

    #[test]

    fn test_check_extension() {
        let test_dir = "/Users/leleiwi1/Desktop/LLNL_postdoc/kmer_tool_16S/test_input";
    }
}
