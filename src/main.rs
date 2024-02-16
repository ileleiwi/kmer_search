use clap::{Command, Arg};
use std::f32::consts::E;
use std::fs;
use std::fs::File;
use std::hash::Hash;
use std::str::FromStr;
use std::collections::HashSet;
use std::collections::HashMap;
use std::io::Write;
use std::path::PathBuf;

pub mod sequence_utils;

const KMER_SIZES: [usize; 5] = [1, 2, 3, 4, 5];
const RMER_SIZES: [usize; 5] = [6, 7, 8, 9, 10];
const EXTENSIONS: [&str; 3] = [".fa", ".fasta", ".fna"];

fn main() {
    // Get command-line arguments
    let matches: clap::ArgMatches = Command::new("kmer_search")
    .arg(Arg::new("input_dir")
            .short('i')
            .long("input_dir")
            .value_name("DIRECTORY")
            .help("Path to directory of input FASTA files")
            .required(true))
    .arg(Arg::new("output_dir")
        .short('o')
        .long("output_dir")
        .value_name("DIRECTORY")
        .help("Output directory")
        .required(false))
    .get_matches();

// Get input directory
let input_dir = matches.get_one::<String>("input_dir").unwrap();

// Use default output file names if not provided
let default_output_dir = String::from("./kmer_search_output");
let output_dir = matches
    .get_one::<String>("output_dir")
    .map_or(PathBuf::from(&default_output_dir), |dir| PathBuf::from(dir));

// Create output directory if it doesn't exist
match fs::create_dir_all(&output_dir) {
    Ok(_) => println!("Output directory created"),
    Err(e) => {println!("Error creating output directory: {}", e);},
}

// Get list of FASTA files in input directory
process_files_with_extensions(&input_dir);

// Make FastaSeq structs from FASTA files

let fast2: sequence_utils::FastaSeq = sequence_utils::fasta_seq_construct(&fasta_file2);





let kmer_counts2 = sequence_utils::count_kmers(&fast2.sequence, kmer_size);

// Find intersection of k-mers
// Clone the key sets
let kmer_keys1 = kmer_counts1.keys().cloned().collect::<HashSet<_>>();
let kmer_keys2 = kmer_counts2.keys().cloned().collect::<HashSet<_>>();

// Intersect the cloned key sets
let shared_kmers: HashSet<String> = kmer_keys1.intersection(&kmer_keys2).cloned().collect();

let mut shared_file = File::create(shared_out_path).unwrap();
writeln!(shared_file, "kmer\tfasta1_count\tfasta2_count").unwrap();
for kmer in &shared_kmers {
    let count1 = kmer_counts1.get(kmer).unwrap_or(&0);
    let count2 = kmer_counts2.get(kmer).unwrap_or(&0);
    match writeln!(shared_file, "{}\t{}\t{}", kmer, count1, count2) {
        Ok(_) => {
            // Do nothing
        },
        Err(e) => {
            eprintln!("Error writing to shared file: {}", e);
            std::process::exit(1);
        }
    };
}

// Find unique k-mers
let unique1: Vec<String> = kmer_keys1.difference(&kmer_keys2).cloned().collect();
let unique2: Vec<String> = kmer_keys2.difference(&kmer_keys1).cloned().collect();

let mut unique_file = File::create(unique_out_path).unwrap();
writeln!(unique_file, "kmer\tfasta1_count\tfasta2_count").unwrap();
match sequence_utils::find_unique_kmers(&unique1, &kmer_counts1, &mut unique_file, true) {
    Ok(()) => println!("Unique kmers from fasta1 written to file"),
    Err(e) => eprintln!("Error writing unique kmers to file: {}", e),
}
match sequence_utils::find_unique_kmers(&unique2, &kmer_counts2, &mut unique_file, false) {
    Ok(()) => println!("Unique kmers from fasta2 written to file"),
    Err(e) => eprintln!("Error writing unique kmers to file: {}", e),
}

}


fn process_files_with_extensions(folder_path: &str) -> Vec<HashMap<String, usize>> {
    // Read the contents of the directory
    let dir = match fs::read_dir(folder_path){
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
        if check_extension(file_path, &EXTENSIONS){
            let fasta: sequence_utils::FastaSeq = sequence_utils::fasta_seq_construct(&file_path.to_str().unwrap());
            // Count k-mers 
            let mut kmer_counts: Vec<HashMap<String, usize>> = Vec::new();
            let mut rmer_counts: Vec<HashMap<String, usize>> = Vec::new();

            for (index, k_size) in KMER_SIZES.iter().enumerate() {
                let rev_idx = index + 5;
                kmer_counts[index] = sequence_utils::count_kmers(&fasta.sequence, k_size);
                kmer_counts[rev_idx] = sequence_utils::count_kmers(&fasta.rev_comp, k_size);                
            }
            
            for (index, r_size) in RMER_SIZES.iter().enumerate() {
                let rev_idx = index + 5;
                rmer_counts[index] = sequence_utils::count_kmers(&fasta.yr_seq, r_size);
                rmer_counts[rev_idx] = sequence_utils::count_kmers(&fasta.rev_comp_yr, r_size);                 
            }
            let mut kmer_hashes = collapse_hashes(&kmer_counts);
            let rmer_hashes = collapse_hashes(&rmer_counts);
           
            let all_file_hashes = kmer_hashes.extend(rmer_hashes);

        };
    }
    all_file_hashes
}

fn check_extension(file_path: PathBuf, extensions: &[&str]) -> bool {
    
    let ext = file_path.extension().unwrap(); // Get extension
    let ext_str = ext.to_string_lossy().to_lowercase(); // Convert to lowercase for case insensitivity
        
    if extensions.iter().any(|&e| e == ext_str) {
        println!("Processing file: {:?}", file_path);
        true
    }else{
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

fn write_tsv(vec_of_hashmaps: &Vec<HashMap<String, usize>>, filename: &str) -> Result<(), std::io::Error> {
    // Collect all unique keys
    let mut unique_keys: HashSet<&String> = HashSet::new();
    for hashmap in vec_of_hashmaps {
        for key in hashmap.keys() {
            unique_keys.insert(key);
        }
    }

    // Open the file for writing
    let mut file = File::create(filename)?;

    // Write to the file
    // Write header

    //todo write function that saves fasta file names in a vector and use that as headers for the output tsv
    writeln!(file, "{}", unique_keys.iter().map(|k| *k).collect::<Vec<&String>>().join("\t"))?;

    // Write data
    for hashmap in vec_of_hashmaps {
        let mut row_values = Vec::new();
        for key in &unique_keys {
            let value = hashmap.get(*key).unwrap_or(&0); // If key not found, default to 0
            row_values.push(value.to_string());
        }
        writeln!(file, "{}", row_values.join("\t"))?;
    }

    Ok(())
}











