use clap::{Command, Arg};
use bio::io::fasta;
use bio::alignment::sparse::hash_kmers;
use std::fs;
use std::str::FromStr;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fs::File;
use std::fs::create_dir_all;
use std::io::Write;
use std::io::Error;
use std::path::PathBuf;


fn main() {
    // Get command-line arguments
    let matches: clap::ArgMatches = Command::new("kmer_search")
    .arg(Arg::new("fasta1")
            .short('1')
            .long("fasta1")
            .value_name("FILE")
            .help("Path to first input FASTA file")
            .required(true))
    .arg(Arg::new("fasta2")
        .short('2')
        .long("fasta2")
        .value_name("FILE")
        .help("Path to second input FASTA file")
        .required(true))
    .arg(Arg::new("shared_out")
        .short('s')
        .long("shared_out")
        .value_name("FILE")
        .help("Name of shared output file")
        .required(false))
    .arg(Arg::new("unique_out")
        .short('u')
        .long("unique_out")
        .value_name("FILE")
        .help("Name of unique output file")
        .required(false))
    .arg(Arg::new("output_dir")
        .short('o')
        .long("output_dir")
        .value_name("DIRECTORY")
        .help("Output directory")
        .required(false))
    .arg(Arg::new("kmer_size")
        .short('k')
        .long("kmer_size")
        .value_name("SIZE")
        .help("k-mer size")
        .required(true))
    .get_matches();

    
let fasta_file1 = matches.get_one::<String>("fasta1").unwrap_or_else(|| {
    eprintln!("Error: fasta1 file not provided");
    std::process::exit(1);
});
let fasta_file2 = matches.get_one::<String>("fasta2").unwrap_or_else(|| {
    eprintln!("Error: fasta2 file not provided");
    std::process::exit(1);
});
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

let default_shared_out = String::from("shared.tsv");
let shared_out = matches
    .get_one::<String>("shared_out")
    .unwrap_or(&default_shared_out);
let shared_out_path = output_dir.join(shared_out);

let default_unique_out = String::from("unique.tsv");
let unique_out = matches
    .get_one::<String>("unique_out")
    .unwrap_or(&default_unique_out);
let unique_out_path = output_dir.join(unique_out);


// Retrieve the value of the kmer_size argument
let kmer_size_str = matches.get_one::<String>("kmer_size").unwrap();

// Parse the value to usize
let kmer_size = match usize::from_str(kmer_size_str) {
    Ok(size) => size,
    Err(_) => {
        eprintln!("Error: Invalid k-mer size provided.");
        // Handle the error appropriately
        return;
    }
};

// Count k-mers in first FASTA file
let kmer_counts1 = count_kmers(&fasta_file1, kmer_size);
// Count k-mers in second FASTA file  
let kmer_counts2 = count_kmers(&fasta_file2, kmer_size);

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
match find_unique_kmers(&unique1, &kmer_counts1, &mut unique_file) {
    Ok(()) => println!("Unique kmers from fasta1 written to file"),
    Err(e) => eprintln!("Error writing unique kmers to file: {}", e),
}
match find_unique_kmers(&unique2, &kmer_counts2, &mut unique_file) {
    Ok(()) => println!("Unique kmers from fasta2 written to file"),
    Err(e) => eprintln!("Error writing unique kmers to file: {}", e),
}

}

fn count_kmers(fasta_file: &str, k: usize) -> HashMap<String, usize> {
 
    // Create a FASTA reader
    let reader = match fasta::Reader::from_file(fasta_file){
        Ok(f) => f,
        Err(e) => {
            eprintln!("Error opening FASTA file: {}", e);
            std::process::exit(1);
        }

    };

    // Create HashMap to store k-mer counts
    let mut kmer_counts: HashMap<String, usize> = HashMap::new();

    // Iterate over the records in the FASTA file
    for record in reader.records() {
        let seq = record.expect("Error reading sequence from FASTA1 file").seq().to_owned();
        let kmer_iter = hash_kmers(&seq, k);
        for (kmer, count_vec) in kmer_iter {
            let kmer_string= String::from_utf8_lossy(kmer).to_string();
            let entry = kmer_counts.entry(kmer_string).or_insert(0);
            *entry += count_vec.len();
        }
    }

    kmer_counts
}

fn find_unique_kmers(kmer_vec: &Vec<String>, kmer_counts: &HashMap<String, usize>, file: &mut File) -> Result<(), Error> {
    for kmer in kmer_vec {
        let count = kmer_counts[kmer];
        match writeln!(file, "{}\t{}\t{}", kmer, 0, count) {
            Ok(_) => {},
            Err(e) => {
                return Err(e);
            }
        };
    }
    Ok(())
}
