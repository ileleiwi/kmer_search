use bio::io::fasta;
use bio::alignment::sparse::hash_kmers;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;

fn count_kmers(fasta_file: &str, k: usize) -> HashMap<String, usize> {
    // Open the FASTA file
    let file = match File::open(fasta_file) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("Error opening FASTA file: {}", e);
            std::process::exit(1);
        }

    };
    let reader = BufReader::new(file);

    // Create HashMap to stroe k-mer counts
    let mut kmer_counts: std::collections::HashMap<String, usize> = std::collections::HashMap::new();

    // Create a FASTA reader
    let fasta_reader = fasta::Reader::new(reader);

    // Iterate over the records in the FASTA file
    for result in fasta_reader.records() {
        let record = result.expect("Error during fasta record parsing");
        let sequence = record.seq();
        // Create an iterator over the kmers
        let kmer_iter = hash_kmers(&sequence, k);

        // Iterate over kmers and update counts
        for (kmer, _) in &kmer_iter {
            let kmer_str = String::from_utf8_lossy(kmer).to_string();
            let count = kmer_counts.entry(kmer_str).or_insert(0);
            *count += 1;
          }
    }

    kmer_counts
}

fn main() {
    // Get command-line arguments
    let args: Vec<String> = std::env::args().collect();

    // Check if a filename argument and k-mer size argument are provided
    if args.len() != 3 {
        eprintln!("Usage: cargo run <fasta_file> <k>");
        std::process::exit(1);
    }

    // Open the file and read the contents into a String
    let fasta_file: &String = &args[1];
    let k: usize = std::env::args().nth(2).expect("Invalid k-mer size").parse().expect("Invalid k-mer size");

    // Count k-mers in the FASTA file
    let kmer_counts: HashMap<String, usize> = count_kmers(fasta_file, k);

    // Print k-mer counts
    for (kmer, count) in &kmer_counts {
        println!("{}: {}", kmer, count); 
      }
}

