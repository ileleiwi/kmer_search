use bio::alignment::sparse::hash_kmers;
use bio::io::fasta;
use std::collections::HashMap;
use std::fs::File;
use std::io::Error;
use std::io::Write;
use std::path::Path;

//TODO: add an revcomp utility to count kmers that works on odd kmer lengths
pub struct FastaSeq {
    pub file_name: String,
    pub sequence: String,
    pub bases: usize,
    pub gc: usize,
}

impl FastaSeq {
    pub fn revcomp_seq(&self) -> String {
        let mut revcomp = String::new();
        for ch in self.sequence.chars().rev() {
            match ch {
                'A' => revcomp.push('T'),
                'C' => revcomp.push('G'),
                'G' => revcomp.push('C'),
                'T' => revcomp.push('A'),
                _ => revcomp.push(ch),
            }
        }
        revcomp
    }

    pub fn yr_seq(&self) -> String {
        let mut yr_seq = String::new();
        for ch in self.sequence.chars().map(|c| c.to_ascii_uppercase()) {
            match ch {
                'A' | 'G' => yr_seq.push('R'),
                'C' | 'T' => yr_seq.push('Y'),
                _ => yr_seq.push('-'),
            }
        }
        yr_seq
    }

    pub fn rev_comp_yr(&self) -> String {
        let mut yr_seq = String::new();
        for ch in self.sequence.chars().map(|c| c.to_ascii_uppercase()) {
            match ch {
                'A' | 'G' => yr_seq.push('R'),
                'C' | 'T' => yr_seq.push('Y'),
                _ => yr_seq.push('-'),
            }
        }
        let mut revcomp_yr = String::new();
        for ch in yr_seq.chars().rev() {
            match ch {
                'Y' => revcomp_yr.push('R'),
                'R' => revcomp_yr.push('Y'),
                _ => revcomp_yr.push(ch),
            }
        }
        revcomp_yr
    }
}

pub fn fasta_seq_construct(fasta_file: &str) -> FastaSeq {
    let base_name = match Path::new(fasta_file)
        .file_stem()
        .and_then(|stem| stem.to_str())
    {
        Some(stem) if stem.ends_with(".fa") => stem.strip_suffix(".fa").unwrap_or(stem),
        Some(stem) if stem.ends_with(".fna") => stem.strip_suffix(".fna").unwrap_or(stem),
        Some(stem) if stem.ends_with(".fasta") => stem.strip_suffix(".fasta").unwrap_or(stem),
        Some(stem) => stem,
        None => fasta_file, // If file_name is not valid, return the original file_name
    };

    let seq = seq_string_from_fasta(fasta_file);
    let mut bases = 0;
    let mut gc = 0;
    for c in seq.chars() {
        match c {
            'A' | 'C' | 'G' | 'T' => {
                bases += 1;
                if c == 'C' || c == 'G' {
                    gc += 1;
                }
            }
            _ => {}
        }
    }

    FastaSeq {
        file_name: base_name.to_string(),
        sequence: seq,
        bases,
        gc,
    }
}

pub fn seq_string_from_fasta(fasta_file: &str) -> String {
    let reader = match fasta::Reader::from_file(fasta_file) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("Error opening FASTA file: {}", e);
            std::process::exit(1);
        }
    };
    let mut fasta_sequence = String::new();
    for record in reader.records() {
        let seq = record
            .expect("Error reading sequence record")
            .seq()
            .to_owned();
        let seq_string = String::from_utf8_lossy(&seq);
        fasta_sequence.push_str(&seq_string);
    }
    fasta_sequence
}

pub fn count_kmers(fasta_seq: &str, k: &usize) -> HashMap<String, usize> {
    // Create k-mer iterator
    let kmer_iter = hash_kmers(fasta_seq.as_bytes(), *k);

    // Create HashMap to store k-mer counts instead of k-mer positions
    let mut kmer_counts: HashMap<String, usize> = HashMap::new();

    // Iterate over k-mers
    for (kmer, count_vec) in kmer_iter {
        let kmer_string = String::from_utf8_lossy(kmer).to_string();
        let entry = kmer_counts.entry(kmer_string).or_insert(0);
        *entry += count_vec.len();
    }

    kmer_counts
}

pub fn find_unique_kmers(
    kmer_vec: &Vec<String>,
    kmer_counts: &HashMap<String, usize>,
    file: &mut File,
    fasta_1: bool,
) -> Result<(), Error> {
    for kmer in kmer_vec {
        let count = kmer_counts[kmer];
        match writeln!(
            file,
            "{}\t{}\t{}",
            kmer,
            if fasta_1 { count } else { 0 },
            if fasta_1 { 0 } else { count }
        ) {
            Ok(_) => {}
            Err(e) => {
                return Err(e);
            }
        };
    }
    Ok(())
}
