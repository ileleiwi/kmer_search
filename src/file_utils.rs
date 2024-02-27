use crate::sequence_utils;
use crate::sequence_utils::normalize_kmers;
use crate::sequence_utils::FastaSeq;
use crate::EXTENSIONS;
use crate::KMER_SIZES;
use crate::RMER_SIZES;
use std::collections::HashMap;
use std::collections::HashSet;
use std::fs;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::path::PathBuf;

pub fn check_extension(file_path: &Path, extensions: &[&str]) -> bool {
    let ext = file_path.extension().unwrap(); // Get extension

    let ext_str = ext.to_string_lossy().to_lowercase(); // Convert to lowercase for case insensitivity

    extensions.iter().any(|&e| e == ext_str) // Return true if extension matches false otherwise
}
pub fn collect_headers(fastas: &Vec<FastaSeq>) -> Vec<String> {
    let mut file_names: Vec<String> = Vec::new();
    for f in fastas {
        file_names.push(f.file_name.clone());
    }
    file_names
}

pub fn collapse_hashes(mer_hash_vec: &Vec<HashMap<String, usize>>) -> HashMap<String, usize> {
    let mut collapsed_hash: HashMap<String, usize> = HashMap::new();
    for hash in mer_hash_vec {
        for (id, count) in hash {
            let mer = collapsed_hash.entry(id.to_string()).or_insert(0);
            *mer += count;
        }
    }
    collapsed_hash
}
pub fn process_files_with_extensions(fastas: &Vec<FastaSeq>) -> Vec<HashMap<String, f64>> {
    // Iterate over the fastas
    let mut all_file_hashes: Vec<HashMap<String, f64>> = Vec::new();
    for f in fastas {
        // Count k-mers
        let mut kmer_counts: HashMap<String, usize> = HashMap::new();
        let mut rev_comp_kmer_counts: HashMap<String, usize> = HashMap::new();
        let mut rmer_counts: HashMap<String, usize> = HashMap::new();
        let mut rev_comp_rmer_counts: HashMap<String, usize> = HashMap::new();

        for k_size in KMER_SIZES.iter() {
            let add_hash = sequence_utils::count_kmers(&f.sequence, k_size);
            kmer_counts.extend(add_hash);
        }

        for k_size in KMER_SIZES.iter() {
            let add_hash = sequence_utils::count_kmers(&f.revcomp_seq(), k_size);
            rev_comp_kmer_counts.extend(add_hash);
        }

        for r_size in RMER_SIZES.iter() {
            let add_hash = sequence_utils::count_kmers(&f.yr_seq(), r_size);
            rmer_counts.extend(add_hash);
        }

        for r_size in RMER_SIZES.iter() {
            let add_hash = sequence_utils::count_kmers(&f.rev_comp_yr(), r_size);
            rev_comp_rmer_counts.extend(add_hash);
        }

        let combined_kr_hashes: Vec<HashMap<String, usize>> = vec![
            kmer_counts,
            rev_comp_kmer_counts,
            rmer_counts,
            rev_comp_rmer_counts,
        ];

        let collapsed_combined_kr_hashes = collapse_hashes(&combined_kr_hashes);
        let collapsed_combined_kr_hashes_norm =
            normalize_kmers(&collapsed_combined_kr_hashes, &f.bases);

        all_file_hashes.push(collapsed_combined_kr_hashes_norm);
    }
    all_file_hashes
}

pub fn instantiate_fastaseq_objects(folder_path: &str) -> Vec<FastaSeq> {
    // Read the contents of the directory
    let dir = match fs::read_dir(folder_path) {
        Ok(dir) => dir,
        Err(e) => {
            eprintln!("Error reading directory: {}", e);
            std::process::exit(1);
        }
    };

    // Iterate over the files in the directory
    let mut fasta_seqs: Vec<FastaSeq> = Vec::new();
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
            fasta_seqs.push(fasta);
        };
    }
    fasta_seqs
}

// TODO: Write a function to write the statistics to a TSV file
pub fn write_stats_tsv(fastas: &Vec<FastaSeq>, out_file: &PathBuf) -> Result<(), std::io::Error> {
    // Open the file for writing
    let mut file = File::create(out_file)?;

    // Write column names
    let column_names = ["fasta", "seq_length", "gc_content", "sequence"];
    writeln!(
        file,
        "{}\t{}\t{}\t{}",
        column_names[0], column_names[1], column_names[2], column_names[3]
    )?;

    for f in fastas {
        // Write data
        writeln!(
            file,
            "{}\t{}\t{}\t{}",
            f.file_name, f.bases, f.gc, f.sequence
        )?;
    }

    Ok(())
}

pub fn write_mer_tsv(
    vec_of_hashmaps: &Vec<HashMap<String, f64>>,
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
            let value = hashmap.get(*key).unwrap_or(&0.0).to_string(); // If key not found, default to 0
            row_values.push(value);
        }
        writeln!(file, "{}", row_values.join("\t"))?;
    }

    Ok(())
}
