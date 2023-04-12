//! This is a simple program that converts a BAM/CRAM file into lossy CRAM file with quality scores of Illumina 8-bins

// Author: Shohei Kojima @ RIKEN

const VERSION: &str = "version 0.0.1";
extern crate my_rust;
extern crate clap;
extern crate num_cpus;
extern crate rust_htslib;
extern crate lazy_static;

mod data_structs;

use std::collections::HashMap;
use std::path::PathBuf;
use std::process::exit;
use std::fs::File;
use std::io::{BufReader, BufRead};
use std::io::Write;
use clap::{AppSettings, Parser};
use rust_htslib::bam::{IndexedReader, Read};
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::record::{Aux, Cigar, CigarStringView};
use lazy_static::lazy_static;
use regex::Regex;
use linear_map::LinearMap;


const MINSEQLEN: usize = 20;
const MIN_BASE_QUAL: u8 = 20;
const ATGC: [char; 4] = ['A', 'T', 'G', 'C'];


fn main() {
    // initial file checks
    let args = Args::parse();
    file_exist_check(&args.i);
    let is_cram = bam_format_check(&args);
    file_absence_check(&args.o, args.n);
    
    // thread number check
    if args.p == 0 {
        eprintln!("Number of threads must be 1 or more.");
        exit(1);
    }
    if args.p > num_cpus::get() {
        eprintln!(
            "Number of threads ({}) exceeds the number of threads ({}) in your machine.",
            args.p,
            num_cpus::get()
        );
        exit(1);
    } else {
        println!(
            "Number of threads ({}) was set; {} threads found in your machine.",
            args.p,
            num_cpus::get()
        );
    }
    
    // make output file
    let mut fstream = crate::data_structs::Fstream::new(&args.o);
    
    // make structs
    let mut read_stat = crate::data_structs::ReadStat::new();
    
    // load reference fasta
    let mut ref_fasta = my_rust::io::FastaRecords::from_file(&args.r);
    ref_fasta.to_upper();
    
    // load target regions
    let beds: Vec<(String, i64, i64)> = match args.b {
        Some(path) => from_bed_file(&path),
        None => from_bai(&args.i),
    };
    
    // make infile obj
    let mut infile = IndexedReader::from_path(&args.i).unwrap();
    if is_cram {
        infile.set_reference(&args.r).unwrap();
    }
    if args.p >= 2 {
        infile.set_threads(args.p - 1).unwrap();
    }
    
    // read bam, per region/chr analysis
    for (chr, start, end) in beds {
        // load fasta, refseq = &str
        let refseq = match ref_fasta.get_seq(&chr) {
            Some(seq) => seq,
            None => panic!("chromosome '{}' was not found in the input reference file.", chr),
        };
        let ref_length = refseq.len();
        let refseqbytes: &[u8] = refseq.as_bytes();
        let mut snvs = crate::data_structs::SNVs::new();
        let mut depth = crate::data_structs::Depth::new(ref_length);
        infile.fetch((&chr, start, end)).unwrap();
        for read in infile.rc_records().map(|x| x.expect("Failure parsing input file")) {
            let mut skip = false;
            read_stat.total += 1;
            if read.is_duplicate() {
                read_stat.duplicate += 1;
                skip = true;
            }
            if read.is_secondary() {
                read_stat.secondary += 1;
                skip = true;
            }
            if read.seq_len() < MINSEQLEN {
                read_stat.too_short += 1;
                skip = true;
            }
            if let Ok(Aux::U8(nh)) = read.aux(b"NH") {
                if nh >= 2 {
                    read_stat.multimap += 1;
                    skip = true;
                }
            }
            if skip {
                continue;
            }
            let s = read.pos(); // left  position, 0-based, i64, equals to C `int64_t& start = b->core.pos`
            let e = read.reference_end(); // right position, 0-based, i64, equals to C `int64_t end = bam_endpos(b)`
            let cigar = read.cigar();
            let seq = read.seq().as_bytes(); // Seq { encoded: &[u8], len: usize }
            let qual: &[u8] = read.qual();
            let spans = cigar_to_ref_and_read_pos(&cigar, s, e);
            for (ref_s, read_s, span) in spans {
                for n in 0..span as usize {
                    let read_pos: usize = read_s + n;
                    // check base quality
                    if qual[read_pos] < MIN_BASE_QUAL {
                        continue;
                    }
                    let ref_pos: usize = ref_s + n;
                    if refseqbytes[ref_pos] == seq[read_pos] {
                        depth.add(ref_pos);
                    } else {
                        snvs.add(ref_pos, seq[read_pos]);
                    }
                }
            }
        }
        // output
        for pos in &snvs.pos {
            let alts = snvs.map.get_mut(pos).unwrap();
            let refcount = depth.v[*pos];
            let alt_sum: u32 = alts.iter().sum();
            let total = refcount + alt_sum;
            if total < args.d {
                continue;
            }
            let total_f64 = total as f64;
            let mut tmp: Vec<(usize, f64)> = Vec::new();
            for n in 0..4 {
                if (alts[n] + refcount) < args.d {
                    continue;
                }
                let vaf = alts[n] as f64 / total_f64;
                if vaf > args.v {
                    tmp.push((n, vaf));
                }
            }
            let mut is_multiallelic: bool = false;
            if tmp.len() >= 2 {
                is_multiallelic = true;
            }
            let refbase = u8_to_char_base(refseqbytes[*pos]);
            if refbase.eq(&'N') {
                continue;
            }
            let index = u8_to_arr_index(refseqbytes[*pos]);
            alts[index] += refcount;
            for (n, vaf) in tmp {
                fstream.writer
                    .write_all(format!("{}\t{}\t{}\t{}\t{:.4}\t{}\t{}\t{}\t{}\t{}\n", chr, pos, refbase, ATGC[n], vaf, alts[0], alts[1], alts[2], alts[3], is_multiallelic as u8).as_bytes())
                    .unwrap();
            }
        }
    }
    fstream.flush_all(); // dropeed soon automatically
    read_stat.print();  // read stats
}

/// convert CigarStringView to vec of ref range
/// returns spans of alignments: Vec<(ref_pos, read_pos, span)>
#[inline]
fn cigar_to_ref_and_read_pos(cigar: &CigarStringView, refstart: i64, refend: i64) -> Vec<(usize, usize, usize)> {
    let mut tmp: Vec<(usize, bool, bool)> = Vec::new(); // (usize, bool, bool) = (length, consume_ref, consume_read)
    for c in (*cigar).iter() {
        match c {
            Cigar::Match(l) => tmp.push((*l as usize, true, true)),
            Cigar::Ins(l) => tmp.push((*l as usize, false, true)),
            Cigar::Del(l) => tmp.push((*l as usize, true, false)),
            Cigar::RefSkip(l) => tmp.push((*l as usize, true, false)),
            Cigar::SoftClip(l) => tmp.push((*l as usize, false, true)),
            Cigar::HardClip(l) => tmp.push((*l as usize, false, true)),
            Cigar::Equal(l) => tmp.push((*l as usize, true, true)),
            Cigar::Diff(l) => tmp.push((*l as usize, true, true)),
            _ => (),
        }
    }
    let mut v: Vec<(usize, usize, usize)> = Vec::new();
    let mut ref_pos = refstart as usize;
    let mut read_pos: usize = 0;
    for (l, b1, b2) in tmp {
        if b1 && b2 {
            v.push((ref_pos, read_pos, l));
        }
        if b1 {
            ref_pos += l;
        }
        if b2 {
            read_pos += l;
        }
    }
    if !ref_pos == refend as usize {
        panic!("cigar parse failed {}", cigar);
    }
    v
}

/// Converts u8 to A, T, G, or C
fn u8_to_char_base(n: u8) -> char {
    match n {
        65 => 'A',
        84 => 'T',
        71 => 'G',
        67 => 'C',
        _ => 'N',
    }
}

/// Converts u8 to A, T, G, or C
fn u8_to_arr_index(n: u8) -> usize {
    match n {
        65 => 0,
        84 => 1,
        71 => 2,
        67 => 3,
        _ => 4,
    }
}

/// read a bed file and converts to Vec
fn from_bed_file(path: &str) -> Vec<(String, i64, i64)> {
    let mut v = Vec::new();
    let fpath = PathBuf::from(path);
    let f = File::open(fpath).expect("input bed file not found.");
    for line in BufReader::new(f).lines() {
         match line {
             Ok(line) => {
                 let ls = line.trim().split("\t").collect::<Vec<&str>>();
                 if ls.len() < 3 {
                     panic!("Column number in the bed file is less than 3. Please check the file format again.");
                 }
                 let t = (ls[0].to_string(), ls[1].parse::<i64>().unwrap(), ls[2].parse::<i64>().unwrap());
                 v.push(t);
             },
             Err(e) => {
                 panic!("{}", e);
             }
         }
    }
    v
}

fn from_bai(path: &str) -> Vec<(String, i64, i64)> {
    let mut v = Vec::new();
    let bam = rust_htslib::bam::Reader::from_path(path).unwrap();
    let header = rust_htslib::bam::Header::from_template(bam.header());
    let map = bam_header_to_hashmap(&header);
    // for (key, records) in header.to_hashmap() {
    for (key, records) in map {
        if key.as_str() != "SQ" {
            continue;
        }
        // if key is "SQ"
        for record in records {
            let mut chr_name: String = String::new();
            let mut chr_length: i64 = 0;
            for (k, v) in record.iter() {
                if k == "SN" {
                    chr_name = (*v).clone();
                }
                if k == "LN" {
                    chr_length = (*v).parse::<i64>().unwrap();
                }
            }
            v.push((chr_name.clone(), 0, chr_length));
        }
    }
    v
}

/// similar to a fn "to_hashmap" in https://github.com/rust-bio/rust-htslib/blob/master/src/bam/header.rs
/// but avoids error when @CO record_type
fn bam_header_to_hashmap(
    header: &rust_htslib::bam::Header,
) -> HashMap<String, Vec<LinearMap<String, String>>> {
    let mut header_map = HashMap::default();
    lazy_static! {
        static ref REC_TYPE_RE: Regex = Regex::new(r"@([A-Z][A-Z])").unwrap();
        static ref TAG_RE: Regex = Regex::new(r"([A-Za-z][A-Za-z0-9]):([ -~]+)").unwrap();
    }
    let header_string = String::from_utf8(header.to_bytes()).unwrap();
    for line in header_string.split('\n').filter(|x| !x.is_empty()) {
        let parts: Vec<_> = line.split('\t').filter(|x| !x.is_empty()).collect();
        let record_type = REC_TYPE_RE
            .captures(parts[0])
            .unwrap()
            .get(1)
            .unwrap()
            .as_str()
            .to_owned();
        if record_type.eq("CO") {
            continue;
        } // added
        let mut field = LinearMap::default();
        for part in parts.iter().skip(1) {
            let cap = TAG_RE.captures(part).unwrap();
            let tag = cap.get(1).unwrap().as_str().to_owned();
            let value = cap.get(2).unwrap().as_str().to_owned();
            field.insert(tag, value);
        }
        header_map
            .entry(record_type)
            .or_insert_with(Vec::new)
            .push(field);
    }
    header_map
}

/// judge whether the input file is BAM or CRAM
fn bam_format_check(args: &Args) -> bool {
    let mut is_cram: bool = false;
    match PathBuf::from(&args.i)
        .extension()
        .unwrap()
        .to_str()
        .unwrap()
    {
        "cram" => {
            is_cram = true;
            println!("Input file ({}) is considered as CRAM format.", args.i);
        }
        "bam" => println!("Input file ({}) is considered as BAM format.", args.i),
        _ => {
            eprintln!("Input file ({}) does not have extension .bam or .cram", args.i);
            exit(1);
        }
    }
    is_cram
}

/// Check whether a file exists.
/// Exit with 1 if the file does not exist.
fn file_exist_check(file_path: &str) {
    if !my_rust::utils::path_exists(file_path) {
        eprintln!("Input file ({}) does not exist.", file_path);
        exit(1);
    } else {
        println!("Input file ({}) found.", file_path);
    }
}

/// Check whether a file does not exist.
fn file_absence_check(file_path: &str, overwrite: u8) {
    if my_rust::utils::path_exists(file_path) {
        if overwrite == 0 {
            println!("Warn: Output will be overwritten in {}.", file_path);
        } else {
            eprintln!("Output file ({}) already exists.", file_path);
            exit(1);
        }
    } else {
        println!("Output will be written in ({}).", file_path);
    }
}

#[derive(Parser, Debug)]
#[clap(author = "Author: Shohei Kojima @ RIKEN", version = VERSION, about = "Converts a BAM/CRAM file to a lossy compressed CRAM file with 8- or 4-binning for quality scores.", setting = AppSettings::DeriveDisplayOrder)]
#[clap(propagate_version = true)]
struct Args {
    /// Specify an input file [Required]
    #[clap(short, long, value_parser, value_name = "FILE")]
    pub i: String,
    
    /// Specify an output file [Required]
    #[clap(short, long, value_parser, value_name = "FILE")]
    pub o: String,
    
    /// Specify a refernece genome file [Required]
    #[clap(short, long, value_parser, value_name = "FILE")]
    pub r: String,
    
    /// File containing target regions (.bed file) [Optional]
    #[clap(short, long, value_parser, value_name = "FILE")]
    pub b: Option<String>,
    
    /// Minimal seq depth required to call variant [Optional]
    #[clap(short, long, value_parser, value_name = "NUMBER", default_value_t = 10)]
    pub d: u32,
    
    /// Minimal VAF required to call variant [Optional]
    #[clap(short, long, value_parser, value_name = "NUMBER", default_value_t = 0.05)]
    pub v: f64,
    
    /// Number of threads to read BAM/CRAM [Optional]
    #[clap(short, long, value_parser, value_name = "NUMBER", default_value_t = 1)]
    pub p: usize,
    
    /// Specify when you don't want to overwrite [Optional]
    #[clap(short, long, action = clap::ArgAction::Count)]
    pub n: u8,
}
