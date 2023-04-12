//! This contains some data structures

use crate::Args;
use crate::VERSION;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

/// Stores basic read stats.
#[derive(Debug)]
pub struct ReadStat {
    pub multimap: i64,
    pub duplicate: i64,
    pub secondary: i64,
    pub too_short: i64,
    pub low_qual: i64,
    pub total: i64,
}

impl ReadStat {
    pub fn new() -> Self {
        ReadStat {
            multimap: 0,
            duplicate: 0,
            secondary: 0,
            too_short: 0,
            low_qual: 0,
            total: 0,
        }
    }

    pub fn print(&self) {
        println!("Number of read processed: {}", self.total);
        println!("Read stats:");
        println!("  multimap  : {}", self.multimap);
        println!("  duplicate : {}", self.duplicate);
        println!("  secondary : {}", self.secondary);
        println!("  too_short : {}", self.too_short);
        println!("  low_qual  : {}", self.low_qual);
    }
}

/// Stores depth of ref bases
#[derive(Debug)]
pub struct Depth {
    pub v: Vec<u32>,
}

impl Depth {
    pub fn new(length: usize) -> Self {
        Depth { v: vec![0; length] }
    }

    pub fn add(&mut self, pos: usize) {
        self.v[pos] += 1;
    }
}

/// Stores SNVs
#[derive(Debug)]
pub struct SNVs {
    pub map: HashMap<usize, [u32; 4]>,
    pub pos: Vec<usize>,
}

impl SNVs {
    pub fn new() -> Self {
        SNVs {
            map: HashMap::new(),
            pos: Vec::new(),
        }
    }

    pub fn add(&mut self, pos: usize, nucl: u8) {
        let index = convert_from_char_to_arr_index(nucl);
        if index == 4 {
            return;
        }
        match self.map.get_mut(&pos) {
            Some(arr) => {
                arr[index] += 1;
            }
            None => {
                let mut arr = [0; 4];
                arr[index] += 1;
                self.map.insert(pos, arr);
                self.pos.push(pos);
            }
        };
    }
}

/// just convert from nucl char to index of arr
#[inline]
fn convert_from_char_to_arr_index(nucl: u8) -> usize {
    match nucl {
        65 => 0,  // A
        84 => 1,  // T
        71 => 2,  // G
        67 => 3,  // C
        97 => 0,  // a
        116 => 1, // t
        103 => 2, // g
        99 => 3,  // c
        _ => 4,
    }
}

/// Stores IO readers and writers.
pub struct Fstream {
    pub writer: BufWriter<File>,
}

impl Fstream {
    pub fn new(args: &Args) -> Self {
        let filename = PathBuf::from(&args.o);
        let mut writer = BufWriter::new(File::create(filename.as_path()).unwrap());

        // write sample info
        writer
            .write_all(format!("#version:{}\n", VERSION).as_bytes())
            .unwrap();
        writer
            .write_all(format!("#input:{}\n", args.i).as_bytes())
            .unwrap();
        writer
            .write_all(format!("#reference:{}\n", args.r).as_bytes())
            .unwrap();
        match &args.b {
            Some(bed) => writer
                .write_all(format!("#bed:{}\n", bed).as_bytes())
                .unwrap(),
            None => writer.write_all(format!("#bed:None\n").as_bytes()).unwrap(),
        }
        writer
            .write_all(format!("#min VAF (-v flag):{}\n", args.v).as_bytes())
            .unwrap();
        writer
            .write_all(format!("#min depth (-d flag):{}\n", args.d).as_bytes())
            .unwrap();
        writer
            .write_all(format!("chr\tpos\tref\talt\tVAF\tA\tT\tG\tC\tis_multiallelic\n").as_bytes())
            .unwrap();

        Fstream { writer: writer }
    }

    pub fn flush_all(&mut self) {
        self.writer.flush().unwrap();
    }
}
