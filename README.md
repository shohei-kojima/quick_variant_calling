# Cheap but quick variant calling from WGS and RNA-seq

## About
This is a simple program that calls variants from RNA-seq and WGS.  
This is intended to be used for quick checking of sample swapping by seeing SNVs.  
Results from this script is not reliable - should NOT be used as a final dataset.  

## Usage
```
git clone https://github.com/shohei-kojima/quick_variant_calling
cd quick_variant_calling
cargo build --release

# takes < 1 min and ~4GB RAM when analysing human RNA-seq with 25M reads
./target/release/quick_variant_calling \
-i /path/to/input.bam \
-o ./test.out \
-b ./test.bed \
-r /path/to/reference_genome.fa \
-d 30 \
-v 0.05 \
-p 2

# help message
./target/release/quick_variant_calling -h
```

