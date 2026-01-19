# Claude Code Bioinformatics Demo Examples

## Basic

| Question | Command/Approach |
|----------|------------------|
| How many sequences in a FASTA? | `grep -c "^>" example.fasta` |
| How many reads in a FASTQ? | `awk 'END {print NR/4}' example.fastq` |
| What's the mapping quality of a read? | `samtools view example.bam \| grep "READ_ID" \| awk '{print $5}'` |

## Intermediate

| Question | Command/Approach |
|----------|------------------|
| Count reads in a FASTQ (robust method) | Divide total lines by 4 (each read = 4 lines) |
| View BAM header info | `samtools view -H example.bam` |
| Check VCF samples | `bcftools query -l file.vcf` |
| View variants in a region | `bcftools view -H -r chr:start-end file.vcf` |

## Advanced

| Question | Command/Approach |
|----------|------------------|
| Extract heterozygous variants with DP>10 | `bcftools view -i 'GT="het" && FORMAT/DP>10' file.vcf` |
| BLAST a sequence from a file | Playwright browser automation to NCBI |
| Interpret CIGAR strings | Parse column 6 of SAM (e.g., `10M1D10M2I10M18S`) |
| Multi-tool pipelines | Pipe samtools → grep → awk for complex filtering |

## Sample Files Available

```
data/
├── example.fasta      # 3 sequences
├── example.fastq      # 10 reads
├── example.sam        # 11 alignments (various CIGARs)
├── example.bam        # Binary alignment + index
├── ex.vcf             # 5 variants, 3 samples
└── aln.bt.vcf.gz      # Larger compressed VCF
```

## Key Teaching Points

1. **Natural language → commands** - Ask questions conversationally, Claude figures out the tools
2. **Tool selection** - Claude chooses samtools vs bcftools vs grep appropriately
3. **Explanation** - Claude explains *why* commands work, not just *what* they do
4. **Error handling** - Claude adapts when approaches don't work (e.g., robust FASTQ counting)

## Example Demo Flow

### Basic: FASTA counting
```
You: How many sequences are in data/example.fasta?
Claude: [runs grep -c "^>" ...] → 3 sequences
```

### Intermediate: BAM query
```
You: What is the mapping quality of ERR003762.5016205?
Claude: [runs samtools view + grep + awk] → MAPQ: 42
```

### Advanced: VCF filtering
```
You: Extract all heterozygous variants with depth > 10
Claude: [runs bcftools view -i 'GT="het" && FORMAT/DP>10' ...] → 5 variants
```

### Web automation: BLAST
```
You: Can you BLAST the first sequence in the FASTA file?
Claude: [extracts sequence, opens NCBI BLAST via Playwright, submits search]
```
