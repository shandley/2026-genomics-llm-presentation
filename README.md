# 🧬 Genomics LLM Presentation 2026 🤖

> 🎯 Demo files for exploring bioinformatics data with Large Language Models! 🚀

## 📁 What's Inside?

This repository contains example genomics files for demonstrating how LLMs can assist with bioinformatics analysis! 🔬✨

### 📊 Data Files

| File | Format | Description | Count |
|------|--------|-------------|-------|
| 🧬 `example.fasta` | FASTA | DNA sequences | 3 sequences |
| 📖 `example.fastq` | FASTQ | Sequencing reads with quality scores | 10 reads |
| 🗺️ `example.sam` | SAM | Sequence alignments (text) | 11 alignments |
| 📦 `example.bam` | BAM | Sequence alignments (binary) | 11 alignments |
| 📇 `example.bam.bai` | BAI | BAM index file | - |
| 🔴 `ex.vcf` | VCF | Variant calls | 5 variants |
| 🗜️ `aln.bt.vcf.gz` | VCF.GZ | Compressed variant calls | - |

## 🛠️ Tools Used

These are the bioinformatics tools we explored:

| Tool | Purpose | Emoji |
|------|---------|-------|
| 🔧 **samtools** | SAM/BAM manipulation | 🗺️ |
| 🔬 **bcftools** | VCF/BCF manipulation | 🔴 |
| 🎯 **bowtie2** | Short read alignment | 🏹 |
| 🚀 **minimap2** | Fast alignment | ⚡ |
| 💥 **BLAST** | Sequence search | 🔍 |
| 🧰 **seqkit** | FASTA/FASTQ toolkit | 🛠️ |
| ✂️ **seqtk** | FASTA/FASTQ processing | 📐 |
| 🏷️ **picard** | BAM processing | 📋 |

## 🎓 Quick Commands Cheat Sheet

### 📊 Counting Things

```bash
# 🧬 Count sequences in FASTA
grep -c "^>" data/example.fasta

# 📖 Count reads in FASTQ
echo $(( $(wc -l < data/example.fastq) / 4 ))

# 🗺️ Count alignments in BAM
samtools view -c data/example.bam

# 🔴 Count variants in VCF
grep -cv "^#" data/ex.vcf
```

### 📈 Getting Statistics

```bash
# 🧬 FASTA stats
seqkit stats data/example.fasta

# 📖 FASTQ stats
seqkit stats data/example.fastq

# 🗺️ BAM alignment stats
samtools flagstat data/example.bam

# 🔴 VCF variant stats
bcftools stats data/ex.vcf
```

## 🎯 Purpose

This demo shows how LLMs can help with:

- 🤔 **Understanding** file formats and bioinformatics concepts
- 💡 **Suggesting** appropriate tools for analysis tasks
- 📝 **Writing** command-line one-liners
- 🔍 **Interpreting** output and results
- 🎓 **Teaching** bioinformatics workflows

## 🌟 Key Takeaways

1. 🧠 LLMs can explain complex bioinformatics concepts in plain language
2. 🛠️ They know which tools work with which file formats
3. 💻 They can write and explain command-line syntax
4. 📊 They can interpret analysis results
5. 🚀 They accelerate learning and productivity!

## 📚 File Format Quick Reference

### 🧬 FASTA
```
>sequence_name description
ATCGATCGATCGATCG...
```

### 📖 FASTQ
```
@read_name
ATCGATCGATCG
+
IIIIIIIIIII
```

### 🔴 VCF
```
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO
chr1    100  .   A    G    30    PASS    DP=50
```

## 🙏 Acknowledgments

- 🤖 Created with assistance from **Claude** (Anthropic)
- 🧬 Example data for educational purposes only
- 🎓 Part of genomics + AI educational materials

---

⭐ **Star this repo if you found it helpful!** ⭐

🐛 Found an issue? Open a ticket! 🎫

📧 Questions? Let's discuss! 💬

---

*Made with 💜 and lots of ☕ by bioinformaticians who love LLMs! 🤖🧬*

🔬🧪🔭🌡️🧫🦠🧬🔗💻🖥️⌨️🖱️📊📈📉🎯🎓📚✨🌟💫⚡🚀🔥💯
