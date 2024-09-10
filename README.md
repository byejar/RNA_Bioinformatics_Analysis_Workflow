# RNA_pipline
预处理


HISAT2 samtools faturecount 
STAR                          （XenofilteR）DESEQ2     RRA

小鼠 移植瘤

cellline
KR

PDX

### GWCNA

# Bioinformatics Analysis Workflow

1. **RNA-seq Data Processing**:
   - Use **HISAT2** for alignment, followed by **Samtools** and **featureCounts**.
   - Alternatively, use **STAR** for alignment (no additional tools needed).

2. **Sample Type Determination**:
   - If sequencing **mouse-hosted human tumor cells**, apply **XenofilteR** to distinguish human and mouse reads.
   - If it's a **cell line**, skip this step.

3. **Differential Expression Analysis**:
   - Use **DESeq2** for differential expression analysis.

4. **Key Regulator Analysis**:
   - If there are multiple key regulators (KR), apply **RRA** (Robust Rank Aggregation) and **WGCNA** (Weighted Gene Co-expression Network Analysis).

5. **Gene Annotation**:
   - Annotate genes using a **GTF file**.

---

## Simplified Flowchart

```plaintext
       +----------------------------+
       |     RNA-seq Processing      |
       | HISAT2 -> Samtools -> Count |
       |      OR STAR Alignment      |
       +----------------------------+
                   |
                   v
       +-----------------------------+
       |  Tumor Type Determination    |
       |  Mouse-hosted Human Tumor?   |
       +-----------------------------+
           /       \          
          Yes       No
           |         |
+------------------+   v
| Apply XenofilteR |   +------------------------+
| Distinguish Reads|   |  DESeq2 Differential   |
+------------------+   |  Expression Analysis   |
                        +------------------------+
                                 |
                                 v
                      +-----------------------------+
                      | Multiple Key Regulators (KR)?|
                      +-----------------------------+
                          /        \
                        Yes        No
                          |          |
                +-----------------+   v
                | RRA & WGCNA      |   +------------------+
                | Network Analysis |   | Gene Annotation  |
                +-----------------+   | with GTF File    |
                                       +------------------+

