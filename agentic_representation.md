# Agentic Multi-View Exploration of a Large-Scale Metagenomic Reference Space

## 1. Project Summary

Metagenomic sequencing has generated an enormous and rapidly growing collection of microbial community data, but these data are difficult to explore at scale. A researcher interested in a new metagenomic sample may want to know not only which existing samples are similar, but also whether that similarity is consistent across different biological representations: nucleotide sequence, protein sequence, functional composition, and taxonomy.

This project proposes an **agentic multi-view metagenomic exploration system** built on an existing corpus of approximately **1.2 million SRA metagenomic accessions**. For these accessions, DNA-level and protein-level Sourmash sketches and FunProfiler functional profiles have already been generated, with YACHT taxonomic profiles available for a subset.

Rather than training a new large language model, the project will use an existing LLM through an agent/skill architecture. The LLM will serve as a research interface and workflow orchestrator, while deterministic computational tools perform similarity searches and quantitative analyses.

The central scientific question is:

> **Does jointly analyzing DNA, protein, and functional representations reveal meaningful relationships among metagenomes that are not captured by any single representation alone?**

The central computational question is:

> **Can an AI agent translate natural-language biological questions into reproducible multi-view retrieval and analysis workflows?**

The project therefore has two complementary components:

1. **Multi-view metagenomic retrieval and analysis** — the primary scientific contribution.
2. **Agentic research interface** — a practical mechanism for researchers to interrogate and reason over the resulting data.

The project deliberately avoids training a new foundation model or attempting exhaustive all-pairs comparison across 1.2 million samples. Instead, it builds on already-generated representations and uses indexed/neighborhood-based retrieval to make the project feasible within a short development period.

---

# 2. Motivation

A major challenge in modern metagenomics is no longer simply generating data. It is making sense of the enormous amount of existing data.

A researcher may have a metagenomic sample and ask questions such as:

* Which previously observed microbiomes resemble this sample?
* Are the most similar samples similar at the DNA level, protein level, or functional level?
* Can two communities be genetically different but functionally similar?
* Can two communities be genetically similar while differing substantially in predicted function?
* When different representations disagree, are the samples taxonomically different?
* Are these patterns reproducible across large numbers of samples?
* Can an AI system perform these analyses from a natural-language research question?

Existing approaches address portions of this problem, but they tend to emphasize one of three strategies:

1. **Similarity search within one or two representation spaces**
2. **Knowledge graphs containing curated biological relationships**
3. **LLM agents that orchestrate bioinformatics workflows or query biological databases**

This project investigates the intersection of these ideas while focusing on an empirical property of the metagenomic corpus itself: **relationships between multiple representations of the same samples.**

---

# 3. Existing Work and Competitive Landscape

A major component of this project is establishing what has already been done. The proposal therefore explicitly distinguishes the proposed system from several relevant approaches.

## 3.1 Microbiome Search Engine 2

**Microbiome Search Engine 2 (MSE2)** is the most important precedent for the retrieval component.

MSE2 provides whole-microbiome searches against a database containing more than 250,000 samples, using taxonomic or functional similarity. It uses indexed search to perform rapid retrieval and reports sub-second search against its database. Its database includes both shotgun metagenomes and 16S samples from hundreds of studies.

This means that the following idea is **not sufficiently novel by itself**:

> "Given a metagenome, find functionally similar metagenomes in a large database."

MSE2 already demonstrates that capability.

### Difference from MSE2

The proposed project instead investigates **joint relationships between multiple representations of the same metagenomic samples**:

```text
                 Sample
              /     |      \
             /      |       \
           DNA    Protein   Function
             \      |       /
              \     |      /
             Multi-view relationship
```

The question is not simply:

> "Which samples are similar?"

It is:

> "How does similarity change depending on what biological representation we use?"

For example, a pair of samples might exhibit:

```text
DNA similarity       LOW
Protein similarity   LOW
Functional similarity HIGH
```

or:

```text
DNA similarity       HIGH
Protein similarity   MODERATE
Functional similarity LOW
```

These are **discordant similarity relationships** that can be systematically investigated.

Importantly, the project will not assume that such discordance has a particular biological cause. It will treat it as an empirical observation and use available taxonomic information and follow-up analyses to investigate possible explanations.

---

# 4. Knowledge Graph Approaches

## 4.1 MetagenomicKG

**MetagenomicKG**, published in 2026, demonstrates that knowledge graphs can integrate taxonomic, functional, pathogenicity, and biomedical information for metagenomic analysis. It connects resources such as GTDB, KEGG, BV-BRC, and biomedical knowledge graphs, and demonstrates applications including hypothesis generation, sample-specific graph embeddings, and pathogen prediction.

Therefore, simply proposing:

> "Build a knowledge graph connecting microbes, functions, diseases, and pathways"

would not constitute a strong novel contribution.

### Difference from MetagenomicKG

MetagenomicKG primarily represents **known biological entities and relationships obtained from reference databases**.

The proposed project instead focuses on **relationships empirically measured across a very large collection of actual metagenomic samples**.

For example:

```text
MetagenomicKG:

Microbe
   ↓
KEGG function
   ↓
Pathway
   ↓
Disease
```

versus:

```text
Proposed system:

Sample A
   ↓
DNA similarity to Sample B = 0.24
Protein similarity = 0.38
Functional similarity = 0.91
Taxonomic similarity = ...
```

The latter relationship is not simply a predefined biological fact. It is a relationship **measured across the corpus**.

A graph representation may eventually be useful for navigating these relationships, but the graph itself is not the proposed scientific contribution.

---

# 5. KODA

KODA is particularly relevant to the proposed agentic component.

KODA is an agentic framework combining LLMs and a Neo4j microbiome knowledge graph. Its multi-agent architecture translates natural-language questions into graph queries and analyzes the resulting information, with a focus on KEGG orthologies, microbial metabolic pathways, and antimicrobial drug-target discovery.

This establishes an important precedent:

> An LLM can act as a research interface over structured microbiome information.

Therefore, the project will **not claim novelty simply from using ChatGPT, multiple agents, or natural-language querying**.

### Difference from KODA

KODA asks questions of a knowledge graph containing curated microbiological relationships.

The proposed system instead asks an agent to **compose computational analyses across multiple representations of empirical metagenomic observations**.

For example:

> "Find samples that are functionally similar to sample X but genetically dissimilar, then determine whether their taxonomic compositions differ."

The agent would need to perform a sequence such as:

```text
1. Retrieve functional neighbors of X
             ↓
2. Calculate DNA similarity to those candidates
             ↓
3. Filter for functional-high / DNA-low samples
             ↓
4. Retrieve protein similarity
             ↓
5. Retrieve YACHT taxonomy where available
             ↓
6. Compare taxonomic composition
             ↓
7. Summarize evidence and uncertainty
```

This is different from simply translating a natural-language question into a Cypher query.

---

# 6. MetaClaw and Agentic Bioinformatics

A particularly important recent development is **MetaClaw**, an auditable AI agent for end-to-end metagenomic and multi-omics analysis.

MetaClaw separates deterministic upstream workflows from configurable downstream skills and uses an LLM as an orchestration layer. Its 2026 revision evaluates the system on four published cohorts totaling 769 metagenomic profiles and emphasizes reproducibility, provenance, workflow registration, and auditable execution.

MetaClaw therefore establishes that:

> "An LLM agent that orchestrates metagenomic analysis workflows"

is already an active research direction.

The proposed project consequently should **not compete with MetaClaw on generic workflow automation**.

Instead, the proposed agent would be specialized around a particular scientific capability:

> **multi-view exploration of a precomputed metagenomic reference space.**

MetaClaw is primarily concerned with executing and reproducing analytical workflows on user datasets. The proposed system is primarily concerned with **searching and reasoning over a large precomputed population of metagenomic observations**.

This distinction is important.

---

# 7. GenomeOcean and Large-Scale Foundation Models

GenomeOcean demonstrates another very different approach to large-scale metagenomics.

GenomeOcean is a 4-billion-parameter genome foundation model trained on more than 600 Gbp of assembled metagenomic sequence derived from approximately 220 TB of datasets. Its goal is to learn a general representation of microbial genomic sequence and support generation and biological modeling.

The proposed project does **not** attempt to compete with this approach.

It will not:

* train a genome foundation model;
* generate genomic sequences;
* require billions of model parameters;
* reprocess hundreds of terabytes of sequence;
* learn a new biological embedding model.

Instead, the project exploits something already available:

> **1.2 million empirical metagenomic observations represented simultaneously in several existing computational spaces.**

This makes the proposed project substantially more feasible.

---

# 8. The Existing Dataset

The project starts with an unusually large precomputed corpus.

Approximately:

**1.2 million SRA accessions**

have already been processed into:

### DNA representation

Sourmash DNA sketches.

These provide a compact representation that can support sequence-level similarity searching without repeatedly comparing raw sequencing reads.

### Protein representation

Sourmash protein sketches.

These provide a second sequence-derived representation at the protein level.

### Functional representation

FunProfiler functional profiles.

These provide a representation of predicted functional composition.

### Taxonomic representation

YACHT taxonomic profiles are available for a subset of the samples.

These provide an opportunity to investigate whether observed multi-view relationships correspond to differences or similarities in taxonomic composition.

The project therefore has the following structure:

```text
1.2M metagenomic samples
        │
        ├── DNA sketches
        │
        ├── Protein sketches
        │
        ├── Functional profiles
        │
        └── Taxonomic profiles
              (subset)
```

---

# 9. Central Hypothesis

The central hypothesis is:

> **Metagenomic samples that appear similar in one representation may not necessarily be similar in another, and systematic analysis of these cross-representation relationships can reveal information that is lost when metagenomes are searched using only a single representation.**

A secondary hypothesis is:

> **An agent capable of composing retrieval and comparison operations across these representations can support biological questions that are difficult to express using a single similarity search.**

The first hypothesis is the scientific hypothesis.

The second is the computational/AI hypothesis.

Keeping these separate is important because the AI interface should not be mistaken for the biological discovery itself.

---

# 10. Proposed Scientific Analysis

The initial analysis will focus on **multi-view nearest-neighbor retrieval**.

For a query sample \(Q\), retrieve its nearest neighbors independently in:

```text
DNA space
Protein space
Functional space
```

For each representation, obtain a ranked list:

```text
DNA:
B, C, D, E, ...

Protein:
B, F, G, H, ...

Function:
B, F, H, J, ...
```

The overlap and disagreement between these neighborhoods become the primary objects of analysis.

---

# 11. Multi-View Similarity

Similarity must be defined according to the actual structure of each representation.

For Sourmash sketches, the appropriate similarity measure will depend on the sketch type and parameters used. Jaccard- or containment-based measures may be appropriate depending on how the sketches were generated.

For FunProfiler, the exact similarity measure will depend on the output format and whether the functional values represent counts, relative abundances, or another quantity. Candidate metrics may include Bray-Curtis or other appropriate profile distances.

The project will therefore first inspect the existing outputs and establish a justified metric for each representation rather than assuming that one metric is appropriate for everything.

The result for each sample pair can be represented as a multi-dimensional relationship:

```text
Sample A ↔ Sample B

DNA similarity       = ...
Protein similarity   = ...
Functional similarity = ...
Taxonomic similarity = ...
```

This produces a **multi-view similarity profile** rather than a single similarity score.

---

# 12. Discordant Neighborhoods

The most interesting cases are expected to be samples whose similarity rankings differ substantially between representations.

For example:

### Case A — functional convergence candidate

```text
DNA similarity        LOW
Protein similarity    LOW/MODERATE
Functional similarity HIGH
```

Possible interpretation:

Different microbial communities may perform similar functional roles.

However, this is only a hypothesis. The system would investigate whether taxonomic differences support this interpretation.

---

### Case B — sequence similarity but functional divergence candidate

```text
DNA similarity        HIGH
Protein similarity    HIGH
Functional similarity LOW
```

This could motivate investigation of differences in functional composition, gene content, or other properties.

Again, the system would report this as an observation requiring explanation rather than automatically assigning a biological mechanism.

---

### Case C — consistent similarity

```text
DNA        HIGH
Protein    HIGH
Function   HIGH
```

These samples represent conventional nearest neighbors where different views agree.

These cases provide an important baseline.

---

# 13. Taxonomic Validation

YACHT data provide an opportunity to investigate the biological context of these patterns.

For samples with taxonomic profiles, the system can ask:

> Do functionally similar but DNA-dissimilar samples have different taxonomic compositions?

Or:

> Do samples with similar DNA and protein representations but divergent functional profiles have particular taxonomic patterns?

Taxonomy therefore serves as a **validation/context layer**, rather than being assumed to be available for every sample.

---

# 14. The Agentic System

Once the multi-view retrieval functions are established, an LLM-based agent will provide a natural-language interface.

The system will not train a new LLM.

Instead, it will leverage an existing LLM through a skill/tool architecture.

The architecture will resemble:

```text
                       Researcher
                           │
                           ↓
                    Natural-language
                        question
                           │
                           ↓
                    Research Agent
                           │
             ┌─────────────┼─────────────┐
             ↓             ↓             ↓
       Retrieval Skill  Analysis Skill  Biology Skill
             │             │             │
             └─────────────┼─────────────┘
                           ↓
                     Tool layer
             ┌─────────────┼─────────────┐
             ↓             ↓             ↓
          Sourmash      FunProfiler      YACHT
```

---

# 15. Agent Roles

The system can emulate the useful conceptual structure of KODA without requiring several independent LLM systems.

## Research/Planning Agent

Interprets the user's question and determines what analysis is required.

Example:

> "Find samples that are functionally similar but genetically different."

The planner converts this into:

```text
functional retrieval
        ↓
DNA comparison
        ↓
filter
        ↓
optional protein/taxonomy analysis
```

---

## Data Agent

Executes deterministic retrieval functions.

Examples:

```text
find_function_neighbors()
find_dna_neighbors()
find_protein_neighbors()
get_taxonomy()
compare_samples()
```

The data agent should return structured results rather than biological speculation.

---

## Biological Interpretation Agent

Receives the quantitative results and explains what they could mean.

It should explicitly distinguish:

```text
Observation
    ↓
Evidence-supported interpretation
    ↓
Possible hypothesis
    ↓
What the current data cannot establish
```

This prevents the LLM from turning correlation into biological mechanism.

---

## Reporting Agent

Produces a reproducible research report containing:

* question;
* query samples;
* retrieval criteria;
* similarity metrics;
* retrieved samples;
* quantitative results;
* taxonomic evidence;
* interpretation;
* uncertainty;
* limitations;
* provenance of every computational result.

---

# 16. Why Use Skills Rather Than Building an Entire Agent Framework?

The project can borrow the useful conceptual design of KODA and MetaClaw without reproducing their entire infrastructure.

The LLM does not need to perform the biological computation itself.

Instead:

```text
LLM = reasoning/orchestration
Python = deterministic computation
Sourmash = sequence retrieval
FunProfiler = functional data
YACHT = taxonomy
```

This separation has an important scientific advantage:

**the numerical result does not depend on the LLM.**

The LLM decides which analysis to perform, but the similarity calculation itself remains deterministic and reproducible.

---

# 17. Evaluation

The project should evaluate both the scientific and agentic components.

## Evaluation 1 — Single-view versus multi-view retrieval

Compare:

```text
DNA-only
Protein-only
Function-only
```

against:

```text
Multi-view retrieval
```

The question is whether combining representations produces useful information that is missed by individual searches.

---

## Evaluation 2 — Neighborhood overlap

For each query:

```text
DNA top-K
Protein top-K
Function top-K
```

calculate:

* overlap;
* rank correlation;
* disagreement;
* representation-specific neighbors.

This establishes whether the representations are redundant or complementary.

---

## Evaluation 3 — Discordance analysis

Identify samples with strong disagreement between representations.

Determine:

* how common discordance is;
* whether it is reproducible;
* whether taxonomic profiles explain some of it;
* whether discordance differs across sample groups/habitats if suitable metadata are available.

---

## Evaluation 4 — Retrieval ablation

Compare:

```text
DNA only
DNA + protein
DNA + function
DNA + protein + function
```

This directly tests whether each representation adds information.

If adding protein or function does not improve retrieval or interpretation, that is an important negative result.

---

## Evaluation 5 — Agent evaluation

Give the agent a set of predefined biological questions.

For each question evaluate:

1. Did it choose the correct tools?
2. Did it retrieve the correct samples?
3. Did it use the correct similarity measure?
4. Did it correctly report the numerical results?
5. Did it distinguish observation from hypothesis?
6. Can another researcher reproduce the analysis from the generated record?

This prevents the project from becoming a subjective "the chatbot seemed useful" demonstration.

---

# 18. Example Research Questions

The agent could eventually support questions such as:

### Similarity

> "Find samples most similar to sample X."

### Cross-view similarity

> "Find samples that are functionally similar to X but genetically dissimilar."

### Representation disagreement

> "Which samples have the largest disagreement between DNA and functional neighborhoods?"

### Protein-level investigation

> "Are there samples with low DNA similarity but high protein similarity?"

### Taxonomic interpretation

> "For the functionally similar but genetically different samples, how different are their taxonomic profiles?"

### Comparative investigation

> "Show me cases where DNA similarity predicts functional similarity poorly."

### Population-level question

> "How common is functional similarity among genetically dissimilar metagenomes?"

These questions are more sophisticated than ordinary nearest-neighbor retrieval because they require the system to **compose multiple analyses**.

---

# 19. Computational Feasibility

A naive implementation would compare every pair of 1.2 million samples.

That would require approximately:

**720 billion pairwise comparisons.**

This is not an appropriate strategy for a short project.

The system will therefore use indexed retrieval, candidate generation, nearest-neighbor search, sampling, or other approximate strategies rather than exhaustive all-pairs comparison.

This is also consistent with the lesson from MSE2, which uses indexing specifically to avoid exhaustive searches and reports substantial acceleration over exhaustive comparison.

The project will begin on a manageable subset to validate the concept and then scale the successful approach.

---

# 20. One-Month Development Plan

The project is intentionally divided into a **minimum viable scientific result** and an optional agentic layer.

## Week 1 — Data and similarity

* Inspect DNA Sourmash sketch type and parameters.
* Inspect protein sketch type and parameters.
* Inspect FunProfiler output structure.
* Determine appropriate similarity/distance measures.
* Implement retrieval on a small subset.
* Validate that the results make sense.

### Deliverable

A reproducible pipeline:

```text
query sample
      ↓
DNA neighbors
protein neighbors
function neighbors
```

---

## Week 2 — Multi-view analysis

Implement:

* neighborhood overlap;
* cross-view rank comparisons;
* discordance detection;
* multi-view candidate retrieval.

Begin investigating YACHT profiles for selected cases.

### Deliverable

Evidence showing whether the three representations provide redundant or complementary information.

---

## Week 3 — Biological analysis

Focus on a small number of scientifically interpretable patterns.

For example:

```text
high function / low DNA
high DNA / low function
high agreement across all views
```

Use YACHT taxonomy where available to investigate these patterns.

### Deliverable

A preliminary biological result.

---

## Week 4 — Agent

Expose the validated computational functions as tools.

Implement a simple research-agent workflow:

```text
question
 ↓
plan
 ↓
retrieve
 ↓
compare
 ↓
interpret
 ↓
report
```

Add provenance and reproducibility.

### Deliverable

A working natural-language research interface demonstrating several predefined biological questions.

---

# 21. What Will NOT Be Built

To keep the project feasible, the following are explicitly outside the initial scope:

* training a new LLM;
* training a genome foundation model;
* building a new protein language model;
* reprocessing the 1.2 million raw SRA datasets;
* exhaustive all-vs-all comparison;
* constructing a giant external biological knowledge graph;
* reproducing the entire KODA architecture;
* reproducing the entire MetaClaw infrastructure;
* building a full production web application;
* claiming causal biological mechanisms from similarity alone.

These can become future directions if the initial experiment succeeds.

---

# 22. Expected Scientific Contribution

The strongest potential contribution is **not the chatbot**.

It is an empirical characterization of how metagenomic communities relate across multiple representations.

The project asks whether:

> **DNA similarity, protein similarity, and functional similarity provide complementary views of metagenomic relatedness at large scale.**

If they do, the project can characterize the structure of this disagreement and investigate its taxonomic context.

If they do not, that is also scientifically informative: it would suggest that the additional representations provide limited independent information for the particular retrieval task.

Either outcome provides a testable result.

---

# 23. Expected Computational Contribution

The computational contribution is an agentic interface that allows researchers to ask questions requiring multiple retrieval operations.

Rather than:

```text
Question → database query → answer
```

the system supports:

```text
Question
   ↓
Plan
   ↓
Multiple computational operations
   ↓
Evidence integration
   ↓
Biological interpretation
```

This makes the agent closer to a **research assistant** than a simple database chatbot.

---

# 24. How the Project Differs from Existing Systems

| System              | Primary object                            | Main capability                                                          | Difference from proposed system                                                                                         |
| ------------------- | ----------------------------------------- | ------------------------------------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------- |
| **MSE2**            | Large microbiome database                 | Taxonomic/functional similarity search                                   | Proposed system explicitly studies relationships across DNA, protein, and function                                      |
| **MetagenomicKG**   | Curated biological KG                     | Connects taxonomy, function, pathogens, disease and biomedical knowledge | Proposed system focuses on empirical relationships among millions of metagenomic observations                           |
| **KODA**            | Microbiome KG                             | Multi-agent natural-language querying and hypothesis generation          | Proposed agent composes multi-view computational retrieval rather than primarily querying a curated KG                  |
| **MetaClaw**        | Metagenomic workflows                     | Agentic workflow execution and auditable analysis                        | Proposed system focuses on exploration of a precomputed reference population rather than end-to-end raw-data processing |
| **GenomeOcean**     | Large-scale genomic sequence              | Foundation model for genomic representation/generation                   | Proposed system does not train a foundation model; it analyzes existing sample-level representations                    |
| **Proposed system** | 1.2M-sample multi-view metagenomic corpus | Cross-representation retrieval + agentic reasoning                       | Tests whether multiple empirical representations provide complementary information                                      |

The distinction should be stated carefully: the project is **not claiming to be categorically "better" than these systems**. They solve different problems. Its goal is to provide a capability that is not the primary focus of those systems.

---

# 25. Why the Combination Is Potentially Valuable

The individual components are not novel in isolation.

```text
Large metagenomic database     → existing
Similarity search              → existing
Functional profiling           → existing
Taxonomic profiling            → existing
Knowledge graphs                → existing
LLM agents                     → existing
Multi-agent reasoning          → existing
```

The potential contribution lies in their **specific combination and the scientific question it enables**:

```text
1.2M real metagenomic samples
          +
DNA representation
          +
protein representation
          +
functional representation
          +
taxonomic validation
          +
agentic orchestration
          ↓
Cross-representation metagenomic exploration
```

The proposal should therefore avoid claiming:

> "We invented AI agents for metagenomics."

or:

> "We built the first microbiome knowledge graph."

Those claims would be contradicted by existing work.

Instead:

> **We investigate and operationalize multi-representation exploration of a large empirical metagenomic reference space, using an AI agent to compose retrieval and analysis across representations.**

---

# 26. Longer-Term Vision

If the initial experiment demonstrates that multi-view analysis produces useful information, the system could eventually evolve into a much larger platform.

The long-term architecture could be:

```text
                         Researcher
                             │
                             ↓
                     AI Research Agent
                             │
             ┌───────────────┼────────────────┐
             ↓               ↓                ↓
       Sequence tools   Functional tools   Taxonomy
             │               │                │
             └───────────────┼────────────────┘
                             ↓
                  1.2M+ metagenomic corpus
                             │
                             ↓
                  Evidence / relationship layer
                             │
              ┌──────────────┼──────────────┐
              ↓              ↓              ↓
          Retrieval       Comparison     Hypothesis
                             │
                             ↓
                    Reproducible report
```

Additional data types could eventually be incorporated, such as metadata, MAGs, pathway information, metabolomics, transcriptomics, or other omics measurements if they become available.

A graph could also eventually be introduced if the analysis reveals that graph traversal provides a useful representation of the observed relationships.

---

# 27. Success Criteria

The project will be considered successful if it can demonstrate all three of the following:

### Scientific

There are measurable and reproducible differences between DNA-, protein-, and function-based neighborhoods.

### Analytical

Combining multiple representations enables at least one useful analysis that cannot be performed as naturally or effectively using a single representation.

### Agentic

An LLM agent can translate natural-language questions into the appropriate sequence of deterministic analyses and produce a reproducible evidence-backed report.

The agent itself is therefore evaluated as an interface to the scientific system, rather than as the scientific result.

---

# 28. Final Proposed Research Question

The project can ultimately be summarized by one question:

> **Can an agentic system use multiple complementary representations of a million-scale metagenomic corpus to discover and explain relationships that are invisible to single-view metagenomic similarity search?**

This framing preserves the part of KODA that is most exciting—the idea of an AI research agent—while avoiding simply recreating KODA.

It also leverages the major asset already available: the 1.2-million-sample corpus and its precomputed DNA, protein, functional, and taxonomic representations.

Most importantly, the project has a natural **fail-fast structure**:

```text
Do the representations contain complementary information?
                    │
             ┌──────┴──────┐
            NO             YES
             │              │
       useful result     investigate
       either way            │
                         build agent
                              │
                         demonstrate
                         scientific use
```

This makes the project feasible as a short-term prototype while leaving a clear path toward a substantially larger research platform.
