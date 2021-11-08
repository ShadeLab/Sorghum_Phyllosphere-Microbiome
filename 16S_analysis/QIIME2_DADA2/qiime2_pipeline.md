
#!/bin/bash --login
########## Define Resources Needed with SBATCH Lines ##########

SBATCH --time=24:00:00               # limit of wall clock time - how long the job will run (same as -t)
SBATCH --ntasks=4
SBATCH --nodes=1-4
SBATCH --cpus-per-task=10             # number of CPUs (or cores) per task (same as -c)
SBATCH --mem=400G                    # memory required per node - amount of memory (in bytes)
SBATCH --job-name read_preprocess        # you can give your job a name for easier identification (same as -J)
SBATCH --mail-user=mechanll@msu.edu
SBATCH --mail-type=BEGIN,END
SBATCH --output=read_preprocess.log

#### set up/ modify paths as needed

cd /mnt/research/ShadeLab/WorkingSpace/MarcoMechan_WorkingSpace/
mkdir Mucilage_16SV4_resequencing/DADA2/

##### First, prepare the tab-delimited manifest file to import the sequence data.

conda activate qiime2-2021.8

##### Import files
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path Mucilage_16SV4_resequencing/manifest_file.tsv --output-path Mucilage_16SV4_resequencing/paired-end-demux.qza --input-format PairedEndFastqManifestPhred33V2

#### All sequence data is stored compressed in the paired-end-demux.qza file
#### You can create a visualization file with
qiime demux summarize --i-data Mucilage_16SV4_resequencing/paired-end-demux.qza --o-visualization Mucilage_16SV4_resequencing/demux.qzv

#### Next step is to run the DADA2 plug-in. 
#### use Zymo Research’s program FIGARO to find the trim parameters to optimize merging Forw and Rev reads. Remove as much of the lower quality portions of the reads as possible and still leaves enough overlap.

qiime dada2 denoise-paired --i-demultiplexed-seqs Mucilage_16SV4_resequencing/paired-end-demux.qza --p-trunc-len-f 237 --p-trunc-len-r 145 --o-table Mucilage_16SV4_resequencing/DADA2/table.qza --o-representative-sequences Mucilage_16SV4_resequencing/DADA2/rep-seqs.qza --o-denoising-stats Mucilage_16SV4_resequencing/DADA2/denoising-stats.qza

#### You can create a visualization file with
qiime metadata tabulate --m-input-file Mucilage_16SV4_resequencing/DADA2/denoising-stats.qza --o-visualization Mucilage_16SV4_resequencing/DADA2/denoising-stats.qzv

qiime feature-table summarize --i-table Mucilage_16SV4_resequencing/DADA2/table.qza --o-visualization Mucilage_16SV4_resequencing/DADA2/table.qzv --m-sample-metadata-file sample-metadata.tsv

qiime feature-table tabulate-seqs --i-data Mucilage_16SV4_resequencing/DADA2/rep-seqs.qza --o-visualization Mucilage_16SV4_resequencing/DADA2/rep-seqs.qzv

#### Align representative sequences and construct a phylogenetic tree
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences Mucilage_16SV4_resequencing/DADA2/rep-seqs.qza --o-alignment Mucilage_16SV4_resequencing/DADA2/aligned-rep-seqs.qza --o-masked-alignment Mucilage_16SV4_resequencing/DADA2/masked-aligned-rep-seqs.qza --o-tree Mucilage_16SV4_resequencing/DADA2/unrooted-tree.qza --o-rooted-tree Mucilage_16SV4_resequencing/DADA2/rooted-tree.qza

##### Taxonomy classifiers for use with q2-feature-classifier : Pre-trained classifiers can be found https://docs.qiime2.org/2021.8/data-resources/
qiime feature-classifier classify-sklearn --i-classifier silva-138-99-515-806-nb-classifier.qza --i-reads Mucilage_16SV4_resequencing/DADA2/rep-seqs.qza --o-classification Mucilage_16SV4_resequencing/DADA2/taxonomy.qza

qiime metadata tabulate --m-input-file Mucilage_16SV4_resequencing/DADA2/taxonomy.qza --o-visualization Mucilage_16SV4_resequencing/DADA2/taxonomy.qzv

mkdir Mucilage_16SV4_resequencing/DADA2/phyloseq

qiime tools export --input-path Mucilage_16SV4_resequencing/DADA2/rep-seqs.qza --output-path Mucilage_16SV4_resequencing/DADA2/phyloseq

qiime tools export --input-path Mucilage_16SV4_resequencing/DADA2/unrooted-tree.qza --output-path Mucilage_16SV4_resequencing/DADA2/phyloseq

cd Mucilage_16SV4_resequencing/DADA2/phyloseq

mv tree.nwk unrooted_tree.nwk
cd ../../..

qiime tools export --input-path Mucilage_16SV4_resequencing/DADA2/rooted-tree.qza --output-path Mucilage_16SV4_resequencing/DADA2/phyloseq

qiime tools export --input-path Mucilage_16SV4_resequencing/DADA2/table.qza --output-path Mucilage_16SV4_resequencing/DADA2/phyloseq

biom convert -i Mucilage_16SV4_resequencing/DADA2/phyloseq/feature-table.biom -o Mucilage_16SV4_resequencing/DADA2/phyloseq/otu_table.txt --to-tsv
