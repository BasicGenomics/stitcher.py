# stitcher.py

_stitcher.py_ reconstructs molecules from BG data and outputs a .bam file which can be used for further analysis. 

## System Requirements for stitcher.py

_stitcher.py_  is a python3 script with dependencies:
```
pandas
numpy
pysam
joblib
pygtrie
portion
```
No further installation is needed.

## Usage

usage: stitcher.py [-h] [-i input] [-o output] [-g gtf] [-iso iso] [-jun jun]
                   [-t threads] [--single-end] [--cells cells]
                   [--contig contig] [-v]

**arguments:**
```
  -h, --help       show this help message and exit
  --i input        Input .bam file
  --o output       Output .bam file
  --g gtf          gtf file with gene information
  -iso iso, --isoform iso json file with isoform information
  -jun jun, --junction jun json file with exon-exon structure
  --t threads      Number of threads
  --cells cells    List of cell barcodes to stitch molecules (text file, one cell barcode per line).
  --contig contig  Restrict stitching to contig
  -v, --version    show program's version number and exit
```
## Input

_stitcher.py_ takes .bam file processed with the basic genomics pipeline together with a gtf file to reconstruct molecules for the genes in the gtf file.

As optional parameter, the user can specify the number of threads used for parallelization, the cells to process, and restrict the reconstruction to a given contig.

## Output 

_stitcher.py_ writes its results to a .bam file as the reads are being processed. Some of the fields have a slightly different interpretation than usual. The algorithm does not handle insertions at the moment, and remove those before stitching the molecule. The query name is in the format "cell:gene:umi". The D character in the CIGAR string indicates missing coverage. The MAPQ is always 255. 

In some cases the UMI reads contain conflicting information regarding the splicing of the molecule. This could either be due to differences in mapping or due to UMI collisions. In those cases, _stitcher.py_ prefer the unspliced option and writes two additional tags which contain the number of bases which have conflicting information and the intervals themselves. 

The .bam file contain many additional custom tags:

```
NR : Number of reads used to stitch.
ER : Number of reads covering an exon.
IR : Number of reads covering an intron.
TC : Number of 3' reads.
IC : Number of internal reads.
FC : Number of 5' reads.
BC : Cell barcode.
XT : Gene name or id.
UB : UMI.
NC : Conflict in the reconstruction, the number of conflicting bases.
IL : Conflict in the reconstruction, the intervals where there is a conflict. Written as start1,end1,start2,end2,...
```
