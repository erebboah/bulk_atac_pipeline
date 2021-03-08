### fastqc_step1.sh 
Takes in a file of sample names (`prefixes_all.txt`) and launches an array of jobs, one per sample. Adjust number of jobs with `--array=1-5`.

Makes a sample folder in the main directory (`/share/crsp/lab/seyedam/share/newATAC/`) that will store results.

### map_hg38_step2.sh
