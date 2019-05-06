#!/bin/bash


n_reads=$3
echo $n_reads
sambamba_path="/home/dario/bioinf/tools"
chip_lines=$($sambamba_path/sambamba-0.6.9 -q view --subsampling-seed=42 -t 4 -c $1)
control_lines=$($sambamba_path/sambamba-0.6.9 -q view --subsampling-seed=42 -t 4 -c $2)
echo $chip_lines $control_lines

# touch result_h3k4me1.txt

for i in {9..9}
  do
    control_s=$(bc -l <<< $n_reads*$i/$control_lines/10)
    if [[ $control_s = 0 ]]
      then 
        chip_s=$(bc -l <<< $n_reads/$chip_lines)
        $sambamba_path/sambamba-0.6.9 -q view --subsampling-seed=42 -f bam -o merged$i.bam -t 4 -s $chip_s $1
      elif [[ $i = 10 ]]
        then 
          $sambamba_path/sambamba-0.6.9 -q view --subsampling-seed=42 -f bam -o merged$i.bam -t 4 -s $control_s $2
      else
        $sambamba_path/sambamba-0.6.9 -q view --subsampling-seed=42 -f bam -o tmp_control.bam -t 4 -s $control_s $2
        control_lines_new=$($sambamba_path/sambamba-0.6.9 -q view -t 4 -c tmp_control.bam)
        chip_s=$(bc -l <<< $n_reads/$chip_lines-$control_lines_new/$chip_lines)
        $sambamba_path/sambamba-0.6.9 -q view --subsampling-seed=42 -f bam -o tmp_chip.bam -t 4 -s $chip_s $1
        $sambamba_path/sambamba-0.6.9 -q merge -t 4 merged$i.bam tmp_chip.bam tmp_control.bam
    fi
    
    echo $i
    #macs2 callpeak -t merged$i.bam -c $2 -n h3k27ac --broad -q 0.000001
    #python3 signal_to_noise_estimation.py merged$i.bam -d $4 --percentiles 0.1 0.9 0.2 0.9 0.1 0.99 >> result_h3k4me1.txt
    #bamCoverage -b merged$i.bam -o merged$i.bw -of bigwig
    #bedtools bamtobed -i merged$i.bam > merged$i.bed
    java -Xmx6G -jar $sambamba_path/span/span-0.11.0.build.jar analyze -t merged$i.bam -c $2 -f 0.05 --threads 4 --cs $sambamba_path/span/hg38.chrom.sizes -p $4_$i.peak
    #sh ../../tools/SICER_V1.1/SICER/SICER.sh . merged$i.bed gvg.bed . hg38 1 200 150 0.75 400 0.05
    rm -rf *.bam
    rm -rf *.bai
    #rm -rf merged$i.bed
    #rm -rf merged$i-1-removed.bed
    #rm -rf huv-1-removed.bed
    rm -rf logs/ 
    rm -rf cache/ 
    rm -rf fit/
  done 
