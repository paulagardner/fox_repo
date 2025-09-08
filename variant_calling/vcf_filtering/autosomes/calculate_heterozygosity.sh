#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH --job-name=het_calc_all
#SBATCH --output=het_results/het_calc_all_%j.log
#SBATCH --time=2-00:00:00

module load bcftools/1.22

output_dir="het_results"
results_file="${output_dir}/heterozygosity_by_chrom.tsv"

mkdir -p "$output_dir"
echo -e "Sample\tChromosome\tHeterozygosity\tTotal_Sites\tJackknife_Mean\tJackknife_SE" > "$results_file"

# Get unique samples from VCF filenames
#basename strips directory path -> SAMPLE.chr1.vcf.gz
#cut -d'.' -f1 takes the part before the first dot -> SAMPLE
# sort -u: unique samples
samples=$(ls het_analysis_filtered/*.vcf.gz | xargs -n1 basename | cut -d'.' -f1 | sort -u)

for sample in $samples; do
    echo "Processing sample: $sample"

    total_hets=0
    total_sites=0
    jk_input_file=$(mktemp)

    #max_chroms=3 #commenting out testing code
    #chrom_count=0
    for vcf in het_analysis_filtered/${sample}.chr*.vcf.gz; do
        chr=$(basename "$vcf" | sed -E 's/.*\.chr(.*)\.vcf\.gz/\1/')
        nr_hets=$(bcftools view -g het -H "$vcf" | wc -l)
        nr_homs=$(bcftools view -g hom -H "$vcf" | wc -l)
        nr_sites=$((nr_hets + nr_homs))

        het_ratio=$(awk "BEGIN { printf \"%.5f\", $nr_hets / $nr_sites }")

        # Write per-chromosome results (no jackknife columns)
        echo -e "$sample\t$chr\t$het_ratio\t$nr_sites\t\t" >> "$results_file"

        # Collect for genome-wide jackknife
        echo -e "${het_ratio}\t${nr_sites}" >> "$jk_input_file"

        total_hets=$((total_hets + nr_hets))
        total_sites=$((total_sites + nr_sites))
        chrom_count=$((chrom_count + 1))
        #if [[ "$chrom_count" -ge "$max_chroms" ]]; then
        #    echo "Stopping after $max_chroms chromosomes for testing"
        #    break
        #fi
    done

    genome_het=$(awk "BEGIN { printf \"%.6f\", $total_hets / $total_sites }")

    # Genome-wide jackknife

    #always using perl to invoke this kind of script in anders' /sw bin will keep you
    #from permissions issues

    cat "$jk_input_file" | perl /gpfs/home/xrq24scu/bergstromlab/sw/bin/jackknife.pl > jk_tmp.txt
    gw_jk_mean=$(awk 'NR==2{print $1}' jk_tmp.txt)
    gw_jk_se=$(awk 'NR==2{print $2}' jk_tmp.txt)

    echo -e "$sample\tGenomeWide\t$genome_het\t$total_sites\t$gw_jk_mean\t$gw_jk_se" >> "$results_file"

    rm jk_tmp.txt "$jk_input_file"
    echo "Sample $sample processed: Genome-wide heterozygosity = $genome_het"
done

{ head -n 1 "$results_file"; grep -E "GenomeWide" "$results_file"; } > jackknife_het.tsv #pipe these results to another file

echo "All samples processed."