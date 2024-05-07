
#include "SuffixArray.hpp"
#include "omp.h"

int main() {
    //print hello from thread i with openmp
#pragma omp parallel
    {
        printf("Hello World... from thread = %d\n",
               omp_get_thread_num());
    }
    // std::string file_location = "../input/fasta_files/influenza_A.fna";
    std::string file_location = "/home/alex/Downloads/ncbi_dataset(1)/ncbi_dataset/data/GCF_016699485.2/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna";
    // std::string file_location = "/home/alex/Downloads/ncbi_dataset/ncbi_dataset/data/GCA_030704535.1/GCA_030704535.1_ASM3070453v1_genomic.fna";
    // std::string file_location = "/home/alex/Downloads/ncbi_dataset/ncbi_dataset/data/GCF_030704535.1/GCF_030704535.1_ASM3070453v1_genomic.fna";
    // std::string file_location = "../input/fasta_files/salmonella.fa";
    uint64_t block_size = 27;
    SuffixArray* sa = SuffixArray::from_fasta(file_location, block_size);
    sa->print_sa();
}
