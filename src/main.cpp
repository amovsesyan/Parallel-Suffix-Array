
#include "SuffixArray.hpp"

int main() {
    // std::string file_location = "/home/alex/Programming/CMSC858N/Parallel-Suffix-Array/input/fasta_files/influenza_A.fna";
    // std::string file_location = "/home/alex/Downloads/ncbi_dataset(1)/ncbi_dataset/data/GCF_016699485.2/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna";
    // std::string file_location = "/home/alex/Downloads/ncbi_dataset/ncbi_dataset/data/GCA_030704535.1/GCA_030704535.1_ASM3070453v1_genomic.fna";
    // std::string file_location = "/home/alex/Downloads/ncbi_dataset/ncbi_dataset/data/GCF_030704535.1/GCF_030704535.1_ASM3070453v1_genomic.fna";
    std::string file_location = "/home/alex/Programming/CMSC858N/Parallel-Suffix-Array/input/fasta_files/salmonella.fa";
    // std::string file_location = "../human-genome/GRCh38_latest_genomic.fna";
    // for(index_t block_size = 1; block_size < 55; block_size++) {
    //     std::cout << "---------------------------------------- USING block_size: " << block_size << std::endl;
    //     SuffixArray* sa = SuffixArray::from_fasta(file_location, block_size);
    //     // sa->print_sa();
    //     delete sa;
    // }
    index_t block_size = 54;
    SuffixArray* sa = SuffixArray::from_fasta(file_location, block_size);
    delete sa;
    // sa->print_sa();
}
