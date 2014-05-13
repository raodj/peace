#include "SFFReader.h"

int
main() {
    SFFReader sff("/home/dmadhava/WindowsShare/temp/CFCN.sff");
    while (sff.getPendingReads() > 0) {
        sff.getNextRead();
    }
    std::ofstream fasta("CFCN_test.fasta");
    EST::dumpESTList(fasta);
    return 0;
}
