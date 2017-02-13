function hasRNAseq(doc) {
    return (!doc.files.rnaseq || !doc.files.rnaseq.length); 
}
