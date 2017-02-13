function hasPhenotype(doc) {  
    return  (!doc.phenotype || !doc.phenotype.length || !doc.phenotype.pheno_id);
}
