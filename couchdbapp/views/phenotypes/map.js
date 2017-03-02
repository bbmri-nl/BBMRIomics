function(doc) {
    if (doc.hasOwnProperty('phenotype'))
        emit(doc.biobank, doc.phenotype);
}
