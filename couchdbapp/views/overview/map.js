/**
 * Map function - use `emit(key, value)1 to generate rows in the output result.
 * @link http://docs.couchdb.org/en/latest/couchapp/ddocs.html#reduce-and-rereduce-functions
 *
 * @param {object} doc - Document Object.
 */
function(doc) {
    // !code vendor/helpers/phenotypes.js
    // !code vendor/helpers/rnaseq.js
    
    doc = doc.doc;

    emit(null, {
        id: doc._id,
        phenotype: hasPhenotype(doc),
        rnaseq: hasRNAseq(doc)
    });
}
