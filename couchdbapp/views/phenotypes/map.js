/**
 * Map function - use `emit(key, value) to generate rows in the output result.
 * @link http://docs.couchdb.org/en/latest/couchapp/ddocs.html#reduce-and-rereduce-functions
 *
 * @param {object} doc - Document Object.
 */
function(doc) {
    var doc = doc.doc;
    if (doc.hasOwnProperty('phenotype')) {
        emit(null, doc.phenotype);
    }
}
