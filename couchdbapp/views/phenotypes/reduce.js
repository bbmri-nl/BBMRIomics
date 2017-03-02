/**
 * NOTE:
 * Built-in reduce functions should be preferred over writing a custom one.
 * Replace the contents of the file with one of these strings:
 *   _sum
 *   _count
 *   _stats
 **/

/**
 * Reduce function
 * @link http://docs.couchdb.org/en/latest/couchapp/ddocs.html#reduce-and-rereduce-functions
 *
 * @param {(array|null)} keys – Array of pairs docid-key for related map function
 * result. Always null if rereduce is running (has true value).
 * @param {array} values – Array of map function result values.
 * @param {boolean} rereduce – Boolean sign of rereduce run.
 *
 * @returns {object} Reduced values
 **/
function(keys, values, rereduce) {
  if (rereduce) {

  } else {

  }
}
