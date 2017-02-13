/**
 * List function - use `start()` and `send()` to output headers and content.
 * @link http://docs.couchdb.org/en/latest/couchapp/ddocs.html#listfun
 *
 * @param {object} head - View Head Information. http://docs.couchdb.org/en/latest/json-structure.html#view-head-info-object
 * @param {object} req - Request Object. http://docs.couchdb.org/en/latest/json-structure.html#request-object
 **/
function(head, req) {
    //adapted from: https://developer.ibm.com/clouddataservices/2015/09/22/export-cloudant-json-as-csv-rss-or-ical/
    var row, first = true;

    // output HTTP headers
    start({
        headers: {  'Content-Type': 'text/csv'  },
    });

    // iterate through the result set
    while(row = getRow()) {

        // get the doc (include_docs=true)		
        //var doc = row.doc;
        var doc = row.value;

        // if this is the first row
        if (first) {

            // output column headers
            send(Object.keys(doc).join(',') + '\n'); // for use in R newline '\n'
            first = false;
        }

        // build up a line of output
        var line = '';

        // iterate through each row
        for(var i in doc) {

            // comma separator
            if (line.length > 0) {
                line += ',';
            }

            // output the value, ensuring values that themselves
            // contain commas are enclosed in double quotes
            var val = doc[i];
            if (typeof val == 'string' && val.indexOf(',') >  -1) {
                line += '"' + val.replace(/"/g,'""') + '"';
            } else {
                line += val;
            }
        }
        line += '\n'; // for use in R newline '\n'

        // send  the line
        send(line);
    }
}
