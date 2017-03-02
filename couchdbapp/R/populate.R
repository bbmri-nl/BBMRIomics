

curl -X GET -u 'usrpwd' -k https://metadatabase.bbmrirp3-lumc.surf-hosted.nl:6984/bios/_all_docs?include_docs=true > dump.json
sed -i -- 's/{"total_rows":6388,"offset":0,"rows"/{"docs"/g' dump.json

curl -d @dump.json -H "Content-type: application/json" -X POST http://127.0.0.1:5984/mdb_test/_bulk_docs

## find a way to validate e.g. use JSON-schema or use validation functions within couchdb


